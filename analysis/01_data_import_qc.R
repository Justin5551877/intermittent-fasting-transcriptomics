# GSE196335: import and paired-sample QC.
# Run scripts 01--04 in order, from the repository root, in fresh R sessions.
# Install dependencies beforehand with BiocManager::install("DESeq2").
# Local intermediates: data/intermediate/*.rds (exclude *.rds in .gitignore).
suppressPackageStartupMessages(library(DESeq2))

main <- function() {
  metadata_file <- "data/metadata_GSE196335_main.csv"
  if (!file.exists(metadata_file)) stop("Run from the repo root; metadata file is missing.")
  meta <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("sample_id", "participant", "group", "time")
  if (!all(required %in% names(meta))) stop("Metadata needs: ", paste(required, collapse = ", "))
  meta[required] <- lapply(meta[required], function(x) trimws(as.character(x)))
  if (nrow(meta) != 56L || anyNA(meta[required]) || any(as.matrix(meta[required]) == "")) {
    stop("Expected 56 samples with complete sample_id, participant, group and time.")
  }
  if (anyDuplicated(meta$sample_id)) stop("Duplicate metadata sample IDs.")
  if (!all(meta$group %in% c("Control", "IF"))) stop("Groups must be Control and IF.")
  meta$time[meta$time == "6 Month"] <- "M6"
  if (!all(meta$time %in% c("Baseline", "M6"))) stop("Time must be Baseline and 6 Month (or M6).")
  meta$group <- factor(meta$group, levels = c("Control", "IF"))
  meta$time <- factor(meta$time, levels = c("Baseline", "M6"))
  meta$participant <- factor(meta$participant)
  if (nlevels(meta$participant) != 28L || any(table(meta$participant, meta$time) != 1L)) {
    stop("Expected 28 participants, each with exactly one Baseline and one M6 sample.")
  }
  if (any(vapply(split(meta$group, meta$participant), function(x) length(unique(x)), integer(1)) != 1L)) {
    stop("A participant is assigned to more than one treatment group.")
  }
  if (any(table(meta$group, meta$time) != 14L)) stop("Expected 14 participants per arm.")
  meta$IF_M6 <- as.integer(meta$group == "IF" & meta$time == "M6")
  rownames(meta) <- meta$sample_id
  design_matrix <- model.matrix(~ participant + time + IF_M6, data = meta)
  if (qr(design_matrix)$rank != ncol(design_matrix)) stop("Paired design is not full rank.")

  files <- list.files("data/gene_counts", pattern = "\\.txt\\.gz$", full.names = TRUE)
  if (length(files) != 56L) stop("Expected 56 .txt.gz count files; found ", length(files), ".")
  # GEO filenames begin with the GSM accession. Exact basenames are also supported.
  stems <- sub("\\.txt\\.gz$", "", basename(files))
  file_ids <- ifelse(grepl("^GSM[0-9]+(?:_|$)", stems, perl = TRUE),
                     sub("^(GSM[0-9]+).*$", "\\1", stems), stems)
  if (anyDuplicated(file_ids) || !setequal(file_ids, meta$sample_id)) {
    stop("Count filenames must identify metadata sample_id exactly or start with its GSM accession.")
  }
  files <- files[match(meta$sample_id, file_ids)]
  read_counts <- function(file) {
    con <- gzfile(file, open = "rt")
    on.exit(close(con))
    x <- read.delim(con, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    if (!all(c("Feature", "Count") %in% names(x))) stop("Missing Feature/Count in ", file)
    ids <- as.character(x$Feature)
    values <- suppressWarnings(as.numeric(as.character(x$Count)))
    if (!length(ids) || anyNA(ids) || any(ids == "") || anyDuplicated(ids)) stop("Invalid/duplicate feature IDs in ", file)
    if (any(!is.finite(values) | values < 0 | values != floor(values) | values > .Machine$integer.max)) {
      stop("Counts must be finite, nonnegative integers in ", file)
    }
    setNames(as.integer(values), ids)
  }
  count_list <- lapply(files, read_counts)
  ids <- names(count_list[[1]])
  if (!all(vapply(count_list, function(x) setequal(names(x), ids), logical(1)))) {
    stop("Count files contain different feature sets.")
  }
  counts <- vapply(count_list, function(x) unname(x[ids]), integer(length(ids)))
  dimnames(counts) <- list(ids, meta$sample_id)
  if (any(colSums(counts) == 0)) stop("A sample has an empty library.")
  # Strip Ensembl version suffixes; stop rather than silently merge ambiguous rows.
  clean_ids <- sub("\\.[0-9]+$", "", ids)
  if (anyDuplicated(clean_ids)) stop("Duplicate IDs after removing Ensembl version suffixes.")
  rownames(counts) <- clean_ids
  keep <- rowSums(counts >= 10L) >= 4L
  if (sum(keep) < 2L) stop("Too few genes pass filtering.")
  dds <- DESeqDataSetFromMatrix(counts[keep, , drop = FALSE], meta,
                                design = ~ participant + time + IF_M6)
  dds <- estimateSizeFactors(dds)
  # Blind VST for exploratory QC; DESeq2 inference uses the raw filtered counts.
  vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
  pca <- plotPCA(vsd, intgroup = c("group", "time", "participant"), returnData = TRUE)
  pca$sample_id <- rownames(pca)
  qc <- data.frame(sample_id = meta$sample_id, participant = meta$participant,
                   group = meta$group, time = meta$time,
                   library_size = colSums(counts), retained_library_size = colSums(counts[keep, , drop = FALSE]),
                   size_factor = sizeFactors(dds))
  dir.create("data/intermediate", recursive = TRUE, showWarnings = FALSE)
  dir.create("tables", showWarnings = FALSE)
  write.csv(qc, "tables/sample_QC.csv", row.names = FALSE)
  write.csv(pca, "tables/PCA_coordinates.csv", row.names = FALSE)
  write.csv(data.frame(features_imported = nrow(counts), genes_retained = sum(keep), samples = ncol(counts)),
            "tables/filtering_summary.csv", row.names = FALSE)
  saveRDS(list(metadata = meta, dds = dds, vsd = vsd, pca = pca,
               percent_var = attr(pca, "percentVar"), qc = qc,
               input_md5 = tools::md5sum(c(metadata_file, files))), "data/intermediate/01_import_qc.rds")
  writeLines(capture.output(sessionInfo()), "data/intermediate/01_session_info.txt")
  message("Imported 56 paired samples; retained ", sum(keep), " genes. Run script 02 next.")
}
main()
