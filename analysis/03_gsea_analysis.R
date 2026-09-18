# GSE196335: GO Biological Process GSEA and leading-edge overlap.
# Run from the repo root after scripts 01 and 02.
# Dependencies: BiocManager::install(c("clusterProfiler", "org.Hs.eg.db",
#                                    "AnnotationDbi", "BiocParallel"))
# References: https://bioconductor.org/packages/clusterProfiler/
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(BiocParallel)
})

main <- function() {
  input <- "data/intermediate/02_deseq2.rds"
  if (!file.exists(input)) stop("Run analysis/02_deseq2_analysis.R first.")
  fit <- readRDS(input)
  res <- fit$results
  stopifnot(all(c("ENSEMBL", "stat") %in% names(res)))
  eligible <- res[is.finite(res$stat), c("ENSEMBL", "stat")]
  valid_keys <- intersect(eligible$ENSEMBL, keys(org.Hs.eg.db, keytype = "ENSEMBL"))
  if (!length(valid_keys)) stop("No Ensembl IDs map to the human annotation database.")
  mapping <- AnnotationDbi::select(org.Hs.eg.db, keys = valid_keys,
                                   keytype = "ENSEMBL", columns = "ENTREZID")
  mapping <- unique(mapping[!is.na(mapping$ENTREZID), ])
  # Avoid replicating one Ensembl statistic across multiple Entrez genes.
  # Ambiguous Ensembl mappings are excluded and reported in the audit table.
  multiplicity <- table(mapping$ENSEMBL)
  ambiguous <- names(multiplicity[multiplicity > 1L])
  unambiguous <- mapping[!mapping$ENSEMBL %in% ambiguous, ]
  ranked <- merge(eligible, unambiguous, by = "ENSEMBL")
  # Multiple Ensembl IDs per Entrez: retain the largest absolute Wald statistic.
  # Ties are resolved by Ensembl ID, deterministically.
  ranked <- ranked[order(-abs(ranked$stat), ranked$ENSEMBL), ]
  ranked <- ranked[!duplicated(ranked$ENTREZID), ]
  ranked <- ranked[order(-ranked$stat, ranked$ENTREZID), ]
  if (nrow(ranked) < 10L) stop("Too few uniquely mapped genes for GSEA.")
  gene_list <- setNames(ranked$stat, ranked$ENTREZID)
  dir.create("tables", showWarnings = FALSE)
  audit <- data.frame(ENSEMBL = eligible$ENSEMBL,
    status = ifelse(eligible$ENSEMBL %in% ranked$ENSEMBL, "retained",
      ifelse(eligible$ENSEMBL %in% ambiguous, "ambiguous_Ensembl_mapping",
        ifelse(eligible$ENSEMBL %in% unambiguous$ENSEMBL, "duplicate_Entrez", "unmapped"))))
  write.csv(audit, "tables/GSEA_mapping_audit.csv", row.names = FALSE)
  write.csv(ranked, "tables/GSEA_ranked_genes.csv", row.names = FALSE)
  # Rank all finite Wald statistics, without a DE significance filter.
  # Keep all tested GO terms; use global BH-adjusted P values for focused results.
  set.seed(196335)
  BiocParallel::register(BiocParallel::SerialParam())
  gsea <- gseGO(geneList = gene_list, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                ont = "BP", minGSSize = 10, maxGSSize = 500,
                pvalueCutoff = 1, pAdjustMethod = "BH", eps = 0,
                by = "fgsea", seed = TRUE, verbose = FALSE)
  all_results <- as.data.frame(gsea)
  if (!nrow(all_results)) stop("GO GSEA returned no tested terms.")
  focus_terms <- c("insulin secretion", "regulation of insulin secretion",
    "cellular response to glucose stimulus", "glucose homeostasis",
    "regulation of lipid metabolic process", "cholesterol efflux",
    "oxidative phosphorylation", "mitochondrial ATP synthesis coupled electron transport",
    "hormone secretion", "peptide hormone secretion", "neuropeptide signaling pathway")
  focus <- all_results[all_results$Description %in% focus_terms, ]
  focus <- focus[order(focus$NES), ]
  missing <- setdiff(focus_terms, focus$Description)
  if (length(missing)) warning("Focused terms absent from tested GO results: ", paste(missing, collapse = "; "))
  if (!nrow(focus)) stop("None of the specified focused pathways were tested.")
  focus$significant_FDR_0.05 <- !is.na(focus$p.adjust) & focus$p.adjust < 0.05
  # Focus is an exploratory, previously selected subset; retain nonsignificant terms.
  leading <- do.call(rbind, lapply(seq_len(nrow(focus)), function(i) {
    ids <- strsplit(as.character(focus$core_enrichment[i]), "/", fixed = TRUE)[[1]]
    ids <- unique(ids[!is.na(ids) & nzchar(ids)])
    if (!length(ids)) return(NULL)
    data.frame(ID = focus$ID[i], Description = focus$Description[i],
                NES = focus$NES[i], p.adjust = focus$p.adjust[i], ENTREZID = ids)
  }))
  if (is.null(leading) || !nrow(leading)) stop("No focused leading-edge genes.")
  symbols <- mapIds(org.Hs.eg.db, keys = unique(leading$ENTREZID),
                    keytype = "ENTREZID", column = "SYMBOL", multiVals = "first")
  leading$SYMBOL <- unname(symbols[leading$ENTREZID])
  # Match the original symbol-based overlap analysis; missing symbols are excluded.
  if (anyNA(leading$SYMBOL)) warning("Unmapped leading-edge symbols excluded from Jaccard analysis.")
  usable <- leading[!is.na(leading$SYMBOL) & nzchar(leading$SYMBOL), ]
  gene_sets <- setNames(lapply(focus$Description, function(term) unique(usable$SYMBOL[usable$Description == term])), focus$Description)
  modules <- list(
    Endocrine = c("insulin secretion", "regulation of insulin secretion", "hormone secretion", "peptide hormone secretion"),
    Glucose = c("cellular response to glucose stimulus", "glucose homeostasis"),
    Mitochondrial = c("oxidative phosphorylation", "mitochondrial ATP synthesis coupled electron transport"))
  module_genes <- lapply(modules, function(terms) unique(usable$SYMBOL[usable$Description %in% terms]))
  jaccard_matrix <- function(sets) {
    result <- matrix(NA_real_, length(sets), length(sets), dimnames = list(names(sets), names(sets)))
    for (i in seq_along(sets)) for (j in seq_along(sets)) {
      union_size <- length(union(sets[[i]], sets[[j]]))
      if (union_size > 0L) result[i, j] <- length(intersect(sets[[i]], sets[[j]])) / union_size
    }
    result
  }
  pathway_jaccard <- jaccard_matrix(gene_sets)
  module_jaccard <- jaccard_matrix(module_genes)
  write.csv(all_results, "tables/GO_GSEA_all_results.csv", row.names = FALSE)
  write.csv(focus, "tables/focused_metabolic_GSEA_clean.csv", row.names = FALSE)
  write.csv(leading, "tables/focused_pathway_leading_edge_genes.csv", row.names = FALSE)
  write.csv(pathway_jaccard, "tables/pathway_Jaccard_matrix.csv", row.names = TRUE)
  write.csv(module_jaccard, "tables/module_Jaccard_matrix.csv", row.names = TRUE)
  saveRDS(list(gsea = gsea, all_results = all_results, focus = focus, leading_edge = leading,
               gene_sets = gene_sets, module_genes = module_genes,
               pathway_jaccard = pathway_jaccard, module_jaccard = module_jaccard,
               deseq_md5 = unname(tools::md5sum(input))), "data/intermediate/03_gsea.rds")
  writeLines(capture.output(sessionInfo()), "data/intermediate/03_session_info.txt")
  message("GSEA and overlap tables saved. Run script 04 next.")
}
main()
