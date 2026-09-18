# GSE196335: paired difference-in-differences inference on raw counts.
# Run from the repo root after 01_data_import_qc.R; no workspace objects required.
suppressPackageStartupMessages(library(DESeq2))

main <- function() {
  input <- "data/intermediate/01_import_qc.rds"
  if (!file.exists(input)) stop("Run analysis/01_data_import_qc.R first.")
  imported <- readRDS(input)
  dds <- imported$dds
  stopifnot(inherits(dds, "DESeqDataSet"), ncol(dds) == 56L,
            identical(colnames(dds), rownames(imported$metadata)))
  design(dds) <- ~ participant + time + IF_M6
  dds <- DESeq(dds, test = "Wald", parallel = FALSE)
  if (!"IF_M6" %in% resultsNames(dds)) stop("Expected IF_M6 coefficient was not fitted.")
  # Positive = greater M6-minus-Baseline change in IF than Control.
  # The participant term absorbs baseline group differences; do not add group.
  result <- results(dds, name = "IF_M6", alpha = 0.05)
  result_table <- data.frame(ENSEMBL = rownames(result), as.data.frame(result), row.names = NULL)
  result_table <- result_table[order(result_table$padj, na.last = TRUE), ]
  dir.create("tables", showWarnings = FALSE)
  write.csv(result_table, "tables/DESeq2_IF_vs_control_change.csv", row.names = FALSE)
  saveRDS(list(dds = dds, results = result_table, coefficient = "IF_M6",
               import_md5 = unname(tools::md5sum(input))), "data/intermediate/02_deseq2.rds")
  writeLines(capture.output(sessionInfo()), "data/intermediate/02_session_info.txt")
  message(sum(result_table$padj < 0.05, na.rm = TRUE), " genes with BH-adjusted P < 0.05. Run script 03 next.")
}
main()
