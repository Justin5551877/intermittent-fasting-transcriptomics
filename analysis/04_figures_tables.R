# GSE196335: export Figures 2--6 as PNG/PDF and descriptive gene tables.
# Run from the repo root after scripts 01--03. No interactive workspace required.
# Dependencies: DESeq2, AnnotationDbi, org.Hs.eg.db, ggplot2, tidyr, dplyr, pheatmap.
suppressPackageStartupMessages({
  library(DESeq2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(pheatmap)
})

main <- function() {
  paths <- c("data/intermediate/01_import_qc.rds", "data/intermediate/02_deseq2.rds", "data/intermediate/03_gsea.rds")
  if (!all(file.exists(paths))) stop("Run analysis scripts 01--03 first.")
  imported <- readRDS(paths[1])
  fit <- readRDS(paths[2])
  enrichment <- readRDS(paths[3])
  if (!identical(fit$import_md5, unname(tools::md5sum(paths[1]))) ||
      !identical(enrichment$deseq_md5, unname(tools::md5sum(paths[2])))) {
    stop("Intermediates come from different runs. Rerun scripts 02 and 03.")
  }
  meta <- imported$metadata
  expr <- assay(imported$vsd)
  stopifnot(identical(colnames(expr), meta$sample_id), all(is.finite(expr)))
  dir.create("figures", showWarnings = FALSE)
  dir.create("tables", showWarnings = FALSE)
  palette <- c(Control = "#526D82", IF = "#D55E00")
  export <- function(plot, name, width, height) {
    for (ext in c("png", "pdf")) {
      ggsave(file.path("figures", paste0(name, ".", ext)), plot = plot,
             width = width, height = height, units = "in", dpi = 300, bg = "white")
    }
  }
  pca <- imported$pca
  pca <- pca[order(pca$participant, pca$time), ]
  fig2 <- ggplot(pca, aes(PC1, PC2, colour = group)) +
    geom_path(aes(group = participant), alpha = 0.35) +
    geom_point(aes(shape = time), size = 2.8) +
    scale_colour_manual(values = palette) + theme_minimal(base_size = 11) +
    labs(x = sprintf("PC1 (%.1f%%)", 100 * imported$percent_var[1]),
         y = sprintf("PC2 (%.1f%%)", 100 * imported$percent_var[2]),
         colour = "Group", shape = "Time", title = "Paired sample PCA",
         subtitle = "Blind VST; top 500 variable genes; lines connect paired samples")
  export(fig2, "Figure2_PCA_paired", 8, 6)

  focus <- enrichment$focus
  focus$Description <- factor(focus$Description, levels = unique(focus$Description[order(focus$NES)]))
  focus$log_fdr <- -log10(pmax(focus$p.adjust, .Machine$double.xmin))
  fig3 <- ggplot(focus, aes(NES, Description)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
    geom_point(aes(size = log_fdr, colour = NES, shape = significant_FDR_0.05)) +
    scale_colour_gradient2(low = "#2166AC", mid = "grey80", high = "#B2182B", midpoint = 0) +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), drop = FALSE) +
    scale_y_discrete(labels = function(x) vapply(x, function(s) paste(strwrap(s, 45), collapse = "\n"), character(1))) +
    theme_minimal(base_size = 11) +
    labs(x = "Normalized enrichment score (NES)", y = NULL, size = "-log10(FDR)",
         shape = "FDR < 0.05", title = "Focused GO Biological Process GSEA",
         subtitle = "Negative NES: more negative six-month change under IF",
         caption = "Exploratory selected pathways; FDR calculated across all tested GO BP terms.")
  export(fig3, "Figure3_focused_GSEA", 11, 8)

  jac <- enrichment$pathway_jaccard
  if (!nrow(jac)) stop("Pathway Jaccard matrix is empty.")
  cluster <- nrow(jac) > 1L && all(is.finite(jac))
  for (ext in c("png", "pdf")) {
    pheatmap(jac, cluster_rows = cluster, cluster_cols = cluster,
      color = colorRampPalette(c("white", "#2166AC"))(100), breaks = seq(0, 1, length.out = 101),
      display_numbers = TRUE, number_format = "%.2f", fontsize = 8,
      na_col = "grey85", main = "Jaccard similarity of leading-edge gene symbols",
      filename = file.path("figures", paste0("Figure4_pathway_Jaccard.", ext)), width = 12, height = 10)
  }

  genes <- c("ABCC8", "FOXA2", "GCG", "PPARGC1A", "COX4I1", "COQ9")
  gene_map <- AnnotationDbi::select(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", columns = "ENSEMBL")
  gene_map <- unique(gene_map[!is.na(gene_map$ENSEMBL) & gene_map$ENSEMBL %in% rownames(expr), ])
  missing <- setdiff(genes, gene_map$SYMBOL)
  if (length(missing)) stop("Representative genes missing after filtering/mapping: ", paste(missing, collapse = ", "))
  # If a symbol has multiple retained Ensembl IDs, choose highest mean normalized
  # count, then Ensembl ID. Export the choice; never duplicate participant rows.
  means <- rowMeans(counts(fit$dds, normalized = TRUE))
  gene_map$mean_normalized_count <- unname(means[gene_map$ENSEMBL])
  gene_map <- gene_map[order(gene_map$SYMBOL, -gene_map$mean_normalized_count, gene_map$ENSEMBL), ]
  gene_map <- gene_map[!duplicated(gene_map$SYMBOL), ]
  write.csv(gene_map, "tables/representative_gene_mapping.csv", row.names = FALSE)
  selected <- expr[gene_map$ENSEMBL, , drop = FALSE]
  rownames(selected) <- gene_map$SYMBOL
  long <- data.frame(SYMBOL = rownames(selected), selected, check.names = FALSE) |>
    tidyr::pivot_longer(cols = -SYMBOL, names_to = "sample_id", values_to = "expression") |>
    dplyr::left_join(meta[, c("sample_id", "participant", "group", "time")], by = "sample_id")
  stopifnot(nrow(long) == 336L, !anyNA(long), !anyDuplicated(long[, c("SYMBOL", "participant", "time")]))
  long$SYMBOL <- factor(long$SYMBOL, levels = genes)
  fig5 <- ggplot(long, aes(time, expression, group = participant, colour = group)) +
    geom_line(alpha = 0.20, linewidth = 0.5) + geom_point(alpha = 0.45, size = 1.4) +
    stat_summary(aes(group = group), fun = mean, geom = "line", linewidth = 1.3) +
    stat_summary(aes(group = group), fun = mean, geom = "point", size = 3) +
    facet_wrap(~ SYMBOL, scales = "free_y", ncol = 3) + scale_colour_manual(values = palette) +
    theme_minimal(base_size = 11) +
    labs(x = NULL, y = "VST expression", colour = "Group",
         title = "Paired expression changes in representative genes",
         caption = "Thin lines: individual participants. Thick lines: group means. Genes selected after pathway analysis.")
  export(fig5, "Figure5_representative_paired_genes", 10, 7)
  changes <- long |>
    dplyr::select(SYMBOL, participant, group, time, expression) |>
    tidyr::pivot_wider(names_from = time, values_from = expression) |>
    dplyr::mutate(delta = M6 - Baseline)
  summary <- changes |>
    dplyr::group_by(SYMBOL, group) |>
    dplyr::summarise(mean_delta = mean(delta), sd_delta = sd(delta), n = dplyr::n(), .groups = "drop")
  stopifnot(all(summary$n == 14L), all(is.finite(changes$delta)))
  effects <- summary |>
    tidyr::pivot_wider(names_from = group, values_from = c(mean_delta, sd_delta, n)) |>
    dplyr::mutate(differential_change = mean_delta_IF - mean_delta_Control,
      se_difference = sqrt(sd_delta_IF^2 / n_IF + sd_delta_Control^2 / n_Control),
      lower_95 = differential_change - 1.96 * se_difference,
      upper_95 = differential_change + 1.96 * se_difference)
  # Match the original descriptive normal-approximation intervals on VST scale.
  # These are not DESeq2 log2 fold changes, multiplicity-adjusted intervals, or
  # independent confirmation of gene significance after selecting these genes.
  fig6 <- ggplot(effects, aes(differential_change, reorder(SYMBOL, differential_change))) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
    geom_errorbar(aes(xmin = lower_95, xmax = upper_95), orientation = "y", width = 0.15) +
    geom_point(size = 3) + theme_minimal(base_size = 11) +
    labs(x = "Mean change in IF - mean change in Control (VST units)", y = NULL,
         title = "Differential six-month expression changes",
         subtitle = "Change = M6 - Baseline; negative values indicate more negative change under IF",
         caption = "Exploratory selected genes; approximate 95% intervals (estimate +/- 1.96 SE).\nIntervals are unadjusted for selection and multiple testing.")
  export(fig6, "Figure6_gene_differential_change", 9, 5.5)
  write.csv(long, "tables/representative_gene_expression.csv", row.names = FALSE)
  write.csv(changes, "tables/representative_gene_participant_changes.csv", row.names = FALSE)
  write.csv(summary, "tables/representative_gene_changes_main.csv", row.names = FALSE)
  write.csv(effects, "tables/representative_gene_differential_changes.csv", row.names = FALSE)
  writeLines(capture.output(sessionInfo()), "session_info.txt")
  message("Figures 2--6 (PNG/PDF), descriptive tables and session_info.txt exported.")
}
main()
