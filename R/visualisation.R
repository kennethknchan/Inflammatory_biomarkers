
plot_volcano <- function(results, settings) {
  data <- results
  
  data$neg_log10_fdr <- -log10(
    pmax(data$FDR, .Machine$double.xmin)
  )
  
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = logFC,
      y = neg_log10_fdr
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(color = significant),
      alpha = 0.6,
      size = 1.5
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "FALSE" = "grey70",
        "TRUE" = "firebrick"
      )
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(settings$fdr_cutoff),
      linetype = "dashed"
    ) +
    ggplot2::geom_vline(
      xintercept = c(
        -settings$logfc_cutoff,
        settings$logfc_cutoff
      ),
      linetype = "dashed"
    ) +
    ggplot2::labs(
      title = "Differential Expression",
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value",
      color = "Significant"
    ) +
    ggplot2::theme_minimal()
}

plot_go_dotplot <- function(
    enrichment,
    title = "GO Enrichment",
    categories = 20) {
  
  if (is.null(enrichment) ||
      nrow(as.data.frame(enrichment)) == 0L) {
    return(NULL)
  }
  
  enrichplot::dotplot(
    enrichment,
    showCategory = categories
  ) +
    ggplot2::ggtitle(title)
}

plot_gsva_heatmap <- function(
    scores,
    metadata,
    group_column) {
  
  annotations <- data.frame(
    Group = metadata[[group_column]]
  )
  
  rownames(annotations) <- metadata$sample_id
  
  pheatmap::pheatmap(
    scores,
    scale = "row",
    annotation_col = annotations,
    show_rownames = nrow(scores) <= 50,
    main = "GSVA Pathway Activity",
    silent = TRUE
  )
}
