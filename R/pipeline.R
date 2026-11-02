
source("config/Config.R")

source("R/io.R")
source("R/differential_expression.R")
source("R/pathway_analysis.R")
source("R/visualisation.R")

dir.create(
  config$output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

output_path <- function(filename) {
  file.path(config$output_dir, filename)
}

# 1. Import and align data
counts <- read_counts(config$input$counts)
metadata <- read_metadata(config$input$metadata)

aligned <- align_samples(counts, metadata)
counts <- aligned$counts
metadata <- aligned$metadata

if (!config$analysis$group_column %in% names(metadata)) {
  stop("Configured group column is missing.")
}

# 2. Differential expression
de <- run_edger(
  counts,
  metadata,
  config$analysis
)

save_table(
  de$results,
  output_path("differential_expression.csv")
)

ggplot2::ggsave(
  output_path("volcano_plot.pdf"),
  plot = plot_volcano(
    de$results,
    config$analysis
  ),
  width = 8,
  height = 6
)

# 3. Gene Ontology
org_db <- getExportedValue(
  config$annotation$organism_db,
  config$annotation$organism_db
)

go_settings <- c(
  config$enrichment,
  list(
    go_ontology = config$annotation$go_ontology
  )
)

go_ora <- run_go_ora(
  de$results,
  config$annotation$gene_id_type,
  org_db,
  go_settings,
  config$analysis
)

go_gsea <- run_go_gsea(
  de$results,
  config$annotation$gene_id_type,
  org_db,
  go_settings
)

if (!is.null(go_ora)) {
  save_table(
    as.data.frame(go_ora),
    output_path("go_ora.csv")
  )
  
  plot <- plot_go_dotplot(
    go_ora,
    "GO Over-Representation Analysis"
  )
  
  if (!is.null(plot)) {
    ggplot2::ggsave(
      output_path("go_ora.pdf"),
      plot = plot,
      width = 9,
      height = 7
    )
  }
}

if (!is.null(go_gsea)) {
  save_table(
    as.data.frame(go_gsea),
    output_path("go_gsea.csv")
  )
  
  plot <- plot_go_dotplot(
    go_gsea,
    "GO Gene Set Enrichment Analysis"
  )
  
  if (!is.null(plot)) {
    ggplot2::ggsave(
      output_path("go_gsea.pdf"),
      plot = plot,
      width = 9,
      height = 7
    )
  }
}

# 4. GSVA
gene_sets <- read_gene_sets(
  config$input$gene_sets
)

gsva_results <- run_gsva(
  de$dge,
  gene_sets,
  config$gsva
)

save_table(
  data.frame(
    pathway = rownames(gsva_results$scores),
    gsva_results$scores,
    check.names = FALSE
  ),
  output_path("gsva_scores.csv")
)

gsva_comparison <- test_gsva_groups(
  gsva_results$scores,
  metadata,
  config$analysis
)

save_table(
  gsva_comparison,
  output_path("gsva_group_comparison.csv")
)

# 5. Reproducibility
writeLines(
  capture.output(sessionInfo()),
  output_path("session_info.txt")
)

message("Pipeline completed successfully.")
