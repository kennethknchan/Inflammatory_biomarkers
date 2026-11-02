
config <- list(
  input = list(
    counts = "data/counts.csv",
    metadata = "data/metadata.csv",
    gene_sets = "data/gene_sets.gmt"
  ),
  
  output_dir = "results",
  
  analysis = list(
    group_column = "group",
    reference_group = "GROUP_A",
    comparison_group = "GROUP_B",
    fdr_cutoff = 0.05,
    logfc_cutoff = 1
  ),
  
  annotation = list(
    gene_id_type = "SYMBOL",
    go_ontology = "BP",
    organism_db = "org.Hs.eg.db"
  ),
  
  enrichment = list(
    min_gs_size = 10,
    max_gs_size = 500,
    pvalue_cutoff = 0.05
  ),
  
  gsva = list(
    min_size = 5,
    max_size = 500,
    kcdf = "Gaussian"
  )
)
