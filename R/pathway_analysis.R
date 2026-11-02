
map_gene_ids <- function(gene_ids, from_type, org_db) {
  mapping <- clusterProfiler::bitr(
    unique(gene_ids),
    fromType = from_type,
    toType = "ENTREZID",
    OrgDb = org_db
  )
  
  mapping <- mapping[
    !is.na(mapping$ENTREZID),
    ,
    drop = FALSE
  ]
  
  mapping
}

run_go_ora <- function(
    results,
    from_type,
    org_db,
    settings,
    analysis_settings) {
  
  mapping <- map_gene_ids(
    results$gene_id,
    from_type,
    org_db
  )
  
  selected <- results$gene_id[
    !is.na(results$FDR) &
      results$FDR < analysis_settings$fdr_cutoff &
      abs(results$logFC) >=
      analysis_settings$logfc_cutoff
  ]
  
  foreground <- unique(
    mapping$ENTREZID[
      mapping[[from_type]] %in% selected
    ]
  )
  
  background <- unique(mapping$ENTREZID)
  
  if (length(foreground) == 0L) {
    warning("No mapped significant genes for GO ORA.")
    return(NULL)
  }
  
  clusterProfiler::enrichGO(
    gene = foreground,
    universe = background,
    OrgDb = org_db,
    keyType = "ENTREZID",
    ont = settings$go_ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = settings$pvalue_cutoff,
    qvalueCutoff = settings$pvalue_cutoff,
    minGSSize = settings$min_gs_size,
    maxGSSize = settings$max_gs_size,
    readable = TRUE
  )
}

run_go_gsea <- function(
    results,
    from_type,
    org_db,
    settings) {
  
  mapping <- map_gene_ids(
    results$gene_id,
    from_type,
    org_db
  )
  
  ranked <- merge(
    results[, c("gene_id", "F")],
    mapping,
    by.x = "gene_id",
    by.y = from_type
  )
  
  # Signed square root of the quasi-likelihood F statistic.
  # Appropriate here for a single-coefficient comparison.
  fold_changes <- results[
    match(ranked$gene_id, results$gene_id),
    "logFC"
  ]
  
  ranked$score <- sign(fold_changes) *
    sqrt(pmax(ranked$F, 0))
  
  ranked <- ranked[
    is.finite(ranked$score),
    ,
    drop = FALSE
  ]
  
  ranked <- ranked[
    order(abs(ranked$score), decreasing = TRUE),
    ,
    drop = FALSE
  ]
  
  ranked <- ranked[
    !duplicated(ranked$ENTREZID),
    ,
    drop = FALSE
  ]
  
  gene_list <- ranked$score
  names(gene_list) <- ranked$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  if (length(gene_list) < settings$min_gs_size) {
    warning("Insufficient mapped genes for GSEA.")
    return(NULL)
  }
  
  clusterProfiler::gseGO(
    geneList = gene_list,
    OrgDb = org_db,
    keyType = "ENTREZID",
    ont = settings$go_ontology,
    minGSSize = settings$min_gs_size,
    maxGSSize = settings$max_gs_size,
    pAdjustMethod = "BH",
    pvalueCutoff = settings$pvalue_cutoff,
    verbose = FALSE
  )
}


read_gene_sets <- function(path) {
  if (!file.exists(path)) {
    stop("Gene-set GMT file not found.")
  }
  
  gene_sets <- GSEABase::getGmt(path)
  
  sets <- lapply(
    gene_sets,
    GSEABase::geneIds
  )
  
  names(sets) <- vapply(
    gene_sets,
    GSEABase::setName,
    character(1)
  )
  
  sets
}

run_gsva <- function(
    dge,
    gene_sets,
    settings) {
  
  expression <- edgeR::cpm(
    dge,
    log = TRUE,
    prior.count = 2
  )
  
  if (anyDuplicated(rownames(expression))) {
    stop("Expression gene IDs must be unique.")
  }
  
  if (length(intersect(
    rownames(expression),
    unique(unlist(gene_sets))
  )) == 0L) {
    stop(
      "No matching identifiers between ",
      "expression matrix and gene sets."
    )
  }
  
  parameters <- GSVA::gsvaParam(
    exprData = expression,
    geneSets = gene_sets,
    minSize = settings$min_size,
    maxSize = settings$max_size,
    kcdf = settings$kcdf
  )
  
  scores <- GSVA::gsva(
    parameters,
    verbose = FALSE
  )
  
  list(
    scores = scores,
    expression = expression
  )
}

test_gsva_groups <- function(
    scores,
    metadata,
    settings) {
  
  group <- factor(
    metadata[[settings$group_column]],
    levels = c(
      settings$reference_group,
      settings$comparison_group
    )
  )
  
  design <- model.matrix(~ group)
  
  fit <- limma::lmFit(scores, design)
  fit <- limma::eBayes(fit)
  
  results <- limma::topTable(
    fit,
    coef = 2,
    number = Inf,
    adjust.method = "BH"
  )
  
  results$pathway <- rownames(results)
  
  results
}
