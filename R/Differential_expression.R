
run_edger <- function(counts, metadata, settings) {
  group <- factor(
    metadata[[settings$group_column]],
    levels = c(
      settings$reference_group,
      settings$comparison_group
    )
  )
  
  if (anyNA(group) || any(table(group) < 2L)) {
    stop("Both groups require at least two samples.")
  }
  
  design <- model.matrix(~ group)
  
  dge <- edgeR::DGEList(counts = counts)
  
  keep <- edgeR::filterByExpr(
    dge,
    design = design
  )
  
  if (!any(keep)) {
    stop("No genes passed expression filtering.")
  }
  
  dge <- dge[
    keep,
    ,
    keep.lib.sizes = FALSE
  ]
  
  dge <- edgeR::calcNormFactors(dge)
  dge <- edgeR::estimateDisp(dge, design)
  
  fit <- edgeR::glmQLFit(
    dge,
    design,
    robust = TRUE
  )
  
  test <- edgeR::glmQLFTest(
    fit,
    coef = 2
  )
  
  results <- edgeR::topTags(
    test,
    n = Inf,
    sort.by = "PValue"
  )$table
  
  results$gene_id <- rownames(results)
  
  results$significant <- (
    results$FDR < settings$fdr_cutoff &
      abs(results$logFC) >= settings$logfc_cutoff
  )
  
  list(
    dge = dge,
    design = design,
    fit = fit,
    test = test,
    results = results
  )
}
