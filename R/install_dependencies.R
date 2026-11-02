
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_packages <- c(
  "ggplot2",
  "pheatmap",
  "testthat"
)

bioc_packages <- c(
  "edgeR",
  "limma",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot",
  "GSVA",
  "GSEABase"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

BiocManager::install(
  bioc_packages,
  ask = FALSE,
  update = FALSE
)
