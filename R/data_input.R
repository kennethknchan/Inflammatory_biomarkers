
read_counts <- function(path) {
  counts <- read.csv(
    path,
    row.names = 1,
    check.names = FALSE
  )
  
  counts <- as.matrix(counts)
  
  if (anyNA(counts) ||
      !is.numeric(counts) ||
      any(counts < 0) ||
      any(counts != round(counts))) {
    stop("Counts must be non-negative integer values.")
  }
  
  if (anyDuplicated(rownames(counts)) ||
      anyDuplicated(colnames(counts))) {
    stop("Duplicate gene or sample identifiers.")
  }
  
  counts
}

read_metadata <- function(path) {
  metadata <- read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  if (!"sample_id" %in% names(metadata)) {
    stop("Metadata requires a sample_id column.")
  }
  
  if (anyNA(metadata$sample_id) ||
      anyDuplicated(metadata$sample_id)) {
    stop("Sample identifiers must be unique and complete.")
  }
  
  metadata
}

align_samples <- function(counts, metadata) {
  if (!setequal(colnames(counts), metadata$sample_id)) {
    stop("Count and metadata sample IDs do not match.")
  }
  
  metadata <- metadata[
    match(colnames(counts), metadata$sample_id),
    ,
    drop = FALSE
  ]
  
  list(
    counts = counts,
    metadata = metadata
  )
}

save_table <- function(data, path) {
  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  write.csv(data, path, row.names = FALSE)
}
