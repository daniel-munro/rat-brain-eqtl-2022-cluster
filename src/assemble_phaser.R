#!/usr/bin/env Rscript

library(tools)
library(argparser)
library(data.table)

# Parse arguments
parser <- arg_parser(description = "Assemble reference and alternate counts tables for ANEVA-H")
parser <- add_argument(
  parser = parser,
  arg = '--input',
  help = "One or more input phASER Gene AE files; mutually exclusive with '-d|--directory'",
  default = NA_character_,
  nargs = Inf,
  short = '-i'
)
parser <- add_argument(
  parser = parser,
  arg = '--directory',
  help = "Directory with input phASER Gene AE files; mutually exclusive with '-i|--input'",
  default = NA_character_,
  short = '-d'
)
parser <- add_argument(
  parser = parser,
  arg = '--output',
  help = "Basename of output file",
  default = file.path(getwd(), 'phaser_assembled'),
  short = '-o'
)
if (!length(x = commandArgs(trailingOnly = TRUE))) {
  print(x = parser)
  stop("Not enough arguments provided")
}
args <- parse_args(parser = parser)

# Check arguments
if (!is.na(x = args$input) && !is.na(x = args$directory)) {
  stop("'-i|--input' and '-d|--directory' are mutually exclusive")
}
args$output <- file_path_sans_ext(x = args$output)
dir.create(path = dirname(path = args$output), showWarnings = FALSE, recursive = TRUE)

# Get input files
if (!is.na(x = args$directory)) {
  args$input <- list.files(
    path = args$directory,
    pattern = '\\.gene_ae\\.txt$',
    full.names = TRUE
  )
}
args$input <- args$input[file.exists(args$input)]
if (is.na(x = args$input) || !length(x = args$input)) {
  stop("No input phASER Gene AE files provided")
}

# Read in the data
nfiles <- length(x = args$input)
message("Working on ", nfiles, " input files")
message("Reading in data")
alt.counts <- ref.counts <- vector(mode = 'list')
for (i in seq_along(along.with = args$input)) {
  message("Reading file ", i, " of ", nfiles)
  fname <- args$input[[i]]
  sample <- gsub(
    pattern = '\\.gene_ae\\.txt$',
    replacement = '',
    x = basename(path = fname)
  )
  df <- fread(
    file = fname,
    sep = '\t',
    quote = FALSE,
    header = TRUE,
    stringsAsFactors = FALSE,
    showProgress = TRUE,
    verbose = FALSE,
    data.table = FALSE
  )
  if (!all(c('name', 'aCount', 'bCount') %in% colnames(x = df))) {
    warning(fname, " is malformed, skipping", immediate. = TRUE)
    next
  }
  genes <- unique(x = df$name)
  pb <- txtProgressBar(char = '=', style = 3, file = stderr())
  for (j in seq_along(along.with = genes)) {
    g <- genes[j]
    dg <- subset(x = df, subset = name == g)
    ref <- sum(dg[['aCount']])
    alt <- sum(dg[['bCount']])
    names(x = ref) <- names(x = alt) <- sample
    ref.counts[[g]] <- append(x = ref.counts[[g]], values = ref)
    alt.counts[[g]] <- append(x = alt.counts[[g]], values = alt)
    setTxtProgressBar(pb = pb, value = j / length(x = genes))
  }
  close(con = pb)
  # setTxtProgressBar(pb = pb, value = i / length(x = args$input))
}

# Assemble data frames
message("Assembling data frames")
CountsDF <- function(counts, gene) {
  df <- as.data.frame(x = counts[[gene]])
  colnames(x = df) <- gene
  return(df)
}
Merge2 <- function(x, y) {
  return(merge(x = x, y = y, by = 0, all = TRUE))
}

ref.counts <- lapply(
  X = names(x = ref.counts),
  FUN = CountsDF,
  counts = ref.counts
)

alt.counts <- lapply(
  X = names(x = alt.counts),
  FUN = CountsDF,
  counts = alt.counts
)

ref.counts <- t(x = do.call(what = 'cbind', args = ref.counts))
alt.counts <- t(x = do.call(what = 'cbind', args = alt.counts))

# Write the outputs
ref.name <- paste0(args$output, '_ref_counts.tsv')
alt.name <- paste0(args$output, '_alt_counts.tsv')

message("Writing reference counts to ", ref.name)
write.table(
  x = ref.counts,
  file = ref.name,
  quote = FALSE,
  sep = '\t',
  row.names = TRUE,
  col.names = TRUE
)

message("Writing alt counts to ", alt.name)
write.table(
  x = alt.counts,
  file = alt.name,
  quote = FALSE,
  sep = '\t',
  row.names = TRUE,
  col.names = TRUE
)
