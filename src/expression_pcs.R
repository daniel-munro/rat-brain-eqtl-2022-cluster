suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
BED_FILE = args[1]
OUT_FILE = args[2]
N_PCS = as.integer(args[3])

df <- read_tsv(BED_FILE, col_types = cols(gene_id = "c", .default = "d")) %>%
    select(-`#chr`, -start, -end) %>%
    column_to_rownames(var = "gene_id")

pca <- prcomp(t(df), center = TRUE, scale = TRUE)
pcs <- round(pca$x[, 1:N_PCS], 6)
pcs <- as_tibble(t(pcs), rownames = "ID")
write_tsv(pcs, OUT_FILE)
