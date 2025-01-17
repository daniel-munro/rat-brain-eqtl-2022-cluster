suppressPackageStartupMessages(library(tidyverse))
library(qvalue)

args <- commandArgs(trailingOnly = TRUE)
PAIR_FILE = args[1]
PERM_FILE = args[2]
COEF_FILE = args[3]
OUT_FILE = args[4]
# PAIR_FILE = "data/qtl2/Acbc_gene_var_lod.tsv.gz"
# PERM_FILE = "data/qtl2/Acbc_gene_lod_perm.tsv.gz"
# COEF_FILE = "data/qtl2/Acbc_gene_var_coefs.tsv.gz"
# OUT_FILE = "data/qtl2/Acbc_gene_var_pval.tsv.gz"

perm <- read_tsv(PERM_FILE, col_types = cols(.default = "d"))

coefs <- read_tsv(COEF_FILE, col_types = cols(gene = "c", snp = "c", .default = "d"))

pairs <- read_tsv(PAIR_FILE, col_types = "ccd") %>%
    arrange(gene, desc(lod)) %>%
    group_by(gene) %>%
    slice(1) %>%
    ## NA for genes with no expression variance:
    mutate(pval = if (is.finite(lod)) 1 - ecdf(perm[[gene]])(lod) else NA) %>%
    ungroup() %>%
    mutate(qval = qvalue(pval)$qvalue) %>%
    left_join(coefs, by = c("gene", "snp"))

pairs %>%
    mutate(across(pval:intercept, ~signif(.x, 6))) %>%
    write_tsv(OUT_FILE)
