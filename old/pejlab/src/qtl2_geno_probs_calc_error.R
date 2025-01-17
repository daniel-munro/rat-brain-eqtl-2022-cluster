suppressPackageStartupMessages(library(tidyverse))
library(qtl2)

# CONTROL <- "data/qtl2/chr/chr12.yaml"
CONTROL <- "data/qtl2/Acbc.yaml"

error_probs <- function(geno, fgeno, pr) {
    ideal_geno <- matrix(nrow = nrow(geno), ncol = ncol(geno), dimnames = dimnames(geno))
    for (r in 1:nrow(ideal_geno)) {
        top_pairs_i <- apply(pr[r, , ], 2, function(x) which.max(rank(x, ties.method = "random")))
        top_pairs <- dimnames(pr)[[2]][top_pairs_i]
        top_strain1 <- str_sub(top_pairs, 1, 1)
        top_strain2 <- str_sub(top_pairs, 2, 2)
        for (col in 1:ncol(ideal_geno)) {
            ideal_geno[r, col] = fgeno[top_strain1[col], col] + fgeno[top_strain2[col], col]
        }
    }
    tibble(individual = rownames(geno)) %>%
        rowwise() %>%
        mutate(error_rate = mean(geno[individual, ] != ideal_geno[individual, ], na.rm = TRUE))
}

cross <- read_cross2(CONTROL)
pr <- calc_genoprob(cross, error_prob = 0.01, cores = 8)

errors <- tibble(chrom = names(pr)) %>%
    group_by(chrom) %>%
    summarise({
        cat(chrom, "\n")
        geno <- cross$geno[[chrom]] - 1L
        fgeno <- (cross$founder_geno[[chrom]] - 1) / 2
        rownames(fgeno) <- LETTERS[1:8]  # A-H like in pr, rather than actual strain names
        error_probs(geno, fgeno, pr[[chrom]])
    }, .groups = "drop")
