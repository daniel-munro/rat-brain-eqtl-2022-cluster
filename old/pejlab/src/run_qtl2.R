# Based on https://kbroman.org/qtl2/assets/vignettes/user_guide.html
library(argparser)
library(qtl2)
suppressPackageStartupMessages(library(tidyverse))

# args <- commandArgs(trailingOnly = TRUE)
# control_file <- args[1]
# tss_file <- args[2]
# outfile <- args[3]
# run_mode <- args[4]
p <- arg_parser("Run R/qtl2 with options for several types of output.")
p <- add_argument(p, "control", help="name of qtl2 control file")
p <- add_argument(p, "tss", help="name of file containing gene, chr, tss")
p <- add_argument(p, "out", help="name of output file")
p <- add_argument(p, "mode", help="nom, perm, or coef")
p <- add_argument(p, "--n_perms", help="for perm mode, number of permutations per gene", default=1000)
p <- add_argument(p, "--error_prob", help="error probability for HMM", default=0.01)
argv <- parse_args(p)

# control_file <- "data/qtl2/Acbc.yaml"
# tss_file <- "data/qtl2/gene_tss.tsv"
# outfile <- "data/qtl2/Acbc_gene_var_lod.tsv.gz"
# run_mode <- "coef" # nom, perm, or coef.

CIS_WINDOW <- 1e6
N_CORES <- 8

# pr_subset <- function(chrom, ids, perm = FALSE) {
#         pr_gene <- pr
#         # scan1perm allows no SNPs in other chromosomes, but scan1 doesn't.
#         other_chrom_index <- if (perm) c() else 1
#         for (chr in names(pr)) {
#             if (chr == chrom) {
#                 pr_gene[[chr]] <- pr_gene[[chr]][, , ids, drop = FALSE]
#             } else {
#                 # scan1 seems to require nonempty arrays. These are dropped later.
#                 pr_gene[[chr]] <- pr_gene[[chr]][, , other_chrom_index, drop = FALSE]
#             }
#         }
#         pr_gene
# }

scan_gene <- function(gene_id) {
    loc <- filter(tss, gene == gene_id)
    stopifnot(nrow(loc) == 1)
    in_window <- snps %>%
        filter(chr == loc$chr) %>%
        filter(abs(pos - loc$tss) <= CIS_WINDOW)
    if (nrow(in_window) > 0) {
        pr_gene <- pr[, loc$chr]
        pr_gene[[loc$chr]] <- pr_gene[[loc$chr]][, , in_window$id, drop = FALSE]
        pheno_gene <- cross$pheno[, gene_id, drop = FALSE]
        out <- scan1(pr_gene, pheno_gene, kinship[loc$chr], cores = N_CORES)
        out %>%
            as_tibble(rownames = "snp") %>%
            rename(lod = {{ gene_id }}) %>%
            mutate(lod = as.double(lod))
    } else {
        tibble(snp = character(0), lod = double(0))
    }
}

scan_gene_perm <- function(gene_id) {
    loc <- filter(tss, gene == gene_id)
    stopifnot(nrow(loc) == 1)
    in_window <- snps %>%
        filter(chr == loc$chr) %>%
        filter(abs(pos - loc$tss) <= CIS_WINDOW)
    if (nrow(in_window) > 0) {
        pr_gene <- pr[, loc$chr]
        pr_gene[[loc$chr]] <- pr_gene[[loc$chr]][, , in_window$id, drop = FALSE]
        pheno_gene <- cross$pheno[, gene_id, drop = FALSE]
        out <- scan1perm(pr_gene, pheno_gene, kinship[loc$chr], n_perm = argv$n_perms, cores = N_CORES)
        tibble(iter = 1:nrow(out),
               threshold = out[, gene_id])
    } else {
        tibble(iter = integer(0), threshold = double(0))
    }
}

scan_gene_coef <- function(gene_id) {
    loc <- filter(tss, gene == gene_id)
    stopifnot(nrow(loc) == 1)
    in_window <- snps %>%
        filter(chr == loc$chr) %>%
        filter(abs(pos - loc$tss) <= CIS_WINDOW)
    if (nrow(in_window) > 0) {
        apr_gene <- apr[, loc$chr]
        apr_gene[[loc$chr]] <- apr_gene[[loc$chr]][, , in_window$id, drop = FALSE]
        pheno_gene <- cross$pheno[, gene_id, drop = FALSE]
        out <- scan1coef(apr_gene, pheno_gene, kinship[[loc$chr]], cores = N_CORES)
        out %>%
            as_tibble(rownames = "snp") %>%
            mutate(across(-snp, as.double))
    } else {
        # set_names(rep(0, 37), c("AA" , "AB" , "BB", "AC", "BC", "CC", "AD", "BD", "CD",
        # "DD", "AE", "BE", "CE", "DE", "EE", "AF", "BF", "CF", "DF", "EF", "FF", "AG", "BG", "CG",
        # "DG", "EG", "FG", "GG", "AH", "BH", "CH", "DH", "EH", "FH", "GH", "HH", "intercept")) %>%
        # bind_rows() %>%
        # add_column(snp = "", .before = 1) %>%
        tibble(snp = "", A = 0, B = 0, C = 0, D = 0, E = 0, F = 0, G = 0, H = 0, intercept = 0) %>%
        slice(NA)
    }
}

cross <- read_cross2(argv$control, quiet = FALSE)
# map <- insert_pseudomarkers(cross$gmap, step = 1)
pr <- calc_genoprob(cross, error_prob = argv$error_prob, cores = 8)
kinship <- calc_kinship(pr, "loco", cores = 8)

tss <- read_tsv(argv$tss, col_types = "cci")
snps <- tibble(id = cross$pmap %>% setNames(NULL) %>% unlist() %>% names()) %>%
    separate(id, c("chr", "pos"), sep = ":", remove = FALSE) %>%
    mutate(chr = str_replace(chr, "chr", ""),
           pos = as.integer(pos))

if (argv$mode == "nom") {

    cat("Nominal mode.\n")
    prog <- 0
    pairs <- tibble(gene = colnames(cross$pheno)) %>%
        group_by(gene) %>%
        summarise({
            prog <- prog + 1
            if (prog %% 100 == 0) {
                cat("Testing gene", prog, "of", length(colnames(cross$pheno)), "\n")
            }
            scan_gene(gene)
        }, .groups = "drop")
    pairs %>%
        mutate(lod = format(lod, nsmall = 6, trim = TRUE)) %>%
        write_tsv(argv$out)

} else if (argv$mode == "perm") {

    cat("Permutation mode.\n")
    prog <- 0
    pairs <- tibble(gene = colnames(cross$pheno)) %>%
        group_by(gene) %>%
        summarise({
            prog <- prog + 1
            if (prog %% 100 == 0) {
                cat("Testing gene", prog, "of", length(colnames(cross$pheno)), "\n")
            }
            scan_gene_perm(gene)
        }, .groups = "drop")
    pairs %>%
        mutate(threshold = format(threshold, nsmall = 6, trim = TRUE)) %>%
        pivot_wider(names_from = gene, values_from = threshold) %>%
        select(-iter) %>%
        write_tsv(argv$out)

} else if (argv$mode == "coef") {

    cat("Coefficient mode.\n")
    apr <- genoprob_to_alleleprob(pr)
    prog <- 0
    pairs <- tibble(gene = colnames(cross$pheno)) %>%
        group_by(gene) %>%
        summarise({
            prog <- prog + 1
            if (prog %% 100 == 0) {
                cat("Testing gene", prog, "of", length(colnames(cross$pheno)), "\n")
            }
            scan_gene_coef(gene)
        }, .groups = "drop")
    pairs %>%
        write_tsv(argv$out)

}
