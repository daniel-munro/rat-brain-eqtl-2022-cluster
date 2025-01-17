# Based on https://kbroman.org/qtl2/assets/vignettes/user_guide.html
library(argparser)
library(qtl2)

# args <- commandArgs(trailingOnly = TRUE)
# infile <- args[1]
# outfile <- args[2]
# n_SNPs <- as.integer(args[3])
p <- arg_parser("Calculate and save 3D array of R/qtl2 genotype probabilities")
p <- add_argument(p, "control", help="name of qtl2 control file")
p <- add_argument(p, "out", help="name of output file (*.rds)")
p <- add_argument(p, "n_SNPs", help="number of representative SNPs to save", default=2000L)
p <- add_argument(p, "--error_prob", help="error probability for HMM", default=0.01)
argv <- parse_args(p)
# infile <- "data/qtl2/chr/chr12/00077E67B5.yaml"
# outfile <- "analysis/haplotype_probs_chr12.rds"
# n_SNPs <- 2000

N_CORES = 8

cross <- read_cross2(argv$control)
pr <- calc_genoprob(cross, error_prob = argv$error_prob, cores = N_CORES)
apr <- genoprob_to_alleleprob(pr)[[1]]

# Save probs for only every Nth (e.g. to get 2000 loci per chromosome) for plotting, etc.
SNP_represent <- round(seq(from = 1, to = dim(apr)[3], length.out = argv$n_SNPs))
apr <- apr[, , SNP_represent, drop = FALSE]

saveRDS(apr, argv$out)
