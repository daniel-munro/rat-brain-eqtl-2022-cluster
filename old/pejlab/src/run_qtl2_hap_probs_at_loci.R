# Based on https://kbroman.org/qtl2/assets/vignettes/user_guide.html
library(argparser)
library(qtl2)
library(tidyverse)

p <- arg_parser("Calculate and save 3D array of R/qtl2 genotype probabilities")
p <- add_argument(p, "control", help="name of qtl2 control file")
# p <- add_argument(p, "loci", help="name of tab-delimited file giving 'chr' and 'pos' for loci to measure at.")
# p <- add_argument(p, "columns", help="column names in loci file corresponding to chromosome and position, separated by comma.")
p <- add_argument(p, "loci", help="name of tab-delimited file giving 'gene_id', 'chr', and 'pos' for TSS to measure at.")
p <- add_argument(p, "out", help="name of output file (*.rds)")
argv <- parse_args(p)

ERROR_PROB <- 0.01
N_CORES <- 8

cross <- read_cross2(argv$control)

# colnames <- str_split(argv$columns, ",")[[1]]
loci <- read_tsv(argv$loci, col_types = cols(gene_id = "c", chr = "i", pos = "i", .default = "-")) %>%
    mutate(pos = pos / 1e6)
loc <- map(unique(loci$chr), ~ loci %>% filter(chr == .x) %>% select(gene_id, pos) %>% deframe())
names(loc) <- unique(loci$chr)
loc <- interp_map(loc, cross$pmap, cross$gmap)

# tol = 0 is needed so that genes close to existing markers are still included:
map <- insert_pseudomarkers(cross$gmap, pseudomarker_map = loc, tol = 0)
pr <- calc_genoprob(cross, map, error_prob = ERROR_PROB, cores = N_CORES)
apr <- genoprob_to_alleleprob(pr)
apr <- clean_genoprob(apr, value_threshold = 1e-3)
print(dim(apr[[1]]))
for (chrom in names(apr)) {
    x <- apr[[chrom]]
    apr[[chrom]] <- x[, , names(loc[[chrom]]), drop = FALSE]
}

saveRDS(apr, argv$out)
