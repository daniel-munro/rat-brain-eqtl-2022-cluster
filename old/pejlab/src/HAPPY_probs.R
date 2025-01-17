library(argparser)
suppressPackageStartupMessages(library(tidyverse))

p <- arg_parser("Convert HAPPY data archive to 3D array RDS file.")
p <- add_argument(p, "prefix", help="path for directory containing .RData files")
p <- add_argument(p, "out", help="name of output file (*.rds)")
argv <- parse_args(p)

# prefix <- "data/baud2014/sdir/full/chr12/"
# load(str_c(prefix, "bp.RData"))
# load(str_c(prefix, "chromosome.RData"))
# load(str_c(prefix, "haploid.RData"))
# load(str_c(prefix, "map.RData"))
# load(str_c(prefix, "markers.RData"))
# load(str_c(prefix, "markers.safe.RData"))
# load(str_c(prefix, "strains.RData"))
# load(str_c(prefix, "subjects.RData"))

load_happy <- function(prefix, rat_indices) {
    env <- new.env()
    load(str_c(prefix, "markers.RData"), env)
    # tibble(marker = markers) %>%
    #     group_by(marker) %>%
    #     summarise({
    #         env <- new.env()
    #         load(str_c(prefix, "@", marker, ".RData"), env)
    #         mat <- env[[marker]][rat_indices, ]
            
    #     })
    probs <- sapply(env$markers, function(marker) {
        load(str_c(prefix, "subjects.RData"), env)
        load(str_c(prefix, "strains.RData"), env)
        load(str_c(prefix, "@", marker, ".RData"), env)
        mat <- env[[marker]]
        mat <- mat[, c(1:6, 8, 7)]  # Swap WN and WKY to match our VCF/qtl2 order.
        # dimnames(mat) <- list(individual = env$subjects, strain = env$strains)
        rownames(mat) <- env$subjects
        # colnames(mat) <- env$strains
        colnames(mat) <- 1:8
        mat[rat_indices, ]
    }, simplify = "array")
    rm(env)
    probs
}

rat_indices <- seq(from = 30, to = 1500, by = 30)  # Choose 50 to vizualize
x <- load_happy(argv$prefix, rat_indices)
saveRDS(x, argv$out)
