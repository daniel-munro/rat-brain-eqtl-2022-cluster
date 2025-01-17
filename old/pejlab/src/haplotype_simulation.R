library(tidyverse)

next_generation <- function(genos) {
    genos %>%
        mutate(father = c(male[2:length(male)], male[1])) %>%
        select(cage, mother = female, father) %>%
        rowwise() %>%
        mutate(female = list(c(sample(mother, 1), sample(father, 1))),
               male = list(c(sample(mother, 1), sample(father, 1)))) %>%
        ungroup() %>%
        select(cage, female, male)
}

maps <- tibble(chrom = 1:20) %>%
    group_by(chrom) %>%
    summarise(read_delim(str_c("data/qtl2/genetic_map/MAP4chr", chrom, ".txt"),
                         " ", col_names = c("pos", "ratio", "cm"), col_types = "idd"),
              .groups = "drop")

maps %>%
    group_by(chrom) %>%
    summarise(gen_length = max(cm) / 100, .groups = "drop") %>%
    mutate(recombs = gen_length * 90)

# genos <- lapply(1:64, function(pair) list(female = sample(LETTERS[1:8], 2, replace = TRUE),
#                                           male = sample(LETTERS[1:8], 2, replace = TRUE)))
genos1 <- tibble(cage = 1:64) %>%
    group_by(cage) %>%
    mutate(female = list(sample(LETTERS[1:8], 2, replace = TRUE)),
           male = list(sample(LETTERS[1:8], 2, replace = TRUE))) %>%
    ungroup()
genos <- list(genos1)
for (gen in 2:90) {
    genos <- c(genos, list(next_generation(genos[[length(genos)]])))
}
genos <- bind_rows(genos, .id = "generation") %>%
    mutate(generation = as.integer(generation))

genos %>%
    group_by(generation) %>%
    