rule qtl2_control_chr:
    output:
        "data/qtl2/chr/{chr}.yaml"
    shell:
        """
        Rscript -e 'qtl2::write_control_file( \
        output_file = "{output}", \
        crosstype = "hs", \
        geno_file = "geno_{wildcards.chr}.csv", \
        founder_geno_file = "founder_geno_{wildcards.chr}.csv", \
        gmap_file = "gmap_{wildcards.chr}.csv", \
        pmap_file = "pmap_{wildcards.chr}.csv", \
        covar_file = "../covar.csv", \
        crossinfo_covar = "generations", \
        geno_codes = c(A = 1L, H = 2L, B = 3L), \
        na.strings = "-", \
        geno_transposed = TRUE, \
        founder_geno_transposed = TRUE)'
        """

# rule qtl2_control_chr_ind:
#     output:
#         "data/qtl2/chr/{chr}/{rat_id}.yaml"
#     shell:
#         """
#         Rscript -e 'qtl2::write_control_file( \
#         output_file = "{output}", \
#         crosstype = "hs", \
#         geno_file = "geno_{wildcards.rat_id}.csv", \
#         founder_geno_file = "../founder_geno_{wildcards.chr}.csv", \
#         gmap_file = "../gmap_{wildcards.chr}.csv", \
#         pmap_file = "../pmap_{wildcards.chr}.csv", \
#         geno_codes = c(A = 1L, H = 2L, B = 3L), \
#         na.strings = "-", \
#         geno_transposed = TRUE, \
#         founder_geno_transposed = TRUE)'
#         """

rule qtl2_subset_chr:
    input:
        "data/qtl2/{filetype}.csv"
    output:
        "data/qtl2/chr/{filetype}_{chr}.csv"
    shell:
        """
        head -1 {input} > {output}
        grep "{wildcards.chr}:" {input} >> {output}
        """

# rule qtl2_geno_chr_ind:
#     input:
#         "data/qtl2/chr/geno_{chr}.csv"
#     output:
#         "data/qtl2/chr/{chr}/geno_{rat_id}.csv"
#     shell:
#         "csvcut -c id,{wildcards.rat_id} {input} > {output}"

rule qtl2_fake_pheno:
    input:
        "data/qtl2/covar.csv"
    output:
        "data/qtl2/chr/pheno.csv"
    shell:
        "csvcut -c id,generations {input} > {output}"

rule qtl2_haplotype_probs_chr:
    input:
        geno = "data/qtl2/chr/geno_{chr}.csv",
        founder_geno = "data/qtl2/chr/founder_geno_{chr}.csv",
        covar = "data/qtl2/covar.csv",
        gmap = "data/qtl2/chr/gmap_{chr}.csv",
        pmap = "data/qtl2/chr/pmap_{chr}.csv",
        control = "data/qtl2/chr/{chr}.yaml"
    output:
        "data/analysis/haplotype_probs_{chr}.rds"
    threads: 8
    shell:
        "Rscript src/run_qtl2_haplotype_probs.R {input.control} {output} 2000"

rule qtl2_haplotype_probs_chr_error_prob:
    input:
        geno = "data/qtl2/chr/geno_{chr}.csv",
        founder_geno = "data/qtl2/chr/founder_geno_{chr}.csv",
        covar = "data/qtl2/covar.csv",
        gmap = "data/qtl2/chr/gmap_{chr}.csv",
        pmap = "data/qtl2/chr/pmap_{chr}.csv",
        control = "data/qtl2/chr/{chr}.yaml"
    output:
        "data/analysis/haplotype_probs_{chr}.error_{error}.rds"
    threads: 8
    shell:
        """
        Rscript src/run_qtl2_haplotype_probs.R {input.control} {output} 2000 \
        --error_prob {wildcards.error}
        """

rule haplotype_sim_chr:
    input:
        geno = "data/qtl2/chr/geno_{chr}.csv",
        founder_geno = "data/qtl2/chr/founder_geno_{chr}.csv",
    output:
        "data/analysis/haplotype_sim_{chr}.npy"
    shell:
        "python3 src/haps_by_identity.py {input.geno} {input.founder_geno} {output} 2000"

rule qtl2_haplotype_probs_baud2014_chr:
    input:
        "data/baud2014/sdir/additive/{chr}/"
    output:
        "data/analysis/baud2014_{chr}.rds"
    shell:
        "Rscript src/HAPPY_probs.R {input} {output}"

# rule qtl2_haplotype_probs_chr_ind:
#     input:
#         geno = "data/qtl2/chr/{chr}/geno_{rat_id}.csv",
#         founder_geno = "data/qtl2/chr/founder_geno_{chr}.csv",
#         covar = "data/qtl2/covar.csv",
#         gmap = "data/qtl2/chr/gmap_{chr}.csv",
#         pmap = "data/qtl2/chr/pmap_{chr}.csv",
#         control = "data/qtl2/chr/{chr}/{rat_id}.yaml"
#     output:
#         "data/qtl2/chr/{chr}/haplotype_probs_{rat_id}.rds"
#     threads: 8
#     shell:
#         "Rscript src/run_qtl2_haplotype_probs.R {input.control} {output} 2000"

# rule qtl2_haplotype_probs_chr:
#     input:
#         expand("data/qtl2/chr/{{chr}}/haplotype_probs_{rat_id}.rds", rat_id=rat_ids)
#     output:
#         "data/analysis/haplotype_probs_{chr}.rds"
#     threads: 8
#     shell:
#         "Rscript src/qtl2_combine_probs.R {input} {output}"
