# Use qtl2 to get haplotype probabilities separately for both chromosome copies. Then combine
# these with ASE to get ASE per gene in relation to founder haplotypes rather than individual SNPs.

rule qtl2_control_file_phase:
    output:
        "data/qtl2/phase/phase{phase}.yaml"
    conda:
        "envs/qtl2.yaml"
    shell:
        """
        Rscript -e 'qtl2::write_control_file( \
        output_file = "{output}", \
        crosstype = "hs", \
        geno_file = "geno_phase{wildcards.phase}.csv", \
        geno_transposed = TRUE, \
        founder_geno_file = "founder_geno.csv", \
        founder_geno_transposed = TRUE, \
        gmap_file = "gmap.csv", \
        pmap_file = "pmap.csv", \
        covar_file = "covar.csv", \
        crossinfo_covar = "generations", \
        geno_codes = c(A = 1L, H = 2L, B = 3L), \
        na.strings = "-")'
        """

rule qtl2_geno_files_phase:
    input:
        pop_vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
        pop_vcfi = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz.tbi",
        founder_vcf = "data/genotype/founders.vcf.gz",
        founder_vcfi = "data/genotype/founders.vcf.gz.tbi",
        snps = "data/genotype/imputing/observed.snplist.txt",
        src = "src/qtl2_geno_phase.py",
    output:
        geno = expand("data/qtl2/phase/geno_phase{phase}.csv", phase=[1, 2]),
        founder_geno = "data/qtl2/phase/founder_geno.csv",
        pmap = "data/qtl2/phase/pmap.csv",
        gmap = "data/qtl2/phase/gmap.csv"
    params:
        gen_maps_dir = "data/genotype/genetic_map"
    conda:
        "envs/biopython.yaml"
    shell:
        """
        python3 src/qtl2_geno_phase.py \
        {input.pop_vcf} {input.founder_vcf} {input.snps} {params.gen_maps_dir} \
        {output.geno} {output.founder_geno} {output.pmap} {output.gmap}
        """

rule tss_for_qtl2_phase:
    input:
        "data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.genes.gtf"
    output:
        "data/qtl2/phase/tss.tsv"
    conda:
        "envs/bioinfo.yaml"
    shell:
        "python3 src/tss_from_gtf.py {input} {output}"

rule qtl2_phased_haplotype_probs:
    input:
        geno = "data/qtl2/phase/geno_phase{phase}.csv",
        founder_geno = "data/qtl2/phase/founder_geno.csv",
        covar = "data/qtl2/phase/covar.csv",
        gmap = "data/qtl2/phase/gmap.csv",
        pmap = "data/qtl2/phase/pmap.csv",
        control = "data/qtl2/phase/phase{phase}.yaml",
        tss = "data/qtl2/phase/tss.tsv"
    output:
        "data/qtl2/phase/hap_probs_phase{phase}.rds"
    conda:
        "envs/qtl2.yaml"
    threads: 8
    shell:
        "Rscript src/run_qtl2_hap_probs_at_loci.R {input.control} {input.tss} {output}"
