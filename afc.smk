localrules:
    phaser_sample_map,
    phaser_gene_var_pairs,


def afc_expr_bed(wildcards):
    """Always use count values."""
    tissue = dict(I="IL", L="LHB", N="Acbc", O="VoLo", P="PL")[wildcards.code[0]]
    return f"data/expression/ensembl-gene_log2_{tissue}.bed.gz"


rule aFC_from_eQTL_model:
    input:
        vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
        vcfi = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz.tbi",
        bed = afc_expr_bed,
        bedi = lambda w: afc_expr_bed(w) + ".tbi",
        qtl = "data/tensorqtl/{code}.cis_qtl.txt.gz",
        qtl_indep = "data/tensorqtl/{code}.cis_independent_qtl.txt.gz",
        covar = covar_file,
    output:
        "data/afc/{code}.aFC.txt"
    params:
        covar_arg = lambda w, input: dict(N="", C=f"--cov {input.covar}")[w.code[2]]
    conda:
        "../envs/afc.yaml"
    resources:
        walltime = 12
    shell:
        """
        python3 ~/tools/aFC/aFC.py \
        --vcf {input.vcf} \
        --pheno {input.bed} \
        --qtl <(python3 src/prepare_qtl_for_afc.py {input.qtl} {input.qtl_indep}) \
        {params.covar_arg} \
        --log_xform 1 \
        --output {output}
        """


###################
## ASE-based aFC ##
###################

rule phaser_sample_map:
    output:
        "data/afc/{tissue}.sample_map.txt"
    run:
        with open(output[0], "w") as out:
            out.write("vcf_sample\tbed_sample\n")
            for sample_id in sample_ids[wildcards.tissue]:
                rat_id = sample_id.split("_")[0]
                out.write(f"{rat_id}\t{sample_id}\n")


rule phaser_gene_var_pairs:
    input:
        "data/afc/{code}.aFC.txt"
    output:
        "data/afc/{code}.pairs.txt"
    run:
        inlines = open(input[0], "r").read().splitlines()
        gene_var = [x.split("\t")[:2] for x in inlines[1:]]
        with open(output[0], "w") as out:
            out.write("gene_id\tvar_id\tvar_contig\tvar_pos\tvar_ref\tvar_alt\n")
            for x in gene_var:
                chrom, pos = tuple(x[1].split(":"))
                chrom = chrom.replace("chr", "")
                out.write(f"{x[0]}\t{x[1]}\t{chrom}\t{pos}\t\t\n")


def phaser_bed_file(wildcards):
    tissue = dict(I="IL", L="LHb", N="NAcc", O="OFC", P="PL")[wildcards.code[0]]
    return f"data/phaser_pop_out/{tissue}/{tissue}.expr_matrix.gw_phased.bed.gz"


def phaser_sample_map_file(wildcards):
    tissue = dict(I="IL", L="LHb", N="NAcc", O="OFC", P="PL")[wildcards.code[0]]
    return f"data/afc/{tissue}.sample_map.txt"


rule phaser_cis_var:
    input:
        bed = phaser_bed_file,
        vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
        gene_var_pairs = "data/afc/{code}.pairs.txt",
        sample_map = phaser_sample_map_file,
    output:
        "data/afc/{code}.ASE_aFC.txt"
    conda:
        "../envs/phaser.yaml"
    resources:
        walltime = 12,
        cpus = 16
    # threads: 16
    shell:
        """
        python ~/tools/phaser/phaser_pop/phaser_cis_var.py \
        --bed {input.bed} \
        --vcf {input.vcf} \
        --pair {input.gene_var_pairs} \
        --map {input.sample_map} \
        --t 16 \
        --o {output}
        """

###################################
## aFC for all significant pairs ##
###################################

rule aFC_all_signif:
    """For data portal visualizations"""
    input:
        vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
        vcfi = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz.tbi",
        bed = afc_expr_bed,
        bedi = lambda w: afc_expr_bed(w) + ".tbi",
        qtl = "data/tensorqtl/{code}.cis_qtl_signif.txt.gz",
        covar = covar_file,
    output:
        "data/afc/all_signif/{code}.aFC_all_signif.txt"
    params:
        covar_arg = lambda w, input: dict(N="", C=f"--cov {input.covar}")[w.code[2]]
    conda:
        "../envs/afc.yaml"
    resources:
        walltime = 12
    shell:
        """
        python3 ~/tools/aFC/aFC.py \
        --vcf {input.vcf} \
        --pheno {input.bed} \
        --qtl <(python3 src/prepare_qtl_for_afc.py {input.qtl}) \
        {params.covar_arg} \
        --log_xform 1 \
        --output {output}
        """
