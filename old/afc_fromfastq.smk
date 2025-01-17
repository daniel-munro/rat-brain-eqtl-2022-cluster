def afc_expr_bed(wildcards):
    return {
        "basic": "data/expression/{region}.rsem_expected_count.bed.gz",
        "basic2": "data/expression/{region}.ensemblGeneV2_Norm-counts-log2+1.bed.gz",
        "basic3": "data/expression/{region}.rsem_TPM.bed.gz",
        "basic4": "data/expression/{region}.rsem_expected_count.bed.gz",
        "basic5": "data/expression/{region}.rsem_expected_count_qnorm.bed.gz",
        "main": "data/expression/{region}.rsem_expected_count.bed.gz",
        "main2": "data/expression/{region}.ensemblGeneV2_Norm-counts-log2+1.bed.gz",
        "main3": "data/expression/{region}.rsem_TPM.bed.gz",
        "main4": "data/expression/{region}.rsem_expected_count.bed.gz",
        "main5": "data/expression/{region}.rsem_expected_count_qnorm.bed.gz",
        "qtl2": "data/expression/{region}.rsem_expected_count.bed.gz",
    }[wildcards.method].format(region=wildcards.region)

def afc_covar_file(wildcards):
    """covar_empty files aren't used, but are here as placeholders."""
    return {
        "basic": "data/tensorqtl/{region}/{region}.basic.covar_empty.txt",
        "basic2": "data/tensorqtl/{region}/{region}.basic2.covar_empty.txt",
        "basic3": "data/tensorqtl/{region}/{region}.basic3.covar_empty.txt",
        "basic4": "data/tensorqtl/{region}/{region}.basic4.covar_empty.txt",
        "basic5": "data/tensorqtl/{region}/{region}.basic5.covar_empty.txt",
        "main": "data/tensorqtl/{region}/{region}.main.combined_covariates.txt",
        "main2": "data/tensorqtl/{region}/{region}.main2.combined_covariates.txt",
        "main3": "data/tensorqtl/{region}/{region}.main3.combined_covariates.txt",
        "main4": "data/tensorqtl/{region}/{region}.main4.combined_covariates.txt",
        "main5": "data/tensorqtl/{region}/{region}.main5.combined_covariates.txt",
        "qtl2": "data/tensorqtl/other_covariates.tsv",
    }[wildcards.method].format(region=wildcards.region)

def afc_qtl_file(wildcards):
    if wildcards.method == "qtl2":
        qtl = "data/qtl2/{region}.gene_var_pval.tsv.gz"
    else:
        qtl = "data/tensorqtl/{region}/{region}.{method}.cis_qtl.txt.gz"
    return qtl.format(region=wildcards.region, method=wildcards.method)

def aFC_qtl_arg(wildcards, input):
    if wildcards.method == "qtl2":
        qtl_call = "<(python3 src/prepare_qtl_for_afc.py {qtl} 1 2)"
    else:
        qtl_call = "<(python3 src/prepare_qtl_for_afc.py {qtl})"
    return qtl_call.format(qtl=input.qtl)


def aFC_covar_arg(wildcards, input):
    if wildcards.method[:5] == "basic":
        return ""
    else:
        return "--cov <(sed 's/-N//g' {covar})".format(covar=input.covar)


rule aFC_from_eQTL_model:
    input:
        # unpack(aFC_inputs)
        vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
        vcfi = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz.tbi",
        bed = afc_expr_bed,
        bedi = lambda w: afc_expr_bed(w) + ".tbi",
        qtl = afc_qtl_file,
        covar = afc_covar_file,
    output:
        "data/afc/{region}.{method}.aFC.txt"
    params:
        qtl_arg = aFC_qtl_arg,
        covar_arg = aFC_covar_arg
    conda:
        "envs/afc.yaml"
    resources:
        walltime = 12
    shell:
        """
        python3 tools/aFC/aFC.py \
        --vcf {input.vcf} \
        --pheno {input.bed} \
        --qtl {params.qtl_arg} \
        {params.covar_arg} \
        --log_xform 1 \
        --output {output}
        """

###################
## ASE-based aFC ##
###################

rule phaser_sample_map:
    output:
        "data/afc/{region}.sample_map.txt"
    run:
        with open(output[0], "w") as out:
            out.write("vcf_sample\tbed_sample\n")
            for sample_id in sample_ids[wildcards.region]:
                rat_id = sample_id.split("_")[0]
                bam_base = "{}.Aligned.sortedByCoord.out".format(sample_id)
                out.write("{}\t{}\n".format(rat_id, bam_base))

rule phaser_gene_var_pairs:
    input:
        "data/afc/{region}.{method}.aFC.txt"
    output:
        "data/afc/{region}.{method}.pairs.txt"
    run:
        inlines = open(input[0], "r").read().splitlines()
        gene_var = [x.split("\t")[:2] for x in inlines[1:]]
        with open(output[0], "w") as out:
            out.write("gene_id\tvar_id\tvar_contig\tvar_pos\tvar_ref\tvar_alt\n")
            for x in gene_var:
                chrom, pos = tuple(x[1].split(":"))
                chrom = chrom.replace("chr", "")
                out.write("{}\t{}\t{}\t{}\t\t\n".format(x[0], x[1], chrom, pos))

rule phaser_cis_var:
    input:
        bed = "data/phaser_pop_out/{region}.expr_matrix.gw_phased.bed.gz",
        # bedi = "data/expression/{region}.expr_matrix.gw_phased.bed.gz.tbi",
        vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
        gene_var_pairs = "data/afc/{region}.{method}.pairs.txt",
        sample_map = "data/afc/{region}.sample_map.txt"
    output:
        "data/afc/{region}.{method}.ASE_aFC.txt"
    conda:
        "envs/phaser.yaml"
    resources:
        walltime = 12
    threads: 8
    shell:
        """
        python tools/phaser/phaser_pop/phaser_cis_var.py \
        --bed {input.bed} \
        --vcf {input.vcf} \
        --pair {input.gene_var_pairs} \
        --map {input.sample_map} \
        --t 8 \
        --o {output}
        """
