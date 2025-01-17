localrules:
    phaser_sample_map,
    phaser_gene_var_pairs,

def afc_expr_bed(wildcards):
    """Always use count-like values."""
    # expr = dict(L="log2", Q="log2", R="rlog")[wildcards.code[1]]
    expr = dict(
        L="log2", Q="log2", R="rlog",
        C="log2_ComBat", D="log2_ComBat", E="rlog_ComBat",
    )[wildcards.code[1]]
    region = dict(A="Acbc", I="IL", L="LHB", P="PL", V="VoLo")[wildcards.code[0]]
    return "data/expression/ensembl-gene_{}_{}.bed.gz".format(expr, region)

# def afc_covar_file(wildcards):
#     """covar_empty files aren't used, but are here as placeholders."""
#     return {
#         "basic": "data/tensorqtl/{region}/{region}.basic.covar_empty.txt",
#         "basic2": "data/tensorqtl/{region}/{region}.basic2.covar_empty.txt",
#         "basic3": "data/tensorqtl/{region}/{region}.basic3.covar_empty.txt",
#         "basic4": "data/tensorqtl/{region}/{region}.basic4.covar_empty.txt",
#         "basic5": "data/tensorqtl/{region}/{region}.basic5.covar_empty.txt",
#         "main": "data/tensorqtl/{region}/{region}.main.combined_covariates.txt",
#         "main2": "data/tensorqtl/{region}/{region}.main2.combined_covariates.txt",
#         "main3": "data/tensorqtl/{region}/{region}.main3.combined_covariates.txt",
#         "main4": "data/tensorqtl/{region}/{region}.main4.combined_covariates.txt",
#         "main5": "data/tensorqtl/{region}/{region}.main5.combined_covariates.txt",
#         "qtl2": "data/tensorqtl/other_covariates.tsv",
#     }[wildcards.method].format(region=wildcards.region)

def afc_qtl_file(wildcards):
    if wildcards.code[3] == "Q":
        qtl = "data/qtl2/{code}.gene_var_pval.tsv.gz"
    else:
        qtl = "data/tensorqtl/{code}.cis_qtl.txt.gz"
    return qtl.format(code=wildcards.code)

def aFC_qtl_arg(wildcards, input):
    if wildcards.code[3] == "Q":
        qtl_call = "<(python3 src/prepare_qtl_for_afc.py {qtl} 1 2)"
    else:
        qtl_call = "<(python3 src/prepare_qtl_for_afc.py {qtl})"
    return qtl_call.format(qtl=input.qtl)

def aFC_covar_arg(wildcards, input):
    if wildcards.code[2] == "N":
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
        covar = covar_file,
    output:
        "data/afc/{code}.aFC.txt"
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
                out.write("{}\t{}\n".format(rat_id, sample_id))

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
                out.write("{}\t{}\t{}\t{}\t\t\n".format(x[0], x[1], chrom, pos))

def phaser_bed_file(wildcards):
    region = dict(A="Acbc", I="IL", L="LHB", P="PL", V="VoLo")[wildcards.code[0]]
    return "data/phaser_pop_out/{}.expr_matrix.gw_phased.bed.gz".format(region)

def phaser_sample_map_file(wildcards):
    region = dict(A="Acbc", I="IL", L="LHB", P="PL", V="VoLo")[wildcards.code[0]]
    return "data/afc/{}.sample_map.txt".format(region)

rule phaser_cis_var:
    input:
        bed = phaser_bed_file,
        vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
        gene_var_pairs = "data/afc/{code}.pairs.txt",
        sample_map = phaser_sample_map_file,
    output:
        "data/afc/{code}.ASE_aFC.txt"
    conda:
        "envs/phaser.yaml"
    resources:
        walltime = 12,
        cpus = 16
    # threads: 16
    shell:
        """
        python tools/phaser/phaser_pop/phaser_cis_var.py \
        --bed {input.bed} \
        --vcf {input.vcf} \
        --pair {input.gene_var_pairs} \
        --map {input.sample_map} \
        --t 16 \
        --o {output}
        """
