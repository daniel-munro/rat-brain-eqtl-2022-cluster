localrules:
    plink_sample_list,
    vcf_to_plink,
    # sample_participant_lookup,
    # prepare_expression,
    # expression_gct_to_bed,
    # expression_gct_to_bed_qnorm,
    expression_tsv_to_bed,
    # expression_gct_to_bed,
    expression_pcs,
    other_covariates,
    combat_covariates,
    combine_covariates,
    empty_covariates,
    tensorqtl_all_signif,

rule plink_sample_list:
    """Assumes any expression file has same individuals as this one for that region."""
    input:
        "data/expression/{base}.bed.gz"
    output:
        "data/tensorqtl/genotypes/{base}.samples.txt"
    shell:
        "zcat {input} | head -1 | cut -f5- | sed 's/\\t/\\n/g' > {output} || true"

def geno_prefix(wildcards):
    # expr = dict(L="log2", Q="inv-quant", R="rlog")[wildcards.code[1]]
    expr = dict(
        L="log2", Q="inv-quant", R="rlog",
        C="log2_ComBat", D="inv-quant_ComBat", E="rlog_ComBat",
    )[wildcards.code[1]]
    region = dict(A="Acbc", I="IL", L="LHB", P="PL", V="VoLo")[wildcards.code[0]]
    return "data/tensorqtl/genotypes/ensembl-gene_{}_{}".format(expr, region)

rule vcf_to_plink:
    input:
        vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
        samples = "data/tensorqtl/genotypes/{base}.samples.txt"
    output:
        # lambda w: multiext(geno_prefix(w), ".bed", ".bim", ".fam")
        multiext("data/tensorqtl/genotypes/{base}", ".bed", ".bim", ".fam")
    params:
        prefix = "data/tensorqtl/genotypes/{base}"
    shell:
        """
        plink2 --make-bed \
        --vcf {input.vcf} \
        --keep {input.samples} \
        --maf 0.0001 \
        --max-maf 0.9999 \
        --out {params.prefix}
        """

# rule vcf_to_plink:
#     input:
#         "data/genotype/P50.rnaseq.88.unpruned.vcf.gz"
#     output:
#         multiext("data/genotype/P50.rnaseq.88.unpruned", ".bed", ".bim", ".fam")
#     params:
#         prefix = "data/genotype/P50.rnaseq.88.unpruned"
#     shell:
#         """
#         plink2 --make-bed \
#         --vcf {input} \
#         --maf 0.0001 \
#         --max-maf 0.9999 \
#         --out {params.prefix}
#         """

rule expression_tsv_to_bed:
    input:
        tsv = "data/expression/ensembl-gene_{expr}_{region}.txt",
        gtf = "data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.genes.gtf"
    output:
        bed = "data/expression/ensembl-gene_{expr}_{region}.bed.gz",
        bedi = "data/expression/ensembl-gene_{expr}_{region}.bed.gz.tbi"
    params:
        interm = "data/expression/ensembl-gene_{expr}_{region}.bed"
    conda:
        "envs/bioinfo.yaml"
    shell:
        """
        python3 src/tsv_to_bed.py {input.tsv} {input.gtf} {params.interm}
        bgzip {params.interm}
        tabix {output.bed}
        """

def expr_bed(wildcards):
    # expr = dict(L="log2", Q="inv-quant", R="rlog")[wildcards.code[1]]
    expr = dict(
        L="log2", Q="inv-quant", R="rlog",
        C="log2_ComBat", D="inv-quant_ComBat", E="rlog_ComBat",
    )[wildcards.code[1]]
    region = dict(A="Acbc", I="IL", L="LHB", P="PL", V="VoLo")[wildcards.code[0]]
    return "data/expression/ensembl-gene_{}_{}.bed.gz".format(expr, region)

def peer_file(wildcards):
    # expr = dict(L="log2", Q="inv-quant", R="rlog")[wildcards.code[1]]
    expr = dict(
        L="log2", Q="inv-quant", R="rlog",
        C="log2_ComBat", D="inv-quant_ComBat", E="rlog_ComBat",
    )[wildcards.code[1]]
    region = dict(A="Acbc", I="IL", L="LHB", P="PL", V="VoLo")[wildcards.code[0]]
    return f"data/expression/peer/ensembl-gene_{expr}_{region}.PEER_covariates.txt"

def expr_pca_file(wildcards):
    # expr = dict(L="log2", Q="inv-quant", R="rlog")[wildcards.code[1]]
    expr = dict(
        L="log2", Q="inv-quant", R="rlog",
        C="log2_ComBat", D="inv-quant_ComBat", E="rlog_ComBat",
    )[wildcards.code[1]]
    region = dict(A="Acbc", I="IL", L="LHB", P="PL", V="VoLo")[wildcards.code[0]]
    return f"data/expression/pca/ensembl-gene_{expr}_{region}.pca.txt"

def covar_file(wildcards):
    # region = dict(A="Acbc", I="IL", L="LHB", P="PL", V="VoLo"}[wildcards.code[0]]
    if wildcards.code[1] in ["C", "D", "E"]:
        covar = dict(N="covar_empty", C="combat_covariates")[wildcards.code[2]]
    else:
        covar = dict(N="covar_empty", C="combined_covariates")[wildcards.code[2]]
    # return "data/tensorqtl/{}/{}.{}.txt".format(region, wildcards.code, covar)
    return "data/tensorqtl/{}.{}.txt".format(wildcards.code, covar)

rule run_peer:
    input:
        # expr_bed
        "data/expression/{base}.bed.gz"
    output:
        # "data/tensorqtl/{code}.PEER_residuals.txt",
        # "data/tensorqtl/{code}.PEER_alpha.txt",
        # "data/tensorqtl/{code}.PEER_covariates.txt"
        "data/expression/peer/{base}.PEER_residuals.txt",
        "data/expression/peer/{base}.PEER_alpha.txt",
        "data/expression/peer/{base}.PEER_covariates.txt"
    params:
        # prefix = "data/tensorqtl/{code}",
        prefix = "data/expression/peer/{base}",
        num_peer = 20
    conda:
        "envs/peer.yaml"
    resources:
        # walltime = "6:00:00"
        walltime = 6
    shell:
        "Rscript src/run_PEER.R {input} {params.prefix} {params.num_peer}"

rule expression_pcs:
    input:
        "data/expression/{base}.bed.gz"
    output:
        "data/expression/pca/{base}.pca.txt"
    params:
        n_pcs = 10
    conda:
        "envs/bioinfo.yaml"
    shell:
        "Rscript src/expression_pcs.R {input} {output} {params.n_pcs}"

rule other_covariates:
    input:
        # "data/DemographicInfo_TN_RNASeq_88_LauraSaba.csv"
        "data/rat_info.txt"
    output:
        "data/tensorqtl/other_covariates.tsv"
    shell:
        "python3 src/tensorqtl_covar.py {input} {output}"

rule similarity_to_founders:
    input:
        pop_vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
        pop_vcfi = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz.tbi",
        founder_vcf = "data/genotype/founders.vcf.gz",
        founder_vcfi = "data/genotype/founders.vcf.gz.tbi"
    output:
        "data/tensorqtl/sim_to_founders.txt"
    conda:
        "envs/biopython.yaml"
    shell:
        "python3 src/sim_to_each_founder.py {input.pop_vcf} {input.founder_vcf} {output}"

rule combine_covariates:
    input:
        # peer = "data/tensorqtl/{code}.PEER_covariates.txt",
        # peer = peer_file,
        pca = expr_pca_file,
        # header = "data/tensorqtl/{code}.covar_empty.txt",
        founder_sim = "data/tensorqtl/sim_to_founders.txt",
        other = "data/tensorqtl/other_covariates.tsv"
    output:
        "data/tensorqtl/{code}.combined_covariates.txt"
    params:
        prefix = "data/tensorqtl/{code}"
    shell:
        """
        python3 src/combine_covariates.py {input.pca} {params.prefix} \
        --add_covariates {input.founder_sim} {input.other}
        """

rule combat_covariates:
    input:
        # header = "data/tensorqtl/{code}.covar_empty.txt",
        pca = expr_pca_file,
        founder_sim = "data/tensorqtl/sim_to_founders.txt",
        other = "data/tensorqtl/other_covariates.tsv"
    output:
        "data/tensorqtl/{code}.combat_covariates.txt"
    params:
        prefix = "data/tensorqtl/{code}",
        tmp = "data/tensorqtl/{code}.combined_covariates.txt"
    shell:
        """
        python3 src/combine_covariates.py {input.pca} {params.prefix} \
        --add_covariates {input.founder_sim} \
        <(grep -vE "^sequencing_batch" {input.other})
        mv {params.tmp} {output}
        """

rule empty_covariates:
    """tensorQTL currently requires a covariates file."""
    input:
        expr_bed
    output:
        "data/tensorqtl/{code}.covar_empty.txt"
    shell:
        "zcat {input} | head -1 | sed 's/gene_id/ID/' | cut -f4- > {output} || true"

rule tensorqtl_perm:
    input:
        # geno = multiext("data/genotype/P50.rnaseq.88.unpruned", ".bed", ".bim", ".fam"),
        geno = lambda w: multiext(geno_prefix(w), ".bed", ".bim", ".fam"),
        bed = expr_bed,
        bedi = lambda w: expr_bed(w) + ".tbi",
        covar = covar_file
    output:
        "data/tensorqtl/{code}.cis_qtl.txt.gz"
    params:
        geno_prefix = "data/genotype/P50.rnaseq.88.unpruned",
        outdir = "data/tensorqtl"
    conda:
        "envs/tensorqtl.yaml"
    resources:
        walltime = 12,
    #     # gpu = ",nodes=1:gpus=1"
        # gpu = ",nodes=1:ppn=8:gtx1080ti -q gpu"
    #     # gpu = ",nodes=1:ppn=20:gtx1080 -q gpu"
    shell:
        """
        # module load cuda
        # pip install -e tools/tensorqtl
        python -m tensorqtl --mode cis \
            {params.geno_prefix} \
            {input.bed} \
            {wildcards.code} \
            --covariates {input.covar} \
            --output_dir {params.outdir}
        """

rule tensorqtl_nominal:
    input:
        # geno = multiext("data/genotype/P50.rnaseq.88.unpruned", ".bed", ".bim", ".fam"),
        geno = lambda w: multiext(geno_prefix(w), ".bed", ".bim", ".fam"),
        bed = expr_bed,
        bedi = lambda w: expr_bed(w) + ".tbi",
        covar = covar_file
    output:
        # "data/tensorqtl/{code}.cis_qtl.txt.gz"
        expand("data/tensorqtl/{{code}}/{{code}}.cis_qtl_pairs.{chrn}.parquet",
               chrn=range(1, 21))
    params:
        geno_prefix = "data/genotype/P50.rnaseq.88.unpruned",
        outdir = "data/tensorqtl/{code}"
    conda:
        "envs/tensorqtl.yaml"
    resources:
        walltime = 12,
    #     # gpu = ",nodes=1:gpus=1"
        # gpu = ",nodes=1:ppn=8:gtx1080ti -q gpu"
    #     # gpu = ",nodes=1:ppn=20:gtx1080 -q gpu"
    shell:
        """
        # module load cuda
        # pip install -e tools/tensorqtl
        mkdir -p {params.outdir}
        python -m tensorqtl --mode cis_nominal \
            {params.geno_prefix} \
            {input.bed} \
            {wildcards.code} \
            --covariates {input.covar} \
            --output_dir {params.outdir}
        """

rule tensorqtl_all_signif:
    input:
        perm = "data/tensorqtl/{code}.cis_qtl.txt.gz",
        nom = expand("data/tensorqtl/{{code}}/{{code}}.cis_qtl_pairs.{chrn}.parquet", chrn=range(1, 21))
    output:
        "data/tensorqtl/{code}.cis_qtl_signif.txt.gz"
    params:
        nom_prefix = "data/tensorqtl/{code}/{code}"
    conda:
        "envs/tensorqtl.yaml"
    shell:
        "python3 src/tensorqtl_all_sig_eqtls.py {input.perm} {params.nom_prefix} {output}"

rule tensorqtl_trans:
    input:
        # geno = multiext("data/genotype/P50.rnaseq.88.unpruned", ".bed", ".bim", ".fam"),
        geno = lambda w: multiext(geno_prefix(w), ".bed", ".bim", ".fam"),
        bed = expr_bed,
        bedi = lambda w: expr_bed(w) + ".tbi",
        covar = covar_file
    output:
        # "data/tensorqtl/{code}.trans_qtl_pairs.parquet"
        "data/tensorqtl/{code}.trans_qtl_pairs.txt.gz"
    params:
        geno_prefix = "data/genotype/P50.rnaseq.88.unpruned",
        # outdir = "data/tensorqtl"
    # conda:
    #     "envs/tensorqtl.yaml"
    # shell:
    #     """
    #     # module load cuda
    #     # pip install -e tools/tensorqtl
    #     python -m tensorqtl --mode trans \
    #         {params.geno_prefix} \
    #         {input.bed} \
    #         {wildcards.code} \
    #         --covariates {input.covar} \
    #         --output_dir {params.outdir}
    #     """
    run:
        # Run in python script so cis variants aren't filtered out.
        # Based on examples at https://github.com/broadinstitute/tensorqtl
        import tensorqtl
        from tensorqtl import genotypeio, trans
        pheno, pheno_pos = tensorqtl.read_phenotype_bed(input.bed)
        covar = pd.read_csv(input.covar, sep="\t", index_col=0).T
        geno = genotypeio.PlinkReader(params.geno_prefix).load_genotypes()
        # Doesn't work due to bug (get_all_genotypes not defined):
        # geno = tensorqtl.read_genotypes(params.geno_prefix)
        d = trans.map_trans(geno, pheno, covar)
        d.to_csv(output[0], sep="\t", index=False, float_format="%.6g")

# rule tensorqtl_perm_chr:
#     """For debugging."""
#     input:
#         # geno = multiext("data/genotype/P50.rnaseq.88.unpruned", ".bed", ".bim", ".fam"),
#         geno = lambda w: multiext(geno_prefix(w), ".bed", ".bim", ".fam"),
#         bed = expr_bed,
#         bedi = lambda w: expr_bed(w) + ".tbi",
#         covar = covar_file
#     output:
#         "data/tensorqtl/{code}.{chrn}.cis_qtl.txt.gz"
#     params:
#         geno_prefix = "data/genotype/P50.rnaseq.88.unpruned",
#         outdir = "data/tensorqtl"
#     conda:
#         "envs/tensorqtl.yaml"
#     resources:
#         walltime = 12,
#     #     # gpu = ",nodes=1:gpus=1"
#         # gpu = ",nodes=1:ppn=8:gtx1080ti -q gpu"
#     #     # gpu = ",nodes=1:ppn=20:gtx1080 -q gpu"
#     shell:
#         """
#         zcat {input.bed} | head -1 > {params.outdir}/{code}.{chrn}.bed
#         zcat {input.bed} | ......
#         python -m tensorqtl --mode cis \
#             {params.geno_prefix} \
#             {input.bed} \
#             {wildcards.code}.{wildcards.chrn} \
#             --covariates {input.covar} \
#             --output_dir {params.outdir}
#         """
