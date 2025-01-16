localrules:
    exon_table,
    combine_splice_covariates,


rat_to_sample = {}
for tissue in tissues:
    rat_to_sample[tissue] = {}
    for sample in sample_ids[tissue]:
        rat = sample.split('_')[0]
        rat_to_sample[tissue][rat] = sample
rat_ids = {}
for tissue in tissues:
    rat_ids[tissue] = list(rat_to_sample[tissue].keys())


def geno_prefix_splice(wildcards):
    """Genotypes are expression dataset-specific to remove monomorphic SNPs."""
    tissue = dict(IL="IL", LHb="LHB", NAcc="Acbc", OFC="VoLo", PL="PL")[wildcards.tissue]
    return f"data/tensorqtl/genotypes/ensembl-gene_inv-quant_{tissue}"


rule junctions:
    """Much faster than LeafCutter's version, so even they now recommend using regtools.
    Use rat IDs instead of sample IDs so IDs in bed file match with genotype IDs.
    """
    input:
        bam = lambda w: f"data/markdup_out/{w.tissue}/{rat_to_sample[w.tissue][w.rat_id]}.bam",
        bai = lambda w: f"data/markdup_out/{w.tissue}/{rat_to_sample[w.tissue][w.rat_id]}.bam.bai"
    output:
        "data/splice/{tissue}/junc/{rat_id}.junc.gz"
    params:
        interm = "data/splice/{tissue}/junc/{rat_id}.junc"
    shell:
        """
        regtools junctions extract -a 8 -m 50 -M 500000 -s 0 {input.bam} -o {params.interm}
        gzip {params.interm}
        """
        # I was using leafcutter_cluster.py before, but now using leafcutter_cluster_regtools_py3.py I don't need this:
        # cut -f1-6 {params.interm} | gzip -c > {output}
        # rm {params.interm}


rule exon_table:
    input:
        "data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.genes.gtf"
    output:
        "data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.genes.exons.tsv"
    shell:
        "Rscript src/splice/exon_table.R {input} {output}"


rule splicing_bed:
    """Note: fails if intermediate files are already present."""
    input:
        junc = lambda w: expand("data/splice/{tissue}/junc/{rat_id}.junc.gz", tissue = w.tissue, rat_id=rat_ids[w.tissue]),
        exons = "data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.genes.exons.tsv",
        gtf = "data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.genes.gtf",
    output:
        bed = "data/splice/{tissue}/{tissue}.leafcutter.bed.gz",
        bedi = "data/splice/{tissue}/{tissue}.leafcutter.bed.gz.tbi",
        clust = "data/splice/{tissue}/{tissue}.leafcutter.clusters_to_genes.txt",
        pcs = "data/splice/{tissue}/{tissue}.leafcutter.PCs.txt",
        groups = "data/splice/{tissue}/{tissue}.leafcutter.phenotype_groups.txt",
    params:
        # file_list = "data/splice/{tissue}/juncfiles.txt",
        junc_list = lambda w: expand("../junc/{rat_id}.junc.gz", rat_id=rat_ids[w.tissue]),
        outdir = "data/splice/{tissue}/clust",
        path_from_outdir = "../../../../",
    conda:
        "../envs/splice.yaml"
    shell:
        # paths must be relative to output_dir
        # echo {input.junc} | tr ' ' '\n' > {params.file_list}
        # {params.path_from_outdir}{params.file_list} \
        # --output_dir data/splice/{wildcards.tissue}
        """
        mkdir -p {params.outdir}
        echo {params.junc_list} | tr ' ' '\n' > {params.outdir}/juncfiles.txt
        cd {params.outdir}
        python2 {params.path_from_outdir}src/splice/cluster_prepare_fastqtl.py \
        juncfiles.txt \
        {params.path_from_outdir}{input.exons} \
        {params.path_from_outdir}{input.gtf} \
        {wildcards.tissue} \
        --leafcutter_dir {params.path_from_outdir}../tools/leafcutter
        mv -i {wildcards.tissue}.leafcutter* ..
        """
        # mv -i {wildcards.tissue}* data/splice/{wildcards.tissue}/
        # mv -i e}* data/splice/{wildcards.tissue}/


rule combine_splice_covariates:
    input:
        pca = "data/splice/{tissue}/{tissue}.leafcutter.PCs.txt",
        founder_sim = "data/tensorqtl/sim_to_founders.txt",
    output:
        "data/splice/{tissue}/{tissue}.combined_covariates.txt"
    params:
        prefix = "data/splice/{tissue}/{tissue}"
    shell:
        """
        python3 src/combine_covariates.py {input.pca} {params.prefix} \
        --add_covariates {input.founder_sim}
        """

rule tensorqtl_perm_splice:
    input:
        # geno = multiext("data/tensorqtl/genotypes/ensembl-gene_inv-quant_{tissue}", ".bed", ".bim", ".fam"),
        geno = lambda w: multiext(geno_prefix_splice(w), ".bed", ".bim", ".fam"),
        bed = "data/splice/{tissue}/{tissue}.leafcutter.bed.gz",
        bedi = "data/splice/{tissue}/{tissue}.leafcutter.bed.gz.tbi",
        covar = "data/splice/{tissue}/{tissue}.combined_covariates.txt",
        groups = "data/splice/{tissue}/{tissue}.leafcutter.phenotype_groups.txt",
    output:
        "data/splice/{tissue}/{tissue}_splice.cis_qtl.txt.gz"
    params:
        # geno_prefix = "data/tensorqtl/genotypes/ensembl-gene_inv-quant_{tissue}",
        geno_prefix = geno_prefix_splice,
        outdir = "data/splice/{tissue}"
    # conda:
    #     "../envs/tensorqtl.yaml"
    resources:
        walltime = 12,
        # partition = "--partition=gpu",
    shell:
        # module load cuda
        """
        python -m tensorqtl \
            --mode cis \
            {params.geno_prefix} \
            {input.bed} \
            {wildcards.tissue}_splice \
            --covariates {input.covar} \
            --phenotype_groups {input.groups} \
            --output_dir {params.outdir}
        """


rule tensorqtl_independent_splice:
    input:
        # geno = lambda w: multiext("data/tensorqtl/genotypes/ensembl-gene_inv-quant_{tissue}", ".bed", ".bim", ".fam"),
        geno = lambda w: multiext(geno_prefix_splice(w), ".bed", ".bim", ".fam"),
        bed = "data/splice/{tissue}/{tissue}.leafcutter.bed.gz",
        bedi = "data/splice/{tissue}/{tissue}.leafcutter.bed.gz.tbi",
        covar = "data/splice/{tissue}/{tissue}.combined_covariates.txt",
        groups = "data/splice/{tissue}/{tissue}.leafcutter.phenotype_groups.txt",
        cis = "data/splice/{tissue}/{tissue}_splice.cis_qtl.txt.gz",
    output:
        "data/splice/{tissue}/{tissue}_splice.cis_independent_qtl.txt.gz"
    params:
        # geno_prefix = "data/tensorqtl/genotypes/ensembl-gene_inv-quant_{tissue}",
        geno_prefix = geno_prefix_splice,
        outdir = "data/splice/{tissue}"
    # conda:
    #     "../envs/tensorqtl.yaml"
    resources:
        walltime = 12,
        # partition = "--partition=gpu",
    shell:
        # module load cuda
        """
        python -m tensorqtl \
            --mode cis_independent \
            {params.geno_prefix} \
            {input.bed} \
            {wildcards.tissue}_splice \
            --covariates {input.covar} \
            --phenotype_groups {input.groups} \
            --cis_output {input.cis} \
            --output_dir {params.outdir}
        """
