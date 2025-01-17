rule rsem_index:
    input:
        fasta = "data/Rnor_6.0/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa",
        gtf = "data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.gtf"
    output:
        "data/rsem_reference/rsem_reference.transcripts.fa"
    conda:
        "envs/rsem.yaml"
    threads: 16
    shell:
        """
        rsem-prepare-reference \
        {input.fasta} \
        data/rsem_reference/rsem_reference \
        --gtf {input.gtf} \
        --num-threads 16
        """

rule rsem:
    input:
        ref = "data/rsem_reference/rsem_reference.transcripts.fa",
        bam = "data/star_out/{region}/{sample_id}.Aligned.toTranscriptome.out.bam"
    output:
        "data/rsem_out/{region}/{sample_id}.genes.results.gz"
    params:
        prefix = "data/rsem_out/{region}/{sample_id}"
    conda:
        "envs/rsem.yaml"
    threads: 16
    shell:
        """
        rsem-calculate-expression \
        --num-threads 16 \
        --quiet \
        --estimate-rspd \
        --no-bam-output \
        --alignments {input.bam} \
        data/rsem_reference/rsem_reference \
        {params.prefix}
        gzip {params.prefix}.genes.results
        rm {params.prefix}.isoforms.results
        rm -r {params.prefix}.stat
        """

def rsem_out_files(wildcards):
    return expand("data/rsem_out/{{region}}/{sample_id}.genes.results.gz",
                  sample_id=sample_ids[wildcards.region])


rule combine_RSEM:
    input:
        rsem_out_files
    output:
        "data/expression/{region}.rsem_{field}.gct.gz"
    shell:
        "python3 src/combine_RSEM.py {input} {wildcards.field} {output}"
