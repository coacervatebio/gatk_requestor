ref = "resources/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"
samples = glob_wildcards("alignments/minified/2_sample/{samp}_minified.cram").samp
container: "docker://broadinstitute/gatk"

rule all:
    input:
        # expand("alignments/minified/2_sample/{sample}_minified.cram.crai", sample=samples)
        # expand("alignments/chr21/{sample}_minified_chr21.cram", sample=samples)
        # expand("alignments/chr21/{sample}_minified_chr21.cram.crai", sample=samples)
        # expand("hc_out/chr21/{sample}_minified_chr21.g.vcf.gz", sample=samples)
        # "combi_out/my_database"
        "geno_out/2_sample_minified_chr21.vcf.gz"

rule index_cram:
    input:
        "alignments/minified/2_sample/{sample}_minified.cram"
    output:
        "alignments/minified/2_sample/{sample}_minified.cram.crai"
    shell:
        "samtools index {input}"

rule split_cram:
    input:
        "alignments/minified/2_sample/{sample}.cram",
    output:
        "alignments/chr21/{sample}_chr21.cram"
    shell:
        "samtools view {input} chr21 -T {ref} "
        "-O cram -o {output}"

rule index_split_cram:
    input:
        "alignments/chr21/{sample}.cram"
    output:
        "alignments/chr21/{sample}.cram.crai"
    shell:
        "samtools index {input}"

rule call_variants:
    input:
        "alignments/chr21/{sample}.cram"
    output:
        "hc_out/chr21/{sample}.g.vcf.gz"
    shell:
        "gatk HaplotypeCaller -I {input} "
        "-O {output} -R {ref} -L chr21 -ERC GVCF"

rule db_combine:
    input:
        expand("hc_out/chr21/{sample}_minified_chr21.g.vcf.gz", sample=samples)
    output:
        directory("combi_out/my_database/")
    params:
        lambda wildcards, input: ' '.join([f'-V {file}' for file in input])
    shell:
        "gatk GenomicsDBImport {params} --genomicsdb-workspace-path {output} -L chr21"

rule db_genotype:
    input:
        "combi_out/my_database/"
    output:
        "geno_out/2_sample_minified_chr21.vcf.gz"
    shell:
        "gatk GenotypeGVCFs -R {ref} -V gendb://{input} -O {output}"


