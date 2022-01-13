import os

container: "docker://broadinstitute/gatk"
SAMP_COUNT = 1
REF = "resources/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"
SAMPLES = [fn.strip("_minified.cram") for fn in os.listdir(f"alignments/minified/{SAMP_COUNT}_sample/")]
REGIONS = [
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    # "chr5",
    # "chr6",
    # "chr7",
    # "chr8",
    # "chr9",
    # "chr10",
    # "chr11",
    # "chr12",
    # "chr13",
    # "chr14",
    # "chr15",
    # "chr16",
    # "chr17",
    # "chr18",
    # "chr19",
    # "chr20",
    # "chr21",
    # "chr22",
    # "chr23"
]


rule all:
    input:
        # expand("alignments/minified/{samp_c}_sample/{sample}_minified.cram.crai", sample=SAMPLES, samp_c=SAMP_COUNT)
        # expand("alignments/{reg}/{sample}_minified_{reg}.cram", sample=SAMPLES, reg=REGIONS)
        # expand("alignments/{reg}/{sample}_minified_{reg}.cram.crai", sample=SAMPLES, reg=REGIONS),
        expand("hc_out/{reg}/{sample}_minified_{reg}.g.vcf.gz", sample=SAMPLES, reg=REGIONS)
        # expand("combi_out/{reg}_database/", reg=REGIONS)
        # expand("geno_out/{samps}_sample_minified_{reg}.vcf.gz", reg=REGIONS, samps=SAMP_COUNT)

rule index_cram:
    input:
        expand("alignments/minified/{samps}_sample/{{sample}}_minified.cram", samps=SAMP_COUNT)
    output:
        expand("alignments/minified/{samps}_sample/{{sample}}_minified.cram.crai", samps=SAMP_COUNT)
    benchmark:
        "benchmarks/index_original/{sample}.tsv"
    shell:
        "samtools index {input}"

rule split_cram:
    input:
        alignments=expand("alignments/minified/{samps}_sample/{{sample}}.cram", samps=SAMP_COUNT),
        indexes=expand("alignments/minified/{samps}_sample/{{sample}}.cram.crai", samps=SAMP_COUNT)
    output:
        "alignments/{reg}/{sample}_{reg}.cram"
    benchmark:
        "benchmarks/split_cram/{sample}_{reg}.tsv"
    shell:
        "samtools view {input} {wildcards.reg} -T {REF} -O cram -o {output}"

rule index_split_cram:
    input:
        "alignments/{reg}/{sample}_minified_{reg}.cram"
    output:
        "alignments/{reg}/{sample}_minified_{reg}.cram.crai"
    benchmark:
        "benchmarks/index_split_cram/{sample}_{reg}.tsv"
    shell:
        "samtools index {input}"

rule call_variants:
    input:
        alignments="alignments/{reg}/{sample}.cram",
        indexes="alignments/{reg}/{sample}.cram.crai"
    output:
        "hc_out/{reg}/{sample}.g.vcf.gz"
    benchmark:
        "benchmarks/call_variants/{sample}_{reg}.tsv"
    shell:
        "gatk --java-options '-Xmx4g' HaplotypeCaller -I {input.alignments} -O {output} -R {REF} -L {wildcards.reg} -ERC GVCF"

# rule db_combine:
#     input:
#         set(expand("hc_out/{reg}/{sample}_minified_{reg}.g.vcf.gz", sample=SAMPLES, reg=REGIONS))
#     output:
#         directory("combi_out/{reg}_database/")
#     params:
#         lambda wildcards, input: ' '.join([f'-V {f_outer}' for f_outer in [f_inner for f_inner in input if wildcards.reg in f_inner]])
#     benchmark:
#         "benchmarks/db_combine/{reg}.tsv"
#     shell:
#        "gatk --java-options '-Xmx4g' GenomicsDBImport {params} -L {wildcards.reg} --genomicsdb-workspace-path {output}"

# rule db_genotype:
#     input:
#         "combi_out/{reg}_database/"
#     output:
#         "geno_out/{SAMP_COUNT}_sample_minified_{reg}.vcf.gz"
#     benchmark:
#         "benchmarks/db_genotype/{SAMP_COUNT}_{reg}.tsv"
#     shell:
#         "gatk --java-options '-Xmx4g' GenotypeGVCFs -R {REF} -V gendb://{input} -O {output}"


