import os

container: "docker://broadinstitute/gatk"
SAMP_COUNT = 2
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
        # expand("alignments/minified/2_sample/{sample}_minified.cram.crai", sample=samples)
        # expand("alignments/{reg}/{sample}_minified_{reg}.cram", sample=SAMPLES, reg=REGIONS)
        expand("alignments/{reg}/{sample}_minified_{reg}.cram.crai", sample=SAMPLES, reg=REGIONS),
        # expand("hc_out/{reg}/{sample}_minified_{reg}.g.vcf.gz", sample=SAMPLES, reg=REGIONS)
        "combi_out/my_database" # gendb combine
        # "geno_out/2_sample_minified_chr21.vcf.gz"

rule index_cram:
    input:
        expand("alignments/minified/{samps}_sample/{{sample}}_minified.cram", samps=SAMP_COUNT)
    output:
        expand("alignments/minified/{samps}_sample/{{sample}}_minified.cram.crai", samps=SAMP_COUNT)
    shell:
        "samtools index {input}"

rule split_cram:
    input:
        expand("alignments/minified/{samps}_sample/{{sample}}.cram", samps=SAMP_COUNT)
    output:
        "alignments/{reg}/{sample}_{reg}.cram"
    shell:
        "samtools view {input} {wildcards.reg} -T {REF} -O cram -o {output}"

rule index_split_cram:
    input:
        "alignments/{reg}/{sample}.cram"
    output:
        "alignments/{reg}/{sample}.cram.crai"
    shell:
        "samtools index {input}"

rule call_variants:
    input:
        "alignments/{reg}/{sample}.cram"
    output:
        "hc_out/{reg}/{sample}.g.vcf.gz"
    shell:
        "gatk HaplotypeCaller -I {input} "
        "-O {output} -R {REF} -L {wildcards.reg} -ERC GVCF"

rule db_combine: # this ones tricky
    input:
        expand("hc_out/{reg}/{sample}_minified_{reg}.g.vcf.gz", sample=SAMPLES, reg=REGIONS)
    output:
        directory("combi_out/{reg}_database/")
    params:
        lambda wildcards, input: ' '.join([f'-V {file}' for file in input])
    shell:
       "gatk GenomicsDBImport {params} --genomicsdb-workspace-path {output} -L chr21"

# rule db_genotype:
#     input:
#         "combi_out/my_database/"
#     output:
#         "geno_out/2_sample_minified_chr21.vcf.gz"
#     shell:
#         "gatk GenotypeGVCFs -R {ref} -V gendb://{input} -O {output}"


