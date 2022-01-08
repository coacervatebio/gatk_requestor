import os

container: "docker://broadinstitute/gatk"
SAMP_COUNT = 6
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
        expand("alignments/{reg}/{sample}_minified_{reg}.cram", sample=SAMPLES, reg=REGIONS)
        # expand("alignments/chr21/{sample}_minified_chr21.cram.crai", sample=samples)
        # expand("hc_out/chr21/{sample}_minified_chr21.g.vcf.gz", sample=samples)
        # "combi_out/my_database" # gendb combine
        # "geno_out/2_sample_minified_chr21.vcf.gz"

# rule index_cram:
#     input:
#         "alignments/minified/2_sample/{sample}_minified.cram"
#     output:
#         "alignments/minified/2_sample/{sample}_minified.cram.crai"
#     shell:
#         "samtools index {input}"

rule split_cram:
    input:
        expand("alignments/minified/{samps}_sample/{{sample}}.cram", samps=SAMP_COUNT)
    output:
        "alignments/{reg}/{sample}_{reg}.cram"
    shell:
        "samtools view {input} {wildcards.reg} -T {REF} -O cram -o {output}"

# rule index_split_cram:
#     input:
#         "alignments/chr21/{sample}.cram"
#     output:
#         "alignments/chr21/{sample}.cram.crai"
#     shell:
#         "samtools index {input}"

# rule call_variants:
#     input:
#         "alignments/chr21/{sample}.cram"
#     output:
#         "hc_out/chr21/{sample}.g.vcf.gz"
#     shell:
#         "gatk HaplotypeCaller -I {input} "
#         "-O {output} -R {ref} -L chr21 -ERC GVCF"

# rule db_combine:
#     input:
#         expand("hc_out/chr21/{sample}_minified_chr21.g.vcf.gz", sample=samples)
#     output:
#         directory("combi_out/my_database/")
#     params:
#         lambda wildcards, input: ' '.join([f'-V {file}' for file in input])
#     shell:
#        "gatk GenomicsDBImport {params} --genomicsdb-workspace-path {output} -L chr21"

# rule db_genotype:
#     input:
#         "combi_out/my_database/"
#     output:
#         "geno_out/2_sample_minified_chr21.vcf.gz"
#     shell:
#         "gatk GenotypeGVCFs -R {ref} -V gendb://{input} -O {output}"


