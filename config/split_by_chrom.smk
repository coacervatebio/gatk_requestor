import os
from utils import file_to_sample

REF = "/config/reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"
SAMPLES = set([file_to_sample(fn) for fn in os.listdir("/datafiles/alignments/full/")])
REGIONS = [
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20",
    "chr21",
    "chr22",
]


rule all:
    input:
        # expand("/datafiles/alignments/full/{sample}.cram.crai", sample=SAMPLES)
        # expand("/datafiles/alignments/{reg}/{sample}_{reg}.cram", sample=SAMPLES, reg=REGIONS)
        # expand("/datafiles/alignments/{reg}/{sample}_{reg}.cram.crai", sample=SAMPLES, reg=REGIONS),
        # expand("/datafiles/hc_out/{reg}/{sample}_{reg}.g.vcf.gz", sample=SAMPLES, reg=REGIONS)
        # expand("/datafiles/combi_out/{reg}_database/", reg=REGIONS)
        # expand("/datafiles/geno_out/combined_{reg}.vcf.gz", reg=REGIONS)
        "/datafiles/gather_out/combined_gathered.vcf.gz"

rule index_cram:
    input:
        "/datafiles/alignments/full/{sample}.cram"
    output:
        temp("/datafiles/alignments/full/{sample}.cram.crai")
    shell:
        "samtools index {input}"

rule split_cram:
    input:
        alignments="/datafiles/alignments/full/{sample}.cram",
        indexes="/datafiles/alignments/full/{sample}.cram.crai"
    output:
        temp("/datafiles/alignments/{reg}/{sample}_{reg}.cram")
    benchmark:
        "/config/benchmarks/split_cram_{sample}_{reg}.tsv"
    shell:
        "samtools view {input.alignments} {wildcards.reg} -T {REF} -O cram -o {output}"

rule index_split_cram:
    input:
        "/datafiles/alignments/{reg}/{sample}_{reg}.cram"
    output:
        temp("/datafiles/alignments/{reg}/{sample}_{reg}.cram.crai")
    shell:
        "samtools index {input}"

rule call_variants:
    input:
        alignments="/datafiles/alignments/{reg}/{sample}.cram",
        indexes="/datafiles/alignments/{reg}/{sample}.cram.crai"
    output:
        called_vcf=temp("/datafiles/hc_out/{reg}/{sample}.g.vcf.gz"),
        index=temp("/datafiles/hc_out/{reg}/{sample}.g.vcf.gz.tbi")
    benchmark:
        "/config/benchmarks/call_vars_{sample}_{reg}.tsv"
    shell:
        "gatk --java-options '-Xmx4g' HaplotypeCaller -I {input.alignments} -O {output.called_vcf} -R {REF} -L {wildcards.reg} -ERC GVCF"

rule combine_region:
    input:
        set(expand("/datafiles/hc_out/{reg}/{sample}_{reg}.g.vcf.gz", sample=SAMPLES, reg=REGIONS))
    output:
        temp(directory("/datafiles/combi_out/{reg}_database/"))
    benchmark:
        "/config/benchmarks/combine_region_{reg}.tsv"
    params:
        lambda wildcards, input: ' '.join([f'-V {fn}' for fn in input if f'{wildcards.reg}.' in fn])
    shell:
       "gatk --java-options '-Xmx4g' GenomicsDBImport {params} -L {wildcards.reg} --genomicsdb-workspace-path {output}"

rule genotype:
    input:
        "/datafiles/combi_out/{reg}_database/"
    output:
        joint_vcf=temp("/datafiles/geno_out/combined_{reg}.vcf.gz"),
        index=temp("/datafiles/geno_out/combined_{reg}.vcf.gz.tbi")
    benchmark:
        "/config/benchmarks/genotype_{reg}.tsv"
    shell:
        "gatk --java-options '-Xmx4g' GenotypeGVCFs -R {REF} -V gendb://{input} -O {output.joint_vcf}"

rule gather_vcfs:
    input:
        expand("/datafiles/geno_out/combined_{reg}.vcf.gz", reg=REGIONS)
    output:
        "/datafiles/gather_out/combined_gathered.vcf.gz" # Needs a better name
    benchmark:
        "/config/benchmarks/gather_vcfs.tsv"
    params:
        lambda wildcards, input: ' '.join([f'-I {fn}' for fn in input])
    shell:
        "gatk --java-options '-Xmx4g' GatherVcfs {params} -O {output}"