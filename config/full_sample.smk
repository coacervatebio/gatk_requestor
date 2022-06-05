import os
from utils import file_to_sample

REF = "/config/reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"
SAMPLES = set([file_to_sample(fn) for fn in os.listdir("/datafiles/alignments/full/")])

rule all:
    input:
        # expand("/datafiles/alignments/full/{sample}.cram.crai", sample=SAMPLES)
        # expand("/datafiles/hc_out/{sample}.g.vcf.gz", sample=SAMPLES)
        # ("/datafiles/combi_out/geno_db/"
        "/datafiles/geno_out/combined.vcf.gz"

rule index_cram:
    input:
        "/datafiles/alignments/full/{sample}.cram"
    output:
        temp("/datafiles/alignments/full/{sample}.cram.crai")
    shell:
        "samtools index {input}"

rule call_variants:
    input:
        alignments="/datafiles/alignments/full/{sample}.cram",
        indexes="/datafiles/alignments/full/{sample}.cram.crai"
    output:
        called_vcf=temp("/datafiles/hc_out/{sample}.g.vcf.gz"),
        index=temp("/datafiles/hc_out/{sample}.g.vcf.gz.tbi")
    benchmark:
        "/config/benchmarks/call_vars_{sample}.tsv"
    shell:
        "gatk --java-options '-Xmx4g' HaplotypeCaller -I {input.alignments} -O {output.called_vcf} -R {REF} -ERC GVCF"

rule combine_samples:
    input:
        expand("/datafiles/hc_out/{sample}.g.vcf.gz", sample=SAMPLES),
        expand("/datafiles/hc_out/{sample}.g.vcf.gz.tbi", sample=SAMPLES)
    output:
        temp(directory("/datafiles/combi_out/geno_db/"))
    benchmark:
        "/config/benchmarks/combine_samples.tsv"
    params:
        lambda wildcards, input: ' '.join([f'-V {fn}' for fn in input])
    shell:
       "gatk --java-options '-Xmx4g' GenomicsDBImport {params} --genomicsdb-workspace-path {output}"

rule genotype:
    input:
        "/datafiles/combi_out/geno_db/"
    output:
        joint_vcf=temp("/datafiles/geno_out/combined.vcf.gz"),
        index=temp("/datafiles/geno_out/combined.vcf.gz.tbi")
    benchmark:
        "/config/benchmarks/genotype.tsv"
    shell:
        "gatk --java-options '-Xmx4g' GenotypeGVCFs -R {REF} -V gendb://{input} -O {output.joint_vcf}"