import os
from utils import dir_to_samples

configfile: "/data/config/config.yml"
samples = dir_to_samples(f"/data/results/alignments/full")

rule all:
    input:
        # expand("/data/results/alignments/{sample}.cram.crai", sample=SAMPLES),
        # expand("/data/results/alignments/{reg}/{sample}_{reg}.cram", sample=samples, reg=config['regs']),
        # expand("/data/results/alignments/{reg}/{sample}_{reg}.cram.crai", sample=samples, reg=config['regs']),
        # expand("/data/results/hc_out/{reg}/{sample}_{reg}.g.vcf.gz", sample=samples, reg=config['regs']),
        # expand("/data/results/combi_out/{reg}_database/", reg=config['regs']),
        # expand("/data/results/geno_out/combined_{reg}.vcf.gz", reg=config['regs']),
        f"/data/results/gather_out/project_output.vcf.gz",

rule index_cram:
    input:
        f"/data/results/alignments/full/{{sample}}.cram"
    output:
        temp(f"/data/results/alignments/full/{{sample}}.cram.crai")
    shell:
        "samtools index {input}"

rule split_cram:
    input:
        alignments=f"/data/results/alignments/full/{{sample}}.cram",
        indexes=f"/data/results/alignments/full/{{sample}}.cram.crai"
    output:
        temp(f"/data/results/alignments/{{reg}}/{{sample}}_{{reg}}.cram")
    shell:
        "samtools view {input.alignments} {wildcards.reg} -T {config[ref]} -O cram -o {output}"

rule index_split_cram:
    input:
        f"/data/results/alignments/{{reg}}/{{sample}}_{{reg}}.cram"
    output:
        temp(f"/data/results/alignments/{{reg}}/{{sample}}_{{reg}}.cram.crai")
    shell:
        "samtools index {input}"

rule call_variants:
    input:
        alignments=f"/data/results/alignments/{{reg}}/{{sample}}_{{reg}}.cram",
        indexes=f"/data/results/alignments/{{reg}}/{{sample}}_{{reg}}.cram.crai"
    output:
        called_vcf=temp(f"/data/results/hc_out/{{reg}}/{{sample}}_{{reg}}.g.vcf.gz"),
        index=temp(f"/data/results/hc_out/{{reg}}/{{sample}}_{{reg}}.g.vcf.gz.tbi")
    shell:
        "gatk --java-options '-Xmx4g' HaplotypeCaller -I {input.alignments} -O {output.called_vcf} -R {config[ref]} -L {wildcards.reg} -ERC GVCF"

rule combine_region:
    input:
        gvcfs=set(expand(f"/data/results/hc_out/{{reg}}/{{sample}}_{{reg}}.g.vcf.gz", sample=samples, reg=config['regs'])),
        indexes=set(expand(f"/data/results/hc_out/{{reg}}/{{sample}}_{{reg}}.g.vcf.gz.tbi", sample=samples, reg=config['regs']))
    output:
        temp(directory(f"/data/results/combi_out/{{reg}}_database/"))
    params:
        lambda wildcards, input: ' '.join([f'-V {fn}' for fn in input.gvcfs if f'{wildcards.reg}.' in fn])
    shell:
       "gatk --java-options '-Xmx4g' GenomicsDBImport {params} -L {wildcards.reg} --genomicsdb-workspace-path {output}"

rule genotype:
    input:
        f"/data/results/combi_out/{{reg}}_database/"
    output:
        joint_vcf=temp(f"/data/results/geno_out/combined_{{reg}}.vcf.gz"),
        index=temp(f"/data/results/geno_out/combined_{{reg}}.vcf.gz.tbi")
    shell:
        "gatk --java-options '-Xmx4g' GenotypeGVCFs -R {config[ref]} -V gendb://{input} -O {output.joint_vcf}"

rule gather_vcfs:
    input:
        expand("/data/results/geno_out/combined_{reg}.vcf.gz", reg=config['regs'])
    output:
        f"/data/results/gather_out/project_output.vcf.gz"
    params:
        lambda wildcards, input: ' '.join([f'-I {fn}' for fn in input])
    shell:
        "gatk --java-options '-Xmx4g' GatherVcfs {params} -O {output}"
