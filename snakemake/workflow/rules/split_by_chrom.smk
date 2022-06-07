import os
from utils import dir_to_samples

configfile: "/mnt/config/config.yml"
samples = dir_to_samples(f"/mnt/results/{config['in_dir']}/alignments/full")

rule all:
    input:
        # expand("/mnt/results/{in_dir}/alignments/{sample}.cram.crai", in_dir=config['in_dir'], sample=SAMPLES)
        # expand("/mnt/results/{in_dir}/alignments/{reg}/{sample}_{reg}.cram", in_dir=config['in_dir'], sample=samples, reg=config['regs'])
        # expand("/mnt/results/{in_dir}/alignments/{reg}/{sample}_{reg}.cram.crai", in_dir=config['in_dir'], sample=samples, reg=config['regs'])
        expand("/mnt/results/{in_dir}/hc_out/{reg}/{sample}_{reg}.g.vcf.gz", in_dir=config['in_dir'], sample=samples, reg=config['regs'])
        # expand("/mnt/results/{in_dir}/combi_out/{reg}_database/", in_dir=config['in_dir'], reg=config['regs'])
        # expand("/mnt/results/{in_dir}/geno_out/combined_{reg}.vcf.gz", in_dir=config['in_dir'], reg=config['regs'])
        # f"/mnt/results/{config['in_dir']}/gather_out/project_{config['in_dir']}_output.vcf.gz"

rule index_cram:
    input:
        f"/mnt/results/{config['in_dir']}/alignments/full/{{sample}}.cram"
    output:
        temp(f"/mnt/results/{config['in_dir']}/alignments/full/{{sample}}.cram.crai")
    shell:
        "samtools index {input}"

rule split_cram:
    input:
        alignments=f"/mnt/results/{config['in_dir']}/alignments/full/{{sample}}.cram",
        indexes=f"/mnt/results/{config['in_dir']}/alignments/full/{{sample}}.cram.crai"
    output:
        f"/mnt/results/{config['in_dir']}/alignments/{{reg}}/{{sample}}_{{reg}}.cram"
    shell:
        "samtools view {input.alignments} {wildcards.reg} -T {config[ref]} -O cram -o {output}"

rule index_split_cram:
    input:
        f"/mnt/results/{config['in_dir']}/alignments/{{reg}}/{{sample}}_{{reg}}.cram"
    output:
        temp(f"/mnt/results/{config['in_dir']}/alignments/{{reg}}/{{sample}}_{{reg}}.cram.crai")
    shell:
        "samtools index {input}"

rule call_variants:
    input:
        alignments=f"/mnt/results/{config['in_dir']}/alignments/{{reg}}/{{sample}}_{{reg}}.cram",
        indexes=f"/mnt/results/{config['in_dir']}/alignments/{{reg}}/{{sample}}_{{reg}}.cram.crai"
    output:
        called_vcf=f"/mnt/results/{config['in_dir']}/hc_out/{{reg}}/{{sample}}_{{reg}}.g.vcf.gz",
        index=temp(f"/mnt/results/{config['in_dir']}/hc_out/{{reg}}/{{sample}}_{{reg}}.g.vcf.gz.tbi")
    shell:
        "gatk --java-options '-Xmx4g' HaplotypeCaller -I {input.alignments} -O {output.called_vcf} -R {config[ref]} -L {wildcards.reg} -ERC GVCF"

rule combine_region:
    input:
        gvcfs=set(expand(f"/mnt/results/{config['in_dir']}/hc_out/{{reg}}/{{sample}}_{{reg}}.g.vcf.gz", sample=samples, reg=config['regs'])),
        indexes=set(expand(f"/mnt/results/{config['in_dir']}/hc_out/{{reg}}/{{sample}}_{{reg}}.g.vcf.gz.tbi", sample=samples, reg=config['regs']))
    output:
        temp(directory(f"/mnt/results/{config['in_dir']}/combi_out/{{reg}}_database/"))
    params:
        lambda wildcards, input: ' '.join([f'-V {fn}' for fn in input.gvcfs if f'{wildcards.reg}.' in fn])
    shell:
       "gatk --java-options '-Xmx4g' GenomicsDBImport {params} -L {wildcards.reg} --genomicsdb-workspace-path {output}"

rule genotype:
    input:
        f"/mnt/results/{config['in_dir']}/combi_out/{{reg}}_database/"
    output:
        joint_vcf=temp(f"/mnt/results/{config['in_dir']}/geno_out/combined_{{reg}}.vcf.gz"),
        index=temp(f"/mnt/results/{config['in_dir']}/geno_out/combined_{{reg}}.vcf.gz.tbi")
    shell:
        "gatk --java-options '-Xmx4g' GenotypeGVCFs -R {config[ref]} -V gendb://{input} -O {output.joint_vcf}"

rule gather_vcfs:
    input:
        expand("/mnt/results/{in_dir}/geno_out/combined_{reg}.vcf.gz", in_dir=config['in_dir'], reg=config['regs'])
    output:
        f"/mnt/results/{config['in_dir']}/gather_out/project_{config['in_dir']}_output.vcf.gz"
    params:
        lambda wildcards, input: ' '.join([f'-I {fn}' for fn in input])
    shell:
        "gatk --java-options '-Xmx4g' GatherVcfs {params} -O {output}"