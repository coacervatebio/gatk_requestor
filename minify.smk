import os

samples = glob_wildcards("alignments/full_exomes/{samp}.cram").samp
REF = "resources/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"
container: "docker://broadinstitute/gatk"

rule all:
    input: 
        expand("alignments/minified/{samp}_minified.cram", samp=samples)

rule minify:
    input: 
        "alignments/full_exomes/{id}.cram"
    output:
        temp("alignments/minified/{id}_minified.sam")
    shell:
        "samtools view -h {input} -T {REF} | python3 minify.py {wildcards.id}"
    
rule to_cram:
    input:
        "alignments/minified/{id}_minified.sam"
    output:
        "alignments/minified/{id}_minified.cram"
    shell:
        "samtools view -h {input} -T {REF} -O cram -o {output}"