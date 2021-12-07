REF = "resources/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"
SAMPLES = ['176', '183']
container: "docker://broadinstitute/gatk"


# rule docker_test:
#     output:
#         "debug/hello.txt"
#     shell:
#         'echo hello > debug/hello.txt'

rule target:
    input:
        expand("alignments/chr21/{sample}_chr21.cram", sample=SAMPLES) #debug split

rule split_cram:
    input:
        "alignments/exomes/{sample}.cram"
    output:
        "alignments/chr21/{sample}_chr21.cram"
    shell:
        "samtools view {input} chr21 -T {REF} "
        "-O cram -o {output}"