rule index_cram:
    input:
        f"/data/results/alignments/full/{{sample}}.cram"
    output:
        temp(f"/data/results/alignments/full/{{sample}}.cram.crai")
    shell:
        "samtools index {input}"

rule split_cram:
    input:
        alignment=f"/data/results/alignments/full/{{sample}}.cram",
        index=f"/data/results/alignments/full/{{sample}}.cram.crai"
    output:
        temp(f"/data/results/alignments/{{reg}}/{{sample}}_{{reg}}.cram")
    params:
        reference=config['ref']
    shell:
        "samtools view {input.alignment} {wildcards.reg} -T {params.reference} -O cram -o {output}"

rule index_split_cram:
    input:
        f"/data/results/alignments/{{reg}}/{{sample}}_{{reg}}.cram"
    output:
        temp(f"/data/results/alignments/{{reg}}/{{sample}}_{{reg}}.cram.crai")
    shell:
        "samtools index {input}"