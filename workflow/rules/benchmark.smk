module to_benchmark:
    snakefile: "/config/split_by_chrom.smk"
    # config: config["to_benchmark"]

use rule * from to_benchmark as *_benchmarked

use rule index_cram from to_benchmark as index_cram_benchmarked with:
    benchmark:
        "/config/benchmarks/index_cram_{sample}.tsv"
use rule split_cram from to_benchmark as split_cram_benchmarked with:
    benchmark:
        "/config/benchmarks/split_cram_{sample}_{reg}.tsv"
use rule call_variants from to_benchmark as call_variants_benchmarked with:
    benchmark:
        "/config/benchmarks/call_vars_{sample}_{reg}.tsv"
use rule combine_region from to_benchmark as combine_region_benchmarked with:
    benchmark:
        "/config/benchmarks/combine_region_{reg}.tsv"
use rule genotype from to_benchmark as genotype_benchmarked with:
    benchmark:
        "/config/benchmarks/genotype_{reg}.tsv"
use rule gather_vcfs from to_benchmark as gather_vcfs_benchmarked with:
    benchmark:
        "/config/benchmarks/gather_vcfs.tsv"

rule all:
    input:
        rules.all_benchmarked.input