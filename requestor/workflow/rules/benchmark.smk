from datetime import datetime

configfile: "/mnt/config/config.yml"

module to_benchmark:
    snakefile: 
        "/mnt/workflow/rules/split_by_chrom.smk"
    config:
        config

use rule * from to_benchmark as *_benchmarked

# Unique dir for each benchmarking run
benchmark_dir = f"/mnt/benchmarks/{config['in_dir']}-{datetime.now().strftime('%Y%m%d-%H%M')}"

# Actual filenames are 3 parts:
# operation - sample (or all) - region (or full)
# This provides better column data via gather_benchmarks.py
use rule index_cram from to_benchmark as index_cram_benchmarked with:
    benchmark:
        f"{benchmark_dir}/index_cram-{{sample}}-full.tsv"
use rule split_cram from to_benchmark as split_cram_benchmarked with:
    benchmark:
        f"{benchmark_dir}/split_cram-{{sample}}-{{reg}}.tsv"
use rule call_variants from to_benchmark as call_variants_benchmarked with:
    benchmark:
        f"{benchmark_dir}/call_vars-{{sample}}-{{reg}}.tsv"
use rule combine_region from to_benchmark as combine_region_benchmarked with:
    benchmark:
        f"{benchmark_dir}/combine-all-{{reg}}.tsv"
use rule genotype from to_benchmark as genotype_benchmarked with:
    benchmark:
        f"{benchmark_dir}/genotype-all-{{reg}}.tsv"
use rule gather_vcfs from to_benchmark as gather_vcfs_benchmarked with:
    benchmark:
        f"{benchmark_dir}/gather-all-full.tsv"

rule all:
    input:
        rules.all_benchmarked.input
