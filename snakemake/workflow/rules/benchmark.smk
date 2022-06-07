from datetime import datetime

configfile: "/mnt/config/config.yml"

module to_benchmark:
    snakefile: 
        "/mnt/workflow/rules/split_by_chrom.smk"
    config:
        config

use rule * from to_benchmark as *_benchmarked

benchmark_dir = f"/mnt/benchmarks/{config['in_dir']}-{datetime.now().strftime('%Y%m%d-%H%M')}"

use rule index_cram from to_benchmark as index_cram_benchmarked with:
    benchmark:
        f"{benchmark_dir}/index_cram_{{sample}}.tsv"
use rule split_cram from to_benchmark as split_cram_benchmarked with:
    benchmark:
        f"{benchmark_dir}/split_cram_{{sample}}_{{reg}}.tsv"
use rule call_variants from to_benchmark as call_variants_benchmarked with:
    benchmark:
        f"{benchmark_dir}/call_vars_{{sample}}_{{reg}}.tsv"
use rule combine_region from to_benchmark as combine_region_benchmarked with:
    benchmark:
        f"{benchmark_dir}/combine_region_{{reg}}.tsv"
use rule genotype from to_benchmark as genotype_benchmarked with:
    benchmark:
        f"{benchmark_dir}/genotype_{{reg}}.tsv"
use rule gather_vcfs from to_benchmark as gather_vcfs_benchmarked with:
    benchmark:
        f"{benchmark_dir}/gather_vcfs.tsv"

rule all:
    input:
        rules.all_benchmarked.input