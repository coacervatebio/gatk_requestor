import os
import sys

# Takes in auto-generated dir from benchmark.smk
# Parses all benchmark files in dir and collects them for analysis
def gather(dir):
    out_ = open(f"{dir}/BENCHMARK_REPORT.tsv", "w")
    outlines = [
        "op\tsample\tregion\ts\th:m:s\tmax_rss\tmax_vms\tmax_uss\tmax_pss\tio_in\tio_out\tmean_load\tcpu_time\n"
    ]
    for bm in os.listdir(dir):
        stats = open(f"{dir}/{bm}", "r").readlines()[1:]
        op_details = bm.rstrip(".tsv").replace("-", "\t")
        outlines.extend([f"{op_details}\t{line}" for line in stats])
    out_.writelines(outlines)
    out_.close()


# Call script with python gather_benchmarks.py [benchmark_dir]
gather(sys.argv[1])
