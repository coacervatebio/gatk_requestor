import os

def gather(dir):
    out_ = open(f'{dir}/BENCHMARK_REPORT.tsv', 'w')
    outlines = ['sample\ts\th:m:s\tmax_rss\tmax_vms\tmax_uss\tmax_pss\tio_in\tio_out\tmean_load\tcpu_time\n']
    for bm in os.listdir(dir):
        stats = open(f'{dir}/{bm}', 'r').readlines()[1:]
        outlines.extend([f'{bm}\t{line}' for line in stats])
        out_.writelines(outlines)
    out_.close()

gather('2_sample_all_chrom')