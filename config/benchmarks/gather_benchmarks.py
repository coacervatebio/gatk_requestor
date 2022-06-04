import os

def gather(dir):
    outlines = ['sample\ts\th:m:s\tmax_rss\tmax_vms\tmax_uss\tmax_pss\tio_in\tio_out\tmean_load\tcpu_time\n']
    for bm in os.listdir(dir):
        stats = open(f'{dir}/{bm}', 'r').readlines()[1:]
        outlines.extend([f'{bm}\t{line}' for line in stats])
    
    with open(f'{dir}/BENCHMARK_REPORT.tsv', 'w') as out_:
        out_.writelines(outlines)

gather('2_sample_all_chrom')