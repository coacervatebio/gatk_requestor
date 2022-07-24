from pathlib import Path, PurePath
from requestor import data

als_path = Path('/home/vagrant/host_shared/snakemake/results/2_sample/alignments')

for task in data(als_path):
    for k,v in task.data.items():
        print(k, ':', type(v))
    break