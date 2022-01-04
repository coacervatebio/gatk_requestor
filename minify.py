import sys
import argparse

parser = argparse.ArgumentParser(description="Sample ID")
parser.add_argument("sample", type=str)
args = parser.parse_args()

counts = {
    "chr1" : 24724,
    "chr2" : 24295,
    "chr3" : 19950,
    "chr4" : 19127,
    "chr5" : 18085,
    "chr6" : 17089,
    "chr7" : 15882,
    "chr8" : 14627,
    "chr9" : 14027,
    "chr10" : 13537,
    "chr11" : 13445,
    "chr12" : 13234,
    "chr13" : 11414,
    "chr14" : 10636,
    "chr15" : 10033,
    "chr16" : 8882,
    "chr17" : 7877,
    "chr18" : 7611,
    "chr19" : 6381,
    "chr20" : 6243,
    "chr21" : 4694,
    "chr22" : 4969,
    "chrX" : 15491,
    "chrY" : 5777,
}

out = []

for line in sys.stdin:
    line = line.rstrip()
    if line.startswith('@'):
        out.append(line)
    else:
        chrom = line.split('\t')[2]

        try:
            if counts[chrom] > 0:
                out.append(line)
                counts[chrom] -= 1
        except KeyError:
            continue

with open(f"alignments/minified/{args.sample}_minified.sam", 'w') as out_:
    out_.write('\n'.join(out))