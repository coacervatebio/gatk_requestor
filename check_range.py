import sys

cur_chrom = ''
end_pos = ''
all_out = []
for line in sys.stdin:
    line = line.rstrip()
    if line.startswith('@'):
        continue
    else:
        line_list = line.split('\t')
        chrom = line_list[2]
        if '_' in chrom or '-' in chrom: continue
        pos = line_list[3]
        if cur_chrom != chrom:
            all_out.append(end_pos)
            cur_chrom = chrom
            all_out.append(chrom)
            all_out.append(pos)
        else:
            end_pos = pos

print('/'.join(all_out[1:]))