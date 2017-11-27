from __future__ import print_function
import sys
import gzip
import subprocess

# read bed files from stdin
bed_files = {int(x.rstrip().split('.')[-3].replace('bin', '')): x.rstrip() for x in sys.stdin}
ordered_bins = sorted(bed_files.keys())

# read summary data from sys.argv[1]
bin_data = {}

for line in open(sys.argv[1]):
    if line.startswith('bin'):
        continue

    bin_id, count, var_type = int(line.split('\t')[0]), int(line.split('\t')[1]), line.rstrip().split('\t')[2]

    if bin_id not in bin_data.keys():
        bin_data[bin_id] = {'ins': 0, 'del': 0}

    bin_data[bin_id][var_type] = count

# process files
merge_point = 0
merge_name = ''
for _bin in ordered_bins:

    no_indels = bin_data[_bin]['ins'] + bin_data[_bin]['del']

    if no_indels < 500:
        merge_point = _bin
        merge_name = bed_files[_bin].split('.bin')[0] + '.bin{}-{}.bed.gz'.format(_bin, ordered_bins[-1])
        print(merge_name)
        break

    print(bed_files[_bin])

# create merged file
merge_bin_bed = merge_name.replace('.gz', '')
with open(merge_bin_bed, 'w') as merge_out:
    for i in range(merge_point, ordered_bins[-1]+1):
        for coord_set in gzip.open(bed_files[i]):
            print(coord_set.rstrip(), file=merge_out)

subprocess.call('cat {} | sort -k1,1 -k2,2n | bgzip -c > {}'.format(merge_bin_bed, merge_name), shell=True)
subprocess.call('tabix -pbed {}'.format(merge_name), shell=True)