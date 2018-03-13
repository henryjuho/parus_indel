from __future__ import print_function
import sys
import subprocess
import pysam
import os
sys.path.insert(0, os.getenv('HOME') + '/parus_indel/anavar_analyses')
from sel_vs_neu_anavar import sfs2counts
from window_sel_vs_neu_anavar import window_call_sites

sfs_cmd = 'tabix -h {vcf} {chromo}:{start}-{end} | ~/sfs_utils/vcf2raw_sfs.py -mode {mode}'

vcf = sys.argv[1]
call_fasta = pysam.FastaFile(sys.argv[2])

for line in sys.stdin:

    out_line = line.rstrip()

    coords = line.split()

    if coords[0] == 'chrZ':
        continue

    for indel_type in ['ins', 'del']:
        sfs = subprocess.Popen(sfs_cmd.format(chromo=coords[0], start=coords[1], end=coords[2],
                                              vcf=vcf, mode=indel_type), shell=True, stdout=subprocess.PIPE
                               ).communicate()[0].split('\n')[:-1]
        sfs_counted = ','.join([str(x) for x in sfs2counts(sfs, 20)])
        out_line += '\t' + sfs_counted

    # get sel call sites
    call = window_call_sites(call_fasta, None, (coords[0], int(coords[1]), int(coords[2])))
    out_line += '\t' + str(call)

    print(out_line)
