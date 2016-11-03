#!/usr/bin/env python

import argparse
from qsub import *

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-t1_r', help='Theta insertions range, format start,stop,step', required=True)
parser.add_argument('-t2_r', help='Theta deletions range, format start,stop,step', required=True)
parser.add_argument('-g1_r', help='Gamma insertions range, format start,stop,step', required=True)
parser.add_argument('-g2_r', help='Gamma deletions range, format start,stop,step', required=True)
parser.add_argument('-e1_r', help='Insertion polarisation error range, format start,stop,step', required=True)
parser.add_argument('-e2_r', help='Deletion polarisation error range, format start,stop,step', required=True)
parser.add_argument('-nrep', help='Number of replicates', default=500, type=int)
parser.add_argument('-out_dir', help='Directory for output', required=True)
args = parser.parse_args()

# variables
t1_r = [int(x) for x in args.t1_r.split(',')]
t2_r = [int(x) for x in args.t2_r.split(',')]
g1_r = [int(x) for x in args.g1_r.split(',')]
g2_r = [int(x) for x in args.g2_r.split(',')]
e1_r = [int(float(x) * 10) for x in args.e1_r.split(',')]
e2_r = [int(float(x) * 10) for x in args.e2_r.split(',')]
nrep = args.nrep
out_dir = args.out_dir
shell_dir = out_dir + 'shell_scripts/'
os.makedirs(shell_dir)
sub_out_dir = ''

# run all combos
output_list = []
jids = []
for t1 in range(t1_r[0], t1_r[1], t1_r[2]):
    for t2 in range(t2_r[0], t2_r[1], t2_r[2]):
        for g1 in range(g1_r[0], g1_r[1], g1_r[2]):
            for g2 in range(g2_r[0], g2_r[1], g2_r[2]):
                command_list = []
                for e1 in range(e1_r[0], e1_r[1], e1_r[2]):
                    e1 /= 10.0
                    for e2 in range(e2_r[0], e2_r[1], e2_r[2]):
                        e2 /= 10.0
                        sub_out_dir = out_dir + '_'.join([str(y) for y in [t1, t2, g1, g2, e1, e2]]) + '/'
                        k_wrapper_cmd = ('./kai_glemin_model_test.py '
                                         '-t1 ' + str(t1) + ' -t2 ' + str(t2) + ' '
                                         '-g1 ' + str(g1) + ' -g2 ' + str(g2) + ' '
                                         '-e1 ' + str(e1) + ' -e2 ' + str(e2) + ' '
                                         '-nrep 500 '
                                         '-out_dir ' + sub_out_dir + ' ')
                        #command_list.append(k_wrapper_cmd)
                        output_list.append(sub_out_dir + 'collated_run_outputs.txt')
                        job_id = 'glemin_' + '_'.join([str(y) for y in [t1, t2, g1, g2, e1, e2]]) + '.sh'
                        jids.append(job_id)
                        q_sub([k_wrapper_cmd], shell_dir + 'glemin_sim', jid=job_id, t=24)

# merge output
counter = 0
merged_out = out_dir + 'merged_glemin_outputs.txt'
cmds = []
for out in output_list:
    counter += 1
    if counter == 1:
        cmds.append('head -n 1 ' + out + ' >> ' + merged_out)
    cmds.append('grep -v ^t ' + out + ' >> ' + merged_out)

q_sub(cmds, out_dir + 'out_merging', hold=jids)

