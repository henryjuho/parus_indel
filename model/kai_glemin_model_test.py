#!/usr/bin/env python

import argparse
from qsub import *
import sys

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-t1', help='Theta insertions', required=True, type=float)
parser.add_argument('-t2', help='Theta deletions', required=True, type=float)
parser.add_argument('-g1', help='Gamma insertions', required=True, type=float)
parser.add_argument('-g2', help='Gamma deletions', required=True, type=float)
parser.add_argument('-e1', help='Insertion polarisation error', required=True, type=float)
parser.add_argument('-e2', help='Deletion polarisation error', required=True, type=float)
parser.add_argument('-nrep', help='Number of replicates', default=500, type=int)
parser.add_argument('-out_dir', help='Directory for output and log files', required=True)
parser.add_argument('-sub', help='If specified will submit script to cluster', action='store_true', default=False)
parser.add_argument('-evolgen', help='If specified will run on lab queue', action='store_true', default=False)
args = parser.parse_args()

# submission loop
if args.sub is True:
    command_line = [' '.join([x for x in sys.argv if x != '-sub'])]
    q_sub(command_line, args.out_dir + 'glemin_sim', evolgen=args.evolgen)
    sys.exit('Script submitted')


# functions
def convergence(v1, v2):
    if abs(v1 - v2) / (0.5 * (v1 + v2)) < 10**-5:
        return True
    else:
        return False

# variables
top_out_dir = args.out_dir
out_dir = top_out_dir + 'out/'
log_dir = top_out_dir + 'log/'
nrep = args.nrep
n = 20

try:
    os.makedirs(out_dir)
    os.makedirs(log_dir)
except OSError:
    pass

t1 = args.t1
t2 = args.t2
g1 = args.g1
g2 = args.g2
e1 = args.e1
e2 = args.e2
tmin = t1 / 100.0
tmax = t2 * 100.0
gmin = -500
gmax = 100
emin = 0
emax = 1
nrun = 50
imprftol = 10**-10
nnoimp = 3
maximp = 50

# call kai's script

# Usage: model_sim result_folder log_folder n nrep #
# "theta_1, theta_2, gamma_1, gamma_2, e_1, e_2" #
# "theta_min, theta_max, gamma_min, gamma_max, e_min, e_max" #
# nrun imprftol nnoimp maximp #

k_command = ('model_sim ' + out_dir + ' ' + log_dir + ' ' + str(n) + ' ' + str(nrep) + ' "' +
             str(t1) + ', ' + str(t2) + ', ' +
             str(g1) + ', ' + str(g2) + ', ' +
             str(e1) + ', ' + str(e2) + '" "' +
             str(tmin) + ', ' + str(tmax) + ', ' +
             str(gmin) + ', ' + str(gmax) + ', ' +
             str(emin) + ', ' + str(emax) + '" ' +
             str(nrun) + ' ' + str(imprftol) + ' ' + str(nnoimp) + ' ' + str(maximp))

subprocess.call(k_command, shell=True)

# process output
with open(top_out_dir + 'collated_run_outputs.txt', 'w') as gathered_out:
    gathered_out.write('t1\tt2\tg1\tg2\te1\te2\tt1_p\tt2_p\tg1_p\tg2_p\te1_p\te2_p\n')

    for k_glem_out in os.listdir(out_dir):
        if k_glem_out.startswith('r'):
            data = [line.split('\t') for line in open(out_dir + k_glem_out).readlines()[1:3]]

            # check exit code
            if not data[0][2] == data[1][2] == '3':
                continue

            # check convergence
            if convergence(float(data[0][9]), float(data[1][9])) is False:
                continue

            # todo check closeness to limits

            # extract
            t1_p, t2_p, g1_p, g2_p, e1_p, e2_p = data[0][3], data[0][4], data[0][5], data[0][6], data[0][7], data[0][8]

            # write output
            gathered_out.write('\t'.join([str(x) for x in [t1, t2, g1, g2, e1, e2]] +
                                         [t1_p, t2_p, g1_p, g2_p, e1_p, e2_p]) + '\n')
