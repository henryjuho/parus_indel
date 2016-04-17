#!/usr/bin/env python

import argparse
import subprocess

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-maf_list', help='List file of maf files to calculate distance matrix from', required=True)
args = parser.parse_args()

# variables
maf_list = [maf_path.rstrip('\n') for maf_path in open(args.maf_list)]
out_dir = maf_list[0][:maf_list[0].rfind('/')+1]

# Commence looping
job_list = []
tree_list = []
for maf in maf_list:
    maf_prefix = maf.rstrip('.maf')
    option_file = maf_prefix + '.maffilter.divergence.options'

    # construct commandline
    with open(option_file, 'w') as options:
        options.write('maf.filter='
                      'DistanceEstimation('
                      'method=ml,'
                      'model=K80(kappa=2),'
                      'rate=Gamma(n=4, alpha=0.5),'
                      'parameter_estimation=initial,'
                      'max_freq_gaps=0.33,'
                      'gap_option=no_gap,'
                      'profiler=std,'
                      'message_handler=std),'
                      'DistanceBasedPhylogeny('
                      'method=bionj,'
                      'dist_mat=MLDistance),'
                      'OutputTrees('
                      'tree=BioNJ,'
                      'file=' + maf_prefix + '.trees.dnd,'
                      'compression=none)')

    tree_list.append(maf_prefix + '.trees.dnd')

    maffilter_cmd = ('"maffilter '
                     'input.file=' + maf_prefix + '.maf '  # Input maf file
                     'input.file.compression=none '
                     'output.log=' + maf_prefix + '.maffilter.log '  # Output log file
                     'param=' + option_file + '"')

    # write shell script and submit
    jid = maf_prefix[maf_prefix.rfind('/')+1:] + '.divergence.sh'
    job_list.append(jid)
    qsub_cmd = ('python qsub_gen.py '
                '-cmd ' + maffilter_cmd + ' '
                '-o ' + maf_prefix + '_divergence '
                '-jid ' + jid + ' '
                '-OM q')

    subprocess.call(qsub_cmd, shell=True)

# Calculate mean divergences
dnd2div = '"~/dnd2div.py '
for dnd in tree_list:
    dnd2div += '-dnd ' + dnd + ' '
dnd2div += '"'

hold = ','.join(job_list)

dnd_qsub = ('python qsub_gen.py '
            '-cmd ' + dnd2div + ' '
            '-o ' + out_dir + 'mean_divergence '
            '-jid mean_divergence.sh '
            '-hold ' + hold + ' '
            '-OM q')

subprocess.call(dnd_qsub, shell=True)
