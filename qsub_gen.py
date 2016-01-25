import argparse
import subprocess

'''Python utility to enable quick generation of job files for qsub'''

# Command line arguments

parser = argparse.ArgumentParser()
parser.add_argument('-cmd', '--command_line', help='Command lines to run in job, require quotation marks',
                    action='append', required=True)
parser.add_argument('-o', '--out', help='Output directory and file prefix for error and output files',
                    required=True)
parser.add_argument('-mo', '--modules', help='Modules to load eg java, python, gatk',
                    action='append', choices=['java', 'python', 'gatk', 'R'], default=[])
parser.add_argument('-t', '--time', help='Time required for job in hours, default 8, max 168', default=8)
parser.add_argument('-rmem', '--real_memory', help='Real memory required for job', default=2)
parser.add_argument('-mem', '--virtual_memory', help='Virtual memory required for job', default=6)
parser.add_argument('-OM', '--output_mode', choices=['p', 'w', 'q'], default='p',
                    help='Determines if output is printed to screen or written to file')
parser.add_argument('-hold', '--hold_jobs', help='Comma separated list of jobs to wait for before executing script',
                    default=False)
parser.add_argument('-jid', '--job_ID', help='Name for job, if not specified follows default of <qsub_filename_job.sh>',
                    required=False, default='DEFAULT')
parser.add_argument('-tr', '--threads', help='Number of threads required for job', default=1)
args = parser.parse_args()

cmd_line = args.command_line  # list of commands to submit in batch job
mods = args.modules
run_time = '#$-l h_rt='+str(args.time)+':00:00\n'
memory = '#$-l mem='+str(args.virtual_memory)+'G\n#$-l rmem='+str(args.real_memory)+'G\n'
threads = args.threads
file_pos = args.out.rfind('/')+1
if args.job_ID == 'DEFAULT':
    output_name = args.out[0:file_pos] + 'qsub_' + args.out[file_pos:] + '_job.sh'
else:
    output_name = args.out[0:file_pos] + args.job_ID
outs = '\n#$-N ' + output_name[output_name.rfind('/')+1:] + '\n#$-o '+args.out+'.out\n#$-e '+args.out+'.error\n'
out_mode = args.output_mode
hold = args.hold_jobs

# module dictionary

available_modules = {'python': 'apps/python/2.7', 'java': 'apps/java/1.7', 'gatk': 'apps/binapps/GATK', 'R': 'apps/R'}

# Write or print qsub job

if out_mode == 'w' or out_mode == 'q':
    output = open(output_name, 'w')
    output.write('#!/bin/bash\n')
    for m in mods:
        output.write('#!module load  '+available_modules[m]+'\n')
    output.write('\n#$-l arch=intel*\n' + run_time + memory + '\n')
    if threads != 1:
        output.write('#$-pe openmp ' + str(threads) + '\n')
    output.write(outs + '\n')
    if hold is not False:
        output.write('#$-hold_jid ' + args.hold_jobs + '\n\n')
    for cmd in cmd_line:
        output.write(cmd + '\n')
    output.close()
    if out_mode == 'q':
        qsub_cmd = 'qsub '+output_name
        subprocess.call(qsub_cmd, shell=True)
elif out_mode == 'p':
    print('#!/bin/bash')
    for m in mods:
        print('#!module load  '+available_modules[m])
    print('\n#$-l arch=intel*\n' + run_time + memory)
    if threads != 1:
        print('#$-pe openmp ' + str(threads))
    print(outs)
    if hold is not False:
        print('#$-hold_jid ' + args.hold_jobs + '\n')
    for cmd in cmd_line:
        print(cmd)
