import argparse
import subprocess

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-bam_list', help='List of bam files to genotype', required=True)
parser.add_argument('-ref', help='Reference genome location', required=True)
parser.add_argument('-out', help='Output_directory and prefix', required=True)
args = parser.parse_args()

# variables
bams = args.bam_list
ref_genome = args.ref
output = args.out
index = ref_genome+'.fai'

# generate chromosome list
chromosome_scaffold_list = [line.split()[0] for line in open(index)]
chromosome_list = [entry for entry in chromosome_scaffold_list if entry.startswith('c')]
chrZ = chromosome_list[len(chromosome_list)-1]
chromosome_list.pop()  # removes Z chromosome
chromosome_list.append('scaffolds')  # adds scaffolds (to come before chrz)
chromosome_list.append(chrZ)  # ensures interval order matches reference

# create scaffold bed file
scaffold_file_name = output+'.scaffolds.bed'
scaffold_file = open(scaffold_file_name, 'w')
for scaffold in open(index):
    if scaffold.startswith('S'):
        scaffold = scaffold.split('\t')
        scaffold_file.write(scaffold[0] + '\t0\t' + scaffold[1] + '\n')
scaffold_file.close()

# generate samtools job for each chromosome
job_list = []
vcf_list = []
for position in chromosome_list:
    new_output = output + '.' + position + '.vcf'
    vcf_list.append(new_output)
    job_name = position + '.jobfile.sh'
    job_list.append(job_name)
    if not position == 'scaffolds':
        # SAMtools
        SAM_commandline = ('"samtools mpileup '
                           '-b ' + bams + ' '
                           '-C 50 '
                           '-f ' + ref_genome + ' '
                           '-r ' + position + ' '
                           '-u '
                           '| bcftools_new call '
                           '-O v '
                           '-m '
                           '-o ' + new_output + '"')

        # submit
        qsub_commandline = ('python qsub_gen.py '
                            '-cmd ' + SAM_commandline + ' '
                            '-o ' + new_output.rstrip('.vcf') + ' '
                            '-mo python '
                            '-t 168 '
                            '-OM q '
                            '-jid ' + job_name)
        subprocess.call(qsub_commandline, shell=True)
        print('Writing VCF to: ' + new_output)
    else:
        # SAMtools
        SAM_commandline = ('"samtools mpileup '
                           '-b ' + bams + ' '
                           '-C 50 '
                           '-f ' + ref_genome + ' '
                           '-l ' + scaffold_file_name + ' '
                           '-u '
                           '| bcftools_new call '
                           '-O v '
                           '-m '
                           '-o ' + new_output + '"')
        # submit
        qsub_commandline = ('python qsub_gen.py '
                            '-cmd ' + SAM_commandline + ' '
                            '-o ' + new_output.rstrip('.vcf') + ' '
                            '-mo python '
                            '-t 168 '
                            '-OM q '
                            '-jid ' + job_name)
        subprocess.call(qsub_commandline, shell=True)
        print('Writing VCF to: ' + new_output)

# merge output vcfs into one vcf using GATK CatVariants
CatVar_commandline = ('"python CatVariants.py '
                      '-ref ' + ref_genome + ' ')
for vcf in vcf_list:
    CatVar_commandline += '-VCF ' + vcf + ' '
CatVar_commandline += '-S"'

# submission
CatVar_submission = ('python qsub_gen.py '
                     '-cmd ' + CatVar_commandline + ' '
                     '-o ' + output + '.vcfmerging '
                     '-mo python '
                     '-mo java '
                     '-mo gatk '
                     '-t 48 '
                     '-rmem 20 '
                     '-mem 20 '
                     '-OM q '
                     '-hold ')
for job in job_list:
    CatVar_submission += job + ','
CatVar_submission = CatVar_submission.rstrip(',')
subprocess.call(CatVar_submission, shell=True)
