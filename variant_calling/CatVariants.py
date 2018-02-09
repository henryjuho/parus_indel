import argparse
import subprocess
import re
import os

parser = argparse.ArgumentParser()
parser.add_argument('-VCF', '--VCF', help='VCF file to add to list to concatenate, or VCF list file',
                    required=True, action='append')
parser.add_argument('-ref', '--reference', help='Reference genome location', required=True)
parser.add_argument('-ncl', '--noclean', help='If specified, will not remove intermediate files generated by GATK',
                    action='store_false', default=True)
parser.add_argument('-S', '--Sorted', help='If specified runs CatVariants with --assumeSorted',
                    action='store_true', default=False)
args = parser.parse_args()

# prepares list of VCF files from either list file or individually specified files
vcf_file_data = args.VCF
if len(vcf_file_data) == 1 and not vcf_file_data[0].endswith('.vcf'):
    input_vcfs = [vcf.rstrip('\n') for vcf in open(vcf_file_data[0])]
else:
    input_vcfs = vcf_file_data

# set remaining variables
ref_genome = args.reference
chr_identifier = re.search(r'\.(chr\d{0,2}[ABLGE]{0,3}\d{0,2})\.vcf', input_vcfs[0]).group(1)  # needs testing
output = re.sub(chr_identifier, 'allsites', input_vcfs[0])
tocleanornottoclean = args.noclean
if args.Sorted is True:
    sort_mode = ' -assumeSorted'
else:
    sort_mode = ''

# CatVariants commandline
cat_var = ('java -Xmx15g -cp $GATKHOME/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants '
           '-R ' + ref_genome + ' ')
for vcf in input_vcfs:
    cat_var += '-V ' + vcf + ' '
cat_var += '-out ' + output + sort_mode

# combine chromosome vcfs into wholegenome vcf
subprocess.call(cat_var, shell=True)

# clean up intermediate vcfs and files
if tocleanornottoclean is True:
    for vcf in input_vcfs:
        os.remove(vcf)  # add some exceptions
        os.remove(vcf+'.idx')
        os.remove(vcf.replace('.vcf', '.out'))
