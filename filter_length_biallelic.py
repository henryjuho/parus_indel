import argparse
import subprocess

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-vcf', help='VCF to filter', required=True)
parser.add_argument('-ref', help='Reference genome location', required=True)
args = parser.parse_args()

# variables
input_vcf = args.vcf
output_vcf_prefix = input_vcf.rstrip('.vcf')+'.maxlength50.biallelic'
ref_genome = args.ref

# gatk select variants
select_var_cmdline = ('"java -Xmx6g -jar $GATKHOME/GenomeAnalysisTK.jar '
                      '-T SelectVariants '
                      '-R ' + ref_genome + ' '
                      '-V ' + input_vcf + ' '
                      '-o ' + output_vcf_prefix + '.vcf '
                      '-restrictAllelesTo BIALLELIC '
                      '--maxIndelSize 50"')

# submit job
qsub_cmdline = ('python qsub_gen.py '
                '-cmd ' + select_var_cmdline + ' '
                '-o ' + output_vcf_prefix + ' '
                '-mo java '
                '-mo gatk '
                '-rmem 10 '
                '-mem 10 '
                '-OM q')
subprocess.call(qsub_cmdline, shell=True)
