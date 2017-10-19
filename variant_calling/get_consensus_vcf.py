import argparse
import subprocess

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-vcf_I', help='First VCF to compare, and to use for naming output', required=True)
parser.add_argument('-vcf_II', help='Second VCF to compare', required=True)
parser.add_argument('-ref', help='Reference genome location', required=True)
parser.add_argument('-out', help='Output_directory', default='input')
args = parser.parse_args()

# variables
vcf_1 = args.vcf_I
vcf_2 = args.vcf_II
ref_genome = args.ref
if args.out == 'input':
    output_prefix = vcf_1.rstrip('.vcf')
else:
    output_prefix = args.out + vcf_1[vcf_1.rfind('/')+1:].rstrip('.vcf')

# GATK select variants to extract INDELs
vcf_2_indels = output_prefix + '.rawindels.vcf'
extract_indels_commandline = ('"java -Xmx6g -jar $GATKHOME/GenomeAnalysisTK.jar '
                              '-T SelectVariants '
                              '-R ' + ref_genome + ' '
                              '-V ' + vcf_2 + ' '
                              '-selectType INDEL  ' 
                              '-trimAlternates '
                              '-env '
                              '-o ' + vcf_2_indels + '"')

# GATK select variants with concordance set
concordance_commandline = ('"java -Xmx6g -jar \$GATKHOME/GenomeAnalysisTK.jar '
                           '-T SelectVariants '
                           '-R ' + ref_genome + ' '
                           '-V ' + vcf_1 + ' '
                           '--concordance ' + vcf_2_indels + ' '
                           '-o ' + output_prefix + '.consensus.vcf"')

# submit to cluster
qsub_commandline = ('python qsub_gen.py '
                    '-cmd ' + extract_indels_commandline + ' '
                    '-cmd ' + concordance_commandline + ' '
                    '-mo java '
                    '-mo gatk '
                    '-mem 10 '
                    '-rmem 10 '
                    '-o ' + output_prefix + '.consensus '
                    '-OM q')
subprocess.call(qsub_commandline, shell=True)
