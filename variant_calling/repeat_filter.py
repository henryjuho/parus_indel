import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('-vcf', '--vcf', help='VCF with variants to filter', required=True)
parser.add_argument('-bed', '--bed_repeats', help='BED file with repeat regions listed', required=True)
parser.add_argument('-ref', '--reference', help='Reference genome', required=True)
args = parser.parse_args()

# files
variants = args.vcf
repeats = args.bed_repeats
reference = args.reference
output_prefix = variants.rstrip('.vcf')

# VariantFiltration
VarFil_cmdline = ('"java -Xmx6g -jar \$GATKHOME/GenomeAnalysisTK.jar '
                  '-T VariantFiltration '
                  '-R ' + reference + ' '
                  '-V ' + variants + ' '
                  '-o ' + output_prefix + '.repeatfilter.vcf '
                  '--mask ' + repeats + ' '
                  '--maskName Repeats"')

# SelectVariants
SelVar_cmdline = ('"java -Xmx6g -jar \$GATKHOME/GenomeAnalysisTK.jar '
                  '-T SelectVariants '
                  '-R ' + reference + ' '
                  '-V ' + output_prefix + '.repeatfilter.vcf '
                  '-o ' + output_prefix + '.repeatfilter.pass.vcf '
                  '--excludeFiltered"')

# remove intermediate file
rm_cmdline = '"rm ' + output_prefix + '.repeatfilter.vcf"'

# submit job
qsub_cmdline = ('python qsub_gen.py '
                '-cmd ' + VarFil_cmdline + ' '
                '-cmd ' + SelVar_cmdline + ' '
                '-cmd ' + rm_cmdline + ' '
                '-o ' + output_prefix + '.repeatfiltering '
                '-mo gatk -mo java '
                '-OM q '
                '-mem 15 -rmem 10')
subprocess.call(qsub_cmdline, shell=True)
