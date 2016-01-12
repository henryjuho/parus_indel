import argparse
import subprocess

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-vcf', help='VCF file for VQSR to be conducted on', required=True)
parser.add_argument('-ref', help='Reference genome location', required=True)
parser.add_argument('-truth_set', help='Truth set of known variants for VQSR', required=True)
parser.add_argument('-out', help='Output_directory', required=True)
args = parser.parse_args()

# variables
target_vcf = args.vcf
ref_genome = args.ref
truth_set = args.truth_set
output_prefix = args.out + target_vcf[target_vcf.rfind('/')+1:].rstrip('.vcf')
recal_file = output_prefix + '.recal'
tranche_file = output_prefix + '.tranches'
tranche_list = ['100.0', '99.9', '99.5', '99.0', '98.0', '90.0']

# variant recalibration
recalibration_commandline = ('"java -Xmx6g -jar \$GATKHOME/GenomeAnalysisTK.jar '
                             '-T VariantRecalibrator '
                             '-R ' + ref_genome + ' '
                             '-input ' + target_vcf + ' '
                             '-resource:hardfilters,known=true,training=true,truth=true,prior=12.0 ' + truth_set + ' '
                             '-an DP '
                             '-an QD '
                             '-an MQ '
                             '-an MQRankSum '
                             '-an ReadPosRankSum '
                             '-an FS '
                             '-an SOR '
                             '-an InbreedingCoeff '
                             '-mode INDEL '
                             '-recalFile ' + recal_file + ' ')
for tranche in tranche_list:
    recalibration_commandline += '-tranche ' + tranche + ' '
recalibration_commandline += ('-tranchesFile ' + tranche_file + ' '
                              '-rscriptFile ' + output_prefix + '.plots.R"')

# apply recalibration and extract passed
apply_commandlines_list = []
for tranche in tranche_list:
    tranche_out_prefix = output_prefix + '.recalibrated.filtered_t' + tranche
    apply_commandline = ('"java -Xmx3g -jar \$GATKHOME/GenomeAnalysisTK.jar '
                         '-T ApplyRecalibration '
                         '-R ' + ref_genome + ' '
                         '-input ' + target_vcf + ' '
                         '--ts_filter_level ' + tranche + ' '
                         '-tranchesFile ' + tranche_file + ' '
                         '-recalFile ' + recal_file + ' '
                         '-mode INDEL '
                         '-o ' + tranche_out_prefix + '.vcf"')
    extract_commandline = ('"java -Xmx6g -jar \$GATKHOME/GenomeAnalysisTK.jar '
                           '-T SelectVariants '
                           '-R ' + ref_genome + ' '
                           '-V ' + tranche_out_prefix + '.vcf '
                           '-o ' + tranche_out_prefix + '.pass.vcf '
                           '--excludeFiltered"')
    rm_commandline = ('"rm ' + tranche_out_prefix + '.vcf"')
    apply_commandlines_list.append(apply_commandline)
    apply_commandlines_list.append(extract_commandline)
    apply_commandlines_list.append(rm_commandline)

# submit jobs
qsub_commandline = ('python qsub_gen.py '
                    '-cmd ' + recalibration_commandline + ' ')
for commandline in apply_commandlines_list:
    qsub_commandline += '-cmd ' + commandline + ' '
qsub_commandline += ('-o ' + output_prefix + '.recal '
                     '-mo java '
                     '-mo gatk '
                     '-OM q '
                     '-mem 10 '
                     '-rmem 10')
subprocess.call(qsub_commandline, shell=True)
