import argparse
import subprocess

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-vcf', help='VCF file to be annotated', required=True)
parser.add_argument('-ref', help='Reference genome location', required=True)
parser.add_argument('-out', help='Output_directory', required=True)
args = parser.parse_args()

# variables
input_vcf = args.vcf
ref_genome = args.ref
output_prefix = args.out + input_vcf[input_vcf.rfind('/')+1:].rstrip('.vcf')

# construct command lines

gatk_var_anno = ('"java -jar \$GATKHOME/GenomeAnalysisTK.jar '
                 '-R ' + ref_genome + ' '
                 '-T VariantAnnotator '
                 '-o ' + output_prefix + '.hr.tr.vcf '
                 '-V ' + input_vcf + ' '
                 '-A TandemRepeatAnnotator '
                 '-A HomopolymerRun"')

qsub_cmd = ('python qsub_gen.py '
            '-cmd ' + gatk_var_anno + ' '
            '-o ' + output_prefix + ' '
            '-mo java '
            '-mo gatk '
            '-OM q')

subprocess.call(qsub_cmd, shell=True)
