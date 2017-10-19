import argparse
import subprocess
import re

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-vcf', help='List of vcfs to look calculate stats for', required=True)
parser.add_argument('-out', help='Output directory and file name', required=True)
args = parser.parse_args()

# variables
input_vcfs = [vcf for vcf in open(args.vcf)]
output = open(args.out, 'w')
output.write('Tranche\tAll sites\tPositive training sites\tNegative training sites\tNovel sites\n')  # file header


# function to get tranche level from file name
def get_tranche(file_name):
    tranche = str(re.search(r't([\d]{2,3}\.\d)', file_name).group(1))
    return tranche

# use grep to count train sites and novel variants
for vcf in input_vcfs:
    vcf = vcf.rstrip('\n')
    negative_train_sites = int(subprocess.check_output('grep -v ^# ' + vcf + ' | grep -c NEGATIVE_TRAIN_SITE',
                                                       shell=True))
    positive_train_sites = int(subprocess.check_output('grep -v ^# ' + vcf + ' | grep -c POSITIVE_TRAIN_SITE',
                                                       shell=True))
    allsites = int(subprocess.check_output(('grep -v ^# -c ' + vcf), shell=True))
    novel_sites = allsites - positive_train_sites - negative_train_sites
    tranche_level = get_tranche(vcf)
    output.write(tranche_level + '\t' +
                 str(allsites) + '\t' +
                 str(positive_train_sites) + '\t' +
                 str(negative_train_sites) + '\t' +
                 str(novel_sites) + '\n')

print('Done')
