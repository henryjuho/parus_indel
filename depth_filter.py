import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('-vcf', '--vcf', help='Location of vcf file to filter', required=True)
parser.add_argument('-DF', '--DepthFilter',
                    help='Defines abnormal depth eg) 2 means abnormal depth is twice and half the mean depth',
                    default=2.0, type=float)
parser.add_argument('-mean_depth', '--mean_depth', help='Mean coverage depth of samples', default=44.0)
parser.add_argument('-N', '--no_individuals', help='Number of individuals in allsites VCF', type=float, default=10.0)
args = parser.parse_args()

# variables
vcf = args.vcf
destination = vcf.rstrip('vcf')+'coveragefiltered.pass.vcf'
output_vcf = open(destination, 'w')
filter_factor = args.DepthFilter
all_data_mean_depth = float(args.mean_depth)
no_indiv = args.no_individuals

# calculate depth cutoffs
lower_depth_limit = all_data_mean_depth / filter_factor
upper_depth_limit = all_data_mean_depth * filter_factor

# filter vcf
number_failed = 0
number_passed = 0
for line in open(vcf):
    if line.startswith('#'):
        output_vcf.write(line)
    else:
        cumulative_depth = float(re.search(r';DP=([\d]*)', line).group(1))
        locus_mean_depth = cumulative_depth / no_indiv
        if locus_mean_depth > upper_depth_limit or locus_mean_depth < lower_depth_limit:
            number_failed += 1
        else:
            number_passed += 1
            output_vcf.write(line)

# descriptive stats printed out at end, ie depth range, sites passed etc
print('Depth range: ' + str(lower_depth_limit) + ' - ' + str(upper_depth_limit))
print(str(number_passed) + ' variants passed, ' + str(number_failed) + ' variants failed')
print('Output written to ' + destination)
