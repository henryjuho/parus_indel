import argparse
import re

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-vcf', help='VCF file to count numbers of TRs and HRs in', required=True)
args = parser.parse_args()

# variables
input_vcf = open(args.vcf)
hr_count = 0
hr_dict = {}
tr_count = 0
tr_dict = {}
tr_hr_overlap = 0
total_indels = 0

# loop through vcf
for line in input_vcf:
    if not line.startswith('#'):
        # count number of indels
        total_indels += 1

        # summarise homopolymer runs
        if re.search(r'HRun=[123456789][\d]{0,3};', line):
            hr_count += 1
            # hr = re.search(r'HRun=([123456789][\d]{0,3});', line).group(1)
            # if hr in hr_dict.keys():
            #     hr_dict[hr] += 1
            # else:
            #     hr_dict[hr] = 1

        # summarise tandem repeats
        if re.search(r'STR;', line):
            tr_count += 1
            # tr_motif = re.search(r'RU=([ATCG]+);', line).group(1)
            # if tr_motif in tr_dict.keys():
            #     tr_dict[tr_motif] += 1
            # else:
            #     tr_dict[tr_motif] = 1

        # count overlap of hr_count and tr_count
        if re.search(r'HRun=[123456789][\d]{0,3};', line) and re.search(r'STR;', line):
            tr_hr_overlap += 1


print('Summary of INDELs in tandem repeats and homopolymer runs\n'
      '--------------------------------------------------------\n'
      'Class\tNumber of INDELs\n'
      '--------------------------------------------------------\n'
      'Hompolymer runs\t' + str(hr_count) + '\n'
      'Tandem repeats\t' + str(tr_count) + '\n'
      'Overlap\t' + str(tr_hr_overlap) + '\n'
      'Total repetative\t' + str((hr_count+tr_count)-tr_hr_overlap) + '\n'
      'Total no. INDELs\t' + str(total_indels) + '\n'
      '% INDEL in repeat\t' + str(int((float(((hr_count+tr_count)-tr_hr_overlap))/float(total_indels))*100.0)) + '\n'
      '--------------------------------------------------------')
