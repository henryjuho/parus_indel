import argparse

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-vcf', help='/input/directory/vcf', required=True)
parser.add_argument('-recomb', help='List file of INDEL positions and recomb rates', required=True)
parser.add_argument('-out', help='/output/directory/', required=False, default=False)
parser.add_argument('-nbin', help='Number of recombination bins required', default=5.0, type=float, required=False)
args = parser.parse_args()

# variables
input_vcf = args.vcf
recomb_info = open(args.recomb)
if args.out is False:
    output_prefix = input_vcf.rstrip('.vcf')
else:
    output_prefix = args.out + input_vcf[input_vcf.rfind('/')+1:].rstrip('.vcf')
bin_number = args.nbin

# read in recomb data
indel_recombs = [(line.rstrip('\n').split('\t')[0],
                  int(line.rstrip('\n').split('\t')[1]),
                  float(line.rstrip('\n').split('\t')[2]))
                 for line in recomb_info]

indel_recombs = sorted(indel_recombs, key=lambda x: x[2])  # list of INDEL positions sorted by recomb rate

# calculate bin size
no_indels = len(indel_recombs)
bin_size = int(no_indels/bin_number)
print('Variant bins will be of size n = ' + str(bin_size) +
      '\nWARNING! If number of variants not exactly divisible by bin size, remainder will be added to final bin.')

# generate bins
bins = []
bin_tracker = 0
for i in range(0, no_indels - bin_size - 1, bin_size):
    bin_tracker += 1
    if bin_tracker == int(bin_number):
        end_coord = no_indels
    else:
        end_coord = i + bin_size
    bin_ID = 'bin_' + str(i+1) + '-' + str(end_coord)
    bins.append([bin_ID, indel_recombs[i:end_coord]])

# write new .vcf file for each bin
for recomb_bin in bins:
    output_path = output_prefix + '.' + recomb_bin[0] + '.vcf'
    output_vcf = open(output_path, 'w')
    bin_set = set([x[0]+'_'+str(x[1]) for x in recomb_bin[1]])
    with open(input_vcf) as vcf_file:
        for line in vcf_file:
            if line.startswith('#'):
                output_vcf.write(line)
            else:
                variant_data = line.split('\t')
                variant_ID = variant_data[0]+'_'+str(variant_data[1])
                if variant_ID in bin_set:
                    output_vcf.write(line)
    print('Varaints in ' + recomb_bin[0] + ' written to ' + output_path)
