import argparse

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-vcf', help='Input vcf to get INDEL positions from', required=True)
parser.add_argument('-out', help='Output directory and file', required=True)
parser.add_argument('-poly', help='File containing variable values for polynomial equations for each chromosome',
                    required=True)
args = parser.parse_args()

# variables
vcf = open(args.vcf)
out = open(args.out, 'w')


# Toni's polynomial function to predict recomb
def predict_recomb(position, a, b, c):
    return 3*a*position**2+2*b*position+c

# store poly variables for each chromosome
polynomial_values = {'chr' + x.split()[0]: (float(x.split()[1]), float(x.split()[2]), float(x.split()[3]))
                     for x in open(args.poly)}

# calculate recomb rate and write to file
for line in vcf:
    if not line.startswith('#'):
        line = line.rstrip('\n').split('\t')
        chromo = line[0]
        position = float(line[1])
        if chromo in polynomial_values.keys():
            a = polynomial_values[chromo][0]
            b = polynomial_values[chromo][1]
            c = polynomial_values[chromo][2]
            recomb_rate = str(predict_recomb((position/1000000.0), a, b, c))
            out.write(chromo + '\t' + str(int(position)) + '\t' + recomb_rate + '\n')

print('Predicted recombination rates written to ' + args.out)
