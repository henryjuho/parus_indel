library(ape)

args = commandArgs(TRUE)

fasta = args[1]
region = args[2]

data = read.FASTA(fasta)

divergance = dist.dna(data, model = "K80")

cat(paste(attr(divergance, "Labels"), sep='_'))
cat('\n')

gt_fc = divergance[1]
zf_fc = divergance[2]
gt_zf = divergance[3]

gt = (gt_fc - zf_fc + gt_zf) / 2

cat(region, gt, sep='\t')
cat('\n')
