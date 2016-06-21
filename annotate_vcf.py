#!/usr/bin/env python

import argparse
import gffutils
import os
import gzip
import math

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-gff', help='GFF file to read annotations from, note squence names must match VCF', required=True)
parser.add_argument('-vcf', help='VCF file to annotate variants in', required=True)
parser.add_argument('-chr', help='Chromosome to annotate, must be consistent with GFF and VCF', required=True)
args = parser.parse_args()

# variables
gff = args.gff
chromo = args.chr
gff_db = gff + '.' + chromo + '.db'

gff_string = ''
for line in gzip.open(gff):
    gff_chromo = line.split('\t')[0]
    if chromo == gff_chromo:
        gff_string += line

vcf_file = args.vcf
annotated_vcf = open(vcf_file.replace('.vcf', '.annotated.' + chromo + '.vcf'), 'w')


# functions
def is_in_region(interval, tool):
    in_region = False
    for position in tool:
        if position[0] == interval[0]:
            if int(position[1]) <= interval[1] < interval[2] <= int(position[2]):
                in_region = True
                return in_region
    return in_region

# make gff databas for chromo
if not os.path.isfile(gff_db):
    annotation_db = gffutils.create_db(gff_string, gff_db, merge_strategy='create_unique', from_string=True)
else:
    annotation_db = gffutils.FeatureDB(gff_db, keep_order=True)

# get introns
introns = annotation_db.create_introns()
intron_intervals = [(x.seqid, x.start, x.end) for x in introns]

# get intergenic
genes = annotation_db.merge(annotation_db.features_of_type('gene', order_by=('seqid', 'start', 'end')),
                            ignore_strand=True)
intergenic = annotation_db.interfeatures(genes, new_featuretype='intergenic')
intergenic_intervals = [(x.seqid, x.start, x.end) for x in intergenic]

# get cds
cds = annotation_db.features_of_type('CDS', order_by=('seqid', 'start', 'end'))
cds_intervals = [(x.seqid, x.start, x.end) for x in cds]

# annotation counters
all_variants = 0
cds_frame = 0
cds_nonframe = 0
inter = 0
intron_count = 0

# loop through vcf and identify category for each variant
previous_line = ''
for line in open(vcf_file):
    if line.startswith('#'):
        if line.startswith('##contig') and previous_line.startswith('##INFO'):
            new_info = '##INFO=<ID=ANNO,Number=1,Type=String,Description="Annotation of genomic region">\n'
            annotated_vcf.write(new_info)
        previous_line = line
        annotated_vcf.write(line)
    elif line.split('\t')[0] == chromo:
        all_variants += 1
        split_line = line.split('\t')
        chromo, start, end = split_line[0], int(split_line[1])-1, int(split_line[1]) + len(split_line[3])
        var_interval = [chromo, start, end]
        in_cds = ['CDS', is_in_region(var_interval, cds_intervals)]
        in_intergenic = ['intergenic', is_in_region(var_interval, intergenic_intervals)]
        in_intron = ['intron', is_in_region(var_interval, intron_intervals)]
        region_true = [x[0] for x in [in_cds, in_intergenic, in_intron] if x[1] is True]
        number_true = len(region_true)

        if number_true == 1:
            region = region_true[0]
            if region == 'CDS':
                variant_length = math.sqrt((len(split_line[3]) - len(split_line[4]))**2)
                if variant_length % 3.0 == 0.0:
                    region = 'CDS_non_frameshift'
                    cds_nonframe += 1
                else:
                    region = 'CDS_frameshift'
                    cds_frame += 1
            if region == 'intron':
                intron_count += 1
            if region == 'intergenic':
                inter += 1

            # write newly annotated line
            anno_line = ('\t'.join(split_line[0:7]) + '\t' + split_line[7] +
                         ';ANNO=' + region + '\t' + '\t'.join(split_line[8:]))
            annotated_vcf.write(anno_line)
        else:
            annotated_vcf.write(line)

    else:
        continue

print('\n'
      '|Category          | Number INDELs   |\n'
      '|:-----------------|:---------------:|\n'
      '|All               | ' + str(all_variants) + ' |\n'
      '|CDS_frameshift    | ' + str(cds_frame) + ' |\n'
      '|CDS_non_frameshift| ' + str(cds_nonframe) + ' |\n'
      '|Intron            | ' + str(intron_count) + ' |\n'
      '|Intergenic        | ' + str(inter) + ' |\n'
      '|Not annotated     | ' + str(all_variants - (cds_frame + cds_nonframe + intron_count + inter)) + ' |')

annotated_vcf.close()
