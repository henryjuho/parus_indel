# Annotating the data 

## RefSeq contig key file generation

```
$ zgrep -A 2 ^##sequence GCF_001522545.1_Parus_major1.0.3_genomic.gff.gz | grep ^NC | cut -f 1,9 | cut -d '=' -f 1,5 | cut -d ';' -f 1 | sed -e 's/ID=/chr/g' > contigs_key.txt
$ cp contigs_key.txt ~/parus_indel/annotation/
```

## Region coordinates

Bed files with coordinates for different genomic contexts were created as follows:

```
$ zgrep -v ^# GCF_001522545.1_Parus_major1.0.3_genomic.rename.gff.gz | cut -f 1-5 | grep region | grep -v chrU | grep -v MT | gff2bed.py > gt_chromosomes.bed
$ zgrep -v ^# GCF_001522545.1_Parus_major1.0.3_genomic.rename.gff.gz | cut -f 1-5 | grep gene | grep -v chrU | grep -v MT | gff2bed.py | sort -k1,1 -k2,2n | bedtools merge > gt_genes.bed
$ zgrep -v ^# GCF_001522545.1_Parus_major1.0.3_genomic.rename.gff.gz | cut -f 1-5 | grep exon | grep -v chrU | grep -v MT | gff2bed.py | sort -k1,1 -k2,2n | bedtools merge > gt_exons.bed
$ zgrep -v ^# GCF_001522545.1_Parus_major1.0.3_genomic.rename.gff.gz | cut -f 1-5 | grep CDS | grep -v chrU | grep -v MT | gff2bed.py | sort -k1,1 -k2,2n | bedtools merge > gt_cds.bed
$ bedtools subtract -a gt_genes.bed -b gt_exons.bed > gt_introns.bed
$ bedtools subtract -a gt_chromosomes.bed -b gt_genes.bed > gt_intergenic.bed
$ cat gt_intergenic.bed gt_introns.bed | sort -k1,1 -k2,2n | bedtools merge > gt_noncoding.bed
$ bedtools subtract -a gt_noncoding.bed.gz -b UCNEs_gtit1.0.4.sorted.bed.gz | bgzip -c > gt_noncoding_withoutUCNE.bed.gz
$ ~/parus_indel/annotation/degen_to_bed_gt.py -cds_fa cds_fasta/GCF_001522545.1_Parus_major1.0.3_cds_from_genomic.fna.gz -degen 0 | sort -k1,1 -k2,2n | bedtools merge > gt_0fold.bed
$ ~/parus_indel/annotation/degen_to_bed_gt.py -cds_fa cds_fasta/GCF_001522545.1_Parus_major1.0.3_cds_from_genomic.fna.gz -degen 2 | sort -k1,1 -k2,2n | bedtools merge > gt_2fold.bed
$ ~/parus_indel/annotation/degen_to_bed_gt.py -cds_fa cds_fasta/GCF_001522545.1_Parus_major1.0.3_cds_from_genomic.fna.gz -degen 3 | sort -k1,1 -k2,2n | bedtools merge > gt_3fold.bed
$ ~/parus_indel/annotation/degen_to_bed_gt.py -cds_fa cds_fasta/GCF_001522545.1_Parus_major1.0.3_cds_from_genomic.fna.gz -degen 4 | sort -k1,1 -k2,2n | bedtools merge > gt_4fold.bed
$ ls *fold.bed | while read i; do bgzip $i; done
$ ls *fold.bed.gz | while read i; do tabix -pbed $i; done
```

## Genomic regions

Firstly INDELs were annotated in the vcf file as belonging to either 'CDS_non_frameshift', 'CDS_frameshift', 'intron' or 'intergenic' as follows:

```
$ ./annotate_all_vcf_chr.py -gff /data/bop15hjb/annotating_gt/GCF_001522545.1_Parus_major1.0.3_genomic.rename.gff.gz-vcf /data/bop15hjb/annotating_gt/bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.repeatfilter.pass.polarised.vcf -evolgen
```

A summary of this annotation is shown below:

|Category             | Number INDELs|
|:--------------------|:------------:|
|All                  | 1221124      |
|CDS                  | 2027         |
|Frameshift           | 973          |
|Non shift            | 1054         |
|Intron               | 646507       |
|Intergenic           | 537246       |
|Unannotated          | 35344        |

## Recombination regions

Secondly the recombination category of each INDEL was annotated using linkage map data to estimate recombination rates. First, 3rd order polynomials were fitted to plots of physical position versus map length. Second, the derivative of each chromosome's polynomial was used to estimate recombination rate for each INDEL start position. This predicition was implemented in the following python script:

```
$ ./annotate_recomb.py -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Analysis_ready_data/bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.repeatfilter.pass.polarised.annotated.vcf -ref /fastdata/bop15hjb/GT_ref/Parus_major_1.04.rename.fa -poly /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Recomb_data/3rd_polynom.txt -out /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Recomb_data/recomb_anno/ -evolgen
```

The file specified by ```-poly``` is a list of variables for each chromosome's polynomial.

# Annotating the SNP data

Firstly SNPs were annotated in the vcf file as belonging to either 'CDS_non_frameshift', 'intron' or 'intergenic' as follows:

```
$ ./annotate_all_vcf_chr.py -gff /data/bop15hjb/annotating_gt_snps/GCF_001522545.1_Parus_major1.0.3_genomic.rename.gff.gz -vcf /data/bop15hjb/annotating_gt_snps/gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.coveragefiltered.pass.repeatfilter.pass.polarised.vcf -evolgen
```

A breakdown of this annotation is shown below:

|Category             | Number SNPs  |
|:--------------------|:------------:|
|All                  | 10643865     |
|CDS                  | 155463       |
|Intron               | 5489556      |
|Intergenic           | 4718020      |
|Unannotated          | 280826       |

Secondly the degeneracy of SNPs in coding seqences was annotated as follows:

```
$ ./annotate_degeneracy.py -gff /fastdata/bop15hjb/GT_ref/GCF_001522545.1_Parus_major1.0.3_genomic.rename.gff.gz -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/SNP_data/gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.coveragefiltered.pass.repeatfilter.pass.polarised.annotated.vcf -ref /fastdata/bop15hjb/GT_ref/Parus_major_1.04.rename.fa -db_dir /data/bop15hjb/databases/greattit/ -out /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/SNP_data/ degeneracy_annotation/ -evolgen
```

Annotation summary below:

|Category          | Number SNPs     |
|:-----------------|:---------------:|
|All               | 10643865        |
|CDS               | 155463          |
|0fold             | 46792           |
|2fold             | 46786           |
|3fold             | 4240            |
|4fold             | 57169           |

## LINEs

Variants in ancestral LINEs were extracted from the post VQSR vcf file, prior to the repeat filtering step as follows:

```
$ bedtools intersect -a bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.vcf.gz -b ../repeat_coordinates/LINE_intersect/Zebrafinch.Flycatcher.Greattit.ancLINEs.sorted.bed.gz -header | bgzip -c > bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.vcf.gz
$ cd /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/polymorphicLINEs
$ zgrep -v ^Scaffold /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/SNP_data/gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.coveragefiltered.pass.vcf.gz | bedtools intersect -a stdin -b /fastdata/bop15hjb/GT_ref/Zebrafinch.Flycatcher.Greattit.ancLINEs.sorted.bed.gz -header | bgzip -c > gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.vcf.gz
$ gunzip gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.vcf.gz 
```

These INDELs were then polarised with the same pipeline as other INDELs as follows:

```
$ ./VARfromMAF.py -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/divergence_data/polymorphic_non_repeat_filtered_vcfs/bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.vcf -maf /fastdata/bop15hjb/bird_alignments/UCSC_pipeline/multiple_zhang_param/Zebrafinch.Flycatcher.Greattit.maf -target_spp Greattit -out /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/divergence_data/polymorphic_non_repeat_filtered_vcfs/ -no_jobs 100
$ polarise_vcf_indels.py -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/divergence_data/polymorphic_non_repeat_filtered_vcfs/bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.vcf -align_data /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/divergence_data/polymorphic_non_repeat_filtered_vcfs/all_variants.alignment_states.txt -target_spp Greattit

$ ~/parus_indel/alignment_and_polarisation/VARfromMAF.py -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/polymorphicLINEs/gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.vcf -maf /fastdata/bop15hjb/bird_alignments/UCSC_pipeline/multiple_zhang_param/Zebrafinch.Flycatcher.Greattit.maf -target_spp Greattit -out /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/polymorphicLINEs/snp_pol/ -no_jobs 100
$ ~/parus_indel/alignment_and_polarisation/polarise_vcf_snps.py -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/polymorphicLINEs/gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.vcf -align_data /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/polymorphicLINEs/snp_pol/all_variants.alignment_states.txt -target_spp Greattit
```

|                 | INDELs | SNPs   |
|:----------------|:------:|:------:|
|Total no variants| 19849  | 241191 |
|INDELs polarised | 15340  | 174838 |
|Hotspots         | 3897   |        |
|Low spp coverage | 10     | 312    |
|Ambiguous        | 486    | 58908  |
|Not in alignment | 116    | 1004   |
|Total unpolarised| 4509   | 66353  |

These variants were then annotated by genomic region as with the main INDEL dataset, and any coding region INDELs found were removed (5 INDELs).

```
$ annotate_all_vcf_chr.py -gff /data/bop15hjb/databases/greattit/GCF_001522545.1_Parus_major1.0.3_genomic.rename.gff.gz -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/divergence_data/polymorphic_non_repeat_filtered_vcfs/bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.polarised.vcf -evolgen
$ zgrep -v CDS bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.polarised.annotated.vcf.gz | bgzip -c > bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.polarised.annotated.noCDS.vcf.gz

$ cd /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/polymorphicLINEs/
$ ~/parus_indel/annotation/annotate_all_vcf_chr.py -gff /data/bop15hjb/databases/greattit/GCF_001522545.1_Parus_major1.0.3_genomic.rename.gff.gz -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/polymorphicLINEs/gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.polarised.vcf
$ zgrep -v CDS gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.polarised.annotated.vcf | bgzip -c > gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.polarised.annotated.noCDS.vcf.gz
```

The annotation was upated for LINEs from 'intergenic' and 'intron' to 'intergenic_ar' and 'intron_ar'.

```
$ zcat bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.polarised.annotated.noCDS.vcf.gz | python ~/parus_indel/annotation/update_line_anno.py | bgzip -c > bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.polarised.annotated.noCDS.retagged.vcf.gz 
$ tabix -pvcf bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.polarised.annotated.noCDS.retagged.vcf.gz

$ zcat gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.polarised.annotated.noCDS.vcf.gz | python ~/parus_indel/annotation/update_line_anno.py | bgzip -c > gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.polarised.annotated.noCDS.retagged.vcf.gz 
$ tabix -pvcf gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.polarised.annotated.noCDS.retagged.vcf.gz 
```

The INDELs and SNPs identified within LINEs were then added to the main vcfs:

```
$ cd /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Analysis_ready_data/final
$ cp ../bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.repeatfilter.pass.polarised.annotated.recomb.vcf.gz* ./
$ cp ../../polymorphicLINEs/bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.polarised.annotated.noCDS.retagged.vcf.gz* ./
$ zgrep ^# bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.repeatfilter.pass.polarised.annotated.recomb.vcf.gz > bgi_10birds.filtered_indels.pol.anno.recomb.line.vcf
$ zcat *vcf.gz | grep -v ^# | sort -k1,1 -k2,2n >> bgi_10birds.filtered_indels.pol.anno.recomb.line.vcf
$ bgzip bgi_10birds.filtered_indels.pol.anno.recomb.line.vcf
$ tabix -pvcf bgi_10birds.filtered_indels.pol.anno.recomb.line.vcf.gz 
$ rm bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.*

$ cd /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Analysis_ready_data/final
$ cp ../../SNP_data/gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.coveragefiltered.pass.repeatfilter.pass.polarised.annotated.degen.vcf.gz bgi_10birds.filtered_snps.pol.anno.degen.vcf.gz
$ zgrep ^# bgi_10birds.filtered_snps.pol.anno.degen.vcf.gz > bgi_10birds.filtered_snps.pol.anno.degen.line.vcf
$ bgzip ../../polymorphicLINEs/gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.polarised.annotated.vcf
$ zcat bgi_10birds.filtered_snps.pol.anno.degen.vcf.gz ../../polymorphicLINEs/gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.polarised.annotated.vcf.gz | grep -v ^# | sort -k1,1 -k2,2n >> bgi_10birds.filtered_snps.pol.anno.degen.line.vcf 
$ bgzip bgi_10birds.filtered_snps.pol.anno.degen.line.vcf 
$ tabix -pvcf bgi_10birds.filtered_snps.pol.anno.degen.line.vcf.gz 
```