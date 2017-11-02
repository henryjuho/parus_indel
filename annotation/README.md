# Annotating the data 
## Genomic regions

Firstly INDELs were annotated in the vcf file as belonging to either 'CDS_non_frameshift', 'CDS_frameshift', 'intron' or 'intergenic' as follows:

```
./annotate_all_vcf_chr.py -gff /data/bop15hjb/annotating_gt/GCF_001522545.1_Parus_major1.0.3_genomic.rename.gff.gz-vcf /data/bop15hjb/annotating_gt/bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.repeatfilter.pass.polarised.vcf -evolgen
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
./annotate_recomb.py -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Analysis_ready_data/bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.repeatfilter.pass.polarised.annotated.vcf -ref /fastdata/bop15hjb/GT_ref/Parus_major_1.04.rename.fa -poly /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Recomb_data/3rd_polynom.txt -out /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Recomb_data/recomb_anno/ -evolgen
```

The file specified by ```-poly``` is a list of variables for each chromosome's polynomial.

## LINEs

INDELs in ancestral LINEs were extracted from the post VQSR vcf file, prior to the repeat filtering step as follows:

```
bedtools intersect -a bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.vcf.gz -b ../repeat_coordinates/LINE_intersect/Zebrafinch.Flycatcher.Greattit.ancLINEs.sorted.bed.gz -header | bgzip -c > bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.vcf.gz
```

These INDELs were then polarised with the same pipeline as other INDELs as follows:

```
./VARfromMAF.py -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/divergence_data/polymorphic_non_repeat_filtered_vcfs/bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.vcf -maf /fastdata/bop15hjb/bird_alignments/UCSC_pipeline/multiple_zhang_param/Zebrafinch.Flycatcher.Greattit.maf -target_spp Greattit -out /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/divergence_data/polymorphic_non_repeat_filtered_vcfs/ -no_jobs 100
polarise_vcf_indels.py -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/divergence_data/polymorphic_non_repeat_filtered_vcfs/bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.vcf -align_data /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/divergence_data/polymorphic_non_repeat_filtered_vcfs/all_variants.alignment_states.txt -target_spp Greattit
```

|                 |       |
|:----------------|:-----:|
|Total no INDELs  | 19849 |
|INDELs polarised | 15340 |
|Hotspots         | 3897  |
|Low spp coverage | 10    |
|Ambiguous        | 486   |
|Not in alignment | 116   |
|Total unpolarised| 4509  |

These variants were then annotated by genomic region as with the main INDEL dataset, and any coding region INDELs found were removed (5 INDELs).

```
annotate_all_vcf_chr.py -gff /data/bop15hjb/databases/greattit/GCF_001522545.1_Parus_major1.0.3_genomic.rename.gff.gz -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/divergence_data/polymorphic_non_repeat_filtered_vcfs/bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.polarised.vcf -evolgen
zgrep -v CDS bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.polarised.annotated.vcf.gz | bgzip -c > bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.polarised.annotated.noCDS.vcf.gz
```

The annotation was upated for LINEs from 'intergenic' and 'intron' to 'intergenic_ar' and 'intron_ar'.

```
$ cd /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/polymorphicLINEs/
$ zcat bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.polarised.annotated.noCDS.vcf.gz | python ~/parus_indel/annotation/update_line_anno.py | bgzip -c > bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.polarised.annotated.noCDS.retagged.vcf.gz 
$ tabix -pvcf bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.polarised.annotated.noCDS.retagged.vcf.gz
```

The INDELs identified within LINEs were then added to the main INDEL vcf:

```
$ cd /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Analysis_ready_data/final
$ cp ../bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.repeatfilter.pass.polarised.annotated.recomb.vcf.gz* ./
$ cp ../../polymorphicLINEs/bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.LINEintersect.polarised.annotated.noCDS.retagged.vcf.gz* ./
$ zgrep ^# bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.repeatfilter.pass.polarised.annotated.recomb.vcf.gz > bgi_10birds.filtered_indels.pol.anno.recomb.line.vcf
$ zcat *vcf.gz | grep -v ^# | sort -k1,1 -k2,2n >> bgi_10birds.filtered_indels.pol.anno.recomb.line.vcf
$ bgzip bgi_10birds.filtered_indels.pol.anno.recomb.line.vcf
$ tabix -pvcf bgi_10birds.filtered_indels.pol.anno.recomb.line.vcf.gz 
$ rm bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.*
```

The final SNP vcf was moved into the same directory and the name shortened to save my keyboard:

```
$ cp ../../SNP_data/gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.coveragefiltered.pass.repeatfilter.pass.polarised.annotated.degen.vcf.gz bgi_10birds.filtered_snps.pol.anno.degen.vcf.gz
$ tabix -pvcf bgi_10birds.filtered_snps.pol.anno.degen.vcf.gz
```