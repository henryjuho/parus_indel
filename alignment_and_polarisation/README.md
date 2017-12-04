# Multispecies alignment and INDEL polarisation

## Alignment

Multispecies alignment was performed between zebra finch, flycatcher and great tit as described here: <https://github.com/henryjuho/bird_alignments/tree/master/Zebrafinch_Flycatcher_Greattit>.

## Getting ancestral repeats using the WGA

Ancestral LINEs coordinates were obtained for the great tit genome using the mutlispecies alignment. Firstly the repeat masker output files were obtained for great tit and zebra finch, and the masking coordinates for the flycatcher were downloaded from <ftp://ftp.ncbi.nlm.nih.gov/genomes/Ficedula_albicollis/masking_coordinates.gz>. LINE coordinates were then extracted from these files as follows.

For great tit and zebra finch:

```
zgrep LINE GCF_001522545.1_Parus_major1.0.3_rm.out.gz | cut -d '.' -f 4- | cut -d ' ' -f 3- | cut -d '(' -f 1 | bgzip -c > GCF_001522545.1_Parus_major1.0.3_rm.out.LINEs.bed.gz
zgrep LINE Taeniopygia_guttata.RepeatMasker.out.gz | cut -d '.' -f 4- | cut -d ' ' -f 3- | cut -d '(' -f 1 | while read i; do echo chr$i; done | bgzip -c > Taeniopygia_guttata.RepeatMasker.out.LINEs.rename.bed.gz
```

For flycatcher LINEs were extracted using a list of LINE repeat types found in zebrafinch and great tit:

```
zcat FC_masking_coordinates.gz | maskingcoord2LINEcoord.py -LINE_type ZF_GT_LINE_types.txt | bgzip -c > FicAlb1.5.masking_coordinates.LINEs.bed.gz
```

Flycatcher and great tit chromosomes were then renamed to format ```chrXX``` from NCBI ref-seq IDs as follows:

```
grep -v ^Loc FC_replicon_info.csv | cut -d ',' -f 3-4 | grep NC | sed 's:,:       :g' | while read i; do echo chr$i; done > FicAlb1.5.chromosome_IDs.txt 
zcat FicAlb1.5.masking_coordinates.LINEs.bed.gz | rename_bedchromsomes.py -chr_IDs FicAlb1.5.chromosome_IDs.txt | bgzip -c > FicAlb1.5.masking_coordinates.LINEs.rename.bed.gz 
```

```
zcat /fastdata/bop15hjb/GT_ref/GCF_001522545.1_Parus_major1.0.3_genomic.gff.gz | GFF2chromoIDs.py | grep NC > GCF_001522545.1_Parus_major1.0.3_chromosome_IDs.txt
zcat GCF_001522545.1_Parus_major1.0.3_rm.out.LINEs.bed.gz | rename_bedchromsomes.py -chr_IDs GCF_001522545.1_Parus_major1.0.3_chromosome_IDs.txt | bgzip -c > GCF_001522545.1_Parus_major1.0.3_rm.out.LINEs.rename.bed.gz 
```

The overlap of all three files was then taken using the whole genome alignment as follows (note this script wraps some utilities from <https://github.com/padraicc/WGAbed>):

```
get_conserved_LINEs.py -wga_bed /fastdata/bop15hjb/bird_alignments/UCSC_pipeline/multiple_zhang_param/Zebrafinch.Flycatcher.Greattit.wga.bed.gz -r_bed Tit_data/BGI_BWA_GATK/divergence_data/repeat_coordinates/GCF_001522545.1_Parus_major1.0.3_rm.out.LINEs.rename.bed.gz -q_bed Zebrafinch,Tit_data/BGI_BWA_GATK/divergence_data/repeat_coordinates/Taeniopygia_guttata.RepeatMasker.out.LINEs.rename.tab.sorted.bed.gz -q_bed Flycatcher,Tit_data/BGI_BWA_GATK/divergence_data/repeat_coordinates/FicAlb1.5.masking_coordinates.LINEs.rename.bed.gz -chr_list /fastdata/bop15hjb/GT_ref/chromosome_list.bed -out Tit_data/BGI_BWA_GATK/divergence_data/repeat_coordinates/LINE_intersect/Zebrafinch.Flycatcher.Greattit -evolgen
zcat Zebrafinch.Flycatcher.Greattit.LINEs.wga.bed.gz | sort -k1,1 -k2,2n | bgzip -c > Zebrafinch.Flycatcher.Greattit.LINEs.sorted.wga.bed.gz
zcat Zebrafinch.Flycatcher.Greattit.LINEs.sorted.wga.bed.gz | cut -f 1-3 | bedtools merge -i stdin | bgzip -c > Zebrafinch.Flycatcher.Greattit.ancLINEs.sorted.bed.gz
```

## INDEL polarisation

Firstly the sequence alignments across species [in addition to 1bp up and down stream of INDEL] for each INDEL in the dataset were pulled out of the multiple alignment file with the following script:

```
~/VARfromMAF.py -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Analysis_ready_data/bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.repeatfilter.pass.vcf -maf /fastdata/bop15hjb/bird_alignments/UCSC_pipeline/multiple_zhang_param/Zebrafinch.Flycatcher.Greattit.maf  -target_spp Greattit -out /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Polarisation/ -no_jobs 100
```

Secondly the output of the above script was used to annotate the ancestral state for each variant in the INDEL vcf as follows:

```
~/polarise_vcf_indels.py -vcf ../Analysis_ready_data/bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.repeatfilter.pass.vcf -align_data all_indels.alignment_states.txt -target_spp Greattit
```

|                 |           |
|:----------------|:---------:|  
|Total no INDELs  | 1240366   |
|INDELs polarised | 593030    |
|Hotspots         | 287659    |
|Low spp coverage | 42732     |
|Ambiguous        | 21464     |
|Not in alignment | 295481    |
|Total unpolarised| 647336    |

## SNP polarisation

Firstly the sequence alignments across species for each SNP in the dataset were pulled out of the multiple alignment file with the following script:

```
./VARfromMAF.py  -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/SNP_data/gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.coveragefiltered.pass.repeatfilter.pass.vcf -maf /fastdata/bop15hjb/bird_alignments/UCSC_pipeline/multiple_zhang_param/Zebrafinch.Flycatcher.Greattit.maf -target_spp Greattit -out /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/SNP_data/polarisation/ -no_jobs 500
```

Secondly the output of the above script was used to annotate the ancestral state for each variant in the INDEL vcf as follows:

```
./polarise_vcf_snps.py -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/SNP_data/gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.coveragefiltered.pass.repeatfilter.pass.vcf -align_data /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/SNP_data/polarisation/all_variants.alignment_states.txt -target_spp Greattit -sub
```

|Category         | Number of SNPs |
|:----------------|:--------------:|
|Total no SNPs    | 10772087       |
|SNPs polarised   | 4430706        |
|Low spp coverage | 446296         |
|Ambiguous        | 2399524        |
|Not in alignment | 2515988        |
|Total unpolarised| 6341381        |
