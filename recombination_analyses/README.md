## Recombination window analysis

2Mb window coordinates were generated and mean recombination rates were calculated per window using scripts from the biased gene conversion project. Windows without recombination rate estimates and with less than 500 polarisable INDELs were excluded.

```
$ qsub prepare_windows.sh
$ cd ~/parus_indel/recombination_analyses/
$ python consolidate_window_info.py > all_2Mb_windows.txt
$ python filter_windows.py all_2Mb_windows.txt
$ ./add_window_indel_subs.py -windows filtered_2Mb_windows.txt -wga_bed /fastdata/bop15hjb/bird_alignments/UCSC_pipeline/multiple_zhang_param/Zebrafinch.Flycatcher.Greattit.wga.bed.gz > filtered_2Mb_windows_with_subs.txt
```

The relationship between theta and recombination rate and Tajima's D and recombination was explored for both insertions and deletions.

```
$ Rscript window_summary.R
```

![sum_stats](window_summary.png)

As well as for INDEL divergence and recombination.

<img src="div_recomb_plot.png" width=600 height=300>


### Runs with sfs corrected by estimated polarisation error

```bash
python consolidate_window_info.py -correct_sfs > all_2Mb_windows_corrected_sfs.txt
python filter_windows.py all_2Mb_windows_corrected_sfs.txt
Rscript window_summary.R filtered_all_2Mb_windows_corrected_sfs.txt window_summary_corrected_sfs.png 

[1] "insertions"
   estimate     p.value statistic   n gp   Method
1 0.1457439 0.008710691  2.639399 324  1 spearman
[1] "deletions"
    estimate   p.value statistic   n gp   Method
1 0.07891814 0.1570563  1.418358 324  1 spearman
```

![](window_summary_corrected_sfs.png)