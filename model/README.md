# Testing an extension of the Glémin et al. (2015) model
Henry Juho Barton  
Department of Animal and Plant Sciences, The University of Sheffield  

# Introduction

This document describes the pipeline for testing an extension of the Glémin et al. (2015)  model, for use on INDEL site frequency data. This novel extension estimates the mutation and selection parameters for insertions and deletions. The model also accounts for demography and polarisation error, thus overcoming many pitfuls of previous INDEL work. We tested this novel maximum likelihood approach with simulated data. For an in-depth description of the model see the note here: <https://github.com/henryjuho/parus_indel/blob/master/model/indel_model_v2.pdf>. 

# Simulation and parameter estimation

Site frequency data was simulated under a range of parameters and fed into the model using a python script as follows:

```
./kai_model_test_range.py -t1_r 500,1501,500 -t2_r 500,1501,500 -g1_r " -50,11,30" -g2_r " -50,11,30" -e1_r 0.1,0.4,0.1 -e2_r 0.1,0.4.0.1 -out_dir /data/bop15hjb/glemin_sim_results/run3
```

# Simulation results

Simulation results can be seen for gamma and theta at different simulated polarisation erros:

  - gamma <https://github.com/henryjuho/parus_indel/blob/master/model/run3_gamma.jpg>
  - theta <https://github.com/henryjuho/parus_indel/blob/master/model/run3_theta.jpg>
