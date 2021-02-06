DDRS: Detection of Drug Response SNPs Specifically in Patients Receiving Drug Treatment
================================
This repository includes an approach to identify drug-response associated SNPs from clinical patients’ follow-up data by integrating cox proportional hazards model and Kaplan-Meier survival analysis. 

![pipeline](https://github.com/ew314/DDRS/blob/main/pipeline/4.figure.1.pipeline.github.jpg)

# PREREQUISITE
The approach were conducted by using Python 2.7.16. 
Following Python packages should be installed:
<ul>
<li><p>pandas</p></li>
<li><p>lifelines</p></li>
</ul>

---
use :
  smartpca.perl: run PCA on input genotype data (calls smartpca)
  smarteigenstrat.perl: run EIGENSTRAT stratification correction.  This program 
    supports all 5 file formats, and supports quantitative phenotypes.
  gc.perl: apply Genomic Control (Devlin and Roeder, 1999) to the
    association statistics computed by EIGENSTRAT.
    
We note that the programs eigenstrat and eigenstratQTL of EIGENSOFT version 2.0
have been replaced by smarteigenstrat.perl.  However, we have retained the old
programs for backwards compatibility (see below).

See ./example.perl and ./exampleQTL.perl for toy examples using our programs.

  smartpca.perl: run PCA on input genotype data (calls smartpca)
  smarteigenstrat.perl: run EIGENSTRAT stratification correction.  This program 
    supports all 5 file formats, and supports quantitative phenotypes.
  gc.perl: apply Genomic Control (Devlin and Roeder, 1999) to the
    association statistics computed by EIGENSTRAT.
    
We note that the programs eigenstrat and eigenstratQTL of EIGENSOFT version 2.0
have been replaced by smarteigenstrat.perl.  However, we have retained the old
programs for backwards compatibility (see below).

See ./example.perl and ./exampleQTL.perl for toy examples using our programs.
