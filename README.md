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

Getting Started
---------------

Usage Examples
--------------

    python DDRS.py -clincal_file ./example_clinical_infor.txt -snp_file ./example_snp.inf -min_g 20 -min_e 5

-clincal_file    : Patients clinical file, contain patients ID, follow-up information (time and event), with drug treatment(1) or not(0)<br>
-snp_file        : Genotype data<br>
-min_g(optional) : (Default is 20) Minimum patients size for each subgroup of pairwise KM survival analysis<br>
-min_e(optional) : (Default is 5)  Minimum patients with events for each subgroup<br>

