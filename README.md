DDRS: Detection of Drug Response SNPs Specifically in Patients Receiving Drug Treatment
================================
This repository includes an approach to identify drug-response associated SNPs from clinical patients’ follow-up data by integrating cox proportional hazards model and Kaplan-Meier survival analysis. By using those patients didn’t receive any drug treatment as control, those prognosis correlated SNPs even without drug treatment were filtered out, and only focus on SNPs related to drug response.

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
1.SNP identification
--------------

    python DDRS.DRS.calculation.py -clincal_file ./example_clinical_infor.txt -snp_file ./example_snp.inf -min_g 20 -min_e 5

`-clincal_file`    : Patients clinical file, contain patients ID, follow-up information (time and event), with drug treatment(1) or not(0)<br>
`-snp_file`        : Genotype data<br>
`-min_g(optional)` : (Default is 20) Minimum patients size for each subgroup of pairwise KM survival analysis<br>
`-min_e(optional)` : (Default is 5)  Minimum patients with events for each subgroup<br>

Output
--------------
    ./output/snp_drug_interaction_cox.txt
            /snp_drug_interaction_cox_filtered.txt

`snp, snp.inf`                                         : SNP id, reference allele(0), alternative allele(1), patients number of Homozygous(reference allele),Heterozygosity and wild type(alternative allele)<br>
`drug.coef,exp_coef,se_coef,z,pval,lower95,upper95`      : CoxPHFitter[1] result of Drug term.<br>
`snp.coef,exp_coef,se_coef,z,pval,lower95,upper95`       : CoxPHFitter result of SNP term.<br>
`snp×drug.coef,exp_coef,se_coef,z,pval,lower95,upper95`  : CoxPHFitter result of SNP×Drug term.<br>

            /KM.subgroup.txt
            /KM.subgroup.fdr.txt
            /KM.subgroup.fdr.filter.txt
`-snp, snp.inf`                                         : SNP id, reference allele(0), alternative allele(1), patients number of Homozygous(reference allele),Heterozygosity and wild type(alternative allele)<br>
`1.WT.drug.used,2.HET.drug.used,3.HOM.drug.used,4.WT.drug.used.not,5.HET.drug.used.not,6.HOM.drug.used.not`: patients numbers of Homozygous,Heterozygosity and wild type in 6 six subgroups<br>
`1&2.sta,1&2.pval,1&2.fdr.pval`: logrank_test[2] result between Homozygous and Heterozygosity patients with drug treatment<br>
`1&3.sta,1&3.pval,1&3.fdr.pval`: logrank_test result between Homozygous and wild type patients with drug treatment<br>
`2&3.sta,2&3.pval,2&3.fdr.pval`: logrank_test result between Heterozygosity and wild type patients with drug treatment<br>
`4&5.sta,4&5.pval,4&5.fdr.pval`: logrank_test result between Homozygous and Heterozygosity patients without drug treatment<br>
`4&6.sta,4&6.pval,4&6.fdr.pval`: logrank_test result between Homozygous and wild type patients without drug treatment<br>
`5&6.sta,5&6.pval,5&6.fdr.pval`: logrank_test result between Heterozygosity and wild type patients without drug treatment<br>
`snp.coef,exp_coef,se_coef,z,pval,lower95,upper95`       : CoxPHFitter result of univariate Cox proportional hazards analysis about SNP in patients with drug treatment.<br>


2.DRS calculation
--------------

    python DDRS.SNP.identify.py -snp_file ./example_snp.inf -filter_SNP ./KM.subgroup.fdr.filter.txt

`-snp_file`        : Genotype data<br>
`-filter_SNP` : identified drug associated SNPs and cox coefficients of SNPs<br>

Output
--------------
    ./output/DRS.score

[1]https://lifelines.readthedocs.io/en/latest/Survival%20Regression.html#cox-s-proportional-hazard-model
[2]https://lifelines.readthedocs.io/en/latest/Examples.html?highlight=logrank#logrank-test
 
---------------------------------------
Yu Rong

ew314ew@stu.xjtu.edu.cn
