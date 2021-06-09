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

    python DDRS.py -clincal_file ./example_clinical_infor.txt -snp_file ./example_snp.inf -min_g 20 -min_e 5

`-clincal_file`    : Patients clinical file, contain patients ID, follow-up information (time and event), with drug treatment(1) or not(0)<br>
`-snp_file`        : Genotype data<br>
`-min_g(optional)` : (Default is 20) Minimum patients size for each subgroup of pairwise KM survival analysis<br>
`-min_e(optional)` : (Default is 5)  Minimum patients with events for each subgroup<br>

Output
--------------
    ./output/snp_drug_interaction_cox.txt
            /snp_drug_interaction_cox_filtered.txt

`-snp, snp.inf`                                         : SNP id, reference allele(0), alternative allele(1), patients number of Homozygous(reference allele),Heterozygosity and wild type(alternative allele)<br>
`-drug.coef,exp_coef,se_coef,z,pval,lower95,upper95`      : CoxPHFitter[1] result of Drug term.<br>
`-snp.coef,exp_coef,se_coef,z,pval,lower95,upper95`       : CoxPHFitter result of SNP term.<br>
`-snp×drug.coef,exp_coef,se_coef,z,pval,lower95,upper95`  : CoxPHFitter result of SNP×Drug term.<br>

            /KM.subgroup.txt
            /KM.subgroup.fdr.txt
            /KM.subgroup.fdr.filter.txt
`-snp, snp.inf`                                         : SNP id, reference allele(0), alternative allele(1), patients number of Homozygous(reference allele),Heterozygosity and wild type(alternative allele)<br>
`1.ale0.drug.used`
`2.ale1.drug.used`
`3.ale2.drug.used`
`4.ale0.drug.used.not`
`5.ale1.drug.used.not`
`6.ale2.drug.used.not`


2.DRS calculation
--------------

    python DDRS.py -clincal_file ./example_clinical_infor.txt -snp_file ./example_snp.inf -min_g 20 -min_e 5

*-clincal_file*    : Patients clinical file, contain patients ID, follow-up information (time and event), with drug treatment(1) or not(0)<br>
*-snp_file*        : Genotype data<br>
*-min_g(optional)* : (Default is 20) Minimum patients size for each subgroup of pairwise KM survival analysis<br>
*-min_e(optional)* : (Default is 5)  Minimum patients with events for each subgroup<br>

Output
--------------
    ./output/snp_drug_interaction_cox.txt
            /snp_drug_interaction_cox_filtered.txt
            /KM.subgroup.txt
            /KM.subgroup.fdr.txt
            /KM.subgroup.fdr.filter.txt
 
---------------------------------------
Yu Rong

ew314ew@stu.xjtu.edu.cn
