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
Relative activity predictor
---------------------------------
this model is trained to predict the relative activity of mismatched sgRNA with its original matched sgRNA.<br>
pre-trained model is https://github.com/ew314/ew314/blob/main/sgRNA_designer/seq2seq_attention/NBT_float_model_120.h5<br>
