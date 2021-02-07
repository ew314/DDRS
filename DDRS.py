import sys,os,re
import pandas as pd
import lifelines
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines import CoxPHFitter
kmf = KaplanMeierFitter()
import argparse

parser = argparse.ArgumentParser(description="file info.")
parser.add_argument('-clincal_file',type=str,action='store',dest='clincal_filename',default=True)
parser.add_argument('-snp_file',type=str,action='store',dest='snp_filename',default=True)
parser.add_argument('-min_g',type=str,action='store',dest='ming',default=False)
parser.add_argument('-min_e',type=str,action='store',dest='mine',default=False)
args = parser.parse_args()

try:
	clincal_file=args.clincal_filename
except:
	print('need clincal file')

try:
	snp_file=args.snp_filename
except:
	print('need snp file')

try:
	ming=int(args.ming)
	if ming==0:
		ming=20
except:
	ming=20

try:
	mine=int(args.mine)
	if mine==0:
		mine=5
except:
	mine=5


globel_loc=os.getcwd()
os.system('mkdir %s/output'%globel_loc)

############## part.1 interaction analysis ###########
fr=open('%s/output/snp_drug_interaction_cox.txt'%globel_loc,'w')
fr.write('snp')
fr.write('\t')
fr.write('snp.inf')
fr.write('\t')

fr.write('drug.coef')
fr.write('\t')
fr.write('drug.exp_coef')
fr.write('\t')
fr.write('drug.se_coef')
fr.write('\t')
fr.write('drug.z')
fr.write('\t')
fr.write('drug.pval')
fr.write('\t')
fr.write('drug.lower95')
fr.write('\t')
fr.write('drug.upper95')
fr.write('\t')

fr.write('snp.coef')
fr.write('\t')
fr.write('snp.exp_coef')
fr.write('\t')
fr.write('snp.se_coef')
fr.write('\t')
fr.write('snp.z')
fr.write('\t')
fr.write('snp.pval')
fr.write('\t')
fr.write('snp.lower95')
fr.write('\t')
fr.write('snp.upper95')
fr.write('\t')

fr.write('snp*drug.coef')
fr.write('\t')
fr.write('snp*drug.exp_coef')
fr.write('\t')
fr.write('snp*drug.se_coef')
fr.write('\t')
fr.write('snp*drug.z')
fr.write('\t')
fr.write('snp*drug.pval')
fr.write('\t')
fr.write('snp*drug.lower95')
fr.write('\t')
fr.write('snp*drug.upper95')
fr.write('\n')
#########################################
f1=open(clincal_file,'r')
m1=f1.readlines()
f1.close()
barc=[]
Tinf={}
Einf={}
alldrugused=[]
alldrugusednot=[]
for i in range(1,len(m1)):
	p1=m1[i].strip().split('\t')
	barc.append(p1[0])
	Tinf[p1[0]]=int(p1[1])
	Einf[p1[0]]=int(p1[2])
	if p1[3] =='1':
		alldrugused.append(p1[0])
	if p1[3] =='0':
		alldrugusednot.append(p1[0])
########################################
fs=open(snp_file,'r')
ms=fs.readlines()
fs.close()
barm=ms[0].strip().split('\t')
barmt=[]
for p in barm:
	p=p.split(';')
	barmt.append(p[0])			
barmind={}
for i in range(0,len(barm)):
	tp=barm[i].split(';')
	barmind[tp[1]]=tp[0]
	
#####################################
allbar=list(set(barmt)&set(barc))	
group_num=20

filtered_SNP_ind=[]
for z in range(2,len(ms)):	
	p2=ms[z].strip().split('\t')
	mut=p2[0]				
	ale0=p2[1].split(';')
	ale1=p2[2].split(';')
	ale2=p2[3].split(';')
	ale0b=[]
	for ii in ale0:
		try:
			ale0b.append(barmind[ii])
		except:
			no=1
			
	ale1b=[]
	for ii in ale1:
		try:
			ale1b.append(barmind[ii])
		except:
			no=1 
	
	ale2b=[]
	for ii in ale2:
		try:
			ale2b.append(barmind[ii])
		except:
			no=1	
	
	ale0us=list(set(alldrugused)&set(ale0b)&set(allbar))
	ale1us=list(set(alldrugused)&set(ale1b)&set(allbar))
	ale2us=list(set(alldrugused)&set(ale2b)&set(allbar))
	ale0un=list(set(alldrugusednot)&set(ale0b)&set(allbar))
	ale1un=list(set(alldrugusednot)&set(ale1b)&set(allbar))
	ale2un=list(set(alldrugusednot)&set(ale2b)&set(allbar))			
	
	panduan=0
	if len(ale0us) >= group_num:
		 panduan+=1
	if len(ale1us) >= group_num:
		 panduan+=1
	if len(ale2us) >= group_num:
		 panduan+=1		 		 

	if len(ale0un) >= group_num:
		 panduan+=1
	if len(ale1un) >= group_num:
		 panduan+=1
	if len(ale2un) >= group_num:
		 panduan+=1		 		 
						
	if  panduan >= 6:
							
		T0us=[]
		E0us=[]
		T1us=[]
		E1us=[]
		T2us=[]
		E2us=[]
	
		T0un=[]
		E0un=[]
		T1un=[]
		E1un=[]
		T2un=[]
		E2un=[]
								
		for i in range(0,len(ale0us)):	
			T0us.append(Tinf[ale0us[i]])
			E0us.append(Einf[ale0us[i]])
	
		for i in range(0,len(ale1us)):	
			T1us.append(Tinf[ale1us[i]])
			E1us.append(Einf[ale1us[i]])
	
		for i in range(0,len(ale2us)):	
			T2us.append(Tinf[ale2us[i]])
			E2us.append(Einf[ale2us[i]])												
	
		for i in range(0,len(ale0un)):	
			T0un.append(Tinf[ale0un[i]])
			E0un.append(Einf[ale0un[i]])
	
		for i in range(0,len(ale1un)):	
			T1un.append(Tinf[ale1un[i]])
			E1un.append(Einf[ale1un[i]])
	
		for i in range(0,len(ale2un)):	
			T2un.append(Tinf[ale2un[i]])
			E2un.append(Einf[ale2un[i]])		
	
		numif=0				
		num0us=E0us.count(1)	
		num1us=E1us.count(1)
		num2us=E2us.count(1)
		num0un=E0un.count(1)
		num1un=E1un.count(1)	
		num2un=E2un.count(1)
		if num0us >= 5:
			numif=numif+1
		if num1us >= 5:
			numif=numif+1														
		if num2us >= 5:
			numif=numif+1			
		
		if num0un >= 5:
			numif=numif+1
		if num1un >= 5:
			numif=numif+1														
		if num2un >= 5:
			numif=numif+1						
			
		
		if numif >=6:
			data=[]
			for i in range(0,len(T0us)):								
				inf=[T0us[i],E0us[i],1,0,0]
				data.append(inf)
	
			for i in range(0,len(T1us)):								
				inf=[T1us[i],E1us[i],1,1,1]
				data.append(inf)
									
			for i in range(0,len(T2us)):								
				inf=[T2us[i],E2us[i],1,2,2]
				data.append(inf)
	
			for i in range(0,len(T0un)):								
				inf=[T0un[i],E0un[i],0,0,0]
				data.append(inf)			
	
			for i in range(0,len(T1un)):								
				inf=[T1un[i],E1un[i],0,1,0]
				data.append(inf)		
									
			for i in range(0,len(T2un)):								
				inf=[T2un[i],E2un[i],0,2,0]
				data.append(inf)									
								
	
																		
			try:
				inde=['date','status','drug','snp','interaction']
				allinfor = pd.DataFrame(data,columns=inde)
				cph = CoxPHFitter()
				cph.fit(allinfor, duration_col='date', event_col='status', show_progress=False)
				result=cph.summary
				name=list(result.axes[0])
				coef=list(result['coef'])
				exp_coef=list(result['exp(coef)'])
				se_coef=list(result['se(coef)'])
				z=list(result['z'])
				pval=list(result['p'])
				lower95=list(result['lower 0.95'])
				upper95=list(result['upper 0.95'])			
		
				fr.write(mut)
				fr.write('\t')
				fr.write(str(len(ale0))+';'+str(len(ale1))+';'+str(len(ale2)))
				fr.write('\t')
								
				fr.write(str(coef[0]))
				fr.write('\t')
				fr.write(str(exp_coef[0]))
				fr.write('\t')
				fr.write(str(se_coef[0]))
				fr.write('\t')
				fr.write(str(z[0]))
				fr.write('\t')
				fr.write(str(pval[0]))
				fr.write('\t')
				fr.write(str(lower95[0]))
				fr.write('\t')
				fr.write(str(upper95[0]))
				fr.write('\t')
									
									
				fr.write(str(coef[1]))
				fr.write('\t')
				fr.write(str(exp_coef[1]))
				fr.write('\t')
				fr.write(str(se_coef[1]))
				fr.write('\t')
				fr.write(str(z[1]))
				fr.write('\t')
				fr.write(str(pval[1]))
				fr.write('\t')
				fr.write(str(lower95[1]))
				fr.write('\t')
				fr.write(str(upper95[1]))
				fr.write('\t')
									
									
				fr.write(str(coef[2]))
				fr.write('\t')
				fr.write(str(exp_coef[2]))
				fr.write('\t')
				fr.write(str(se_coef[2]))
				fr.write('\t')
				fr.write(str(z[2]))
				fr.write('\t')
				fr.write(str(pval[2]))
				fr.write('\t')
				fr.write(str(lower95[2]))
				fr.write('\t')
				fr.write(str(upper95[2]))								
				fr.write('\n')
				filtered_SNP_ind.append(z)
			except:
				nosce=1



############## part.2 subgroup KM analysis ###########
focus_snp=[]
ft=open('%s/output/snp_drug_interaction_cox.txt'%globel_loc,'r')
mt=ft.readlines()
ft.close()

f2=open('%s/output/snp_drug_interaction_cox_filtered.txt'%globel_loc,'w')
f2.write(mt[0])
for i in range(1,len(mt)):
	pt=mt[i].strip().split('\t')
	try:
		if pt[6] != 'nan' and  pt[13] != 'nan' and pt[20] != 'nan' :
			if float(pt[20])<0.05 and float(pt[6]) < 0.05:
				focus_snp.append(pt[0])	
				f2.write(mt[i])
	except:
		no=1	
f2.close()

fr=open('%s/output/KM.subgroup.txt'%globel_loc,'w')
fr.write('snp')
fr.write('\t')
fr.write('snp.inf')
fr.write('\t')
fr.write('1.WT.drug.used')
fr.write('\t')
fr.write('2.HET.drug.used')
fr.write('\t')
fr.write('3.HOM.drug.used')
fr.write('\t')
fr.write('4.WT.drug.used.not')
fr.write('\t')
fr.write('5.HET.drug.used.not')
fr.write('\t')
fr.write('6.HOM.drug.used.not')
fr.write('\t')

fr.write('1&2.sta')
fr.write('\t')
fr.write('1&2.pval')
fr.write('\t')

fr.write('1&3.sta')
fr.write('\t')
fr.write('1&3.pval')
fr.write('\t')

fr.write('2&3.sta')
fr.write('\t')
fr.write('2&3.pval')
fr.write('\t')

fr.write('4&5.sta')
fr.write('\t')
fr.write('4&5.pval')
fr.write('\t')

fr.write('4&6.sta')
fr.write('\t')
fr.write('4&6.pval')
fr.write('\t')

fr.write('5&6.sta')
fr.write('\t')
fr.write('5&6.pval')
fr.write('\t')

fr.write('snp.coef')
fr.write('\t')
fr.write('snp.exp_coef')
fr.write('\t')
fr.write('snp.se_coef')
fr.write('\t')
fr.write('snp.z')
fr.write('\t')
fr.write('snp.pval')
fr.write('\t')
fr.write('snp.lower95')
fr.write('\t')
fr.write('snp.upper95')
fr.write('\n')

	
for z in range(2,len(ms)):
	p2=ms[z].strip().split('\t')
	if len(p2)==4 and p2[0] in focus_snp:
		mut=p2[0]				
		ale0=p2[1].split(';')
		ale1=p2[2].split(';')
		ale2=p2[3].split(';')
		ale0b=[]
		for ii in ale0:
			try:
				ale0b.append(barmind[ii])
			except:
				no=1
			
		ale1b=[]
		for ii in ale1:
			try:
				ale1b.append(barmind[ii])
			except:
				no=1 
	
		ale2b=[]
		for ii in ale2:
			try:
				ale2b.append(barmind[ii])
			except:
				no=1	
	
		ale0us=list(set(alldrugused)&set(ale0b)&set(allbar))
		ale1us=list(set(alldrugused)&set(ale1b)&set(allbar))
		ale2us=list(set(alldrugused)&set(ale2b)&set(allbar))
		ale0un=list(set(alldrugusednot)&set(ale0b)&set(allbar))
		ale1un=list(set(alldrugusednot)&set(ale1b)&set(allbar))
		ale2un=list(set(alldrugusednot)&set(ale2b)&set(allbar))			
		

		panduan=0
		if len(ale0us) >= group_num:
			 panduan+=1
		if len(ale1us) >= group_num:
			 panduan+=1
		if len(ale2us) >= group_num:
			 panduan+=1		 		 
	
		if len(ale0un) >= group_num:
			 panduan+=1
		if len(ale1un) >= group_num:
			 panduan+=1
		if len(ale2un) >= group_num:
			 panduan+=1		 		 
									
		if  panduan >= 6:
									
			T0us=[]
			E0us=[]
			T1us=[]
			E1us=[]
			T2us=[]
			E2us=[]
		
			T0un=[]
			E0un=[]
			T1un=[]
			E1un=[]
			T2un=[]
			E2un=[]
									
			for i in range(0,len(ale0us)):	
				T0us.append(Tinf[ale0us[i]])
				E0us.append(Einf[ale0us[i]])
		
			for i in range(0,len(ale1us)):	
				T1us.append(Tinf[ale1us[i]])
				E1us.append(Einf[ale1us[i]])
		
			for i in range(0,len(ale2us)):	
				T2us.append(Tinf[ale2us[i]])
				E2us.append(Einf[ale2us[i]])												
		
			for i in range(0,len(ale0un)):	
				T0un.append(Tinf[ale0un[i]])
				E0un.append(Einf[ale0un[i]])
		
			for i in range(0,len(ale1un)):	
				T1un.append(Tinf[ale1un[i]])
				E1un.append(Einf[ale1un[i]])
		
			for i in range(0,len(ale2un)):	
				T2un.append(Tinf[ale2un[i]])
				E2un.append(Einf[ale2un[i]])		
		
			numif=0				
			num0us=E0us.count(1)	
			num1us=E1us.count(1)
			num2us=E2us.count(1)
			num0un=E0un.count(1)
			num1un=E1un.count(1)	
			num2un=E2un.count(1)
			if num0us >= 5:
				numif=numif+1
			if num1us >= 5:
				numif=numif+1														
			if num2us >= 5:
				numif=numif+1			
			
			if num0un >= 5:
				numif=numif+1
			if num1un >= 5:
				numif=numif+1														
			if num2un >= 5:
				numif=numif+1						
		
			if numif >=6:
				
				pot0us=[]
				pot1us=[]
				pot2us=[]
				pot0un=[]
				pot1un=[]
				pot2un=[]
		
				for i in range(0,len(E0us)):
					if E0us[i] == 1:
						pot0us.append(T0us[i])
		
				for i in range(0,len(E1us)):
					if E1us[i] == 1:
						pot1us.append(T1us[i])
		
				for i in range(0,len(E2us)):
					if E2us[i] == 1:
						pot2us.append(T2us[i])									
		
				for i in range(0,len(E0un)):
					if E0un[i] == 1:
						pot0un.append(T0un[i])
		
				for i in range(0,len(E1un)):
					if E1un[i] == 1:
						pot1un.append(T1un[i])
		
				for i in range(0,len(E2un)):
					if E2un[i] == 1:
						pot2un.append(T2un[i])	
								
								
				results = logrank_test(T0us, T1us, event_observed_A=E0us, event_observed_B=E1us)
				pval12=results.p_value
				sta12=results.test_statistic
			
				results = logrank_test(T0us, T2us, event_observed_A=E0us, event_observed_B=E2us)
				pval13=results.p_value
				sta13=results.test_statistic							

				results = logrank_test(T1us, T2us, event_observed_A=E1us, event_observed_B=E2us)
				pval23=results.p_value
				sta23=results.test_statistic

				results = logrank_test(T0un, T1un, event_observed_A=E0un, event_observed_B=E1un)
				pval45=results.p_value
				sta45=results.test_statistic
	
				results = logrank_test(T0un, T2un, event_observed_A=E0un, event_observed_B=E2un)
				pval46=results.p_value
				sta46=results.test_statistic							
	
				results = logrank_test(T1un, T2un, event_observed_A=E1un, event_observed_B=E2un)
				pval56=results.p_value
				sta56=results.test_statistic


				data=[]
				for i in range(0,len(T0us)):								
					inf=[T0us[i],E0us[i],0]
					data.append(inf)
	
				for i in range(0,len(T1us)):								
					inf=[T1us[i],E1us[i],1]
					data.append(inf)
									
				for i in range(0,len(T2us)):								
					inf=[T2us[i],E2us[i],2]
					data.append(inf)
	
	
				inde=['date','status','snp']
				allinfor = pd.DataFrame(data,columns=inde)
				cph = CoxPHFitter()
				cph.fit(allinfor, duration_col='date', event_col='status', show_progress=False)
				result=cph.summary
				name=list(result.axes[0])
				coef=list(result['coef'])
				exp_coef=list(result['exp(coef)'])
				se_coef=list(result['se(coef)'])
				z=list(result['z'])
				pval=list(result['p'])
				lower95=list(result['lower 0.95'])
				upper95=list(result['upper 0.95'])			
									
				fr.write(mut)
				fr.write('\t')
				fr.write(str(len(ale0))+';'+str(len(ale1))+';'+str(len(ale2)))
				fr.write('\t')
				fr.write(str(len(ale0us)))
				fr.write('\t')
				fr.write(str(len(ale1us)))
				fr.write('\t')
				fr.write(str(len(ale2us)))
				fr.write('\t')
				fr.write(str(len(ale0un)))
				fr.write('\t')
				fr.write(str(len(ale1un)))
				fr.write('\t')
				fr.write(str(len(ale2un)))
				fr.write('\t')
								
				fr.write(str(sta12))
				fr.write('\t')
				fr.write(str(pval12))
				fr.write('\t')
								
				fr.write(str(sta13))
				fr.write('\t')
				fr.write(str(pval13))
				fr.write('\t')
								
				fr.write(str(sta23))
				fr.write('\t')
				fr.write(str(pval23))
				fr.write('\t')
								
				fr.write(str(sta45))
				fr.write('\t')
				fr.write(str(pval45))
				fr.write('\t')
				
				fr.write(str(sta46))
				fr.write('\t')
				fr.write(str(pval46))
				fr.write('\t')
								
				fr.write(str(sta56))
				fr.write('\t')
				fr.write(str(pval56))
				fr.write('\t')


				fr.write(str(coef[0]))
				fr.write('\t')
				fr.write(str(exp_coef[0]))
				fr.write('\t')
				fr.write(str(se_coef[0]))
				fr.write('\t')
				fr.write(str(z[0]))
				fr.write('\t')
				fr.write(str(pval[0]))
				fr.write('\t')
				fr.write(str(lower95[0]))
				fr.write('\t')
				fr.write(str(upper95[0]))
				fr.write('\n')
fr.close()
					
		
def groupfdr(pset):
	psort=pset[:]
	psort.sort()	
	pg={}
	for i in range(0,len(psort)):
		pg[psort[i]]=i+1
	np={}
	num=len(pset)
	for i in range(0,len(pset)):
		ind=pg[pset[i]]
		np[pset[i]]=float(pset[i])*float(num)/float(ind)
	return np
	
pval=0.05


############## part.3 subgroup fdr filter ###########
pgroup=[]
f1=open('%s/output/KM.subgroup.txt'%globel_loc,'r')
m1=f1.readlines()
f1.close()
for i in range(1,len(m1)):
	p1=m1[i].strip().split('\t')
	pgroup.append(float(p1[9]))
	pgroup.append(float(p1[11]))
	pgroup.append(float(p1[13]))
	pgroup.append(float(p1[15]))
	pgroup.append(float(p1[17]))
	pgroup.append(float(p1[19]))					
npgroup=groupfdr(pgroup)	


f1=open('%s/output/KM.subgroup.fdr.txt'%globel_loc,'w')
bar=m1[0].strip().split('\t')
for i in range(0,10):
	f1.write(bar[i])
	f1.write('\t')
f1.write('1&2.fdr.pval')
f1.write('\t')

for i in [10,11]:
	f1.write(bar[i])
	f1.write('\t')
f1.write('1&3.fdr.pval')
f1.write('\t')

for i in [12,13]:
	f1.write(bar[i])
	f1.write('\t')
f1.write('2&3.fdr.pval')
f1.write('\t')

for i in [14,15]:
	f1.write(bar[i])
	f1.write('\t')
f1.write('4&5.fdr.pval')
f1.write('\t')

for i in [16,17]:
	f1.write(bar[i])
	f1.write('\t')
f1.write('4&6.fdr.pval')
f1.write('\t')

for i in [18,19]:
	f1.write(bar[i])
	f1.write('\t')
f1.write('5&6.fdr.pval')
f1.write('\t')

for i in range(20,len(bar)):
	f1.write(bar[i])
	f1.write('\t')
f1.write('\n')
	
	
for i in range(1,len(m1)):
	p1=m1[i].strip().split('\t')	
	np1_2fdr=npgroup[float(p1[9])]
	np1_3fdr=npgroup[float(p1[11])]
	np2_3fdr=npgroup[float(p1[13])]
	np4_5fdr=npgroup[float(p1[15])]
	np4_6fdr=npgroup[float(p1[17])]
	np5_6fdr=npgroup[float(p1[19])]

	for i in range(0,10):
		f1.write(p1[i])
		f1.write('\t')
	f1.write(str(np1_2fdr))
	f1.write('\t')
	for i in [10,11]:
		f1.write(p1[i])
		f1.write('\t')
	f1.write(str(np1_3fdr))
	f1.write('\t')
	for i in [12,13]:
		f1.write(p1[i])
		f1.write('\t')
	f1.write(str(np2_3fdr))
	f1.write('\t')
	for i in [14,15]:
		f1.write(p1[i])
		f1.write('\t')
	f1.write(str(np4_5fdr))
	f1.write('\t')
	for i in [16,17]:
		f1.write(p1[i])
		f1.write('\t')
	f1.write(str(np4_6fdr))
	f1.write('\t')
	for i in [18,19]:
		f1.write(p1[i])
		f1.write('\t')
	f1.write(str(np4_6fdr))
	f1.write('\t')
	for i in range(20,len(p1)):
		f1.write(p1[i])
		f1.write('\t')
	f1.write('\n')
f1.close
	

f1=open('%s/output/KM.subgroup.fdr.txt'%globel_loc,'r')
m1=f1.readlines()
f1.close()
f1=open('%s/output/KM.subgroup.fdr.filter.txt'%globel_loc,'w')
bar=m1[0].strip().split('\t')
for p in bar:
	f1.write(p)
	f1.write('\t')
f1.write('sig_infor')
f1.write('\n')
for i in range(1,len(m1)):
	p1=m1[i].strip().split('\t')	
	np1_2fdr=float(p1[10])
	np1_3fdr=float(p1[13])
	np2_3fdr=float(p1[16])
	np4_5=float(p1[18])
	np4_6=float(p1[21])
	np5_6=float(p1[24])
	
	sig_12='nan'
	sig_13='nan'
	sig_23='nan'
	num=0	
	if np4_5 > pval and np4_6 > pval and np5_6 > pval:		
		if np1_2fdr < pval:
			sig_12='ref_ref.vs.ref_mut'
			num+=1
		if np1_3fdr < pval:
			sig_13='ref_ref.vs.mut_mut'
			num+=1		
		if np2_3fdr < pval:
			sig_23='ref_mut.vs.mut_mut'
			num+=1	
	if num > 0:
		for p in p1:
			f1.write(p)
			f1.write('\t')
		f1.write(sig_12+';'+sig_13+';'+sig_23)
		f1.write('\n')
f1.close()			