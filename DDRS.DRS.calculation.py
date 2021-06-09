import sys,os,re
import pandas as pd
import lifelines
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines import CoxPHFitter
kmf = KaplanMeierFitter()
import argparse

parser = argparse.ArgumentParser(description="file info.")
parser.add_argument('-snp_file',type=str,action='store',dest='snp_filename',default=True)
parser.add_argument('-filter_SNP',type=str,action='store',dest='filter_file',default=True)
args = parser.parse_args()

try:
	snp_file=args.snp_filename
except:
	print('need snp file')

try:
	filter_file=args.filter_file
except:
	print('need snp score file')



globel_loc=os.getcwd()
os.system('mkdir %s/output'%globel_loc)

############## part.4 DRS calculation ###########
	

f1=open(filter_file,'r')
m1=f1.readlines()
f1.close()
SNP_set=[]
SNP_score={}
for i in range(1,len(m1)):
	p1=m1[i].strip().split('\t')
	snp_inf=p1[0].split(';')
	SNP_set.append(snp_inf[0])
	SNP_score[snp_inf[0]]=float(p1[26])

fs=open(snp_file,'r')
ms=fs.readlines()
fs.close()
barm=ms[0].strip().split('\t')
barmt=[]
barscore={}
	
for p in barm:
	p=p.split(';')
	barmt.append(p[0])		
	barscore[p[0]]=0
	
barmind={}
for i in range(0,len(barm)):
	tp=barm[i].split(';')
	barmind[tp[1]]=tp[0]

for z in range(1,len(ms)):		
	p2=ms[z].strip().split('\t')
	snp=p2[0].split(';')[0]
	if snp in SNP_set:
		ale0=p2[1].split(';')
		ale1=p2[2].split(';')
		ale2=p2[3].split(';')
		for ii in ale1:
			try:
				barscore[barmind[ii]]=barscore[barmind[ii]]+SNP_score[snp]
			except:
				no=1 
	
		for ii in ale2:
			try:
				barscore[barmind[ii]]=barscore[barmind[ii]]+2*SNP_score[snp]
			except:
				no=1 
				
f1=open('%s/output/DRS.score'%globel_loc,'w')
for p in barmt:
	f1.write(p)
	f1.write('\t')
	f1.write(str(barscore[p]))
	f1.write('\n')
f1.close()
	
	


