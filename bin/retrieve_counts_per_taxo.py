import sys
import re 
import os
from collections import defaultdict

abondances = open(sys.argv[1],'r').readlines()

design_exp = open(sys.argv[2],'r').readlines()

samples = []
for li in design_exp[1:]:
	li = li.strip().split('\t')
	samples.append(li[0])

taxo_to_count = defaultdict(int)
total = 0
columns = []
for i in range(len(abondances[0].strip().split('\t'))):
	if abondances[0].strip().split('\t')[i] in samples:
		columns.append(i)


for li in abondances[1:]:
	li = li.strip().split('\t')
	taxo = li[1].split(';')[4]

	for i in columns:
		total+= int(li[i])
		taxo_to_count[taxo] += int(li[i])

for a,b in taxo_to_count.items():
	print(a+'\t'+str(b)+'\t'+str(b/total*100))





	
