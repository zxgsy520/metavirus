#!/usr/bin/env python
import sys
GOinput=sys.argv[1]  ## GO.annot
GO_level4=sys.argv[2] ## GO_level4.deal.txt
f=open(GOinput,"r")
d={}
for line in f.readlines():
        k = line.strip().split('\t')
        for i in k[1:]:
                d.setdefault(i,[]).append(k[0])
f.close()

print('GO:level1\tGO:level2\tGO:level3\tGO:level4\tGene_Number\tGenes_id')
f=open(GO_level4,"r")
dl={}
for line in f.readlines():
	k = line.strip().split('\t')
	if len(k)==4:
		level4=k[3].split('(')[0]
		if level4 in d:
			print('\t'.join(k)+'\t'+str(len(d[level4]))+'\t'+','.join(d[level4]))
	if len(k)>4:
		n=0
		gene=''
		for i in k[3:]:
			level = i.split('(')[0]
			if level in d:
				n += len(d[level])
				gene+=','.join(d[level])+','
		if n >0:
			gene0 = list(set(gene.strip(',').split(',')))
			n0 = len(gene0)
			print('\t'.join(k[0:4])+'\t'+str(n0)+'\t'+','.join(gene0))
f.close()
				
