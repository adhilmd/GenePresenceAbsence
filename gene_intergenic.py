# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 12:17:47 2016

@author: adhil
"""

import re
import pandas as pd


#gff = args.gff
#allchr = args.achr
#flankinter = args.fint

gff = "/home/adhil/temp/Homo_sapiens.GRCh37.75.gtf"
allchr = "/home/adhil/temp/allchr.bed"
flankinter = 100

f1 = open(gff, 'r')
exon = []
gene = []
required_features = ["IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "polymorphic_pseudogene", "protein_coding", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene"]
i = 0
for line in f1:
    if not line.startswith('#'):
        temp = line.split('\t')
        if (str(temp[1]) in required_features):
            if ((str(temp[0]).find("_") == -1) and (str(temp[0]).find(".") == -1)):
                if (str(temp[2]).startswith('exon')):
                    ex = [str('chr')+ str(temp[0]), str(temp[3]), str(temp[4]), str(temp[8]).strip('\n')]
                    exon.append(ex)
                if (str(temp[2]).startswith('gene')):
                    i = i+1                
                    b = [str(str('chr') + str(temp[0])), str(int(temp[3])-flankinter), str(int(temp[4])+flankinter), str(temp[8]).strip('\n')]
                    gene.append(b)        

te = pd.DataFrame(gene, columns=['chr','start','end','feature'])

f2 = open(allchr, 'r')
fin1 = []
for line in f2:
    fin1.append(line.strip('\n').split('\t'))

te1 = pd.DataFrame(fin1, columns=['chr','start','end','feature'])
te1['start'] = te1['end']
frames = [te,te1]
result = pd.concat(frames)
result[['start', 'end']] = result[['start', 'end']].astype(int)
uniquechr = list(set(result['chr'].tolist()))
intergenic = pd.DataFrame(data=None, columns=['chr','start','end','feature'])

for item in uniquechr:
    temp = result.loc[result['chr'] == item]
    temp = temp.sort(['start'])
    end = temp['end'].tolist()
    start = temp['start'].tolist()
    if (len(start)>1):
        startfi = start[0]
        del start[0]
        del end[len(end)-1]
        if (startfi != 0):
            start.insert(0,startfi)
            end.insert(0,1)
        chrv = re.sub(',$','',len(start)*str(str(item)+str(','))).split(',')
        feat = [str(item)+str('_')+str('intergenic')+str(k) for k in range(1,len(start)+1)] 
        tempd = pd.DataFrame({'chr':chrv, 'start':end, 'end':start, 'feature':feat})
        frames = [intergenic,tempd]
        intergenic = pd.concat(frames)
        

cols = ['chr','start','end','feature']
intergenicpd = intergenic[cols]
intergenicpd = intergenicpd[intergenicpd['start'] < intergenicpd['end']]
exonpd = pd.DataFrame(exon, columns=['chr','start','end','feature'])
genepd = te



print intergenicpd
print exonpd
