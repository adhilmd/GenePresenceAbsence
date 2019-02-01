import pysam, glob,numpy, time, sys
import multiprocessing
from os import getpid
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy, datetime
import math
from math import log as log
import re
import scipy.stats
import pandas as pd
import argparse
import pdb
import os
suffix = ''

def argparser():
    parser = argparse.ArgumentParser(prog = '''\nRNAPA_call.py''', description='''\n---------------RNA PA call----------------- \n
    \n[Date: 25 February 2015], \n[help: python PAcall -h], \n[Author: Mohamood Adhil]\n''', usage = 'PAcall.py [-h] [-i <Folder with bam and bai files> -g <GFF/GTF format file> -od <Output directory path> -ch <Chromosomes bed file> -t <TissueName> -te <TissueLevelExpressionFile> -c <P-value threshold> -n <Number Core Process>]')
    parser.add_argument('-i','--InputBam', type=str, dest='i', help="Path for bam (Mandatory)", action = 'store', required = True)
    parser.add_argument('-od','--Output', type=str, dest='o', help="Output directory path(Mandatory)", action = 'store', required = True)
    parser.add_argument('-gf','--gfffile', type=str, dest='gff', help="GFF file path (Mandatory)", action = 'store', required = True, default ='')
    parser.add_argument('-cf','--chrlen', type=str, dest='chr', help="Chromosome length file (Mandatory)", action = 'store', required = True, default ='')
    parser.add_argument('-tf','--tissuefile', type=str, dest='tfile', help="Base line tissue expression file path (Tab seperated where rownames are genes and column names are the tissue name)", action = 'store', required = False, default ='')
    parser.add_argument('-t','--TissueName', type=str, dest='t', help="Tissue for expression annotation (Column name in the Tissue File)", action = 'store', required = False, default ='')
    parser.add_argument('-c','--Pval', type=float, dest='c', help="P-value for call (Mandatory)", action = 'store', required = True)
    parser.add_argument('-n','--numProc', type=int, dest='np', help="Number of process for parallel computing (Mandatory)", action = 'store', required = True)
    parser.add_argument('-sampleN','--sample', type=str, dest='s', help="sample name (Mandatory)", action = 'store', required = True)
    args = parser.parse_args()
    return(args)

def gffparser(gff,allchr,flankinter):
    f1 = open(gff, 'r')
    exon = []
    gene = []
    required_features = ["IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "polymorphic_pseudogene", "protein_coding", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene"]
    i = 0
    for line in f1:
        if not line.startswith('#'):
            temp = line.split('\t')
            if (str(temp[1]) in required_features):
                if ((str(temp[0]).find("_") == -1) and (str(temp[0]).find(".") == -1) and (str(temp[0]).find("MT") == -1)):
                    if (str(temp[2]).startswith('exon')):
                        ex = [str(temp[0]), str(temp[3]), str(temp[4]), str(temp[8]).strip('\n')]
                        exon.append(ex)
                    if (str(temp[2]).startswith('gene')):
                        i = i+1                
                        b = [str(temp[0]), str(int(temp[3])-flankinter), str(int(temp[4])+flankinter), str(temp[8]).strip('\n')]
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
    return(exonpd,intergenicpd)

def process(x):
    path,reg,bn = x
    bamfile = pysam.AlignmentFile(path,'rb')  
    counts = dict()
    for i in reg:
        total, mapped = 0, 0
        for rec in bamfile.fetch(i[0],int(i[1]),int(i[2])):
            total += 1
            
            if rec.is_unmapped:
                continue
        
            mapped += 1
        counts[i[3]] = [mapped,total,int(i[2])-int(i[1])+1]
        bamfile.reset()
    bamfile.close()
    return counts
    
def counts (bam,reftype,outdir,reg):
    batchsize = 1000
    if reftype == 'gen':
        suffix = 'gen'
    elif reftype == 'int':
        suffix = 'int'
    out = '.counts.'+suffix
    start = time.time()    
    params = [(bam,reg[this:this+batchsize],i) for i,this in enumerate(range(0,len(reg),batchsize))]
    pool = multiprocessing.Pool(processes=nprocess)
    q = pool.map(process,params)
    counts = dict((k,v) for d in q for (k,v) in d.items())
    print bam,time.time()-start
    return(counts,out,time.time()-start)

def main():
    start_timestamp = str(datetime.datetime.fromtimestamp(time.mktime(datetime.datetime.now().timetuple())))
    gfffile = args.gff
    allbams = [args.i]
    print (allbams)
    allchr = args.chr
    outdir = args.o
    tfile = args.tfile
    tissue = args.t
    thresh = args.c
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    print outdir
    exon,intergenic=gffparser(gfffile,allchr,200)
    logf = open(outdir+'/'+args.s + '_' + 'PAcall.log', 'wb')
    tm_trk = str(datetime.datetime.fromtimestamp(time.mktime(datetime.datetime.now().timetuple())))
    logf.writelines(tm_trk+'\n')
    logf.writelines(tm_trk+': Started running PA-call at ' + start_timestamp + '\n')
    #Exon:
    reg = []
    tm_trk = str(datetime.datetime.fromtimestamp(time.mktime(datetime.datetime.now().timetuple())))
    logf.writelines(tm_trk+': Identiying exonic regions\n')
    for row in exon.iterrows():
        line=list(list(row)[1])
        reg.append(line)
    #Intergenic:
    tm_trk = str(datetime.datetime.fromtimestamp(time.mktime(datetime.datetime.now().timetuple())))
    logf.writelines(tm_trk+ ': Identiying intergenic regions\n')
    reg2 = []
    for row in intergenic.iterrows():
        line=list(list(row)[1])
        reg2.append(line)
    for bam in allbams:
        print (bam)
        tm_trk = str(datetime.datetime.fromtimestamp(time.mktime(datetime.datetime.now().timetuple())))
        logf.writelines(tm_trk+': Beginning analysis for file: ' + bam + '\n')
        counts1,out,time_tk = counts(bam,'gen',outdir,reg)
        logf.writelines('\t Time taken for exon region mapping: ' + str(time_tk) + '\n')
        out1 = args.o+'/'+args.s+out
        f = open(out1,'w')
        for i in counts1:
           f.write('\t'.join([i]+map(str,counts1[i]))+'\n')
        f.close()
        counts2,out,time_tk = counts(bam,'int',outdir,reg2)
        logf.writelines('\t Time taken for intergenic region mapping: ' + str(time_tk) + '\n')
        out1 = args.o+'/'+args.s+out
        f = open(out1,'w')
        for i in counts2:
           f.write('\t'.join([i]+map(str,counts2[i]))+'\n')
        f.close()
        tm_trk = str(datetime.datetime.fromtimestamp(time.mktime(datetime.datetime.now().timetuple())))
        logf.writelines('\t '+tm_trk+': Calculating RPN values'+ '\n')
        genic_reads = numpy.array([float(counts1[gcounts][0]) for gcounts in counts1])
        genic_lengths = numpy.array([float(counts1[gcounts][2]) for gcounts in counts1])
        genic_ratios = genic_reads/genic_lengths
        genes = [feature.split(';') for feature in counts1 ]
        gen_names = [f[3].strip() for f in genes]
        genic_X = [[gen_names[i], genic_ratios[i]] for i in xrange(0,len(genes))]
        genic_df = pd.DataFrame(genic_X)
        genic_df.columns = ['gene_name', 'ratio']
        g2 = genic_df.groupby('gene_name')
        g3 = g2.aggregate(numpy.mean)
        g4=g3.values.tolist()
        genes = g3.index
        gene_names= genes.values.tolist()
        x = numpy.log([10**-20 if d[0]==0 else d[0] for d in g4])
        intergenic_reads = numpy.array([float(counts2[intcounts][0]) for intcounts in counts2])
        intergenic_lengths = numpy.array([float(counts2[intcounts][2]) for intcounts in counts2])
        y = numpy.log([10**-20 if d==0 else d for d in list(intergenic_reads/intergenic_lengths)])
        #Plot distribution
        tm_trk = str(datetime.datetime.fromtimestamp(time.mktime(datetime.datetime.now().timetuple())))
        logf.writelines('\t '+ tm_trk+': Calculating RPN frequencies' + '\n')
        x_freq, x_bins = numpy.histogram(x, bins = 100, range = (-15,10))
        x_mid = [(x_bins[i]+x_bins[i+1])/2 for i in xrange(0,len(x_freq))]
        xhist = [[x_mid[i],x_freq[i]] for i in xrange(0,len(x_freq))]
        xhistdf = pd.DataFrame(xhist, columns = ['Mid', 'Frequency'])
        histfile1 = outdir+'/'+args.s+'_genic.hist'
        xhistdf.to_csv(histfile1, sep  ='\t', index = False)
        y_freq, y_bins = numpy.histogram(y, bins = 200, range = (-15,10))
        y_mid = [(y_bins[i]+y_bins[i+1])/2 for i in xrange(0,len(y_freq))]
        yhist = [[y_mid[i],y_freq[i]] for i in xrange(0,len(y_freq))]
        yhistdf = pd.DataFrame(yhist, columns = ['Mid', 'Frequency'])
        histfile2 = outdir+'/' + args.s+'_intergenic.hist'
        yhistdf.to_csv(histfile2, sep  ='\t', index = False)
        #Make PA calls
        #Define noise distribution:
        tm_trk = str(datetime.datetime.fromtimestamp(time.mktime(datetime.datetime.now().timetuple())))
        logf.writelines('\t '+ tm_trk+': Identifying distributions' + '\n')
        noise_dist = [val for val in y if val != min(y)]
        n_mu,n_sd = scipy.stats.norm.fit(noise_dist)
        #EM decomposition
        from sklearn import mixture
        data = numpy.log([i[0] for i in g4 if i[0]!=0])
        data = data.reshape((len(data),1))
        g = mixture.GMM(n_components=2)
        g.fit(data)
        wts, means, cov = g.weights_, g.means_, g.covars_ 
        mu1 = means[0]
        sd1 = math.sqrt(cov[0])
        mu2 = means[1]
        sd2 = math.sqrt(cov[1])
        if mu1 > mu2:
            le_mu = mu2
            le_sd = sd2
            he_mu = mu1
            he_sd = sd1
        elif mu2 > mu1:
            le_mu =mu1
            le_sd = sd1
            he_mu = mu2
            he_sd = sd2
        pval_A = [1- scipy.stats.norm.cdf((val-n_mu)/n_sd) for val in x]
        pval_LE = [2*(1- scipy.stats.norm.cdf(math.fabs(val-le_mu)/le_sd)) for val in x]
        pval_HE = [1*(scipy.stats.norm.cdf((val-he_mu)[0]/he_sd)) for val in x]
        tm_trk = str(datetime.datetime.fromtimestamp(time.mktime(datetime.datetime.now().timetuple())))
        logf.writelines('\t ' + tm_trk +': Making PA calls' + '\n')
        #Append with tissue data
        if ((tfile != '') and (tissue != '')):
            tissue_exprs = pd.DataFrame.from_csv(tfile,sep='\t')
            baseline_expr = tissue_exprs[tissue]
            pcalls = []
        for i in xrange(0,len(gene_names)):
            nam = gene_names[i].split()[1].replace('"','').strip()
            basl = 'NA'
            if ((tfile != '') and (tissue != '')):
                if nam in baseline_expr.index:
                    basl = baseline_expr[nam]
            if (pval_HE[i] > thresh):
                pacall = "P"
            elif (pval_A[i] > thresh):
                pacall = "A"
            elif ((pval_A[i] < thresh) and (pval_HE[i] < thresh)):
                pacall = "M"
            temp = [nam,x[i],pacall,pval_A[i], basl, pval_LE[i], pval_HE[i]]
            pcalls.append(temp)
        pcallsdf = pd.DataFrame(pcalls, columns=['Gene Name','Reads/Length','P/A call','Absence p-val','Tissue baseline','LE p-val', 'HE p-val'])
        #Write file
        resfile = outdir+'/'+ args.s+'.calls'
        pcallsdf.to_csv(resfile, sep ='\t', index = False)
        finresfile = outdir+'/' + args.s+'final.calls'
        finresdf = pcallsdf.loc[:,'Gene Name':'Tissue baseline']
        finresdf.to_csv(finresfile, sep ='\t', index = False)
    end_time = str(datetime.datetime.fromtimestamp(time.mktime(datetime.datetime.now().timetuple())))
    logf.writelines('Completed PA-call run at ' + tm_trk +'\n')
    logf.close()
    return 0
    
if (__name__ == "__main__"):
    args = argparser()
    nprocess = args.np
    status = main()
    sys.exit(status)
