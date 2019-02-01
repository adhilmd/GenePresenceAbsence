# GenePresenceAbsence
This tool uses RNA seq data to identify gene presence and absence using off target reads

The RNA-Seq reads are aligned to the exonic regions (coding regions) of genes. 
The number of reads aligning to the exonic region of a gene is proportional to the expression level of the gene and the 
length of the exonic region. However, due to sequencing errors and the statistical nature of reference assembly, errors are 
incorporated  into the read alignment process. This results in some leakage or spillover in the data. There will be a few 
reads which will align to genes purely by chance (noise). Similarly, there will also be some reads that will align to the 
intergenic regions by chance (noise). In order to identify the noise, PA-call uses the feature level 
Reads Per Nucleotide (RPN) statistic. RPN is defined as,

RPN = R/L where R is the number of reads in a region of length L nucleotides

It may be noted that while, the read alignment by chance at a whole genome level is modeled using Poisson distribution, 
the RPN distribution is lognormal. 

The RPN statistic is calculated for the genic as well as intergenic regions to identify the noise. PA-call models 
the exonic (collapsed to genic level) and non-coding/intergenic level RPN, and uses these to generate p-values for 
each gene and to obtain an RPN cut-off. Using this RPN cut-off, the sample level presence/absence of a gene is identified.

![rpn](https://user-images.githubusercontent.com/18418058/52130792-6de28b00-263b-11e9-95fd-d904c33b5502.jpeg)


RPN distribution across the genome for genic and intergenic (noise) along with regions of 
Absence, Marginal, and Presence




