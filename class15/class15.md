---
title: "Bioinformatics Class 15"
author: "Barry Grant"
output: 
  html_document: 
    keep_md: yes
---



## Differential Expression Analysis

Import data files for DESeq analysis (countData and colData).


```r
library(DESeq2)
```

```
## Warning: package 'DESeq2' was built under R version 3.4.2
```

```
## Warning: package 'S4Vectors' was built under R version 3.4.2
```

```
## Warning: package 'BiocGenerics' was built under R version 3.4.2
```

```
## Warning: package 'IRanges' was built under R version 3.4.2
```

```
## Warning: package 'GenomicRanges' was built under R version 3.4.3
```

```
## Warning: package 'GenomeInfoDb' was built under R version 3.4.2
```

```
## Warning: package 'SummarizedExperiment' was built under R version 3.4.3
```

```
## Warning: package 'Biobase' was built under R version 3.4.2
```

```
## Warning: package 'DelayedArray' was built under R version 3.4.2
```


```r
metaFile <- "data/GSE37704_metadata.csv"
countFile <- "data/GSE37704_featurecounts.csv"

# Import metadata and take a peak
colData <-  read.csv(metaFile, row.names=1)
head(colData)
```

```
##               condition
## SRR493366 control_sirna
## SRR493367 control_sirna
## SRR493368 control_sirna
## SRR493369      hoxa1_kd
## SRR493370      hoxa1_kd
## SRR493371      hoxa1_kd
```


```r
# Import countdata
countData <- read.csv(countFile, row.names=1)
head(countData)
```

```
##                 length SRR493366 SRR493367 SRR493368 SRR493369 SRR493370
## ENSG00000186092    918         0         0         0         0         0
## ENSG00000279928    718         0         0         0         0         0
## ENSG00000279457   1982        23        28        29        29        28
## ENSG00000278566    939         0         0         0         0         0
## ENSG00000273547    939         0         0         0         0         0
## ENSG00000187634   3214       124       123       205       207       212
##                 SRR493371
## ENSG00000186092         0
## ENSG00000279928         0
## ENSG00000279457        46
## ENSG00000278566         0
## ENSG00000273547         0
## ENSG00000187634       258
```

Remove the first `length` col from the `countData` data.frame


```r
# Note we need to remove the odd first $length col
countData <- as.matrix(countData[,-1])
head(countData)
```

```
##                 SRR493366 SRR493367 SRR493368 SRR493369 SRR493370
## ENSG00000186092         0         0         0         0         0
## ENSG00000279928         0         0         0         0         0
## ENSG00000279457        23        28        29        29        28
## ENSG00000278566         0         0         0         0         0
## ENSG00000273547         0         0         0         0         0
## ENSG00000187634       124       123       205       207       212
##                 SRR493371
## ENSG00000186092         0
## ENSG00000279928         0
## ENSG00000279457        46
## ENSG00000278566         0
## ENSG00000273547         0
## ENSG00000187634       258
```

This looks better but there are lots of zero entries in there so let’s get rid of them as we have no data for these.


```r
# Filter count data where you have 0 read count across all samples.
countData = countData[rowSums(countData)>1, ]
head(countData)
```

```
##                 SRR493366 SRR493367 SRR493368 SRR493369 SRR493370
## ENSG00000279457        23        28        29        29        28
## ENSG00000187634       124       123       205       207       212
## ENSG00000188976      1637      1831      2383      1226      1326
## ENSG00000187961       120       153       180       236       255
## ENSG00000187583        24        48        65        44        48
## ENSG00000187642         4         9        16        14        16
##                 SRR493371
## ENSG00000279457        46
## ENSG00000187634       258
## ENSG00000188976      1504
## ENSG00000187961       357
## ENSG00000187583        64
## ENSG00000187642        16
```


Setup the DESeqDataSet object required for the DESeq() function and then run the DESeq pipeline. This is again similar to our last days hands-on session.


```r
dds <- DESeqDataSetFromMatrix(countData=countData,
                             colData=colData,
                             design=~condition)
```

Now we can run our DE analysis


```r
dds <- DESeq(dds)
```

```
## estimating size factors
```

```
## estimating dispersions
```

```
## gene-wise dispersion estimates
```

```
## mean-dispersion relationship
```

```
## final dispersion estimates
```

```
## fitting model and testing
```


```r
dim(countData)
```

```
## [1] 15280     6
```



```r
dds
```

```
## class: DESeqDataSet 
## dim: 15280 6 
## metadata(1): version
## assays(3): counts mu cooks
## rownames(15280): ENSG00000279457 ENSG00000187634 ...
##   ENSG00000276345 ENSG00000271254
## rowData names(21): baseMean baseVar ... deviance maxCooks
## colnames(6): SRR493366 SRR493367 ... SRR493370 SRR493371
## colData names(2): condition sizeFactor
```

## Extracting our results table


```r
res <- results(dds)
summary(res)
```

```
## 
## out of 15280 with nonzero total read count
## adjusted p-value < 0.1
## LFC > 0 (up)     : 4352, 28% 
## LFC < 0 (down)   : 4400, 29% 
## outliers [1]     : 0, 0% 
## low counts [2]   : 590, 3.9% 
## (mean count < 1)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

Too many genes here so lets increase our thresholds of log2 fold change and decrease our FDR *p* value.
 

```r
resSig <- results(dds, lfcThreshold=2, alpha=0.05)
summary(resSig)
```

```
## 
## out of 15280 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)     : 99, 0.65% 
## LFC < 0 (down)   : 134, 0.88% 
## outliers [1]     : 0, 0% 
## low counts [2]   : 1482, 9.7% 
## (mean count < 2)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```



Let’s reorder these results by p-value and call summary() on the results object to get a sense of how many genes are up or down-regulated at the default FDR of 0.1.



```r
res <- res[order(res$pvalue),]
summary(res)
```

```
## 
## out of 15280 with nonzero total read count
## adjusted p-value < 0.1
## LFC > 0 (up)     : 4352, 28% 
## LFC < 0 (down)   : 4400, 29% 
## outliers [1]     : 0, 0% 
## low counts [2]   : 590, 3.9% 
## (mean count < 1)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```




### Running DESeq
Once we have our `DESeqDataSet` setup we can run differential expression analtsis with *DESeq*:


```r
#dds <- DESeq(dds)
```

A new `DESeqDataSet` is returned that contains all the results (and the input `DESeqDataSet` parameters) within it.  

Extracting out the tables of results we actually want from this object can be a bit tricky. The next section describes one common way to do this. 

### Extracting our results table
Calling the *DESeq* packages **results()** function on our `DESeqDataSet` without any arguments will extract the estimated log2 fold changes and p values like so:


```r
#res <- results(dds)
res
```

```
## log2 fold change (MLE): condition hoxa1 kd vs control sirna 
## Wald test p-value: condition hoxa1 kd vs control sirna 
## DataFrame with 15280 rows and 6 columns
##                    baseMean log2FoldChange      lfcSE          stat
##                   <numeric>      <numeric>  <numeric>     <numeric>
## ENSG00000117519    4483.627      -2.422719 0.06001850     -40.36620
## ENSG00000183508    2053.881       3.201955 0.07241968      44.21388
## ENSG00000159176    5692.463      -2.313737 0.05757255     -40.18820
## ENSG00000150938    7442.986      -2.059631 0.05386627     -38.23601
## ENSG00000116016    4423.947      -1.888019 0.04318301     -43.72134
## ...                     ...            ...        ...           ...
## ENSG00000166507  13.0403308  -5.824072e-04 0.49235436 -1.182903e-03
## ENSG00000175279  15.1411485  -4.751747e-04 0.45034963 -1.055124e-03
## ENSG00000133997 199.4222834   1.866463e-05 0.13315341  1.401739e-04
## ENSG00000140600   0.6821063  -1.068898e-04 2.23351213 -4.785728e-05
## ENSG00000160948 788.3419630  -2.521036e-06 0.06945116 -3.629941e-05
##                    pvalue      padj
##                 <numeric> <numeric>
## ENSG00000117519         0         0
## ENSG00000183508         0         0
## ENSG00000159176         0         0
## ENSG00000150938         0         0
## ENSG00000116016         0         0
## ...                   ...       ...
## ENSG00000166507 0.9990562 0.9993283
## ENSG00000175279 0.9991581 0.9993622
## ENSG00000133997 0.9998882 0.9999710
## ENSG00000140600 0.9999618 0.9999710
## ENSG00000160948 0.9999710 0.9999710
```

The returned `res` object is not a standard R data.frame but one that carries extra meatadata on the meaning of the columns:


```r
mcols(res, use.names = TRUE)
```

```
## DataFrame with 6 rows and 2 columns
##                        type
##                 <character>
## baseMean       intermediate
## log2FoldChange      results
## lfcSE               results
## stat                results
## pvalue              results
## padj                results
##                                                                description
##                                                                <character>
## baseMean                         mean of normalized counts for all samples
## log2FoldChange log2 fold change (MLE): condition hoxa1 kd vs control sirna
## lfcSE                  standard error: condition hoxa1 kd vs control sirna
## stat                   Wald statistic: condition hoxa1 kd vs control sirna
## pvalue              Wald test p-value: condition hoxa1 kd vs control sirna
## padj                                                  BH adjusted p-values
```

The column `log2FoldChange` is the effect size estimate. It tells us how much the gene's expression seems to have changed due to treatment with dexamethasone in comparison to untreated samples.  This value is reported on a logarithmic scale to base 2: for example, a log2 fold change of 1.5 means that the gene's expression is increased by a multiplicative factor of \(2^{1.5} \approx 2.82\).

*DESeq2* performs for each gene a *hypothesis test* to see whether evidence is sufficient to decide against the *null hypothesis* that there is zero effect of the treatment on the gene and that the observed difference between treatment and
control was merely caused by experimental variability (i.e., the type of variability that you can expect between different
samples in the same treatment group). As usual in statistics, the result of this test is reported as a *p* value, and it is found in the column `pvalue`. Remember that a *p* value indicates the probability that a fold change as strong as the observed one, or even stronger, would be seen under the situation described by the null hypothesis.  

We can also summarize the results with the *DESeq2* specific version of the **summary()** function. This will report some additional useful information: 


```r
summary(res)
```

```
## 
## out of 15280 with nonzero total read count
## adjusted p-value < 0.1
## LFC > 0 (up)     : 4352, 28% 
## LFC < 0 (down)   : 4400, 29% 
## outliers [1]     : 0, 0% 
## low counts [2]   : 590, 3.9% 
## (mean count < 1)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

Note that there are many many genes with differential expression reported above. Let's therefore be more strict about which set of genes are considered 'significant'. There are two main ways we can do this:

* Lower the false discovery rate threshold (i.e. the threshold on the adjusted p-value (`padj`) in the results table)
* Raise the log2 fold change threshold from 0 to a higher value.



> **Q.** In the summary of our results printed above (and by default) the FDR level is set to 10% (i.e. adjusted p-value < 0.1) and the log2 fold change threshold is set to 0. Use the `alpha` and `lfcThreshold` input arguments to the **results()** function to change these to an FDR of 5% and a log2 fold change threshold of 2. Then use the **summary()** function to find out how many genes are up and down at these thresholds.  




```
## 
## out of 15280 with nonzero total read count
## adjusted p-value < 0.1
## LFC > 0 (up)     : 4352, 28% 
## LFC < 0 (down)   : 4400, 29% 
## outliers [1]     : 0, 0% 
## low counts [2]   : 590, 3.9% 
## (mean count < 1)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

You could also use the ever useful **table()** function on your output of the **results()** function like so:


```r
table(resSig$padj < 0.05)
```

```
## 
## FALSE  TRUE 
## 13565   233
```


```r
table(resSig$log2FoldChange > 2)
```

```
## 
## FALSE  TRUE 
## 14723   557
```

Then combining to determine the number of genes that meet both the *p* value and log2 fold change thresholds (*UP* genes: 99; and *DOWN* genes: 134): 


```r
table( resSig$padj < 0.05, resSig$log2FoldChange > 2)
```

```
##        
##         FALSE  TRUE
##   FALSE 13292   273
##   TRUE    134    99
```


> **Side-Note:** In high-throughput biology, we are careful to not use the *p* values directly as evidence against the null, but to correct for **multiple testing**.  
>
> What would happen if we were to simply threshold the *p* values at a low value, say 0.05? There are 8594 genes with a *p* value below 0.05 among the 15280 genes for which the test succeeded in reporting a *p* value:


```r
table(res$pvalue < 0.05)
```

```
## 
## FALSE  TRUE 
##  6686  8594
```

> *DESeq2* uses the Benjamini-Hochberg (BH) adjustment as implemented in the base R **p.adjust()** function; in brief, this method calculates for each gene an adjusted *p* value that answers the following question: if one called significant all genes with an adjusted *p* value less than or equal to this gene's adjusted *p* value threshold, what would be the fraction of false positives (the *false discovery rate*, FDR) among them, in the sense of the calculation outlined above? These values, called the BH-adjusted *p* values, are given in the column `padj` of the `res` object.  
>
> The FDR is a useful statistic for many high-throughput experiments, as we are often interested in reporting or focusing on a set of interesting genes, and we would like to put an upper bound on the percent of false positives in this set.  
>
> Hence, if we consider a fraction of 5% false positives acceptable, we can consider all genes with an adjusted *p* value below 5% = 0.05 as significant. How many such genes are there?  


```r
table(res$padj < 0.05)
```

```
## 
## FALSE  TRUE 
##  6539  8151
```


We can now subset the results table to extract the genes with adjusted *p* value less than 0.05 and then sort them by their log2 fold change estimate to get the significant genes with the strongest down-regulation:


```r
# Make a new results object 'resSig' with only significant genes
resSig <- subset(res, padj < 0.05)

# Print the first 10 strongest DOWN genes
ord.down <- order(resSig$log2FoldChange)
head(resSig[ ord.down, ], 10)
```

```
## log2 fold change (MLE): condition hoxa1 kd vs control sirna 
## Wald test p-value: condition hoxa1 kd vs control sirna 
## DataFrame with 10 rows and 6 columns
##                  baseMean log2FoldChange     lfcSE      stat       pvalue
##                 <numeric>      <numeric> <numeric> <numeric>    <numeric>
## ENSG00000171587  2.561658      -4.902903  1.514295 -3.237746 1.204782e-03
## ENSG00000155011  7.971242      -4.601385  1.147782 -4.008936 6.099288e-05
## ENSG00000101306  1.888372      -4.425978  1.644059 -2.692105 7.100260e-03
## ENSG00000179855  3.760403      -4.424858  1.400644 -3.159161 1.582241e-03
## ENSG00000188581  1.795033      -4.405196  1.637739 -2.689804 7.149406e-03
## ENSG00000164342  3.517655      -4.323085  1.415825 -3.053403 2.262620e-03
## ENSG00000117598  1.680397      -4.261884  1.686659 -2.526821 1.151003e-02
## ENSG00000180332  9.377022      -4.258232  0.966178 -4.407296 1.046691e-05
## ENSG00000136943  3.367218      -4.250652  1.454036 -2.923347 3.462904e-03
## ENSG00000204183  1.542573      -4.168844  1.679468 -2.482241 1.305591e-02
##                         padj
##                    <numeric>
## ENSG00000171587 0.0027291049
## ENSG00000155011 0.0001615844
## ENSG00000101306 0.0142919728
## ENSG00000179855 0.0035179536
## ENSG00000188581 0.0143825518
## ENSG00000164342 0.0049307055
## ENSG00000117598 0.0224039090
## ENSG00000180332 0.0000301340
## ENSG00000136943 0.0073542374
## ENSG00000204183 0.0251331758
```

> **Q.** Do the same as above but print out the top 10 strongest up-regulated genes. HINT: see the help for the **order()** function to see how to return the decreasing ordered indices you will want for accesing your `resSig` result.  


```r
# Print the first 10 strongest YP genes
ord.up <- order(resSig$log2FoldChange, decreasing = TRUE)
head(resSig[ ord.up, ], 10)
```

```
## log2 fold change (MLE): condition hoxa1 kd vs control sirna 
## Wald test p-value: condition hoxa1 kd vs control sirna 
## DataFrame with 10 rows and 6 columns
##                    baseMean log2FoldChange     lfcSE      stat
##                   <numeric>      <numeric> <numeric> <numeric>
## ENSG00000128052  158.249550       8.822086 1.0347214  8.526049
## ENSG00000141668   43.941562       8.820638 1.2040490  7.325813
## ENSG00000004799   97.260234       8.117454 1.0373919  7.824868
## ENSG00000109321 3907.190542       7.479166 0.2880675 25.963240
## ENSG00000162892   25.278929       7.052106 1.2210031  5.775666
## ENSG00000073756 1357.627066       6.337747 0.3399843 18.641296
## ENSG00000239590    7.399064       6.250557 1.2943354  4.829164
## ENSG00000151322    7.309569       6.237074 1.2988709  4.801920
## ENSG00000135363    6.035727       5.955831 1.3350911  4.460992
## ENSG00000175868    5.981592       5.945427 1.3207432  4.501577
##                        pvalue          padj
##                     <numeric>     <numeric>
## ENSG00000128052  1.514306e-17  9.338857e-17
## ENSG00000141668  2.374553e-13  1.195414e-12
## ENSG00000004799  5.081914e-15  2.783494e-14
## ENSG00000109321 1.288826e-148 1.147446e-146
## ENSG00000162892  7.664930e-09  2.907251e-08
## ENSG00000073756  1.485963e-77  4.409858e-76
## ENSG00000239590  1.371077e-06  4.320274e-06
## ENSG00000151322  1.571514e-06  4.919145e-06
## ENSG00000135363  8.158097e-06  2.375470e-05
## ENSG00000175868  6.745119e-06  1.982349e-05
```




## Annotating our genes and mapping to Entrez IDs


```r
#res = res[order(res$pvalue),]
#summary(res)

#resSig = resSig[order(resSig$pvalue),]
#summary(resSig)

#res = (resSig[order(resSig$pvalue),])
```

Since we mapped and counted against the Ensembl annotation, our results only have information about Ensembl gene IDs. However, our pathway analysis downstream will use KEGG pathways, and genes in KEGG pathways are annotated with Entrez gene IDs. So lets add them as we did the last day.


```r
#biocLite("AnnotationDbi")
#biocLite("org.Hs.eg.db")
```


```r
library("AnnotationDbi")
```

```
## Warning: package 'AnnotationDbi' was built under R version 3.4.2
```

```r
library("org.Hs.eg.db")
```

```
## 
```

```r
columns(org.Hs.eg.db)
```

```
##  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT" 
##  [5] "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"    
##  [9] "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"       
## [13] "IPI"          "MAP"          "OMIM"         "ONTOLOGY"    
## [17] "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
## [21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
## [25] "UNIGENE"      "UNIPROT"
```


```r
res$symbol = mapIds(org.Hs.eg.db,
                    keys=row.names(res), 
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")
```

```
## 'select()' returned 1:many mapping between keys and columns
```

```r
res$entrez = mapIds(org.Hs.eg.db,
                    keys=row.names(res), 
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first")
```

```
## 'select()' returned 1:many mapping between keys and columns
```

```r
res$name =   mapIds(org.Hs.eg.db,
                    keys=row.names(res), 
                    column="GENENAME",
                    keytype="ENSEMBL",
                    multiVals="first")
```

```
## 'select()' returned 1:many mapping between keys and columns
```

```r
head(res, 10)
```

```
## log2 fold change (MLE): condition hoxa1 kd vs control sirna 
## Wald test p-value: condition hoxa1 kd vs control sirna 
## DataFrame with 10 rows and 9 columns
##                   baseMean log2FoldChange      lfcSE        stat
##                  <numeric>      <numeric>  <numeric>   <numeric>
## ENSG00000279457   29.91358     0.17927483 0.32459294   0.5523066
## ENSG00000187634  183.22965     0.42644724 0.14017817   3.0421802
## ENSG00000188976 1651.18808    -0.69272061 0.05484412 -12.6307173
## ENSG00000187961  209.63794     0.72975918 0.13178350   5.5375609
## ENSG00000187583   47.25512     0.04055411 0.27169055   0.1492658
## ENSG00000187642   11.97975     0.54275443 0.52117592   1.0414035
## ENSG00000188290  108.92213     2.05704993 0.19679801  10.4525953
## ENSG00000187608  350.71687     0.25737707 0.10267368   2.5067482
## ENSG00000188157 9128.43942     0.38990881 0.04673847   8.3423531
## ENSG00000131591  156.47908     0.19659758 0.14551980   1.3510022
##                       pvalue         padj       symbol      entrez
##                    <numeric>    <numeric>  <character> <character>
## ENSG00000279457 5.807383e-01 6.846746e-01 LOC102723897   102723897
## ENSG00000187634 2.348712e-03 5.109223e-03       SAMD11      148398
## ENSG00000188976 1.429690e-36 1.745815e-35        NOC2L       26155
## ENSG00000187961 3.067131e-08 1.109758e-07       KLHL17      339451
## ENSG00000187583 8.813439e-01 9.191354e-01      PLEKHN1       84069
## ENSG00000187642 2.976883e-01 4.016657e-01        PERM1       84808
## ENSG00000188290 1.425695e-25 1.226198e-24         HES4       57801
## ENSG00000187608 1.218475e-02 2.358287e-02        ISG15        9636
## ENSG00000188157 7.282663e-17 4.347108e-16         AGRN      375790
## ENSG00000131591 1.766947e-01 2.604501e-01     C1orf159       54991
##                                                                     name
##                                                              <character>
## ENSG00000279457                        WAS protein family homolog 2-like
## ENSG00000187634                 sterile alpha motif domain containing 11
## ENSG00000188976 NOC2 like nucleolar associated transcriptional repressor
## ENSG00000187961                              kelch like family member 17
## ENSG00000187583                 pleckstrin homology domain containing N1
## ENSG00000187642             PPARGC1 and ESRR induced regulator, muscle 1
## ENSG00000188290                   hes family bHLH transcription factor 4
## ENSG00000187608                            ISG15 ubiquitin-like modifier
## ENSG00000188157                                                    agrin
## ENSG00000131591                      chromosome 1 open reading frame 159
```

Great, this is looking good so far. Now lets see how pathway analysis can help us make further sense out of this ranked list of differentially expressed genes.


## Pathway Analysis

Here we are going to use the gage package for pathway analysis. Once we have a list of enriched pathways, we’re going to use the pathview package to draw pathway diagrams, shading the molecules in the pathway by their degree of up/down-regulation.

### KEGG


```r
#source("http://bioconductor.org/biocLite.R")
#biocLite( c("pathview", "gage", "gageData") )
```



```r
library(pathview)
```

```
## Warning: package 'pathview' was built under R version 3.4.3
```

```
## ##############################################################################
## Pathview is an open source software package distributed under GNU General
## Public License version 3 (GPLv3). Details of GPLv3 is available at
## http://www.gnu.org/licenses/gpl-3.0.html. Particullary, users are required to
## formally cite the original Pathview paper (not just mention it) in publications
## or products. For details, do citation("pathview") within R.
## 
## The pathview downloads and uses KEGG data. Non-academic uses may require a KEGG
## license agreement (details at http://www.kegg.jp/kegg/legal.html).
## ##############################################################################
```

```r
library(gage)
```

```
## Warning: package 'gage' was built under R version 3.4.3
```

```r
library(gageData)
```



```r
data(kegg.sets.hs)
data(sigmet.idx.hs)

kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs, 3)
```

```
## $`hsa00232 Caffeine metabolism`
## [1] "10"   "1544" "1548" "1549" "1553" "7498" "9"   
## 
## $`hsa00983 Drug metabolism - other enzymes`
##  [1] "10"     "1066"   "10720"  "10941"  "151531" "1548"   "1549"  
##  [8] "1551"   "1553"   "1576"   "1577"   "1806"   "1807"   "1890"  
## [15] "221223" "2990"   "3251"   "3614"   "3615"   "3704"   "51733" 
## [22] "54490"  "54575"  "54576"  "54577"  "54578"  "54579"  "54600" 
## [29] "54657"  "54658"  "54659"  "54963"  "574537" "64816"  "7083"  
## [36] "7084"   "7172"   "7363"   "7364"   "7365"   "7366"   "7367"  
## [43] "7371"   "7372"   "7378"   "7498"   "79799"  "83549"  "8824"  
## [50] "8833"   "9"      "978"   
## 
## $`hsa00230 Purine metabolism`
##   [1] "100"    "10201"  "10606"  "10621"  "10622"  "10623"  "107"   
##   [8] "10714"  "108"    "10846"  "109"    "111"    "11128"  "11164" 
##  [15] "112"    "113"    "114"    "115"    "122481" "122622" "124583"
##  [22] "132"    "158"    "159"    "1633"   "171568" "1716"   "196883"
##  [29] "203"    "204"    "205"    "221823" "2272"   "22978"  "23649" 
##  [36] "246721" "25885"  "2618"   "26289"  "270"    "271"    "27115" 
##  [43] "272"    "2766"   "2977"   "2982"   "2983"   "2984"   "2986"  
##  [50] "2987"   "29922"  "3000"   "30833"  "30834"  "318"    "3251"  
##  [57] "353"    "3614"   "3615"   "3704"   "377841" "471"    "4830"  
##  [64] "4831"   "4832"   "4833"   "4860"   "4881"   "4882"   "4907"  
##  [71] "50484"  "50940"  "51082"  "51251"  "51292"  "5136"   "5137"  
##  [78] "5138"   "5139"   "5140"   "5141"   "5142"   "5143"   "5144"  
##  [85] "5145"   "5146"   "5147"   "5148"   "5149"   "5150"   "5151"  
##  [92] "5152"   "5153"   "5158"   "5167"   "5169"   "51728"  "5198"  
##  [99] "5236"   "5313"   "5315"   "53343"  "54107"  "5422"   "5424"  
## [106] "5425"   "5426"   "5427"   "5430"   "5431"   "5432"   "5433"  
## [113] "5434"   "5435"   "5436"   "5437"   "5438"   "5439"   "5440"  
## [120] "5441"   "5471"   "548644" "55276"  "5557"   "5558"   "55703" 
## [127] "55811"  "55821"  "5631"   "5634"   "56655"  "56953"  "56985" 
## [134] "57804"  "58497"  "6240"   "6241"   "64425"  "646625" "654364"
## [141] "661"    "7498"   "8382"   "84172"  "84265"  "84284"  "84618" 
## [148] "8622"   "8654"   "87178"  "8833"   "9060"   "9061"   "93034" 
## [155] "953"    "9533"   "954"    "955"    "956"    "957"    "9583"  
## [162] "9615"
```

The main `gage()` function requires a named vector of fold changes, where the names of the values are the Entrez gene IDs.


```r
#foldchanges = resSig$log2FoldChange
#names(foldchanges) = resSig$entrez
#head(foldchanges)
```


```r
foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
head(foldchanges)
```

```
##   102723897      148398       26155      339451       84069       84808 
##  0.17927483  0.42644724 -0.69272061  0.72975918  0.04055411  0.54275443
```

Run the pathway analysis...


```r
# Get the results
keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
```

Look at results


```r
attributes(keggres)
```

```
## $names
## [1] "greater" "less"    "stats"
```


```r
head(keggres$less)
```

```
##                                      p.geomean stat.mean        p.val
## hsa04110 Cell cycle               1.004024e-05 -4.353447 1.004024e-05
## hsa03030 DNA replication          8.909718e-05 -3.968605 8.909718e-05
## hsa03013 RNA transport            1.333746e-03 -3.038773 1.333746e-03
## hsa04114 Oocyte meiosis           2.496109e-03 -2.840329 2.496109e-03
## hsa03440 Homologous recombination 2.942017e-03 -2.868137 2.942017e-03
## hsa00240 Pyrimidine metabolism    5.800212e-03 -2.549616 5.800212e-03
##                                         q.val set.size         exp1
## hsa04110 Cell cycle               0.001606438      120 1.004024e-05
## hsa03030 DNA replication          0.007127774       36 8.909718e-05
## hsa03013 RNA transport            0.071133118      143 1.333746e-03
## hsa04114 Oocyte meiosis           0.094144560       98 2.496109e-03
## hsa03440 Homologous recombination 0.094144560       28 2.942017e-03
## hsa00240 Pyrimidine metabolism    0.138500584       95 5.800212e-03
```

Each keggres$greater and keggres$less object is data matrix with gene sets as rows sorted by p-value. Lets look at both up (greater), down (less), and statistics by calling head() with the lapply() function. As always if you want to find out more about a particular function or its return values use the R help system (e.g. ?gage or ?lapply).


```r
lapply(keggres, head)
```

```
## $greater
##                                         p.geomean stat.mean       p.val
## hsa04640 Hematopoietic cell lineage   0.002709366  2.857393 0.002709366
## hsa04630 Jak-STAT signaling pathway   0.005655916  2.557207 0.005655916
## hsa04142 Lysosome                     0.008948808  2.384783 0.008948808
## hsa00140 Steroid hormone biosynthesis 0.009619717  2.432105 0.009619717
## hsa04740 Olfactory transduction       0.010846448  2.357877 0.010846448
## hsa04916 Melanogenesis                0.017613266  2.123151 0.017613266
##                                           q.val set.size        exp1
## hsa04640 Hematopoietic cell lineage   0.3470863       49 0.002709366
## hsa04630 Jak-STAT signaling pathway   0.3470863      103 0.005655916
## hsa04142 Lysosome                     0.3470863      117 0.008948808
## hsa00140 Steroid hormone biosynthesis 0.3470863       26 0.009619717
## hsa04740 Olfactory transduction       0.3470863       39 0.010846448
## hsa04916 Melanogenesis                0.4696871       85 0.017613266
## 
## $less
##                                      p.geomean stat.mean        p.val
## hsa04110 Cell cycle               1.004024e-05 -4.353447 1.004024e-05
## hsa03030 DNA replication          8.909718e-05 -3.968605 8.909718e-05
## hsa03013 RNA transport            1.333746e-03 -3.038773 1.333746e-03
## hsa04114 Oocyte meiosis           2.496109e-03 -2.840329 2.496109e-03
## hsa03440 Homologous recombination 2.942017e-03 -2.868137 2.942017e-03
## hsa00240 Pyrimidine metabolism    5.800212e-03 -2.549616 5.800212e-03
##                                         q.val set.size         exp1
## hsa04110 Cell cycle               0.001606438      120 1.004024e-05
## hsa03030 DNA replication          0.007127774       36 8.909718e-05
## hsa03013 RNA transport            0.071133118      143 1.333746e-03
## hsa04114 Oocyte meiosis           0.094144560       98 2.496109e-03
## hsa03440 Homologous recombination 0.094144560       28 2.942017e-03
## hsa00240 Pyrimidine metabolism    0.138500584       95 5.800212e-03
## 
## $stats
##                                       stat.mean     exp1
## hsa04640 Hematopoietic cell lineage    2.857393 2.857393
## hsa04630 Jak-STAT signaling pathway    2.557207 2.557207
## hsa04142 Lysosome                      2.384783 2.384783
## hsa00140 Steroid hormone biosynthesis  2.432105 2.432105
## hsa04740 Olfactory transduction        2.357877 2.357877
## hsa04916 Melanogenesis                 2.123151 2.123151
```

Now, let’s try out the pathview() function from the pathview package to make a pathway plot with our result shown in color. To begin with lets manually supply a pathway.id (namely the first part of the "hsa04110 Cell cycle") that we could see from the print out above.


```r
pathview(gene.data=foldchanges, pathway.id="hsa04110")
```

```
## 'select()' returned 1:1 mapping between keys and columns
```

```
## Info: Working in directory /Users/barry/Desktop/courses/bggn213_S18/bggn213_github/class15
```

```
## Info: Writing image file hsa04110.pathview.png
```

Here is the default low resolution raster PNG output from the first pathview() call above:



![](hsa04110.pathview.png)



Note how many of the genes in this pathway are pertubed (i.e. colored) in our results.

Now, let's process our results a bit more to automagicaly pull out the top 5 upregulated pathways, then further process that just to get the IDs needed by the **pathview()** function. We'll use these KEGG pathway IDs for plotting below.


```r
## Focus on top 5 upregulated pathways here for demo purposes only
keggrespathways <- rownames(keggres$greater)[1:5]

# Extract the IDs part of each string
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids
```

```
## [1] "hsa04640" "hsa04630" "hsa04142" "hsa00140" "hsa04740"
```

Finally, lets pass these IDs in `keggresids` to the **pathview()** function to draw plots for all the top 5 pathways.


```r
pathview(gene.data=foldchanges, pathway.id=keggresids, species="hsa")
```

Here are the plots:

![](hsa00140.pathview.png)

![](hsa04142.pathview.png)

![](hsa04630.pathview.png)

![](hsa04640.pathview.png)

![](hsa04740.pathview.png)



# Section 3. Gene Ontology (GO)

We can also do a similar procedure with gene ontology. Similar to above, go.sets.hs has all GO terms. go.subs.hs is a named list containing indexes for the BP, CC, and MF ontologies. Let’s only do Biological Process.


```r
data(go.sets.hs)
data(go.subs.hs)
gobpsets = go.sets.hs[go.subs.hs$BP]

gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)

head(gobpres$less, 10)
```

```
##                                                     p.geomean stat.mean
## GO:0000279 M phase                               1.475275e-16 -8.323757
## GO:0048285 organelle fission                     7.498096e-16 -8.160310
## GO:0000280 nuclear division                      2.135081e-15 -8.034815
## GO:0007067 mitosis                               2.135081e-15 -8.034815
## GO:0000087 M phase of mitotic cell cycle         5.927745e-15 -7.891754
## GO:0007059 chromosome segregation                1.055849e-11 -6.988384
## GO:0051301 cell division                         1.088128e-10 -6.424636
## GO:0000236 mitotic prometaphase                  1.535091e-10 -6.717799
## GO:0000226 microtubule cytoskeleton organization 1.097649e-09 -6.088264
## GO:0007017 microtubule-based process             3.605409e-09 -5.849001
##                                                         p.val        q.val
## GO:0000279 M phase                               1.475275e-16 5.867170e-13
## GO:0048285 organelle fission                     7.498096e-16 1.490996e-12
## GO:0000280 nuclear division                      2.135081e-15 2.122805e-12
## GO:0007067 mitosis                               2.135081e-15 2.122805e-12
## GO:0000087 M phase of mitotic cell cycle         5.927745e-15 4.714928e-12
## GO:0007059 chromosome segregation                1.055849e-11 6.998521e-09
## GO:0051301 cell division                         1.088128e-10 6.182122e-08
## GO:0000236 mitotic prometaphase                  1.535091e-10 7.631319e-08
## GO:0000226 microtubule cytoskeleton organization 1.097649e-09 4.850388e-07
## GO:0007017 microtubule-based process             3.605409e-09 1.433871e-06
##                                                  set.size         exp1
## GO:0000279 M phase                                    492 1.475275e-16
## GO:0048285 organelle fission                          373 7.498096e-16
## GO:0000280 nuclear division                           349 2.135081e-15
## GO:0007067 mitosis                                    349 2.135081e-15
## GO:0000087 M phase of mitotic cell cycle              359 5.927745e-15
## GO:0007059 chromosome segregation                     141 1.055849e-11
## GO:0051301 cell division                              452 1.088128e-10
## GO:0000236 mitotic prometaphase                        84 1.535091e-10
## GO:0000226 microtubule cytoskeleton organization      269 1.097649e-09
## GO:0007017 microtubule-based process                  407 3.605409e-09
```

```r
#lapply(gobpres, head)
```

# Reactome Pathway Analysis Online

Reactome, such as many other tools, has an online software available (https://reactome.org/) and R package 

Write out our list for use online


```r
sig_genes <- res[res$padj <= 0.05 & !is.na(res$padj), "symbol"]
print(paste("Total number of significant genes:", length(sig_genes)))
```

```
## [1] "Total number of significant genes: 8151"
```

```r
write.table(sig_genes, file="significant_genes.txt", 
            row.names=FALSE, col.names=FALSE, quote=FALSE)
```



## Bounus: Gene clustering, heatmaps and PCA

Many statistical methods for analysis of multidimensional data, for example *clustering* and *principal components analysis* (PCA), work best for data that generally has the same range of variance at different ranges of the mean values.  

However, for counts from RNA-seq the expected variance grows with the mean. For example, if one performs PCA or clustering directly on a matrix of counts then the results will be heavely influenced by the genes with the highest counts (because they show the largest absolute differences between samples).  To address this problem the *DESeq2* package offers the **vst()** function (that stands for Variance Stabilizing Transformation) that stabilizes the variance of count data across different mean values. We will use **vst()** here as input for our clustering.


```r
vsd <- vst(dds, blind = FALSE)
```

Since gene clustering is only really relevant for genes that actually carry a signal, one usually would only cluster a subset of the most highly variable genes. Here, for demonstration purposes we select the 20 genes with the highest variance across samples.



```r
library("genefilter")
```

```
## Warning: package 'genefilter' was built under R version 3.4.2
```

```
## 
## Attaching package: 'genefilter'
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     rowSds, rowVars
```

```r
#row.variance <- apply(assay(vsd), 1, var)
row.variance <- rowVars(assay(vsd))
ord.variance <- order( row.variance, decreasing = TRUE) 

# Focus on top 20 most variable genes for demo purposes
mat  <- assay(vsd)[ ord.variance[1:20], ]
```

The heatmap becomes more interesting if we do not look at absolute expression strength but rather at the amount by which each gene deviates in a specific sample from the gene’s average across all samples. To do this we center each genes’ values across samples by subtracting their mean values, and then plot the heatmap (figure below). 


```r
library(pheatmap)
mat.center  <- mat - rowMeans(mat)
pheatmap(mat.center)
```

![](class15_files/figure-html/unnamed-chunk-45-1.png)<!-- -->


> **Side-note:** We can do PCA with our `vsd` object.


```r
pcaData <- plotPCA(vsd, intgroup="condition", returnData = TRUE)
pcaData
```

```
##                 PC1        PC2         group     condition      name
## SRR493366 -12.74667  0.3046597 control_sirna control_sirna SRR493366
## SRR493367 -13.04690 -0.0278821 control_sirna control_sirna SRR493367
## SRR493368 -13.28760 -0.3782467 control_sirna control_sirna SRR493368
## SRR493369  13.38297 -0.1131647      hoxa1_kd      hoxa1_kd SRR493369
## SRR493370  11.48330  1.1139575      hoxa1_kd      hoxa1_kd SRR493370
## SRR493371  14.21490 -0.8993238      hoxa1_kd      hoxa1_kd SRR493371
```


```r
library(ggplot2)

ggplot(pcaData, aes(x = PC1, y = PC2, color = condition) ) +
  geom_point(size =3) 
```

![](class15_files/figure-html/unnamed-chunk-46-1.png)<!-- -->


