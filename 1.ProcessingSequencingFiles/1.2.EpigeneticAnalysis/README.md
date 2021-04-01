# Epigenetic analysis with differentially-methylated regions from msGBS library

Killfish Project Pipeline |  by  Waldir M. Berbel-Filho [![Foo](../ORCID-iD.png)](https://orcid.org/0000-0001-6991-4685)

This documentation outlines the pipelines used for epigenetic analysis in the preprint manuscript [Berbel-Filho et al. (2020)](https://www.biorxiv.org/content/10.1101/2020.07.01.182022v1.full).

Last Modified: 21 July 2020

Please, contact waldirmbf@gmail.com should any question arise.
__________________________________________
 For the epigenetic analysis, given the evidence of hybridisation between _K. ocellatus_ and _K. hermaphroditus_ in only only two populatins in southeast Brazil (FUN and GUA), we  subset the samples containing only _K. ocellatus_ and _K. hermaphroditus_ from those populations. We used the R package **msgbsR v. 1.12.0** in R v. 3.50 to run the differential methylation analysis.

## 1. msGBS in R:
##### 1.1) Loading required packages
 ```
 > library(msgbsR)
 > library(GenomicRanges)
 > library(SummarizedExperiment)
```
##### 1.2) Loading input BAM files
```
> my_path <- system.file("extdata", package = "msgbsR")
> se <- rawCounts(bamFilepath = "/data/home/waldir/Desktop/msGBS_data/demultiplexed_samples/Sorted_Indexed_BAM/New_names/Filtered/Hybrids/", threads = 10)
```

##### 1.3) Setting group names and sample sizes

To set the group (either one of the parental species or hybrids), we used the following command:
```
> colData(se) <- DataFrame(Group = c(rep("Hyb", 5), rep("Kher", 21), rep("Koce", 11)),
> row.names = colnames(assay(se)))
```
 The five BAM files with the prefix "Hyb" represented the hybrid F1 individuals identified using the genetic analysis.

 The 21 BAM files with the prefix "Kher" represented pure _K. hermaphroditus_ individuals identified using the genetic analysis.

 The 11 BAM files with the prefix "Koce" represented pure _K. ocellatus_ individuals identified using the genetic analysis.

 ##### 1.4) Extracting cut sites
 To extract the cut sites from the methylation-sensitive enzyme used in the library, in this case HpaII (C*CGG).
 ```
> cutSites <- rowRanges(se)
> cutSites
```
 ##### 1.5) Adjusting the cut sites to overlap recognition site on each strand
```
> start(cutSites) <- ifelse(test = strand(cutSites) == '+',
                          yes = start(cutSites) - 1, no = start(cutSites) - 2)
> end(cutSites) <- ifelse(test = strand(cutSites) == '+',
                          yes = end(cutSites) + 2, no = end(cutSites) + 1)
```

 ##### 1.6) Running checkCuts to filter out incorrect cuts
```
> library(BSgenome)
> available.genomes()
```
As the reference genome of _K. marmoratus_ (NCBI:ASM164957v1) used during aligments was not among the ones available in BSgenome, we have forged the genome locally
into BSgenome using the instructions provided [here](https://www.bioconductor.org/packages/devel/bioc/vignettes/BSgenome/inst/doc/BSgenomeForge.pdf).

The forged genome data pakcage was locally called "BSgenome.KmarmoratusRefSeq.SU2"

Back to R, we have now filtered out the incorrect cut sites with the reference genome as follows:
```
> library(BSgenome.KmarmoratusRefSeq.SU2)
> correctCuts <- checkCuts(cutSites = cutSites, genome = "BSgenome.KmarmoratusRefSeq.SU2", seq = "CCGG")
> se <- subsetByOverlaps(se, correctCuts)
```

##### 1.8) Plotting counts per cut sites
To plot the number of cut sites vs number of reads per sample, we used the following command:

```
> pdf('NofReadsvsNofCutSites.pdf')
> plotCounts(se = se, cateogory = "Group")
> dev.off()
```

##### 1.9) Differential methylation analysis
The dataset is now ready for the differential methylation analysis between groups. We first compared samples of _K. ocellatus_ vs. _K. hermaphroditus_ as follows:
```
#Total DMCs
> KoceVSKher <- diffMeth(se = se, cateogory = "Group",
                condition1 = "Koce", condition2 = "Kher",
                cpmThreshold = 1, thresholdSamples = 11)
#Write file         
> write.csv(KoceVSKher,"diffMeth_KocevsKher_total.csv")

## Total number of DMCs = 56,705 between Koce and Kher.

#Filter DMCs with FDR < 0.01
> KoceVSKher_filtered <- GRanges(KoceVSKher$site[which(KoceVSKher$FDR < 0.01)])
> summary(KoceVSKher_filtered)
> KoceVSKher_filtered

#Write file
write.csv(KoceVSKher_filtered,"diffMeth_KocevsKher_filtered.csv")

## Total number of DMCs with  FDR < 0.01 = 38,473 between Koce and Kher.
```

Similarly, we then compared samples of _K. ocellatus_ vs.  F1 hybrids as follows:
```
#Total DMCs
> KoceVSHyb <- diffMeth(se = se, cateogory = "Group",
                condition1 = "Koce", condition2 = "Hyb",
                cpmThreshold = 1, thresholdSamples = 5)
#Write file         
> write.csv(KoceVSHyb,"diffMeth_KocevsHyb_total.csv")

## Total number of DMCs = 54,296 between Koce and Hyb.

#Filter DMCs with FDR < 0.01
> KoceVSHyb_filtered <- GRanges(KoceVSHyb$site[which(top$FDR < 0.01)])
> summary(KoceVSHyb_filtered)
> KoceVSHyb_filtered
#Write file
write.csv(KoceVSHyb_filtered,"diffMeth_KocevsHyb_filtered.csv")

## Total number of DMCs with FDR < 0.01 = 12,006 between Koce and Kher.

```
Finally, we compared samples of _K. hermaphroditus_ vs. hybrids as follows:
```
#Total DMCs
> KherVSHyb <- diffMeth(se = se, cateogory = "Group",
                condition1 = "Kher", condition2 = "Hyb",
                cpmThreshold = 1, thresholdSamples = 7)
#Write file         
> write.csv(KherVSHyb,"diffMeth_KhervsHyb_total.csv")

## Total number of DMCs = 60,212 between Koce and Hyb.

#Filter DMCs with FDR < 0.01
> KherVSHyb_filtered <- GRanges(KherVSHyb$site[which(KherVSHyb$FDR < 0.01)])
> summary(KherVSHyb_filtered)
> KherVSHyb_filtered
#Write file
> write.csv(KherVSHyb_filtered,"diffMeth_KhervsHyb_filtered.csv")

## Total number of DMCs with FDR < 0.01 = 12,221 between Kher and Hyb.
```

All output '.csv' files containing both total and filtered DMCs are attached in this repository.

## 2.  Coverage for DMCs In shell:

##### 2.1) Extracting coverage information for DMCs
For downstream analysis, we need to have the number of reads for each DMC involved in the comparisons between groups. For that we, first created two '.bed' files with the coordinates for the DMCs between parental species (one with the 38,473 DMMs between _K.ocellatus_ and _K. hermaphroditus_, named as ['diffMeth_KocevsHyb_filtered.bed'](diffMeth_KocevsHyb_filtered.bed), and another one containing the  5,800 DMCs shared in the comparisons between hybrids and both parental species, named as ['diffMeth_HybvsParentals.bed'](diffMeth_HybvsParentals.bed). For donwstream analysis, we got the number of reads across all individuals for the DMCs per positions in the '.bed' files using **samtools v.1.9**, as follows:

First, for the coverage of DMCs between _K.ocellatus_ and _K. hermaphroditus_:
```
samtools depth -b diffMeth_KocevsHyb_filtered.bed -f listBAM.txt > rawCounts_KocevsKher.csv
```
`-b` = Input '.bed' file.

`-f` = List of input BAM fies.

We manually edited this file, including the first two columns (EntrezID and position) and sample names for each column. The edited file was called ['rawCounts_KocevsKher.txt'](rawCounts_KocevsKher.txt).


For the coverage of DMCs between hybrids and both parental species:

```
samtools depth -b diffMeth_HybvsParentals.bed -f listBAM.txt > rawCounts_HybvsPar.csv
```

We manually edited this file, including the first two columns (EntrezID and position) and sample names for each column. The edited file was called ['rawCounts_HybvsPar.txt'](rawCounts_HybvsPar.txt).


For the coverage of all sites covered in the library:
```
samtools depth -a -f listBAM.txt > rawCounts_Allsites.csv
```
`-a` = Output all positions (including those with zero depth).

Output file was too big to be deposited in GitHub. Analysis on this file follow the same pipelines as described below but were locally.

## 3. Normalisation of reads in R:
To normalise the raw counts for downstream analysis, we followed a pipeline commonly used in RNAseq analysis with the R package **edgeR v. 3.11** in R v. 3.5. Following we will only show the file containing the DMCs between hybrids and parental species as a example of the pipeline also used for the normalisation of reads. The same script was also used for the DMCs between parental species and also for the normalisation of all sites. As follows:

````
#Load edgeR
> library(edgeR)

#Input file with samples information (i.e. SampleName, SamplingSite, Species)
> sampleinfo<-read.delim('SampleInfo.txt')

#Input raw counts data
> seqdata<-read.delim('rawCounts_HybvsPar.txt',stringsAsFactors=FALSE)

#Create file with only counts (remove EntrezID and Lenght columns)
> countdata<- seqdata[,-(1:2)]

#Obtain cpms(counts per million 1 in at least 7 samples)
> myCPM<-cpm(countdata)

#Filter out cpm values >0.5
> thresh<-myCPM > 0.5
> keep <- rowSums(thresh) >=2
> counts.keep <- countdata[keep,]

#Convert counts data in a DGELIst object
> dgeObj  <- DGEList(counts.keep)

#Get log2 counts per million
> logcounts <- cpm(dgeObj, log=TRUE)

#Boxplot to check overall distrIbution of log2 counts files
> pdf('log2counts.pdf')
> boxplot(logcounts, xlab="", ylab="Log2 counts per million", las=2) + abline(h=median(logcounts),col="blue")
> dev.off()


#calculate normalisation factors for each sample
> dgeObj <- calcNormFactors(dgeObj)
> write.txt(dgeObj$Samples"normalisationFactors.txt")

#Get normalised counts
> norm<-cpm(dgeObj, normalized.lib.sizes=FALSE)
````

'norm' is the variable with the normalised counts that were used in the overall variation (MDS) and hierarchical clustering analysis.

### Citation

Berbel-Filho, Waldir M., et al. "Additive and non-additive epigenetic signatures of hybridisation between fish species with different mating systems." bioRxiv (2020). https://doi.org/10.1101/2020.07.01.182022
