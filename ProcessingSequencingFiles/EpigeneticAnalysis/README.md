# Epigenetic analysis with differentially-methylated regions from msGBS library

Killfish Project Pipeline |  by  Waldir M. Berbel-Filho [![Foo](../ORCID-iD.png)](https://orcid.org/0000-0001-6991-4685)

This documentation outlines the pipelines used for epigenetic analysis in Berbel-Filho et al. (XXXX)

Last Modified: 17 July 2020

Please, contact waldirmbf@gmail.com should any question arise.
__________________________________________
 For the epigenetic analysis, given the evidence of hybridisation between _K. ocellatus_ and _K. hermaphroditus_ in only only two populatins in southeast Brazil (FUN and GUA), we  subset the samples containing only _K. ocellatus_ and _K. hermaphroditus_ from those populations. We used the R package **msgbsR v. 1.12.0** in R v. 3.50 to run the differential methylation analysis.

## 1. In R:
##### 1.1) Loading required packages
 ```
 library(msgbsR)
 library(GenomicRanges)
 library(SummarizedExperiment)
```
##### 1.2) Loading input BAM files
```
my_path <- system.file("extdata", package = "msgbsR")
se <- rawCounts(bamFilepath = "/data/home/waldir/Desktop/msGBS_data/demultiplexed_samples/Sorted_Indexed_BAM/New_names/Filtered/Hybrids/", threads = 10)
```

##### 1.3) Setting group names and sample sizes

To set the group (either one of the parental species or hybrids), we used the following command:
```
colData(se) <- DataFrame(Group = c(rep("Hyb", 7), rep("Kher", 21), rep("Koce", 11)),
row.names = colnames(assay(se)))
```
 The seven BAM files with the prefix "Hyb" represented the hybrid individuals identified using the genetic analysis.

 The 21 BAM files with the prefix "Kher" represented pure _K. hermaphroditus_ individuals identified using the genetic analysis.

  The 11 BAM files with the prefix "Koce" represented pure _K. ocellatus_ individuals identified using the genetic analysis.

 ##### 1.4) Extracting cut sites
 To extract the cut sites from the methylation-sensitive enzyme used in the library, in this case HpaII (C*CGG).
 ```
cutSites <- rowRanges(se)
cutSites
```
 ##### 1.5) Adjusting the cut sites to overlap recognition site on each strand
```
start(cutSites) <- ifelse(test = strand(cutSites) == '+',
                          yes = start(cutSites) - 1, no = start(cutSites) - 2)
end(cutSites) <- ifelse(test = strand(cutSites) == '+',
                          yes = end(cutSites) + 2, no = end(cutSites) + 1)
```

 ##### 1.6) Running checkCuts to filter out incorrect cuts
```
library(BSgenome)
available.genomes()
```
As the reference genome of _K. marmoratus_ (NCBI:ASM164957v1) used during aligments was not among the ones available in BSgenome, we have forged the genome locally
into BSgenome using the instructions provided [here](https://www.bioconductor.org/packages/devel/bioc/vignettes/BSgenome/inst/doc/BSgenomeForge.pdf).

The forged genome data pakcage was locally called "BSgenome.KmarmoratusRefSeq.SU2"

Back to R, we have now filtered out the incorrect cut sites with the reference genome as follows:
```
library(BSgenome.KmarmoratusRefSeq.SU2)
correctCuts <- checkCuts(cutSites = cutSites, genome = "BSgenome.KmarmoratusRefSeq.SU2", seq = "CCGG")
se <- subsetByOverlaps(se, correctCuts)
```

##### 1.8) Differential methylation analysis
The dataset is now ready for the differential methylation analysis between groups. We first compared samples of _K. ocellatus_ vs. _K. hermaphroditus_ as follows:
```
#Total DMCs
KoceVSKher <- diffMeth(se = se, cateogory = "Group",
                condition1 = "Koce", condition2 = "Kher",
                cpmThreshold = 1, thresholdSamples = 11)
#Write file         
write.csv(KoceVSKher,"diffMeth_KocevsKher_total.csv")

## Total number of DMCs = 56,705 between Koce and Kher.

#Filter DMCs with FDR < 0.05
KoceVSKher_filtered <- GRanges(KoceVSKher$site[which(top$FDR < 0.05)])
summary(KoceVSKher_filtered)
KoceVSKher_filtered

#Write file
write.csv(KoceVSKher_filtered,"diffMeth_KocevsKher_filtered.csv")

## Total number of DMCs with < 0.05 = 37,664 between Koce and Kher.
```

Similarly, we then compared samples of _K. ocellatus_ vs. hybrids as follows:
```
#Total DMCs
KoceVSHyb <- diffMeth(se = se, cateogory = "Group",
                condition1 = "Koce", condition2 = "Hyb",
                cpmThreshold = 1, thresholdSamples = 7)
#Write file         
write.csv(KoceVSHyb,"diffMeth_KocevsHyb_total.csv")

## Total number of DMCs = 50,143 between Koce and Hyb.

#Filter DMCs with FDR < 0.05
KoceVSHyb_filtered <- GRanges(KoceVSHyb$site[which(top$FDR < 0.05)])
summary(KoceVSHyb_filtered)
KoceVSHyb_filtered
#Write file
write.csv(KoceVSHyb_filtered,"diffMeth_KocevsHyb_filtered.csv")

## Total number of DMCs with < 0.05 = 10,620 between Koce and Kher.

```
Finally, we compared samples of _K. hermaphroditus_ vs. hybrids as follows:
```
#Total DMCs
KherVSHyb <- diffMeth(se = se, cateogory = "Group",
                condition1 = "Kher", condition2 = "Hyb",
                cpmThreshold = 1, thresholdSamples = 7)
#Write file         
write.csv(KherVSHyb,"diffMeth_KocevsKher_total.csv")

## Total number of DMCs = 56,390 between Koce and Hyb.

#Filter DMCs with FDR < 0.05
KherVSHyb_filtered <- GRanges(KherVSHyb$site[which(top$FDR < 0.05)])
summary(KherVSHyb_filtered)
KherVSHyb_filtered
#Write file
write.csv(KherVSHyb_filtered,"diffMeth_KhervsHyb_filtered.csv")

## Total number of DMCs with < 0.05 = 13,905 between Koce and Kher.
```
## 2. In shell:
For downstream analysis, we need to have the number of reads for each DMC involved in the comparisons between groups.
