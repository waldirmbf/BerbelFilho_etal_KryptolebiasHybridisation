# **Genetic analysis with SNPs from msGBS library**

Killfish Project Pipeline | Ultra Documentation - by George Pacheco and Waldir M. Berbel-Filho

This documentation outlines the pipelines used for genetic analysis (SNPs genotypes and sites extracted from msGBS library) in Berbel-Filho et al. (XXXX)

Last Modified: 22nd June 2020              

Please, contact george.pacheco@snm.ku.dk or waldirmbf@gmail.com should any question arise.
__________________________________________
The original methylation-sensitive genotype-by-sequencing (msGBS) library contained samples from four _Kryptolebias_ species, on which we were interested to investigate the relationships between genetic and epigenetic variation across mating systems. Given the evidence of hybridisation between _K. ocellatus_ and _K. hermaphroditus_, we decided to subset the samples containing only _K. ocellatus_ and _K. hermaphroditus_ populations. We set up a threshold of only samples with at least 500.000 reads to be included in this dataset. The number of total number of  reads per sample in the library was retrieved from the demultiplexing results in [gbsDemultiplex.stats.txt](gbsDemultiplex.stats.txt).

__________________________________________

### 1) Filtering samples according to threshold of reads

A [file containing a list of samples](KFP--GoodSamplesReads.list)  with number of reads above the threshold (>500k reads) was created was created as follows:

```
tail -n +2 ~/Desktop/msGBS_data/George/GBS_Data/KFP-Demultiplexed_GBSX--v1.3/gbsDemultiplex.stats | awk '$4>500000' | cut -f1 > ~/Desktop/msGBS_data/George/KFP-Mapped/KFP--GoodSamplesReads.list
```

### 2)  Further filtering samples according to species and populations of interest

To make a list of BAM files containing only populations of _K. ocellatus_ and _K. hermaphroditus_, we further filtered the data removing samples from freshwater _Kryptolebias_ species (_K. brasiliensis and K. gracilis_) and another _Kryptolebias_ population from Espírito-Santo in Brazil (ES). As follows:

```
find ~/Desktop/msGBS_data/George/KFP-Mapped/SortedIndexed/*.bam | grep -f ~/Desktop/msGBS_data/George/KFP-Mapped/KFP--GoodSamplesReads.list | grep -f ~/Desktop/msGBS_data/George/KFP-Mapped/KFP--GoodSamplesReadsNoKbraKgra.list | fgrep -v -f ~/Desktop/msGBS_data/George/KFP-Mapped/KFP--GoodSamples_OnlyES1_Kher.list > ~/Desktop/msGBS_data/George/KFP-Mapped/KFP--GoodSamplesReads_NoES1_Kher.BAMlist
```
The FastQC files inclued in this project can be acessed at
[NCBI acession PRJNA563625](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA563625)

### 3)  Number of sites covered in the dataset (Dataset  I)
The software **ANGSD v. 0.92** was used to calculate the total number of reads in the dataset under the following parameters.

```
~/Desktop/msGBS_data/Tools/ngsTools/angsd/angsd -nThreads 2 -ref ~/Desktop/msGBS_data/George/Genome/GCF_001649575.1_ASM164957v1_genomic.Edited.fasta -bam ~/Desktop/msGBS_data/George/KFP-Mapped/KFP--GoodSamplesReads_NoES1_Kher.BAMlist -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 -minMapQ 30 -minQ 20 -minInd $((53*95/100)) -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 1 -doPost 2 -doGeno 3 -doPlink 2 -geno_minDepth 3 -setMaxDepth $((53*500)) -dumpCounts 2 -postCutoff 0.95 -doHaploCall 1 -doVcf 1 -out ~/Desktop/msGBS_data/George/KFP-ANGSDRuns/KFP--GoodSamplesReads_NoES1_Kher_SITES--Article--Ultra
```
`-remove_bads 1` = Discard 'bad' reads, (flag>=256).

`-uniqueOnly 1` = Discard reads that doesn't map uniquely.

`-baq 1`= Adjust qscores around indels (1=normal baq 2= extended(as SAMtools)).

`-C 50` = Adjust mapQ for excessive mismatches (as SAMtools)

`-minMapQ 30` = Discard reads with mapping quality below 30.

`-minQ 20` = Discard bases with base quality below 20.

`-minInd $((53*95/100))` = Only use site if at least 95%  of all 53 individuals has data.

`-doCounts 1` = Count the number A,C,G,T. All sites, All samples.

`-GL 1` = Genotype-likelihood using SAMtools models.

`-doGlf 2`= Output the log genotype likelihoods to a file as beagle genotype likelihood format (use directly for imputation).

`-doMajorMinor 1` = Infer major and minor allele from GL.

`-doMaf 1` = Known major, and Known minor.

`-doPost 2` = Calculate posterior probabiblity based on uniform distribution.

`-doGeno 3` = Write the called genotype directly (eg. AA,AC).

`-doPlink 2` = Create Plink file in a tfam/tped format.

`-geno_minDepth 3` = Minimum depth of 3 reads per site across all individuals.

`-setMaxDepth $((53*500))` = Maximum depth to include site is 500.

`-dumpCounts 2` = Print the depth for each site for each individual.

`-postCutoff 0.95` = Only genotype to missing if below this threshold (5%).

`-doHaploCall 1` = When haploid calling, use a random base.

`-doVcf 1` =  Create VCF file.


The total number of sites in the dataset was __597.333__.


#### 3.1)  Real coverage calculation per sample in Dataset I
The following script was used to calculate the average depth of reads across all sites per individual in the dataset I.

```
zcat ~/Desktop/msGBS_data/George/KFP-ANGSDRuns/KFP--GoodSamplesReads_NoES1_Kher_SITES--Article--Ultra.counts.gz | tail -n +2 | gawk ' {for (i=1;i<=NF;i++){a[i]+=$i;++count[i]}} END{ for(i=1;i<=NF;i++){print a[i]/count[i]}}' | paste ~/Desktop/msGBS_data/George/KFP-Mapped/KFP--GoodSamplesReads_NoES1_Kher.labels - > ~/Desktop/msGBS_data/George/KFP-Analyses/KFP-Miscellaneous/KFP-RealCoverage/KFP--GoodSamplesReads_NoES1_Kher_SITES--Article--Ultra.GL-RealCoverage.txt
```
The results can be found in this [file](KFP--GoodSamplesReads_NoES1_Kher_SITES--Article--Ultra.GL-RealCoverage).

#### 3.2)  Average percentage of missing data per sample in Dataset I
The following script was used to the average percentage of missing data across all sites per individual in the dataset I.

```
zcat ~/Desktop/msGBS_data/George/KFP-ANGSDRuns/KFP--GoodSamplesReads_NoES1_Kher_SITES--Article--Ultra.beagle.gz | tail -n +2 | perl ~/Desktop/msGBS_data/Tools/Scripts/call_geno.pl --skip 3 | cut -f 4- | awk '{ for(i=1;i<=NF; i++){ if($i==-1)x[i]++} } END{ for(i=1;i<=NF; i++) print i"\t"x[i] }' | paste ~/Desktop/msGBS_data/George/KFP-Mapped/KFP--GoodSamplesReads_NoES1_Kher.labels - | awk '{print $1"\t"$3"\t"$3*100/597733}' > ~/Desktop/msGBS_data/George/KFP-Analyses/KFP-Miscellaneous/KFP-MissingData/KFP--GoodSamplesReads_NoES1_Kher_SITES--Article--Ultra.GL-MissingData.txt
```

The results can be found in this [file](KFP--GoodSamplesReads_NoES1_Kher_SITES--Article--Ultra.GL-MissingData.txt).

### 4) Calling SNPs (Dataset  II)
The software **ANGSD v. 0.92** was used to call SNPs genotypes calculate the  in the dataset under the following parameters.

```
~/Desktop/msGBS_data/Tools/ngsTools/angsd/angsd -nThreads 2 -ref ~/Desktop/msGBS_data/George/Genome/GCF_001649575.1_ASM164957v1_genomic.Edited.fasta -bam ~/Desktop/msGBS_data/George/KFP-Mapped/KFP--GoodSamplesReads_NoES1_Kher.BAMlist -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 -minMapQ 30 -minQ 20 -minInd $((53*95/100)) -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 1 -MinMaf 0.02 -SNP_pval 1e-6 -doPost 2 -doGeno 3 -doPlink 2 -geno_minDepth 3 -setMaxDepth $((53*500)) -dumpCounts 2 -postCutoff 0.95 -doHaploCall 1 -doVcf 1 -out ~/Desktop/msGBS_data/George/KFP-ANGSDRuns/KFP--GoodSamplesReads_NoES1_Kher--Article--Ultra
```
`-remove_bads 1` = Discard 'bad' reads, (flag>=256).

`-uniqueOnly 1` = Discard reads that doesn't map uniquely.

`-baq 1`= Adjust qscores around indels (1=normal baq 2= extended(as SAMtools)).

`-C 50` = Adjust mapQ for excessive mismatches (as SAMtools)

`-minMapQ 30` = Discard reads with mapping quality below 30.

`-minQ 20` = Discard bases with base quality below 20.

`-minInd $((53*95/100))` = Only use site if at least 95% of all 53 individuals has data.

`-doCounts 1` = Count the number A,C,G,T. All sites, All samples.

`-GL 1` = Genotype-likelihood using SAMtools models.

`-doGlf 2`= Output the log genotype likelihoods to a file as beagle genotype likelihood format (use directly for imputation).

`-doMajorMinor 1` = Infer major and minor allele from GL.

`-doMaf 1` = Known major, and Known minor.

`-MinMaf 0.02`= Only work with sites with Maf above 0.02.

`-SNP_pval 1e-6` = Only work with sites with a p-value less than 1e-6.

`-doPost 2` = Calculate posterior probabiblity based on uniform distribution.

`-doGeno 3` = Write the called genotype directly (eg. AA,AC).

`-doPlink 2` = Create Plink file in a tfam/tped format.

`-geno_minDepth 3` = Minimum depth of 3 reads per site across all individuals.

`-setMaxDepth $((53*500))` = Maximum depth to include site is 500.

`-dumpCounts 2` = Print the depth for each site for each individual.

`-postCutoff 0.95` = Only genotype to missing if below this threshold (5%).

`-doHaploCall 1` = When haploid calling, use a random base.

`-doVcf 1` =  Create VCF file.

The total number of SNPs in the dataset was __5.477__.


#### 4.1)  Real coverage calculation per sample in Dataset II
The following script was used to calculate the average depth of reads across all sites per individual in the dataset II.
```
zcat ~/Desktop/msGBS_data/George/KFP-ANGSDRuns/KFP--GoodSamplesReads_NoES1_Kher--Article--Ultra.counts.gz | tail -n +2 | gawk ' {for (i=1;i<=NF;i++){a[i]+=$i;++count[i]}} END{ for(i=1;i<=NF;i++){print a[i]/count[i]}}' | paste ~/Desktop/msGBS_data/George/KFP-Mapped/KFP--GoodSamplesReads_NoES1_Kher.labels - > ~/Desktop/msGBS_data/George/KFP-Analyses/KFP-Miscellaneous/KFP-RealCoverage/KFP--GoodSamplesReads_NoES1_Kher--Article--Ultra.GL-RealCoverage.txt
```

The results can be found in this [file](KFP-RealCoverageKFP--GoodSamplesReads_NoES1_Kher--Article--Ultra.GL-RealCoverage.txt).

#### 4.2)  Average percentage of missing data per sample in Dataset II
The following script was used to the average percentage of missing data across all sites per individual in the dataset II.
```
zcat ~/Desktop/msGBS_data/George/KFP-ANGSDRuns/KFP--GoodSamplesReads_NoES1_Kher--Article--Ultra.beagle.gz | tail -n +2 | perl ~/Desktop/msGBS_data/Tools/Scripts/call_geno.pl --skip 3 | cut -f 4- | awk '{ for(i=1;i<=NF; i++){ if($i==-1)x[i]++} } END{ for(i=1;i<=NF; i++) print i"\t"x[i] }' | paste ~/Desktop/msGBS_data/George/KFP-Mapped/KFP--GoodSamplesReads_NoES1_Kher.labels - | awk '{print $1"\t"$3"\t"$3*100/5477}' > ~/Desktop/msGBS_data/George/KFP-Analyses/KFP-Miscellaneous/KFP-MissingData/KFP--GoodSamplesReads_NoES1_Kher--Article--Ultra.GL-MissingData.txt
```

The results can be found in this [file](KFP--GoodSamplesReads_NoES1_Kher--Article--Ultra.GL-MissingData.txt).

### 5) Genomic analysis with dataset I
#### 5.1) Sites distribution
To have an ideia of the sites density across the reference genome in our dataset, the number of scaffolds with at least one site was reriteved as follows:
```
zcat ~/Desktop/msGBS_data/George/KFP-ANGSDRuns/KFP--GoodSamplesReads_NoES1_Kher_SITES--Article--Ultra.mafs.gz | tail -n +2 | sort -u -k 1,1 | wc -l
```
 __929 (30.24% of the total)__ scaffolds with had at least one site.

For sites density per scaffold, we used the following scripts:

```
zcat ~/Desktop/msGBS_data/George/KFP-ANGSDRuns/KFP--GoodSamplesReads_NoES1_Kher_SITES--Article--Ultra.mafs.gz | tail -n +2 | cut -f1 | sort | uniq -c | awk '{print $2"\t"$1}' | sort -n -k 2,2 > ~/Desktop/msGBS_data/George/KFP-ANGSDRuns/KFP--GoodSamplesReads_NoES1_Kher_SITES--Article--Ultra.SITESDensity.txt
```
Number of sites per scaffold which had at least one site is contained in this [file](KFP--GoodSamplesReads_NoES1_Kher_SITES--Article--Ultra.SITESDensity). Now, we extracted how long each scaffold is in the genome reference using:
```
awk 'BEGIN{OFS="\t"} NR==FNR{x[$1]=$2} NR!=FNR && $2>1000{if(!x[$1])x[$1]=0; print $1,$2,x[$1]}' ~/Desktop/msGBS_data/George/KFP-ANGSDRuns/KFP--GoodSamplesReads_NoES1_Kher_SITES--Article--Ultra.SITESDensity.txt ~/Desktop/msGBS_data/George/Genome/GCF_001649575.1_ASM164957v1_genomic.Edited.fasta.fai | sort -n -k 2,2 > ~/Desktop/msGBS_data/George/KFP-ANGSDRuns/KFP--GoodSamplesReads_NoES1_Kher_SITES--Article--Ultra.ScaffoldInfo.txt
```
 Number of base pairs in each scaffold in the reference genome is contained in this [file](KFP--GoodSamplesReads_NoES1_Kher_SITES--Article--Ultra.ScaffoldInfo). We generated a new file now containing lenght information for the scaffolds containting at least one site as follows:
```
awk '{if ($3!=0) print;}' ~/Desktop/msGBS_data/George/KFP-ANGSDRuns/KFP--GoodSamplesReads_NoES1_Kher_SITES--Article--Ultra.ScaffoldInfo.txt > ~/Desktop/msGBS_data/George/KFP-ANGSDRuns/KFP--GoodSamplesReads_NoES1_Kher_SITES--Article--Ultra.ScaffoldInfo_OnlyWithSites.txt
```
Finally to correlate sites density and scaffold lenght, we used this [**R script**](KFP--ToPlot_ScaffoldLength-NumberOfSites.R), and following is the plot:

![Scattter plot with number of sites and scaffold lenght](SitesvsScaffoldLenght.jpg)
