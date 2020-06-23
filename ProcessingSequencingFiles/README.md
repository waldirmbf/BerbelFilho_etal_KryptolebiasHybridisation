# **Processing sequencing files**

Killfish Project Pipeline | Ultra Documentation - by George Pacheco and Waldir M. Berbel-Filho

This documentation outlines the pipelines used for processing sequencing files used in both genetic and epigenetic analysis  in Berbel-Filho et al. (XXXX)

Last Modified: 3rd June 2020              

Please, contact george.pacheco@snm.ku.dk or waldirmbf@gmail.com should any question arise.
__________________________________________

## **Filtering samples and extracting sorted and indexed BAM files**
### 1) Access to raw data and local storage

The GBS raw data was directly downloaded from the Genomics Research Hub - Cardiff University using an ordinary -wget" command and it is now stored on Sonia's Consuegra local server at Swansea University. The MD5SUM numbers were confirmed for all downloaded files.

### 2) Genotype-by-sequencing (GBS) data

A general sequencing quality check of each plate was performed using the software FastQC using default options. We considered that all sequenced lanes passed this general quality check. The results for each lane are stored inside their respective folders.

```
~/Desktop/msGBS_data/George/GBS_Data/S0_R1_001.fastq.gz
~/Desktop/msGBS_data/George/GBS_Data/S0_R2_001.fastq.gz
```
### 3) Demultiplexing GBS reads

The software **GBSX--v1.3** was used to demultiplex the data based on the barcode **(KFP_Barcodes.info)** info in the respective key files. We allowed for one mismatch in the barcode sequences using the parameter`-mb 1`
```
java -jar ~/Desktop/msGBS_data/Tools/GBSX-master/releases/GBSX_v1.3/GBSX_v1.3.jar --Demultiplexer -t 8 -f1 ~/Desktop/msGBS_data/George/GBS_Data/S0_R1_001.fastq.gz -f2 ~/Desktop/msGBS_data/George/GBS_Data/S0_R2_001.fastq.gz -i ~/Desktop/msGBS_data/George/GBS_Data/KFP_Barcodes.info -mb 1 -me 0 -ca false -gzip true -o ~/Desktop/msGBS_data/George/GBS_Data/KFP-Demultiplexed_GBSX--v1.3/
```
### 4) Merging paired reads

The function **_bbmerge_** in **BBMAP--v.38.32**, was used to merge paired reads from the same individual,  as follows:

````
while read i j k l m;
do
bash ~/Desktop/msGBS_data/Tools/bbmap/bbmerge.sh in1=~/Desktop/msGBS_data/George/GBS_Data/KFP-Demultiplexed_GBSX--v1.3/$i in2=~/Desktop/msGBS_data/George/GBS_Data/KFP-Demultiplexed_GBSX--v1.3/$j out=~/Desktop/msGBS_data/George/GBS_Data/KFP-Merged/$k outu1=~/Desktop/msGBS_data/George/GBS_Data/KFP-Merged/$l outu2=~/Desktop/msGBS_data/George/GBS_Data/KFP-Merged/$m qtrim=r ordered=t mininsert=-150 mininsert0=-150 minlength=25
done < ~/Desktop/msGBS_data/George/GBS_Data/KFP-Merged/KFP-MergeFastq.txt
````

### 5) Genome editing before mapping

Genome for _K.marmoratus_ from [Rhee et al. (2015)](https://www.ncbi.nlm.nih.gov/assembly/GCF_001649575.1) was used for mapping. Genome editing was done to remove unecessary information in the scaffolds name. For the genetic analysis, the scaffold containing the mtDNA genome was removed manually for genetic analysis.

````
cat /data/home/waldir/Downloads/GCF_001649575.1_ASM164957v1_genomic.fna | awk '{split($0,a," "); print a[1]'} > ~/Desktop/msGBS_data/George/Genome/GCF_001649575.1_ASM164957v1_genomic.Edited.fasta
````

Edited genome was indexed using **bowtie2**, as follows:

````
bowtie2-build ~/Desktop/msGBS_data/George/Genome/GCF_001649575.1_ASM164957v1_genomic.Edited.fasta KFP-Genome
````

### 6) Filtering reads, mapping, sorting and indexing samples

Filtering of the reads was done using the function **_bbduk_** in  **BBMAP--v.38.32**. Filtered files were mapped to _K.marmoratus_ genome using **bowtie2**. Resulting SAM files were transformed in BAM files, sorted and indexed using **samtools**. The last line of coding was included to keep  a log file (**KFP-ToRunMapping.log**) with the mapping statistics for each samples as follows:

````
while read i j k l;
do
zcat ~/Desktop/msGBS_data/George/GBS_Data/KFP-Merged/$i.unR1.fastq.gz ~/Desktop/msGBS_data/George/GBS_Data/KFP-Merged/$i.merged.fastq.gz | bash ~/Desktop/msGBS_data/Tools/bbmap/bbduk.sh in=stdin.fq out=~/Desktop/msGBS_data/George/KFP-Mapped/Filtered/$i.fq qtrim=r minlength=25 int=t && ~/Desktop/msGBS_data/Tools/bowtie2-2.3.5-linux-x86_64/bowtie2 -q --threads 8 -x ~/Desktop/msGBS_data/George/Genome/KFP-Genome -U ~/Desktop/msGBS_data/George/KFP-Mapped/Filtered/$i.fq -S ~/Desktop/msGBS_data/George/KFP-Mapped/SAMs/$i.sam && samtools view -bS ~/Desktop/msGBS_data/George/KFP-Mapped/SAMs/$i.sam > ~/Desktop/msGBS_data/George/KFP-Mapped/First_BAMs/$i.bam && samtools sort ~/Desktop/msGBS_data/George/KFP-Mapped/First_BAMs/$i.bam -o ~/Desktop/msGBS_data/George/KFP-Mapped/SortedIndexed/$i.bam && samtools index ~/Desktop/msGBS_data/George/KFP-Mapped/SortedIndexed/$i.bam && mv ~/Desktop/msGBS_data/George/KFP-Mapped/SortedIndexed/$i.bam.bai ~/Desktop/msGBS_data/George/KFP-Mapped/SortedIndexed/$i.bai
done < ~/Desktop/msGBS_data/George/KFP-Mapped/KFP-Mapping.txt

~/Desktop/msGBS_data/George/KFP-Mapped/KFP-ToRunMapping.txt >& ~/Desktop/msGBS_data/George/KFP-Mapped/Logs/KFP-ToRunMapping.log
````

With this pipeline, we ended up with sorted and indexed BAM files for every sample in the msGBS library to be used in the genetic and epigenetic analysis.
