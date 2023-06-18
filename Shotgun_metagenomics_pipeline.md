# **Shotgun Metagenomic Sequencing Data Pipeline for Canola Root Metagenome**
## <center>Yunliang Li </center> 
## <center>liyunliang81@gmail.com</center>
## <center>June 17, 2023</center>
---
## **Removal of host and spike DNA**  
### *Download canola genome and phix-spike genome*
```shell
mkdir /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/ref_seq_conta
# 4 canola genomes
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/020/379/485/GCF_020379485.1_Da-Ae/GCF_020379485.1_Da-Ae_genomic.fna.gz -P /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/ref_seq_conta

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/686/985/GCF_000686985.2_Bra_napus_v2.0/GCF_000686985.2_Bra_napus_v2.0_genomic.fna.gz -P /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/ref_seq_conta

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/026/770/255/GCA_026770255.1_ASM2677025v1/GCA_026770255.1_ASM2677025v1_genomic.fna.gz -P /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/ref_seq_conta

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/026/770/265/GCA_026770265.1_ASM2677026v1/GCA_026770265.1_ASM2677026v1_genomic.fna.gz -P /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/ref_seq_conta

# phix-spike genome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/814/215/GCF_002814215.1_Sedor1/GCF_002814215.1_Sedor1_genomic.fna.gz  -P /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/ref_seq_conta
```
### Unzip the compressed .gz file and concatenate the depressed genome data
```
gunzip /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/ref_seq_conta/*.gz
cat /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/ref_seq_conta/*.fna > phix_canola.fasta
```
### *Build reference genomes idex files*
```
bowtie2-build phix_canola.fasta /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/ref_seq_conta/index/ --threads 30
```
### *Removal of contaminated DNA using kneaddata*
``` shell
# my sequencing files are D2000{120:145}, but no 133 and 136； using shell "for" loop to run kneaddata for each sample

#!/bin/sh 
for i in $(seq 120 145)
do
  if [ $i -ne 133 ] && [ $i -ne 136 ]
  then
cp /u1/SqueezeMeta/Categorized_Meta_Genome_Sequences/MetaGenome_Combined/D2000"$i"* .
cat D2000"$i"_S*R1_* > D2000"$i"_R1.fastq.gz
cat D2000"$i"_S*R2_* > D2000"$i"_R2.fastq.gz
mkdir D2000$i
kneaddata -t 40 -i1 /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/D2000"$i"_R1.fastq.gz -i2 /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/D2000"$i"_R2.fastq.gz -o /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/D2000"$i" -db /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/ref_seq_conta/index/phix_canola.index --bowtie2-options "--very-sensitive --dovetail --reorder" --trimmomatic /home/yul713/miniconda3/envs/kneaddata/bin --trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:50" --remove-intermediate-output

rm D2000"$i"_*
done
```
## **Taxonomy annotation at reads level using kraken2**

### *Build kraken2 database and generate the Bracken database files*
```shell
# Download the NCBI taxonomy
kraken2-build --db /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/krakendb --download-taxonomy 

# Download reference libraries
kraken2-build --download-library bacteria --db /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/krakendb
kraken2-build --download-library fungi --db /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/krakendb
kraken2-build --download-library archaea --db /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/krakendb
kraken2-build --download-library viral --db /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/krakendb
kraken2-build --download-library protozoa --db /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/krakendb

# The phix_canola.fasta was also added into the reference database
kraken2-build --add-to-library /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/ref_seq_conta/phix_canola.fasta --db /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/krakendb

# Build the kraken2 reference database
kraken2-build --db /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/krakendb --build --threads 30

# Generate the Bracken database files
bracken-build -d /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/krakendb -t 30 -k 35 -l 100
```
### *Taxonomy assignment* 

```shell
# make two folders for the results
mkdir /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/kraken_outputs
mkdir /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/kreports

# taxonomy assignment for each sample using shell "for" loop

#!/bin/sh 
for i in $(seq 120 145)
do
  if [ $i -ne 133 ] && [ $i -ne 136 ]
  then
    kraken2 --db /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/krakendb --threads 16 --report /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/kreports/D2000"$i".k2report \
    --report-minimizer-data --minimum-hit-groups 3 /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/D2000"$i"/D2000"$i"_R1_kneaddata_paired_1.fastq \
    /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/D2000"$i"/D2000"$i"_R1_kneaddata_paired_2.fastq > /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/kraken_outputs/D2000"$i".kraken2
  fi
done
```
### *Run bracken for abundance estimation of microbiome samples*

```shell
# make two folders for the results
mkdir /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/bracken_outputs
mkdir /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/breports

# abundance estimation for each sample using shell "for" loop

#!/bin/sh 
for i in $(seq 120 145)
do
  if [ $i -ne 133 ] && [ $i -ne 136 ]
  then
 bracken -d /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/krakendb \
-i /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/kreports/D2000"$i".k2report \
-r 100 -l S -t 10 \
-o /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/bracken_outputs/D2000"$i".bracken \
-w /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/breports/D2000"$i".breport

fi
done
```
## **Assembly reads to contigs using MetaSPAdes**
```shell
# Assembly reads of each individual sample to contigs using shell "for" loop 

#!/bin/sh 
for i in $(seq 120 145)
do
  if [ $i -ne 133 ] && [ $i -ne 136 ]
  then
mkdir -p /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/Asse_D2000"$i"

metaspades.py -t 40 --only-assembler -m 1000 -1 /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/D2000"$i"/D2000"$i"_R1_kneaddata_paired_1.fastq -2 /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/D2000"$i"/D2000"$i"_R1_kneaddata_paired_2.fastq -o /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/Asse_D2000"$i"

fi
done
```
## **Reorder reads** 
If you forgot to use bowtie2-options '--reorder' in kneaddata step, you have to run this step to keep the reads in order, which is required for mapping reads to contig using bowtie2. If you already used the bowtie2-options '--reorder', this step can be skipped

### *Reorder reads using fastq-repair*
``` shell
#!/bin/sh 
for i in $(seq 120 145)
do
  if [ $i -ne 133 ] && [ $i -ne 136 ]
  then

fastq_pair  /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/D2000"$i"/D2000"$i"_R1_kneaddata_paired_1.fastq /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/D2000"$i"/D2000"$i"_R1_kneaddata_paired_2.fastq
fi
done
```
## **Mapping reads to contigs**
```shell
#!/bin/sh 
for i in $(seq 120 145)
do
  if [ $i -ne 133 ] && [ $i -ne 136 ]
  then

# Index contig
bowtie2-build /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/Asse_D2000"$i"/contigs.fasta index_contig_D2000"$i" 

# Mapping 
bowtie2 --no-unal -p 40 -x index_contig_D2000"$i" -1 /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/D2000"$i"/D2000"$i"_R1_kneaddata_paired_1.fastq.paired.fq \
-2 /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/D2000"$i"/D2000"$i"_R1_kneaddata_paired_2.fastq.paired.fq -S /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/Asse_D2000"$i"/D2000"$i".sam

fi
done
```
## **Generate Metagenome-Associated Genome (MAG) using MetaBAT2, Binning**

### *Combine contigs.fasta of each individual sample*
```shell
mkdir /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/merged4bin

# copy contigs.fasta to the folder merged4bin

#!/bin/sh 
for i in $(seq 120 145)
do
  if [ $i -ne 133 ] && [ $i -ne 136 ]
  then

cp /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/Asse_D2000"$i"/contigs.fasta  /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/merged4bin/D2000"$i"_contigs.fasta

fi
done

# concatenate all contigs.fasta as merged.fasta
cat /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/merged4bin/D2000*.fasta > /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/merged4bin/merged.fasta
```
### **Combine decontaminated reads (R1 and R2) of each individual sample respectively**
```shell
# copy the decontaminated reads (R1 and R2) to merged4bin folder

#!/bin/sh 
for i in $(seq 120 145)
do
  if [ $i -ne 133 ] && [ $i -ne 136 ]
  then

cp /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/D2000"$i"/D2000"$i"_R1_kneaddata_paired_1.fastq.paired.fq  /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/merged4bin/D2000"$i"_R1_kneaddata_paired_1.fastq.paired.fq

cp /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/D2000"$i"/D2000"$i"_R1_kneaddata_paired_2.fastq.paired.fq  /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/merged4bin/D2000"$i"_R1_kneaddata_paired_2.fastq.paired.fq

fi
done

# combine all R1 reads as merged_R1.fq, and all R2 reads as merged_R2.fq 
cat /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/merged4bin/D2000*_1.fastq.paired.fq > /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/merged4bin/merged_R1.fq

cat /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/merged4bin/D2000*_1.fastq.paired.fq > /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/merged4bin/merged_R2.fq
```
### *Map reads to contigs*
```shell
bowtie2-build /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/merged4bin/merged.fasta \
/u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/merged4bin/index_merged_contig

bowtie2 --no-unal -p 40 -x /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/merged4bin/index_merged_contig -1 /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/merged4bin/merged_R1.fq \
-2 /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/merged4bin/merged_R2.fq -S /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/merged4bin/merged.sam
```
### *Bin contigs using MetaBAT2*
```shell
# SAM file need to be sorte
samtools sort /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/merged4bin/merged.sam -o /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/merged4bin/merged_sorted.bam -@ 30

# binning contigs
jgi_summarize_bam_contig_depths --outputDepth /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/merged4bin/contigs.fasta.depth.txt /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/merged4bin/merged_sorted.bam
 
metabat2 --inFile  /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/merged4bin/merged.fasta --outFile /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/merged4bin/bin --abdFile /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/merged4bin/contigs.fasta.depth.txt --minContig 2000
```
## **Assessment, classification and refinement of MAGs using mdmcleaner**
### *Build mdmcleaner database*
```
mkdir /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/mdm_db
mdmcleaner makedb -o /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/mdm_db --verbose
```
### *Assess contigs*
```shell
mdmcleaner set_configs --db_basedir /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/mdm_db --threads 30

#!/bin/sh 
for i in $(seq 1 18) # the number depends on how many bins were generated
do

mdmcleaner clean -i /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/merged4bin/bin."$i".fa -o /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/merged4bin -t 30 -overview_files_basename "overview"$i""

done
```
### *Combine summary results*
```shell
cat /u2/SqueezeMeta/Canola_shotgun_metagenome_Peter/assembly/merged4bin/overview*tsv > overview_merged.tsv
```