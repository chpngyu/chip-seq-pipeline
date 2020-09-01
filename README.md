# Pipeline of ChIP-seq analysis

<img src="https://github.com/chpngyu/pipeline-of-chip-seq/blob/master/data/pipeline.png">
Computational pipeline of ENCODE ChIP-seq analysis

For runing whole pipene line for a TF, the command line is
perl ....

The detail processes are described followings step-by-step 
## (1) Download raw data from ENCODE and prepare a reference genome sequence

## (2) Preprocessing and Aligning reads to the reference genome
The following code is used to deal with paired-end reads
```Shell
BASE=/home/your/working/directory
Trimmomatic_Path=/where/is/Trimmomatic
Trimmomatic='java -jar /path/Trimmomatic-0.36/trimmomatic-0.36.jar'
cd $BASE/Sample_one/
R1=LGJ20-XQ41_R1_001.fastq.gz
R2=LGJ20-XQ41_R2_001.fastq.gz
mkdir QC fastqc
$Trimmomatic PE -threads 1 $R1 $R2 QC/${R1%.fastq.gz}_trim_paired.fastq.gz QC/${R1%.fastq.gz}_trim_unpaired.fastq.gz QC/${R2%.fastq.gz}_trim_paired.fastq.gz QC/${R2%.fastq.gz}_trim_unpaired.fastq.gz ILLUMINACLIP:$Trimmomatic_Path/adapters/TruSeq3-PE-2.fa:2:40:12:8:true LEADING:10 SLIDINGWINDOW:4:15 MINLEN:50 2> read_processing.log
$FASTQC QC/*trim_paired.fastq.gz -o fastqc
```
