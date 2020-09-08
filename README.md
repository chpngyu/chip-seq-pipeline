# Pipeline of ChIP-seq analysis

<img src="https://github.com/chpngyu/pipeline-of-chip-seq/blob/master/data/flowchart.png">
Computational pipeline of ENCODE ChIP-seq analysis

For runing whole pipene line for a TF, the command line is
perl ....


## The computational pipeline:
```
(A)	Data download 
(B)	Select a control (the DNA control)
(C)	ChIP-seq read quality QC
(D)	Map reads to the human genome (GRCh38) using bowtie2
(E)	Read peaks calling using MACS2 
(F)	Remove peak regions overlapping any blacklisted regions
(G)	Merge reads from replicates? 
  G1.	For each replicate, use MEME-chip to compute the top 5 PWMs from the top 500 peaks (200bp each, ±100bp from the peak summit). If none of the PWMs is supported by >100 peaks, remove that replicate. If < 2 replicate passes this test, abandon that experiment.
  G2.	(using PWM similarity instead of IDR? ) If more than one replicate passed the above test, conduct the Irreproducible Discovery Rate (IDR) test: If the IDR score (-log[IDR]) <1.5, abandon the experiment.
(H)	Single experiment: 
If more than one replicate passed the above two tests, merge the “passed” replicates. Select the top 500 peaks to calculate PWMs by MEME-chip: Report the PWMs supported by > 100 peaks	
(I)	Multiple experiments: 
Treat each experiment as above. If only one experiment passes the above tests, do as the case of one experiment.
If > 1 experiment passed the tests. Use the ranking criterion to select the top 500 peaks from the experiments to compute consensus PWMs and report those PWMs supported by > 100 peaks. 
However, before computing the ranking scores we should compute the PCCs between the PWMs of every 2 experiments in 2 different cell lines or samples. If there are more than 2 experiments available, throw out any experiment that does not show at least one motif having PCC>0.80 with a motif from any other experiment. If there are only two experiments and neither experiment has at least one motif with PCC>0.80, select the experiment with PWMs better supported by larger numbers of peaks.
```

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
