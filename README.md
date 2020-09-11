# Pipeline of ChIP-seq analysis


## <a name="contents"></a>Contents

* [Computational Pipeline](#computational-pipeline)
	* [Download raw data from ENCODE and prepare a reference genome sequence](#download-genome)
	* [Preprocessing and aligning reads to the reference genome](#preprocessing-aligning)
* [Supplementary Perl Scripts](#supplementary-perl-scripts)
* [Supplementary Data](#supplementary-data)
* [Acknowledgments](#acknowledgments)
* [Citation](#citation)
* [References](#references)


## <a name="computational-pipeline"></a>Computational Pipeline:
```
(A)	Data download 
(B)	Select a control (the DNA control)
(C)	ChIP-seq read quality QC
(D)	Map reads to the genome (e.g., GRCh38 for human) using bowtie2
(E)	Call read peaks using MACS2 
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

<img src="https://github.com/chpngyu/pipeline-of-chip-seq/blob/master/images/flowchart.png">


The detail processes are described followings step-by-step 
### (1) <a name="download-genome"></a>Download raw data from ENCODE and prepare a reference genome sequence

### (2) <a name="preprocessing-aligning"></a>Preprocessing and aligning reads to the reference genome
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


## <a name="supplementary-perl-scripts"></a>Supplementary Perl Scripts

Exact commands for running the pipeline for a single TF are provided in the following Perl scripts:

1. `chip_seq_download.pl`. This script downloads a list of ENCODE experiments provided as a TAB-delimited input file with the following six (named) columns:

	* `Target name`. The name of the transcription factor (*e.g.*, **adb-A**).
	* `directory_name`. The name of the directory to create in which to store the transcription factor's data (*e.g.*, **adb_A**)
	* `AC`. Accession ID of the relevant experiment (*e.g.*, **ENCSR609JDR**). Note that one experiment will be associated with more than one replicate and control, i.e., each Accession ID will have more than one row in the file. A new subdirectory within `directory_name` will be created for each unique `AC`.
	* `type`. Must contain the value `replicate` or `control`. 
	* `Library`. Library (assay) accession ID for the specific replicate or control (e.g., **ENCLB367MZH**). Individual assays that are used as controls will have `_control` appended to the name of their `Library` subdirectory.
	* `url`. The ENCODE url for direct download of the FASTQ data associated with this `AC`/`Library` (*e.g.*, **https://www.encodeproject.org/files/ENCFF036RZV/@@download/ENCFF036RZV.fastq.gz**). This can be constructed by appending the FASTQ file ID to the end of `https://www.encodeproject.org/files/ENCFF036RZV/@@download/`.

	For example, a short file for downloading some *Drosophila* FASTQ data is shown below, followed by a figure explaining the subdirectory structure that will be created in the process. This structure will be expected by the subsequent scripts (*e.g.*, `chip_seq_pipeline.pl`). Note that replicate and control subdirectories will contain one FASTQ file for single-end (SE) reads and two FASTQ files for paired-end (PE) reads; this is how the pipeline differentiates between the two.

```
Target name	directory_name	AC	type	Library	url
abd-A	abd_A	ENCSR609JDR	replicate	ENCLB367MZH	https://www.encodeproject.org/files/ENCFF036RZV/@@download/ENCFF036RZV.fastq.gz
abd-A	abd_A	ENCSR609JDR	replicate	ENCLB506QBD	https://www.encodeproject.org/files/ENCFF999PWE/@@download/ENCFF999PWE.fastq.gz
abd-A	abd_A	ENCSR609JDR	replicate	ENCLB589GMA	https://www.encodeproject.org/files/ENCFF646QDZ/@@download/ENCFF646QDZ.fastq.gz
abd-A	abd_A	ENCSR609JDR	control	ENCLB428GIT_control	https://www.encodeproject.org/files/ENCFF120UZS/@@download/ENCFF120UZS.fastq.gz
Abd-B	Abd_B	ENCSR465HPZ	replicate	ENCLB009GTL	https://www.encodeproject.org/files/ENCFF469WCT/@@download/ENCFF469WCT.fastq.gz
Abd-B	Abd_B	ENCSR465HPZ	replicate	ENCLB205LSV	https://www.encodeproject.org/files/ENCFF969SEG/@@download/ENCFF969SEG.fastq.gz
Abd-B	Abd_B	ENCSR465HPZ	control	ENCLB730TLW_control	https://www.encodeproject.org/files/ENCFF730IEL/@@download/ENCFF730IEL.fastq.gz
Abd-B	Abd_B	ENCSR692UBK	replicate	ENCLB526KET	https://www.encodeproject.org/files/ENCFF820QJV/@@download/ENCFF820QJV.fastq.gz
Abd-B	Abd_B	ENCSR692UBK	replicate	ENCLB588NTL	https://www.encodeproject.org/files/ENCFF252IVG/@@download/ENCFF252IVG.fastq.gz
Abd-B	Abd_B	ENCSR692UBK	control	ENCLB728VIS_control	https://www.encodeproject.org/files/ENCFF354CHN/@@download/ENCFF354CHN.fastq.gz
achi	achi	ENCSR959SWC	replicate	ENCLB061SFK	https://www.encodeproject.org/files/ENCFF979YUJ/@@download/ENCFF979YUJ.fastq.gz
achi	achi	ENCSR959SWC	replicate	ENCLB064AEO	https://www.encodeproject.org/files/ENCFF881ITO/@@download/ENCFF881ITO.fastq.gz
achi	achi	ENCSR959SWC	replicate	ENCLB722DLY	https://www.encodeproject.org/files/ENCFF721FQK/@@download/ENCFF721FQK.fastq.gz
achi	achi	ENCSR959SWC	control	ENCLB240LXI_control	https://www.encodeproject.org/files/ENCFF548BRW/@@download/ENCFF548BRW.fastq.gz
```

<img src="https://github.com/chpngyu/pipeline-of-chip-seq/blob/master/images/directory_structure.png">

<img src="https://github.com/chpngyu/pipeline-of-chip-seq/blob/master/images/example.png">

2. `chip_seq_pipeline.pl`. Assuming the experimental data have been downloaded using `chip_seq_download.pl`, this script navigates the directory structure set up during download to carry out the first steps of the pipeline, as shown below (SE=single-end; PE=paired-end reads)...

<img src="https://github.com/chpngyu/pipeline-of-chip-seq/blob/master/images/perl_pipeline.png">

3. `etc.`...


## <a name="supplementary-data"></a>Supplementary Data
Supplementary data is available in the data folder within this repository. This includes the following files:

* File 1. Description.
* File 2. Description.
* ...

## <a name="acknowledgments"></a>Acknowledgments
The authors acknowledgme XYZ...


## <a name="citation"></a>Citation
When using this software, please refer to and cite:

>THE PAPER

and this page:

>https://github.com/chpngyu/pipeline-of-chip-seq/


## <a name="references"></a>References
* EXAMPLE: Yang Z. 2014. <a target="_blank" href="https://www.oxfordscholarship.com/view/10.1093/acprof:oso/9780199602605.001.0001/acprof-9780199602605">*Molecular Evolution: A Statistical Approach*</a>. New York, NY: Oxford University Press.
