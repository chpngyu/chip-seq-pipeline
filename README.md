# Pipeline of ChIP-seq analysis


## <a name="contents"></a>Contents

* [Computational Pipeline](#computational-pipeline)
	* (1) [Data download](#data-download)
	* (2) [Control selection](#control-selection)
	* (3) [Read quality control](#read-quality-control)
	* (4) [Read mapping](#read-mapping)
	* (5) [Peak calling](#peak-calling)
	* (6) [Motif discovery](#motif-discovery)
	* (7) [Replicate quality control and merging](#replicate-quality-control)
	* (8) [Representative motif selection](#representative-motif-selection)
* [Supplementary Perl Scripts](#supplementary-perl-scripts)
* [Supplementary Data](#supplementary-data)
* [Acknowledgments](#acknowledgments)
* [Citation](#citation)
* [References](#references)


## <a name="computational-pipeline"></a>Computational Pipeline

Our pipeline is represented in the following figure and described in detail below.

<img src="https://github.com/chpngyu/pipeline-of-chip-seq/blob/master/images/flowchart.png">


### (1) <a name="data-download"></a>Data download

Download raw data from ENCODE and prepare a reference genome sequence


### (2) <a name="control-selection"></a>Control selection

Select a control (the DNA control)


### (3) <a name="read-quality-control"></a>Read quality control

ChIP-seq read quality QC

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


### (4) <a name="read-mapping"></a>Read mapping

Preprocessing and aligning reads to the reference genome
Map reads to the genome (e.g., GRCh38 for human) using bowtie2


### (5) <a name="peak-calling"></a>Peak calling

Call read peaks using MACS2 
Remove peak regions overlapping any blacklisted regions


### (6) <a name="motif-discovery"></a>Motif discovery


### (7) <a name="replicate-quality-control"></a>Replicate quality control and merging

Merge reads from replicates? 

For each replicate, use MEME-chip to compute the top 5 PWMs from the top 500 peaks (200bp each, ±100bp from the peak summit). If none of the PWMs is supported by >100 peaks, remove that replicate. If < 2 replicate passes this test, abandon that experiment.

(using PWM similarity instead of IDR? ) If more than one replicate passed the above test, conduct the Irreproducible Discovery Rate (IDR) test: If the IDR score (-log[IDR]) <1.5, abandon the experiment.

Single experiment: 
If more than one replicate passed the above two tests, merge the “passed” replicates. Select the top 500 peaks to calculate PWMs by MEME-chip: Report the PWMs supported by > 100 peaks	

Multiple experiments: 
Treat each experiment as above. If only one experiment passes the above tests, do as the case of one experiment.
If > 1 experiment passed the tests. Use the ranking criterion to select the top 500 peaks from the experiments to compute consensus PWMs and report those PWMs supported by > 100 peaks. 


### (8) <a name="representative-motif-selection"></a>Representative motif selection

One per cell line or tissue

However, before computing the ranking scores we should compute the PCCs between the PWMs of every 2 experiments in 2 different cell lines or samples. If there are more than 2 experiments available, throw out any experiment that does not show at least one motif having PCC>0.80 with a motif from any other experiment. If there are only two experiments and neither experiment has at least one motif with PCC>0.80, select the experiment with PWMs better supported by larger numbers of peaks.


## <a name="analysis-scripts"></a>Analysis Scripts

In addition to published tools, our analyses utilized the following custom scripts:

1. `pwm.py`. **!!ADD DESCRIPTION!!**
2. `motif_cluster.py`. **!!ADD DESCRIPTION!!**
3. `consensus_pwm.py`. **!!ADD DESCRIPTION!!**
4. `correlation.py`. **!!ADD DESCRIPTION!!**
5. `peak_motif_ranges.R`. **!!ADD DESCRIPTION!!**
6. `score_quantile.py`. **!!ADD DESCRIPTION!!**
7. `top_peak_sel.py`. **!!ADD DESCRIPTION!!**


## <a name="supplementary-perl-scripts"></a>Supplementary Perl Scripts

Exact commands for running the pipeline for a single transcription factor (TF) are provided in the following Perl scripts:

1. `chip_seq_download.pl`. This script downloads a user-provided list of ENCODE experiments, *i.e.*, the FASTQ data from ENCODE TF ChIP-seq assays. It should be called from the directory you wish to populate with directories and subdirectories containing the data. Input should be provided as a TAB-delimited input file with the following six (named) columns:

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

2. `chip_seq_pipeline.pl`. This script can be used after the experimental data have been downloaded using `chip_seq_download.pl`. Alternatively, the user may download their own data, provided they have placed them in directories precisely matching the structure described above. 

	This script **must be called from within the directory corresponding to a single Target (TF)**. This means the user must first navigate to one of the `directory_name` values specified above, e.g., abd_A (input file example above) or ATF3 (input figure example above). This means that multiple TFs can be run in parallel, but that a single TF cannot. The reason for this choice is that several steps in the pipeline (e.g., BOWTIE2) allow the user to specify multiple CPUs for further parallelism.
	
	This script navigates the directory structure starting at this point, carrying out the first steps of the pipeline, proceeding from quality control (Trimmomatic; FASTQC) to read mapping (BOWTIE2) to peak calling (MACS2) to motif discovery (MEME-chip). **Of utmost importance, the user must alter the code block at the top of the script, `### MANUALLY SET GLOBAL VARIABLES ###`, to set paths to each tool installed on their system.** Specifically, the user must provide paths or values from the following:

	* `FASTQC` (v0.11.8): path to software
	* `TRIMMOMATIC` (v0.39): path to software
	* `BOWTIE2` (v2.3.5 / 64-bit): path to software
	* `SAMTOOLS` (v1.9 using htslib 1.9): path to software
	* `BEDTOOLS` (v2.28.0): path to software
	* `MACS2` (v2.1.2): path to software
	* `MACS2_gsize`: a value, e.g. `hs` for *Homo sapeins* or `dm` for *Drosophila melanogaster*
	* `PEAK_RADIUS`: a value, the radius of read calling peaks, *e.g.*, `100` to extend peaks 100 bp in either direction from the peak summit.
	* `MIN_PEAK_SCORE`: a value, the minimum peak score required, *e.g.*, `13`.
	* `MEME_CHIP` (v5.0.5): path to software
	* `CCUT`: a value, the size (bp) to which peak regions should be trimmed for MEME-chip, e.g., `100`. This allows MEME-chip to examine the central region of the peaks for motifs while comparing to the flanking regions as a control for local sequence content. **!!CHECK WITH CHUN-PING!!**
	* `adaptor_SE`: path to file (FASTA format) containing sequencing adaptors for single-end (SE) experiments, *e.g.*, `TruSeq3-SE.fa`
	* `adaptor_PE`: path to file (FASTA format) containing sequencing adaptors for paired-end experiments, *e.g.*, `TruSeq3-PE-2.fa`
	* `NUM_TOP_PEAKS`: a value, the number of top peaks to consider, *e.g.*, 500
	* `MFOLD_MIN`: a value, the minimum fold depth enrichment required to call a peak for MACS2, *e.g.*, `5`. **!!CHECK WITH CHUN-PING!!**
	* `MFOLD_MAX`: a value, the maximum fold depth enrichment allowed to call a peak for MACS2, *e.g.*, `50`. **!!CHECK WITH CHUN-PING!!**
	* `blacklist`: path to file (BED format) containing a blacklist (excluded genome regions, including repetitive and low-complexity regions), *e.g.*, `ENCFF023CZC_sorted.bed`. See Amemiya *et al.* (2008).
	* `GENOME_FASTA`: path to file (FASTA format) containing the primary genome assembly for the organism of interest, *e.g.*, `Homo_sapiens.GRCh38.dna.primary_assembly.fa`.
	* `GENOME_IDX_PREFIX`: path to file (genome index) created using BOWTIE2, *e.g.*, `Homo_sapiens.GRCh38.dna.primary_assembly`. This can be accomplished using the following command: `bowtie2-build -f Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38.dna.primary_assembly`.
	
	The pipeline is shown below. Note that these tools use slightly different commands for single-end (SE) and paired-end (PE) reads (see script source code).

<img src="https://github.com/chpngyu/pipeline-of-chip-seq/blob/master/images/perl_pipeline.png">

3. `PWM_pipeline.pl`. IDR and PCC. **!!ADD DESCRIPTION!!**


## <a name="supplementary-data"></a>Supplementary Data
Supplementary data is available in the data folder within this repository. This includes the following files:

* File 1. Description.
* File 2. Description.
* ...

## <a name="acknowledgments"></a>Acknowledgments
The authors acknowledgme XYZ...


## <a name="citation"></a>Citation
When using this software, please refer to and cite:

>THE PAPER **!!ADD INFO!!**

and this page:

>https://github.com/chpngyu/pipeline-of-chip-seq/


## <a name="references"></a>References
* Amemiya HM, Kundaje A, Boyle AP. 2019. <a target="_blank" href="https://www.nature.com/articles/s41598-019-45839-z">The ENCODE Blacklist: Identification of Problematic Regions of the Genome</a>. *Scientific Reports* **9**:9354.

