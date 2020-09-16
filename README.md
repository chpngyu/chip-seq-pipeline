# ChIP-seq analysis: supplementary material


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
* [Software and Data Versions Used](#software-and-data)
* [Analysis Scripts](#analysis-scripts)
* [Supplementary Perl Scripts](#supplementary-perl-scripts)
* [Supplementary Data](#supplementary-data)
* [Acknowledgments](#acknowledgments)
* [Citation](#citation)
* [References](#references)


## <a name="computational-pipeline"></a>Computational pipeline

Our pipeline is represented in the following figure and described in detail below.

<img src="https://github.com/chpngyu/pipeline-of-chip-seq/blob/master/images/flowchart.png">


### (1) <a name="data-download"></a>Data download

Raw transcription factor ChIP-seq FASTQ reads were downloaded directly from ENCODE. For example, for *GATA1*, to obtain the SE36nt reads associated with Accession ID ENCFF000YND (experiment ENCSR000EFT, library ENCLB209AJT), the following URL was used:

`https://www.encodeproject.org/files/ENCFF000YND/@@download/ENCFF000YND.fastq.gz`


### (2) <a name="control-selection"></a>Control selection

Controls were selected as described in the manuscript, with preference for type input DNA, as shown in the figure below.

<img src="https://github.com/chpngyu/pipeline-of-chip-seq/blob/master/images/choosing_controls.png">

### (3) <a name="read-quality-control"></a>Read quality control

Read quality control (QC) was performed using Trimmomatic, with command-line argument values as decribed in the manuscript. For paired-end (PE) reads, the following commands were used:

```Shell
#set global variables
BASE=/home/your/working/directory
Trimmomatic_Path=/where/is/Trimmomatic
Trimmomatic='java -jar /path/Trimmomatic-0.36/trimmomatic-0.36.jar'

# replace LGJ20-XQ41_R1_001 and LGJ20-XQ41_R2_001 with your read IDs
cd $BASE/Sample_one/
R1=LGJ20-XQ41_R1_001.fastq.gz
R2=LGJ20-XQ41_R2_001.fastq.gz
mkdir QC fastqc

# QC
$Trimmomatic PE -threads 1 $R1 $R2 QC/${R1%.fastq.gz}_trim_paired.fastq.gz QC/${R1%.fastq.gz}_trim_unpaired.fastq.gz QC/${R2%.fastq.gz}_trim_paired.fastq.gz QC/${R2%.fastq.gz}_trim_unpaired.fastq.gz ILLUMINACLIP:$Trimmomatic_Path/adapters/TruSeq3-PE-2.fa:2:40:12:8:true LEADING:10 SLIDINGWINDOW:4:15 MINLEN:50 2> read_processing.log

# check read quality
$FASTQC QC/*trim_paired.fastq.gz -o fastqc
```

For single-end (SE) reads, the Trimmomatic call was replaced with the following:

```Shell
$Trimmomatic SE -threads 1 $READ QC/${READ%.fastq.gz}_trim_paired.fastq.gz ILLUMINACLIP:/home/cpyu/bin/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:40:12 LEADING:10 SLIDINGWINDOW:4:15 MINLEN:30 2> read_processing.log
```


### (4) <a name="read-mapping"></a>Read mapping

Reads were preprocessed and aligned reads to the appropriate reference genome (e.g., GRCh38 for human) using bowtie2, as shown below.

```Shell
#set global variables
export GENOME=/home/user/human/genome/Homo_sapiens.GRCh38.dna.primary_assemblyexport
BASE_PATH=/home/user/human/ChIP-seq/
export SUFFIX=_trim_paired.fastq.gz
cd $BASE_PATH/GENE_ID/replicate1

# replace _read_ID_ with your IDs
export R1=QC/_read_ID_$SUFFIX
export R2=QC/_read_ID_$SUFFIX

# do alignment
time bowtie2 -p 4 -x $GENOME -1 $R1 -2 $R2 -S alignment.sam 2> log.txt

# if single-read, use
export READ=QC/_read_ID_$SUFFIX
time bowtie2 -p 4 -x $GENOME -U $READ -S alignment.sam 2> log.txt

#remove unmapped reads and duplicated reads (268= Read unmapped (4) or  Mate unmapped (8) or  Not primary alignment (256))
samtools view -h -F 268 -q 5 -bS alignment.sam > unique_alignment.bam
samtools sort unique_alignment.bam -o unique_alignment_sorted.bam
samtools rmdup unique_alignment_sorted.bam unique_alignment_sorted_rd.bam
rm alignment.sam unique_alignment.bam unique_alignment_sorted.bam
```

### (5) <a name="peak-calling"></a>Peak calling

Peak calling was performed using MACS2 (200bp each, ±100bp from the peak summit), followed by removal of peak regions overlapping any blacklisted regions.

<img src="https://github.com/chpngyu/pipeline-of-chip-seq/blob/master/images/motif_discovery_blacklist.png">

```Shell
#set global variables
export CURR=/home/user/human/ChIP-seq/
cd $CURR/GENE_ID/replicate1
export CONTROL="control1.bam control2.bam" # control string may include multiple files

# individual replicate, paired-end (PE) reads
macs2 callpeak -t unique_alignment_sorted_rd.bam -c $CONTROL -f BAMPE --gsize hs --outdir macs2

# individual replicate, single-read (SE) reads
macs2 callpeak -t unique_alignment_sorted_rd.bam -c $CONTROL -f BAM --gsize hs --outdir macs2

## merging two replicates, PE reads
macs2 callpeak -t replicate1.bam replicate2.bam -c $CONTROL $ CONTROL -f BAMPE --gsize hs --outdir macs2

# mergin two replicates, SE reads 
macs2 callpeak -t replicate1.bam replicate2.bam -c $CONTROL $ CONTROL -f BAM --gsize hs --outdir macs2
```


### (6) <a name="motif-discovery"></a>Motif discovery

Motif discovery was performed using MEME-chip as described in the manuscript, using the top 500 peaks to determine five motifs per analysis.

```Shell
#set global variables
export CURR=/home/user/human/ChIP-seq/
export GENOME=/home/user/human/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa
export BlackList=/home/human/blacklist/ENCFF023CZC_sorted.bed
export NoTopPeak=500 # number of top peaks to analyze

# navigate to replicate directory, extract peaks
cd $CURR/GENE_ID/replicate1/
awk 'BEGIN {WIDTH=100} {if($2<=WIDTH) print $1 "\t1\t" $2+WIDTH "\t" $4 "\t" $5; else print $1 "\t" $2-WIDTH "\t" $2+WIDTH "\t" $4 "\t" $5}' macs2/NA_summits.bed > extended_peaks.bed

# remove blacklist regions from peaks
bedtools subtract -a extended_peaks.bed -b $BlackList -A > extended_bk_removal_peaks.bed
sort -r -k5 -n extended_bk_removal_peaks.bed| head -n $NoTopPeak > top_peaks.bed
bedtools getfasta -bed top_peaks.bed -fi $GENOME > top_peaks.fa

# motif discovery
meme-chip -meme-nmotifs 5 top_peaks.fa
```


### (7) <a name="replicate-quality-control"></a>Replicate quality control and merging

Replicates were merged as described in the manuscript. Briefly, replicates were required to yield:

1. at least one position weight matrix (PWM) supported by >100 peaks and an E-value of <0.0001
2. and IDR score (-log[IDR]) ≥1.5 when compared to the other replicate of the same experiment

Experiments were required to have at least 2 passing replicates, and excluded otherwise.

```Shell
# normalizing scores for each peak
python score_quantile.py -s sample1_bk_removal.bed -g sample1 -o sample1_score.bed # repeat for all samples

# concatenating all scoring peaks and grouping peaks into clusters
cat sample1_score.bed > total.bed
cat sample2_score.bed >> total.bed #... for all samples

# sort and cluster
bedtools sort -i total.bed > temp.bed
bedtools cluster -i temp.bed > total_cluter.bed
rm *_score.bed total.bed temp.bed

#select top 500 peaks
python top_peak_sel.py -i total_cluter.bed -o top_combined.bed

#calculate motif position weight matrix (PWM) correlations
python correlation.py -m1 combined.meme -o motif_pcc.txt
<<<<<<< Updated upstream
python motif_cluster.py motif_pcc.txt occurrences.txt
```


### (8) <a name="representative-motif-selection"></a>Representative motif selection

**Single experiment**. 
For transcription factors represented by **one passing experiment** in the ENCODE data, replicates were merged and step 5 ([peak calling](#peak-calling)) was repeated using the merged replicates. Final representative motifs were called by repeating and 6 ([motif discovery](#motif-discovery)) using these merged-replicate peaks.

**!!TODO** *CP or CH please check that the BELOW section is correct*

**Multiple experiments**. 
For transcription factors represented by **more than one passing experiment** from **more than one biosample (*i.e.*, cell line or tissue type)** in the ENCODE data, one experiment was selected to represent each biosample. To do this, we first computed the similarities (KFV cos; Xu and Su 2010) between all pairs of motifs (*i.e.*, position weight matrices; PWMs) for all pairs of experiments utilizing different biosamples. For each biosample, the experiment was selected which shows the highest KFV cos value with another PWM from another experiment using a different biosample. 

For TFs with two biosamples, if neither biosample’s representative experiment had at least one PWM with cos>=0.80, we selected the experiment whose top PWM had the highest number of occurrences in peak regions. For TFs with >2 biosamples, experiments lacking at least one motif with KFV cos>=0.80 were excluded. Finally, using the representative experiments for all available biosamples, the ranking method described in the manuscript was used to select the top 500 peaks across all experiments. Final representative motifs were called by repeating step 6 ([motif discovery](#motif-discovery)) using these top peaks.

For transcription factors represented by **more than one passing experiment** from **only one biosample** in the ENCODE data, different experiments were not compared. Instead, we directly used the ranking method described in the manuscript to select the top 500 peaks across all experiments. Final representative motifs were called by repeating step 6 ([motif discovery](#motif-discovery)) using these top peaks.


## <a name="software-and-data"></a>Software and Data Versions Used
* Genome
	* Human: GRCh38
	* Mouse: GRCm38
* Blacklists
	* Amemiya et al. (2019)
* FASTQC v0.11.8
* TRIMMOMATIC v0.39 (Bolger et al. 2014)
* BOWTIE2 v2.3.5 (64-bit) (Langmead et al. 2012)
* MACS2 v2.1.2 (https://github.com/taoliu/MACS) (Zhang et al. 2008)
* MEME_CHIP v5.0.5 (Bailey et al. 2009)
* SAMTOOLS v1.9 (using htslib 1.9) (Li et al. 2009)
* BEDTOOLS v2.28.0 (Quinlan et al. 2010)
* IDR (https://github.com/nboley/idr) (Li et al. 2011)


## <a name="analysis-scripts"></a>Analysis Scripts

In addition to published tools, our analyses utilized the following custom scripts: 
(Yu: I'll add descriptions, Chen-Hao: Chen-Hao can help!)

1. `pwm.py`. **!!TODO: please check this is correct. Can the input have only one, or multiple PWMs?**
	* **Decription**: Trims PWMs from their ends until the remaining termini have an information content greater than or equal to some threshold value.
	* **Input**: (1) `-i`, PWMs in MEME format; (2) `--trim`, threshold information content; and (3) `-o` path/name of output file.
	* **Output**: A MEME file containing trimmed PWMs.
	* **Example**:

		`python3.6 pwm.py -i myPWM.meme --trim 0.3 -o myPWM_trimmed.meme`

2. `motif_cluster.py`. **!!TODO: please check this is correct.**
	* **Description**: Determines groups (clusters) of similar PWMs based on their correlations.
	* **Input**: (1) the motif_pcc.txt files produced by `correlation.py` (first argument, unnamed); and (2) **!!TODO??**
	* **Output**: Groups (clusters) of PWMs... **!!TODO??**. Printed to STOUT, here redirected to the file `motif_cluster.txt`.

		`python3.6 motif_cluster.py motif_pcc.txt occurrences.txt > motif_cluster.txt`

3. `consensus_pwm.py`. **!!TODO** please add

4. `correlation.py`. **!!TODO: please check this is correct.**
	* **Description**: Calculates similarity between PWMs in MEME format with user-defined correlation methods described in KFV.
	* **Input**: (1) `--method`, method for calculating correlation (cos, pcc, eucl, and kl); (2) `-m1`, PWMs in MEME format (if one PWM file, calculates pairwise correlations within the file; if two PWM files, calculates pairwise correlation between the two files); and (3) `-o` path/name of output file.
	* **Output**: a plain text file giving the similarity between each pair of motifs.
	* **Example**:

		`python3.6 correlation.py --method pcc -m1 PWM_trimmed.meme -o motif_pcc.txt`

5. `peak_motif_ranges.R`. **!!TODO** please add

6. `score_quantile.py`. **!!TODO** please add

7. `top_peak_sel.py`. **!!TODO** please add


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
	* `CCUT`: a value, the size (bp) to which peak regions should be trimmed for MEME-chip, e.g., `100`. This allows MEME-chip to examine the central region of the peaks for motifs while comparing to the flanking regions as a control for local sequence content. **!!TODO** please check this
	* `adaptor_SE`: path to file (FASTA format) containing sequencing adaptors for single-end (SE) experiments, *e.g.*, `TruSeq3-SE.fa`
	* `adaptor_PE`: path to file (FASTA format) containing sequencing adaptors for paired-end experiments, *e.g.*, `TruSeq3-PE-2.fa`
	* `NUM_TOP_PEAKS`: a value, the number of top peaks to consider, *e.g.*, 500
	* `MFOLD_MIN`: a value, the minimum fold depth enrichment required to call a peak for MACS2, *e.g.*, `5`. Not employed in our final analyses.
	* `MFOLD_MAX`: a value, the maximum fold depth enrichment allowed to call a peak for MACS2, *e.g.*, `50`. Not employed in our final analyses.
	* `blacklist`: path to file (BED format) containing a blacklist (excluded genome regions, including repetitive and low-complexity regions), *e.g.*, `ENCFF023CZC_sorted.bed`. See Amemiya *et al.* (2008).
	* `GENOME_FASTA`: path to file (FASTA format) containing the primary genome assembly for the organism of interest, *e.g.*, `Homo_sapiens.GRCh38.dna.primary_assembly.fa`.
	* `GENOME_IDX_PREFIX`: path to file (genome index) created using BOWTIE2, *e.g.*, `Homo_sapiens.GRCh38.dna.primary_assembly`. This can be accomplished using the following command: `bowtie2-build -f Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38.dna.primary_assembly`.
	
	The pipeline is shown below. Note that these tools use slightly different commands for single-end (SE) and paired-end (PE) reads (see script source code).

<img src="https://github.com/chpngyu/pipeline-of-chip-seq/blob/master/images/perl_pipeline.png">

3. `PWM_pipeline.pl`. IDR and PCC. **!!TODO** Chase will add description


## <a name="supplementary-data"></a>Supplementary Data
Supplementary data is available in the data folder within this repository. This includes the following files:

* File 1. Description. **!!TODO**
* File 2. Description. **!!TODO**
* ...

## <a name="acknowledgments"></a>Acknowledgments
The authors acknowledgme XYZ...


## <a name="citation"></a>Citation
When using this software, please refer to and cite:

>THE PAPER **!!TODO**

and this page:

>https://github.com/chpngyu/pipeline-of-chip-seq/


## <a name="references"></a>References
* Amemiya HM, Kundaje A, Boyle AP. 2019. <a target="_blank" href="https://www.nature.com/articles/s41598-019-45839-z">The ENCODE Blacklist: Identification of Problematic Regions of the Genome</a>. *Scientific Reports* **9**: 9354.
* Bailey TL, Boden M, Buske FA, *et al.* 2009. <a target="_blank" href="https://academic.oup.com/nar/article/37/suppl_2/W202/1135092">MEME SUITE: tools for motif discovery and searching</a>. *Nucleic Acids Research* **37**: W202-W208.
* Bolger AM, Lohse M, Usadel B. 2014. <a target="_blank" href="https://academic.oup.com/bioinformatics/article/30/15/2114/2390096">Trimmomatic: a flexible trimmer for Illumina sequence data</a>. *Bioinformatics* **30**(15): 2114-2120.
* Langmead B, Salzberg S. Langmead B, Salzberg SL. 2012. <a target="_blank" href="https://www.nature.com/articles/nmeth.1923">Fast gapped-read alignment with Bowtie 2</a>. *Nature Methods* **9**(4): 357-359.
* Li H, Handsaker B, Wysoker A, *et al.* 2009. <a target="_blank" href="https://pubmed.ncbi.nlm.nih.gov/19505943/">The Sequence Alignment/Map format and SAMtools</a>. *Bioinformatics* **25**(16): 2078-2079.
* Li Q, Brown JB, Huang H, Bickel PJ. 2011. <a target="_blank" href="https://projecteuclid.org/euclid.aoas/1318514284">Measuring reproducibility of high-throughput experiments</a>. *Annals of Applied Statistics* **5**(3): 1752--1779.
* Quinlan AR, Hall IM. 2010. <a target="_blank" href="https://academic.oup.com/bioinformatics/article/26/6/841/244688">BEDTools: a flexible suite of utilities for comparing genomic features</a>. *Bioinformatics* **26**(6): 841-842.
* Xu M, Su Z. 2010. <a target="_blank" href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2808352/">A Novel Alignment-Free Method for Comparing Transcription Factor Binding Site Motifs</a>. *PLoS ONE* **5**:e8797.
* Zhang Y, Liu T, Meyer CA, *et al.* 2008. <a target="_blank" href="https://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-9-r137">Model-based analysis of ChIP-Seq (MACS)</a>. *Genome Biology* **9**(9): R137.
