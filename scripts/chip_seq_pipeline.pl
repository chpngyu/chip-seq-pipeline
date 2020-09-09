#! /usr/bin/env perl

# DESCRIPTION: Perl script to automate the ChIP-seq pipeline through meme-chip

#########################################################################################
# EXAMPLE CALL. Must be called within a "Target" name directory, e.g., ~/GATA1/
#########################################################################################
# EXAMPLE1: chip_seq_pipeline.pl --ncpus=${NCPUS} --working_directory=${WORK} 2>&1 | tee output.txt
# EXAMPLE2: chip_seq_pipeline.pl --ncpus=4 --working_directory=/home/cnelson/ChIP-seq/human/GATA1 2>&1 | tee output.txt
#########################################################################################

# PIPELINE DESIGNERS: Chun-Ping Yu, Chen-Hao Kuo
# SCRIPT AUTHOR: Chase W. Nelson
# Copyright (C) 2019 Chase W. Nelson
# DATE CREATED: September 2019 - present

# CONTACT: cnelson@gate.sinica.edu.tw

# AFFILIATION1: Biodiversity Research Center, Academia Sinica, Taipei 115, Taiwan
# AFFILIATION2: Institute for Comparative Genomics, American Museum of Natural
#     History, New York, NY 10024, USA

use strict;
use Getopt::Long;


#########################################################################################
#########################################################################################
### MANUALLY SET GLOBAL VARIABLES ###
#########################################################################################
#########################################################################################
# Modules already loaded by job shell script that calls this one, run_chip_seq_pipeline_S*.sh
my $FASTQC = '/cnelson/bin/anaconda2/bin/fastqc'; # FastQC v0.11.8
my $TRIMMOMATIC = '/cnelson/bin/anaconda2/bin/trimmomatic'; # 0.39 
my $BOWTIE2 = '/cnelson/bin/anaconda2/bin/bowtie2'; # version 2.3.5 / 64-bit
my $SAMTOOLS = '/cnelson/bin/anaconda2/bin/samtools'; # samtools 1.9 / using htslib 1.9
my $BEDTOOLS = '/cnelson/bin/anaconda2/bin/bedtools'; # bedtools v2.28.0
my $MACS2 = '/cnelson/bin/anaconda2/bin/macs2'; # macs2 2.1.2
my $MACS2_gsize = 'hs'; # 'dm'; # 'hs';
my $PEAK_RADIUS = 100; # default=100 400
my $MIN_PEAK_SCORE = 13; # default, nothing # 20
my $MEME_CHIP = '/cnelson/bin/anaconda2/bin/meme-chip'; # 5.0.5
my $CCUT = 100; # 400; # default=100
my $adaptor_SE = '/home/cnelson/bin/Trimmomatic-0.39/adapters/TruSeq3-SE.fa';
my $adaptor_PE = '/home/cnelson/bin/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa';
my $NUM_TOP_PEAKS = 500; #100; #default=500;
my $MFOLD_MIN = 5; #2; # 5
my $MFOLD_MAX = 50; # 50
#my $blacklist = '/home/cnelson/ChIP-seq/blacklist/dm6-blacklist.v2.REMOVEchrPREFIX.bed'; # '/home/cnelson/ChIP-seq/blacklist/ENCFF023CZC_sorted.bed';
my $blacklist = '/home/cnelson/ChIP-seq/blacklist/ENCFF023CZC_sorted.bed';
#my $GENOME_FASTA = '/home/cnelson/genomes/Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa'; # '/home/cnelson/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa';
my $GENOME_FASTA = '/home/cnelson/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa';
#my $GENOME_IDX_PREFIX = '/home/cnelson/genomes/bowtie_index_dm6/Drosophila_melanogaster.BDGP6.22.dna.toplevel'; # '/home/cnelson/genomes/bowtie_2.3.5.1_index/Homo_sapiens.GRCh38.dna.primary_assembly'; # no extension; bowtie-indexed
my $GENOME_IDX_PREFIX = '/home/cnelson/genomes/bowtie_2.3.5.1_index/Homo_sapiens.GRCh38.dna.primary_assembly'; # no extension; bowtie-indexed
# the above has been accomplished using run_bowtie_indexing_dm6.sh
#########################################################################################
#########################################################################################


#########################################################################################
# INITIALIZE
print "\n################################################################################" .
"\n##                                                                            ##" .
"\n##                         ChIP-Seq Analysis Initiated!                       ##" .
"\n##                                                                            ##" .
"\n################################################################################\n";

#########################################################################################
# Get user input, if given. If a Boolean argument is passed, its value is 1; else undef

my $ncpus;
my $working_directory;

GetOptions( "ncpus=i" => \$ncpus,
			"working_directory=s" => \$working_directory )
			
	or die "\n### WARNING: problem with command-line arguments. Script terminated.\n\n";
	# If an argument is called as a flag, its value is 0; if not called, it's null

# Determine number of CPUs from input, if given; if not, use only 1
unless($ncpus) {
	$ncpus = 1;
}

# Determine working directory using input, if given
# Should be a particular TF directory. Examples:
# /home/cnelson/ChIP-seq/human/ATF2
# /home/cnelson/ChIP-seq/human/ATF3
# /home/cnelson/ChIP-seq/human/GATA1

if($working_directory) { # $working_directory passed as argument
	if(-d $working_directory) { # $working_directory is a real directory
		print "\n### DIRECTORY TO BEGIN ANALYSIS:\n";
		chomp($working_directory);
		print "WORKING_DIRECTORY=$working_directory\n";
	} else { # $working_directory is NOT a real directory
		print "\n### WARNING: DIRECTORY $working_directory DOESN'T EXIST. USING WORKING DIRECTORY TO BEGIN ANALYSIS:\n";
		$working_directory = `pwd`; # /home/cnelson/ChIP-seq/human/ATF2 
		chomp($working_directory);
		print "WORKING_DIRECTORY=$working_directory\n";
	}
} else {
	print "\n### BEGIN ANALYSIS IN WORKING DIRECTORY:\n";
	$working_directory = `pwd`; # /home/cnelson/ChIP-seq/human/ATF2 
	chomp($working_directory);
	print "WORKING_DIRECTORY=$working_directory\n";
}

my $WORK = $working_directory; # to keep nomenclature consistent
chdir("$WORK"); # no longer need to change manually


#########################################################################################
# Loop through subdirectories and run jobs
my %replicate2read_type;
my @wd_contents = glob "*";


print "\n################################################################################\n";
print "POTENTIAL DATASETS TO EXAMINE: @wd_contents\n";

print "\n################################################################################\n";
print "################################################################################\n";
print "TRIMMING, MAPPING, AND SORTING READS\n";

# LOOP ALL DATASETS (TFs and also controls)
DATASET: foreach my $dataset (@wd_contents) { # i.e., experiment ID
	
	if(-d "$WORK\/$dataset" && $dataset =~ /^ENC/) { # e.g., /home/cnelson/ChIP-seq/human/ENCSR000BKC
		
		print "\n################################################################################\n";
		print "################################################################################\n";
		print "Examining dataset $dataset\n";
		
		chdir("$WORK\/$dataset");
		
		my @wd_contents = glob "*"; 
		print "REPLICATES TO EXAMINE: @wd_contents\n";
		
		# LOOP ALL REPLICATES
		REPLICATE: foreach my $replicate (@wd_contents) { # i.e., isogenic replicate ID
	
			if(-d "$WORK\/$dataset\/$replicate" && $replicate =~ /^ENC/) { # e.g., /home/cnelson/ChIP-seq/human/ENCSR000BKC/ENCLB559JTU
				
				print "\n################################################################################\n";
				print "Examining replicate $replicate\n";
				
				chdir("$WORK\/$dataset\/$replicate");
				
				# GATHER AND ANALYZE ALL FASTQ FILES
				
				# Identify FASTQ file names; there will be 1 if SE and 2 if PE
				my @fastq_file_list = glob "*.fastq.gz";
				my $fastq_file_count = scalar(@fastq_file_list);
				print "FASTQ FILES TO EXAMINE: @fastq_file_list\n";
				
				# Make output directories; what does trimmomatic do if files already exist?
				unless(-d "$WORK\/$dataset\/$replicate\/QC") {
					mkdir("QC");
				}
				
				unless(-d "$WORK\/$dataset\/$replicate\/fastqc1") {
					mkdir("fastqc1");
				}
				
				unless(-d "$WORK\/$dataset\/$replicate\/fastqc2") {
					mkdir("fastqc2");
				}
				
				##########################################################################
				# S1: read quality and trimming
				##########################################################################
				
				##########################################################################
				# Run FASTQC BEFORE trimming (raw reads)
				foreach my $fastq_file (@fastq_file_list) {
					print "\n### fastq_file=$fastq_file\n";
					my $fastqc1_command = "$FASTQC --threads $ncpus --outdir $WORK\/$dataset\/$replicate\/fastqc1 $fastq_file 2>&1 | tee FASTQC_$fastq_file\.out";
					# NOTE that, with 'tee', no need to overwrite using >|, because tee overwrites automatically
					print "\n### fastqc1_command:\n$fastqc1_command\n";
					
					 # RUN FASTQC
				    `time $fastqc1_command`;
				}
				
				# RUN PAIRED- OR SINGLE-END ANALYSIS DEPENDING ON NUMBER OF FASTQ FILES
				if($fastq_file_count == 2) { # PAIRED-END ANALYSIS (if two inputs, or read_type =~ /^PE/)
					print "read_type=PE\n";
					$replicate2read_type{$replicate} = 'PE';
				    
				    my $R1 = $fastq_file_list[0];
				    my $R2 = $fastq_file_list[1];
				    
				    # DETERMINE READ LENGTH FOR R1 (SEQUENCE and QUALITIES are on lines 2 and 4 of .fastq)
				    my $first_read_seq_R1 = `zcat $R1 | head -n 2 | tail -n +2`; # `zcat $R1 | sed -n '2p'`;
				    chomp($first_read_seq_R1);
				    my $first_read_seq_len_R1 = length($first_read_seq_R1);
				    my $R1_first_read_qual = `zcat $R1 | head -n 4 | tail -n +4`; # `zcat $R1 | sed -n '4p'`;
				    chomp($R1_first_read_qual);
				    my $first_read_qual_len_R1 = length($first_read_seq_R1);
				    
				    # TRY THIS IN PERL:
				    #open(IN, "gunzip -c $file |") or die "gunzip $file: $!";
				    
				    my $min_len_R1;
				    if($first_read_seq_len_R1 == $first_read_qual_len_R1) { # should match
				    	print "first_read_seq_len_R1=$first_read_seq_len_R1\; ";
				    	
				    	if($first_read_seq_len_R1 < 25) {
				    		die "### WARNING: Read length too short or program bug. TERMINATED. ###\n";
				    	} elsif($first_read_seq_len_R1 <= 30) {
				    		$min_len_R1 = 25;
				    		#$min_len_R1 = $first_read_seq_len_R1 - 3;
				    	} elsif($first_read_seq_len_R1 <= 50) {
				    		$min_len_R1 = 30;
				    	} else {
				    		$min_len_R1 = 50;
				    	}
				    	
				  		print "min_len_R1=$min_len_R1\n";
				    } else {
				    	print "### WARNING: expected read and quality lengths unequal for R1.";
				    	next REPLICATE;
				    }
				    
				    # DETERMINE READ LENGTH FOR R2 (SEQUENCE and QUALITIES are on lines 2 and 4 of .fastq)
				    my $first_read_seq_R2 = `zcat $R2 | head -n 2 | tail -n +2`;
				    chomp($first_read_seq_R2);
				    my $first_read_seq_len_R2 = length($first_read_seq_R2);
				    my $R2_first_read_qual = `zcat $R2 | head -n 4 | tail -n +4`;
				    chomp($R2_first_read_qual);
				    my $first_read_qual_len_R2 = length($first_read_seq_R2);
				    
				    my $min_len_R2;
				    if($first_read_seq_len_R2 == $first_read_qual_len_R2) { # should match
				    	print "first_read_seq_len_R2=$first_read_seq_len_R2\; ";
				    	
				    	if($first_read_seq_len_R2 < 25) {
				    		die "### WARNING: Read length too short or program bug. TERMINATED. ###\n";
				    	} elsif($first_read_seq_len_R2 <= 30) {
				    		$min_len_R2 = 25;
				    		#$min_len_R2 = $first_read_seq_len_R2 - 3;
				    	} elsif($first_read_seq_len_R2 <= 50) {
				    		$min_len_R2 = 30;
				    	} else {
				    		$min_len_R2 = 50;
				    	}
				    	
				  		print "min_len_R2=$min_len_R2\n";
				    } else {
				    	print "### WARNING: expected read and quality lengths unequal for R2.";
				    	next REPLICATE;
				    }
				    
				    # Make sure the files agree
				    my $min_len;
				    if($min_len_R1 != $min_len_R2 || $first_read_seq_len_R1 != $first_read_seq_len_R2) {
				    	print "### WARNING: expected read lengths unequal between replicate.";
				    	next REPLICATE;
				    } else {
				    	$min_len = $min_len_R1;
				    }
				    
				    print "min_len=$min_len\n";
				    
				    ######################################################################
				    # Trimmomatic with paired-end reads
				    my $trimmomatic_command = "$TRIMMOMATIC PE -threads $ncpus $R1 $R2 " .
				        "$WORK\/$dataset\/$replicate\/QC\/$R1\_trimmed_paired.fastq.gz " .
				        "$WORK\/$dataset\/$replicate\/QC\/$R1\_trimmed_unpaired.fastq.gz " .
				        "$WORK\/$dataset\/$replicate\/QC\/$R2\_trimmed_paired.fastq.gz " .
				        "$WORK\/$dataset\/$replicate\/QC\/$R2\_trimmed_unpaired.fastq.gz " .
				        "ILLUMINACLIP:$adaptor_PE:2:40:12:8:true " .
				        "LEADING:10 SLIDINGWINDOW:4:15 MINLEN:$min_len 2>&1 | tee TRIMMOMATIC.out";
				    
				    print "\n### trimmomatic_command:\n$trimmomatic_command\n";
				    
				    # RUN TRIMMOMATIC
				    `time $trimmomatic_command`;
				    
				} elsif ($fastq_file_count == 1) { # SINGLE-END ANALYSIS (if read_type =~ /^SE/)
					print "read_type=SE\n";
					$replicate2read_type{$replicate} = 'SE';
					
				    my $READ = $fastq_file_list[0];
				    
				    # Check read length, SEQUENCE and QUALITIES are on lines 2 and 4
				    my $first_read_seq = `zcat $READ | head -n 2 | tail -n +2`;
				    chomp($first_read_seq);
				    my $first_read_seq_len = length($first_read_seq);
				    my $first_read_qual = `zcat $READ | head -n 4 | tail -n +4`;
				    chomp($first_read_qual);
				    my $first_read_qual_len = length($first_read_seq);
				    
				    my $min_len;
				    if($first_read_seq_len == $first_read_qual_len) { # should match
				    	print "Read length is $first_read_seq_len\; ";
				    	
				    	if($first_read_seq_len < 25) {
				    		die "### WARNING: Read length too short or program bug. TERMINATED. ###\n";
				    	} elsif($first_read_seq_len <= 30) {
				    		$min_len = 25;
				    		#$min_len = $first_read_seq_len - 3;
				    	} elsif($first_read_seq_len <= 50) {
				    		$min_len = 30;
				    	} else {
				    		$min_len = 50;
				    	}
				    	
				  		print "min length selected is $min_len\n";
				    } else {
				    	print "### WARNING: expected sequence and quality lengths unequal.";
				    	next REPLICATE;
				    }
				    
				    ######################################################################
				    # Trimmomatic with single-end reads
				    my $trimmomatic_command = "$TRIMMOMATIC SE -threads $ncpus $READ " .
				    	"$WORK\/$dataset\/$replicate\/QC\/$READ\_trimmed_single.fastq.gz " . # why this paired?
				    	"ILLUMINACLIP:$adaptor_SE:2:40:12 " .
				    	"LEADING:10 SLIDINGWINDOW:4:15 MINLEN:$min_len 2>&1 | tee TRIMMOMATIC.out";
				    
				    print "\n### trimmomatic_command:\n$trimmomatic_command\n";
				    
				    # RUN TRIMMOMATIC
				    `time $trimmomatic_command`;
				    
				} else {
					die "There were 0 or >2 reads files! TERMINATED.";
				}
				
				##########################################################################
				# FASTQC for read quality AFTER trimmomatic
				#my @trimmed_fastq_file_list = glob "./QC\/*trimmed_paired.fastq.gz";
				my @trimmed_fastq_file_list = glob "$WORK\/$dataset\/$replicate\/QC\/*.fastq.gz";
				
				foreach my $trimmed_fastq_file (@trimmed_fastq_file_list) {
					if($trimmed_fastq_file =~ /^.+\/(.+)(single|paired).fastq.gz$/) {
						$trimmed_fastq_file = "$1$2\.fastq\.gz";
						print "\n### trimmed_fastq_file=$trimmed_fastq_file\n"; # e.g., ./QC\/ENCFF000ZOO.fastq_trimmed_paired.fastq.gz
						#my $fastqc_command = "$FASTQC $trimmed_fastq_file -o fastqc";
						my $fastqc2_command = "$FASTQC --threads $ncpus --outdir $WORK\/$dataset\/$replicate\/fastqc2 ./QC/$trimmed_fastq_file 2>&1 | tee FASTQC_$trimmed_fastq_file\.out";
						print "\n### fastqc2_command:\n$fastqc2_command\n";
						
						 # RUN FASTQC
					    `time $fastqc2_command`;
					    #^I CONFIRM that the above command overwrites whatever file already existed by the same name (e.g., here, ENCFF000OMW.fastq.gz_trimmed_single.fastq.gz)
				    }
				}
				
				##########################################################################
				# S2a: read mapping and removal of unmapped and duplicate reads
				##########################################################################
				
				##########################################################################
				# BOWTIE: identify quality-trimmed reads and map them to genome
				if($replicate2read_type{$replicate} eq 'PE') { # $fastq_file_count == 2) { # PAIRED-END (PE)
					my $R1 = $fastq_file_list[0];
				    my $R2 = $fastq_file_list[1];
				    my $R1_trimmed_paired = "$WORK\/$dataset\/$replicate\/QC\/$R1\_trimmed_paired.fastq.gz";
				    my $R2_trimmed_paired = "$WORK\/$dataset\/$replicate\/QC\/$R2\_trimmed_paired.fastq.gz";
				    
				    #`rm $WORK\/$dataset\/$replicate\/alignment.sam`; # redundant if clobber with >|
				    my $bowtie2_command = "$BOWTIE2 --threads $ncpus -x $GENOME_IDX_PREFIX -1 $R1_trimmed_paired -2 $R2_trimmed_paired -S $WORK\/$dataset\/$replicate\/alignment.sam 2>&1 | tee BOWTIE2.out";
					print "\n### bowtie2_command:\n$bowtie2_command\n";
					`time $bowtie2_command`;
					
				} elsif ($replicate2read_type{$replicate} eq 'SE') { #$fastq_file_count == 1) { # SINGLE-END (SE)
					my $READ = $fastq_file_list[0];
				    my $READ_trimmed_single = "$WORK\/$dataset\/$replicate\/QC\/$READ\_trimmed_single.fastq.gz";
				    
				    #`rm $WORK\/$dataset\/$replicate\/alignment.sam`; # redundant if clobber with >|
				    my $bowtie2_command = "$BOWTIE2 --threads $ncpus -x $GENOME_IDX_PREFIX -U $READ_trimmed_single -S $WORK\/$dataset\/$replicate\/alignment.sam 2>&1 | tee BOWTIE2.out";
					print "\n### bowtie2_command:\n$bowtie2_command\n";
					`time $bowtie2_command`;
					
				} else {
					die "There were 0 or >2 reads files! TERMINATED.";
				}
				
				##########################################################################
				# S2b: SAMTOOLS to convert to BAM, SORT, and REMOVE DUPLICATES
				# Same procedure for SE and PE
				##########################################################################
				
				# Convert SAM file to BAM file. Other samtools commands require BAM input
				
				# BELOW, consider that an exhaustive select statement on Huxley would be:
				# #PBS -l select=1:ncpus=64:mem=256gb
				# However, should be fine to leave mem blank UNLESS the job not only HIGH MEMORY but also
				# LOW CPU.
				# Just make sure nothing goes wrong with parallelism using samtools sort. Leave -m blank at start.


				# SAMTOOLS VIEW
				# HELP
				# Usage: samtools view [options] <in.bam>|<in.sam>|<in.cram> [region ...]
				# -@: Number of **ADDITIONAL** threads to use [0]
				# -h: include header in SAM output
				# -F 268: only include reads WITHOUT the FLAGS in INT (=268) present [0]
				#     # 0x10c    268    UNMAP,MUNMAP,SECONDARY # THIS REMOVES NON-UNIQUELY MAPPED*** 
				#     # UNMAP: segment unmapped
				#     # MUNMAP: next segment in the template unmapped (?)
				#     # SECONDARY: secondary alignment
				# -q 5: only include reads with mapping quality >= INT (=5) [0]
				# -bS: output BAM; and input format SAM? "ignored (input format is auto-detected)"
				# -o: output file name [stdout]
				# N.B.: '-m' is NOT USED for memory with view!
				#EX1: samtools view -h -F 268 -q 5 -bS alignment.sam > unique_alignment.bam
				#EX2: samtools view -@ ${thread} -h -F 268 -q 5 -bS ${outdir}/Rep1_trim_map.sam -o ${outdir}/unique_alignment.bam
				my $samtools_view_command = "$SAMTOOLS view -@ " . ($ncpus - 1) . " -h -F 268 -q 5 -bS $WORK\/$dataset\/$replicate\/alignment.sam -o $WORK\/$dataset\/$replicate\/unique_alignment.bam 2>&1 | tee samtools_view.out";
				print "\n### samtools_view_command:\n$samtools_view_command\n";
				`time $samtools_view_command`;
				
				
				# SAMTOOLS SORT
				# HELP
				# Usage: samtools sort [options...] [in.bam]
				# -m: Set maximum memory per thread; suffix K/M/G recognized [768M]
				#     # N.B.: Huxley has 1024 virtual cores and 4G/VCore
				# -@: Number of **ADDITIONAL** threads to use [0]
				# -o: Write final output to FILE rather than standard output
				#EX1: samtools sort unique_alignment.bam -o unique_alignment_sorted.bam
				#EX2: samtools sort -@ ${thread} ${outdir}/unique_alignment.bam -o ${outdir}/unique_alignment_sorted.bam |& tee ${base}/${AC}/Rep1/logs/Rep1_sort.log
				# may add in '-m 4G' is fails due to memory
				my $samtools_sort_command = "$SAMTOOLS sort -@ " . ($ncpus - 1) . " $WORK\/$dataset\/$replicate\/unique_alignment.bam -o $WORK\/$dataset\/$replicate\/unique_alignment_sorted.bam 2>&1 | tee samtools_sort.out";
				print "\n### samtools_sort_command:\n$samtools_sort_command\n";
				`time $samtools_sort_command`;
				
				
				# SAMTOOLS RMDUP
				# HELP # IF MORE THAN ONE READ MAPS TO SAME COORDINATES, THIS KEEPS ONLY THE ONE OF HIGHEST QUALITY: https://www.htslib.org/doc/samtools-rmdup.1.html
				# samtools rmdup [-sS] <input.srt.bam> <output.bam>
				#EX1: samtools rmdup unique_alignment_sorted.bam unique_alignment_sorted_rd.bam
				#EX2: samtools rmdup ${outdir}/unique_alignment_sorted.bam ${outdir}/unique_alignment_sorted_rd.bam |& tee ${base}/${AC}/Rep1/logs/Rep1_redup.log
				my $samtools_rmdup_command = "$SAMTOOLS rmdup $WORK\/$dataset\/$replicate\/unique_alignment_sorted.bam $WORK\/$dataset\/$replicate\/unique_alignment_sorted_rd.bam 2>&1 | tee samtools_rmdup.out";
				print "\n### samtools_rmdup_command:\n$samtools_rmdup_command\n";
				`time $samtools_rmdup_command`;
				#^ I CONFIRM that the above OVERWRITES whatever file already existed by the same name (here, unique_alignment_sorted_rd.bam)
				
				# RM BAM FILES
				# Remove unnecessary files to save space
				#EX1: rm alignment.sam unique_alignment.bam unique_alignment_sorted.bam
				#EX2: rm ${outdir}/unique_alignment.bam ${outdir}/unique_alignment_sorted.bam
				my $rm_bam_files_command = "rm $WORK\/$dataset\/$replicate\/alignment.sam $WORK\/$dataset\/$replicate\/unique_alignment.bam $WORK\/$dataset\/$replicate\/unique_alignment_sorted.bam  2>&1 | tee rm_bam_files.out";
				print "\n### rm_bam_files_command:\n$rm_bam_files_command\n";
				`time $rm_bam_files_command`;
				
				
				### DONE
				print "We reached the end of replicate $replicate.\n";
				
			} # end this replicate (passed test)
		} # end all replicates
		
		print "We reached the end of dataset $dataset.\n";
		
	} # end this dataset (passed test)
} # end all datasets (experiments)


print "\n################################################################################\n";
print "################################################################################\n";
print "PEAK CALLING FOR INDIVIUAL REPLICATES\n";

# LOOP ALL INDIVIDUAL REPLICATES (TFs BUT NOT controls)
DATASET: foreach my $dataset (@wd_contents) { # i.e., experiment ID
	
	if(-d "$WORK\/$dataset" && $dataset =~ /^ENC/) { # e.g., /home/cnelson/ChIP-seq/human/ENCSR000BKC
		
		print "\n################################################################################\n";
		print "################################################################################\n";
		print "PEAK CALLING FOR INDIVIDUAL REPLICATES IN EXPERIMENT=$dataset\n"; # always only one per TF, so far
		
		chdir("$WORK\/$dataset");
		
		my @wd_contents = glob "*"; 
		print "REPLICATES TO EXAMINE: @wd_contents\n";
		
		# LOOP ALL REPLICATES
		REPLICATE: foreach my $replicate (@wd_contents) { # i.e., isogenic replicate ID
	
			if(-d "$WORK\/$dataset\/$replicate" && $replicate =~ /^ENC/) { # e.g., /home/cnelson/ChIP-seq/human/ENCSR000BKC/ENCLB559JTU # && ! $replicate =~ /control$/ didn't work
				unless($replicate =~ /control/) { # NOT controls
					print "\n################################################################################\n";
					print "Examining replicate $replicate (not a control)\n";
					
					chdir("$WORK\/$dataset\/$replicate");
					
					# ANALYZE BAM FILES: 1 replicate with ≥1 controls, which are shared by any other replicate(s)
					
					# There should always be one BAM file, by the name 'unique_alignment_sorted_rd.bam'
					print "BAM FILE TO EXAMINE: unique_alignment_sorted_rd.bam\n";
					
					# Make output directory? What does MACS2 do if files already exist?
					unless(-d "$WORK\/$dataset\/$replicate\/macs2") {
						mkdir("$WORK\/$dataset\/$replicate\/macs2");
					}
					
	
					##########################################################################
					# S3: MACS2 for PEAK CALLING
					# DIFFERENT procedure for SE and PE
					##########################################################################
					
					# This required a CIRCUS of Python 2.7 installation MADNESS
					
					# Gather names of control files, which are .bam files that have been placed in $WORK/$dataset/controls/
					my @control_bam_files = glob "$WORK\/$dataset\/*_control\/unique_alignment_sorted_rd.bam";
					#my $CONTROL_STRING = '';
					#foreach my $control_fastq_file (@control_bam_files) {
					#	if($control_fastq_file =~ /^.+\/(.+).bam$/) {
					#		#$control_fastq_file = "$1\.bam"; # we want the whole path this time
					#		$CONTROL_STRING .= "$control_fastq_file "; # space at the end!
					#	}
					#}
					
					my $CONTROL_STRING = "@control_bam_files"; # this should do it
					$CONTROL_STRING =~ s/\s+$//; # Remove trailing whitespace (hopefully just one space)
					print "\n### CONTROL_STRING:\n$CONTROL_STRING\n";
					
					# MACS2
					# HELP
					# Usage: macs2 callpeak [-h] -t TFILE [TFILE ...] [-c [CFILE [CFILE ...]]]...
					# -t: ChIP-seq treatment file. If multiple files are given as '-t A B C', then they will all be read and pooled together. REQUIRED.
					# -c: Control file. If multiple files are given as '-c A B C', they will be pooled to estimate ChIP-seq background noise.
					# -f: Format of tag file... BAM, ... The default AUTO option will let MACS decide.
					# --gsize: Effective genome size. It can be 1.0e+9 or 1000000000, or shortcuts:'hs' for human (2.7e9), 'mm' for mouse (1.87e9), ... Default:hs
					# --outdir: If specified all output files will be written to that directory. Default: the current working directory
					
					#	-q QVALUE, --qvalue QVALUE
					#	Minimum FDR (q-value) cutoff for peak detection.
					#	DEFAULT: 0.05. -q, and -p are mutually exclusive.
                    
					if($replicate2read_type{$replicate} eq 'PE') { # $fastq_file_count == 2) { # PAIRED-END (PE)
						
						# Call peaks for this replicate ### WILL THE BELOW CREATE A macs2 DIRECTORY??
						
						# INDIVIDUAL REPLICATE
					    my $macs2_command = "$MACS2 callpeak -t $WORK\/$dataset\/$replicate\/unique_alignment_sorted_rd.bam -c $CONTROL_STRING -f BAMPE --gsize $MACS2_gsize --mfold $MFOLD_MIN $MFOLD_MAX --outdir $WORK\/$dataset\/$replicate\/macs2 2>&1 | tee MACS2.out";
						print "\n### macs2_command:\n$macs2_command\n";
						`time $macs2_command`;
						
					} elsif ($replicate2read_type{$replicate} eq 'SE') { # $fastq_file_count == 1) { # SINGLE-END (SE)
						
						# Call peaks for this replicate
						
						# INDIVIDUAL REPLICATE
					    my $macs2_command = "$MACS2 callpeak -t $WORK\/$dataset\/$replicate\/unique_alignment_sorted_rd.bam -c $CONTROL_STRING -f BAM --gsize $MACS2_gsize --mfold $MFOLD_MIN $MFOLD_MAX --outdir $WORK\/$dataset\/$replicate\/macs2 2>&1 | tee MACS2.out";
						print "\n### macs2_command:\n$macs2_command\n";
						`time $macs2_command`;
						
						
					} else {
						die "There were 0 or >2 reads files! TERMINATED.";
					}
					
					
					
					##########################################################################
					# S4: MOTIF DISCOVERY
					# Same procedure for SE and PE
					##########################################################################
					
					### EXPLANATION:
					# awk: remember "AWKward tables"; awk is for tabular data
					#     -awk processes input one record at a time: $0 stores the entire record (line); $1 stores the first field (column); $2 stores the second field; etc.
					#     -Usage: awk pattern { action }
					#     -"pattern" is evaluated like an IF statement; the ACTION is performed only if the PATTERN evaluates to TRUE for a RECORD
					#     -"BEGIN" specifies what to do BEFORE the first record is read, e.g., define a variable.
					# BED FORMAT uses 0-based positions, half-closed, half-open: [start, end). Length is therefore = end - start
					#     -$1 is the sequence name
					#     -$2 is the range start (0-based, closed)
					#     -$3 is the range end (0-based, open)
					#     -($3 - $2) is length
					# THIS BED FILE, produced by macs2 callpeaks, is described:
					#NAME_summits.bed is in BED format, which contains the peak summits locations for every peaks. 
					#The 5th column in this file is the same as NAME_peaks.narrowPeak. 
					#If you want to find the motifs at the binding sites, this file is recommended. 
					#The file can be loaded directly to UCSC genome browser. 
					#Remove the beginning track line if you want to analyze it by other tools.
					# NAME_peaks.narrowPeak 5th column is described:
					#5th: integer score for display. It's calculated as int(-10*log10pvalue) or int(-10*log10qvalue) 
					#depending on whether -p (pvalue) or -q (qvalue) is used as score cutoff. 
					#Please note that currently this value might be out of the [0-1000] range defined in UCSC Encode narrowPeak format. 
					#You can let the value saturated at 1000 (i.e. p/q-value = 10^-100) by using the following 1-liner awk: 
					#awk -v OFS="\t" '{$5=$5>1000?1000:$5} {print}' NAME_peaks.narrowPeak
					# https://github.com/taoliu/MACS/
					
					# TOP OF: ATF3/ENCSR000BKC/ENCLB559JTU/macs2/NA_summits.bed
					#1	629938  629939  NA_peak_1       26.02448
					#1	634028  634029  NA_peak_2       35.47879
					#1	778706  778707  NA_peak_3       46.39153
					#1	1000913 1000914 NA_peak_4       7.42275
					#my $awk_on_macs2_command = "awk 'BEGIN {WIDTH=$PEAK_RADIUS} {if(\$2<=WIDTH) print \$1 \"\\t1\\t\" \$2+WIDTH \"\\t\" \$4 \"\\t\" \$5; else print \$1 \"\\t\" \$2-WIDTH \"\\t\" \$2+WIDTH \"\\t\" \$4 \"\\t\" \$5}' $WORK\/$dataset\/$replicate\/macs2/NA_summits.bed > $WORK\/$dataset\/$replicate\/extended_peaks.bed"; # COMEBACK to figure out how to tee this, just the error
					my $awk_on_macs2_command = "awk 'BEGIN {WIDTH=$PEAK_RADIUS} {if(\$2<=WIDTH) print \$1 \"\\t1\\t\" \$2+WIDTH \"\\t\" \$4 \"\\t\" \$5; else print \$1 \"\\t\" \$2-WIDTH \"\\t\" \$2+WIDTH \"\\t\" \$4 \"\\t\" \$5}' $WORK\/$dataset\/$replicate\/macs2/NA_summits.bed | awk '\$5 \>\= $MIN_PEAK_SCORE' > $WORK\/$dataset\/$replicate\/extended_peaks.bed"; # COMEBACK to figure out how to tee this, just the error
					# Verified the above produces the correct command string in "test_perl_code.pl": 
					# awk 'BEGIN {WIDTH=100} 
					#     {if($2<=WIDTH) print $1 "\t1\t" $2+WIDTH "\t" $4 "\t" $5; 
					#     else print $1 "\t" $2-WIDTH "\t" $2+WIDTH "\t" $4 "\t" $5}' 
					#     some_dir/macs2/NA_summits.bed > extended_peaks.bed
					# -The point of this juts appears to be extended the summits by 100 sites on either side,
					#     making sure not to go before site 1.
					print "\n### awk_on_macs2_command:\n$awk_on_macs2_command\n";
					`time $awk_on_macs2_command`;
					# Outputs: $1=seq_name, $2=extended_start, $3=extended_end, $4=peak_ID, $5=peak_score
					
					### NEW: FILTER IF SCORE ISN'T HIGH ENOUGH
					###
					
					# Remove peaks in blacklist with bedtools
					my $bedtools_subtract_command = "$BEDTOOLS subtract -a $WORK\/$dataset\/$replicate\/extended_peaks.bed -b $blacklist -A > $WORK\/$dataset\/$replicate\/extended_bk_removal_peaks.bed";
					print "\n### bedtools_subtract_command:\n$bedtools_subtract_command\n";
					`time $bedtools_subtract_command`;
					
					# Sort peaks by highest-to-lowest score and keep top 500 peaks only
					# -r reverse the output, sorting highest to lowest 
					# -k specifies which column(s) to sort by (here, 5, which is scores)
					my $sort_peaks_command = "sort -r -k5 -n $WORK\/$dataset\/$replicate\/extended_bk_removal_peaks.bed | head -n $NUM_TOP_PEAKS > $WORK\/$dataset\/$replicate\/top_peaks.bed";
					print "\n### sort_peaks_command:\n$sort_peaks_command\n";
					`time $sort_peaks_command`;
					
					# Get FASTA sequences of 500 top peaks
					my $bedtools_getfasta_command = "$BEDTOOLS getfasta -bed $WORK\/$dataset\/$replicate\/top_peaks.bed -fi $GENOME_FASTA > $WORK\/$dataset\/$replicate\/top_peaks.fa";
					print "\n### bedtools_getfasta_command:\n$bedtools_getfasta_command\n";
					`time $bedtools_getfasta_command`;
					
					# Motif discovery with meme-chip
					my $meme_chip_dir = 'meme_chip';
					#my $meme_chip_dir_num = 1;
					#while(-d "$WORK\/$dataset\/$replicate\/$meme_chip_dir") {
					#	$meme_chip_dir = 'meme_chip_' . $meme_chip_dir_num;
					#	$meme_chip_dir_num++;
					#}
					
					# Remove the directory if it already exists
					if(-d "$WORK\/$dataset\/$replicate\/$meme_chip_dir") {
						`rm -r $WORK\/$dataset\/$replicate\/$meme_chip_dir`;
					}
					
					my $meme_chip_command = "$MEME_CHIP -meme-nmotifs 5 -ccut $CCUT $WORK\/$dataset\/$replicate\/top_peaks.fa -o $WORK\/$dataset\/$replicate\/$meme_chip_dir 2>&1 | tee MEME_CHIP.out"; # needs to create a new, non-existent directory
					print "\n### meme_chip_command:\n$meme_chip_command\n";
					`time $meme_chip_command`;
					
					### DONE with this INDIVIDUAL replicate
					print "We reached the end of replicate $replicate.\n";
					
				} # this replicate wasn't a control	
			} # end this replicate (passed test)
		} # end all replicates
		
		print "We reached the end of dataset $dataset.\n";
		
	} # end this dataset (passed test)
} # end all datasets (experiments)


print "\n################################################################################\n";
print "################################################################################\n";
print "MERGED PEAK CALLING FOR EXPERIMENTS (MERGED)\n";

##########################################################################################
# MERGED PEAK CALLING: MERGED PEAKS FOR EACH EXPERIMENT
MERGED: foreach my $dataset (@wd_contents) { # i.e., experiment ID
	
	if(-d "$WORK\/$dataset" && $dataset =~ /^ENC/) { # e.g., /home/cnelson/ChIP-seq/human/ENCSR000BKC
		
		print "\n################################################################################\n";
					print "Examining experiment $dataset (not a control)\n";
		
		chdir("$WORK\/$dataset");
		
		my @wd_contents = glob "*"; 
		print "REPLICATES TO EXAMINE: @wd_contents\n";
		
		# LOOP ALL REPLICATES, ADD THEIR BAM FILES TO THE REPLICATE OR CONTROL LIST
		my @bam_files = glob "$WORK\/$dataset\/ENC*\/unique_alignment_sorted_rd.bam";
		my @control_bam_files;
		my @replicate_bam_files;
		
		foreach my $bam_file (@bam_files) {
			if($bam_file =~ /_control/) {
				push(@control_bam_files, $bam_file);
			} else {
				push(@replicate_bam_files, $bam_file);
			}
		}
		
		# Sort files
		@control_bam_files = sort(@control_bam_files);
		@replicate_bam_files = sort(@replicate_bam_files);
		
		my $CONTROL_STRING = "@control_bam_files"; # this should do it
		$CONTROL_STRING =~ s/\s+$//; # Remove trailing whitespace (hopefully just one space)
		print "\n### CONTROL_STRING:\n$CONTROL_STRING\n";
		
		# Make output directory for MERGED results
		unless(-d "$WORK\/$dataset\/merged") {
			mkdir("$WORK\/$dataset\/merged");
		}
		
		chdir("$WORK\/$dataset\/merged");
		
		my $bam_file_count = scalar(@replicate_bam_files);
		
		# LOOP ALL BAM FILES
		for(my $i = 0; $i < scalar(@replicate_bam_files); $i++) {
			my $bam_file_i = $replicate_bam_files[$i];
			my @bam_file_group = ($bam_file_i);
			
			for(my $j = $i + 1; $j < scalar(@replicate_bam_files); $j++) {
				my $bam_file_j = $replicate_bam_files[$j];
				push(@bam_file_group, $bam_file_j);
				
				
				##########################################################################
				##########################################################################
				# Do it for this unique pair (i/j, usually i/i+1)
				# The ≥1 controls are shared by any/all other replicate(s)
				##########################################################################
				##########################################################################
				
				print "### BAM FILES TO EXAMINE:\n" .
					"bam1=$bam_file_i\n" . 
					"bam2=$bam_file_j\n";
					
				# Obtain replicate ID's
				my $replicate_i_ID = $bam_file_i; # e.g., $WORK\/$dataset\/ENC*\/unique_alignment_sorted_rd.bam
				$replicate_i_ID =~ s/^$WORK\/$dataset\///; # trim leading
				$replicate_i_ID =~ s/\/unique_alignment_sorted_rd.bam\s*$//; # trim trailing
				
				# Obtain replicate ID's
				my $replicate_j_ID = $bam_file_j; # e.g., $WORK\/$dataset\/ENC*\/unique_alignment_sorted_rd.bam
				$replicate_j_ID =~ s/^$WORK\/$dataset\///; # trim leading
				$replicate_j_ID =~ s/\/unique_alignment_sorted_rd.bam\s*$//; # trim trailing
				
				
				##########################################################################
				# S3: MACS2 for PEAK CALLING
				# DIFFERENT procedure for SE and PE
				##########################################################################
				
				# Python 2.7
				
				# MACS2
				if($replicate2read_type{$replicate_i_ID} eq 'PE' && $replicate2read_type{$replicate_j_ID} eq 'PE') { # $fastq_file_count == 2) { # PAIRED-END (PE)
					
					# REPLICATE PAIR (PE)
				    my $macs2_command = "$MACS2 callpeak -t $bam_file_i $bam_file_j -c $CONTROL_STRING -f BAMPE --gsize $MACS2_gsize --mfold $MFOLD_MIN $MFOLD_MAX --outdir $WORK\/$dataset\/merged\/macs2_$replicate_i_ID\_$replicate_j_ID 2>&1 | tee $WORK\/$dataset\/merged\/macs2_$replicate_i_ID\_$replicate_j_ID\/MACS2.out";
					print "\n### macs2_command:\n$macs2_command\n";
					`time $macs2_command`;
					
				} elsif ($replicate2read_type{$replicate_i_ID} eq 'SE' && $replicate2read_type{$replicate_j_ID}) { # $fastq_file_count == 1) { # SINGLE-END (SE)
					
					# REPLICATE PAIR (SE)
				    my $macs2_command = "$MACS2 callpeak -t $bam_file_i $bam_file_j -c $CONTROL_STRING -f BAM --gsize $MACS2_gsize --mfold $MFOLD_MIN $MFOLD_MAX --outdir $WORK\/$dataset\/merged\/macs2_$replicate_i_ID\_$replicate_j_ID 2>&1 | tee $WORK\/$dataset\/merged\/macs2_$replicate_i_ID\_$replicate_j_ID\/MACS2.out";
					print "\n### macs2_command:\n$macs2_command\n";
					`time $macs2_command`;
					
				} else {
					die "There were 0 or >2 reads files or conflicting read types! TERMINATED.";
				}
				
				
				##########################################################################
				# S4: MOTIF DISCOVERY
				# Same procedure for SE and PE
				##########################################################################
				
				#my $awk_on_macs2_command = "awk 'BEGIN {WIDTH=$PEAK_RADIUS} {if(\$2<=WIDTH) print \$1 \"\\t1\\t\" \$2+WIDTH \"\\t\" \$4 \"\\t\" \$5; else print \$1 \"\\t\" \$2-WIDTH \"\\t\" \$2+WIDTH \"\\t\" \$4 \"\\t\" \$5}' $WORK\/$dataset\/merged\/macs2_$replicate_i_ID\_$replicate_j_ID\/NA_summits.bed > $WORK\/$dataset\/merged\/macs2_$replicate_i_ID\_$replicate_j_ID\/extended_peaks.bed";
				my $awk_on_macs2_command = "awk 'BEGIN {WIDTH=$PEAK_RADIUS} {if(\$2<=WIDTH) print \$1 \"\\t1\\t\" \$2+WIDTH \"\\t\" \$4 \"\\t\" \$5; else print \$1 \"\\t\" \$2-WIDTH \"\\t\" \$2+WIDTH \"\\t\" \$4 \"\\t\" \$5}' $WORK\/$dataset\/merged\/macs2_$replicate_i_ID\_$replicate_j_ID\/NA_summits.bed | awk '\$5 \>\= $MIN_PEAK_SCORE' > $WORK\/$dataset\/merged\/macs2_$replicate_i_ID\_$replicate_j_ID\/extended_peaks.bed";
				
				print "\n### awk_on_macs2_command:\n$awk_on_macs2_command\n";
				`time $awk_on_macs2_command`;
				# Outputs: $1=seq_name, $2=extended_start, $3=extended_end, $4=peak_ID, $5=peak_score
				
				# Remove peaks in blacklist with bedtools
				my $bedtools_subtract_command = "$BEDTOOLS subtract -a $WORK\/$dataset\/merged\/macs2_$replicate_i_ID\_$replicate_j_ID\/extended_peaks.bed -b $blacklist -A > $WORK\/$dataset\/merged\/macs2_$replicate_i_ID\_$replicate_j_ID\/extended_bk_removal_peaks.bed";
				print "\n### bedtools_subtract_command:\n$bedtools_subtract_command\n";
				`time $bedtools_subtract_command`;
				
				# Sort peaks by highest-to-lowest score and keep top 500 peaks only
				my $sort_peaks_command = "sort -r -k5 -n $WORK\/$dataset\/merged\/macs2_$replicate_i_ID\_$replicate_j_ID\/extended_bk_removal_peaks.bed | head -n $NUM_TOP_PEAKS > $WORK\/$dataset\/merged\/macs2_$replicate_i_ID\_$replicate_j_ID\/top_peaks.bed";
				print "\n### sort_peaks_command:\n$sort_peaks_command\n";
				`time $sort_peaks_command`;
				
				# Get FASTA sequences of 500 top peaks
				my $bedtools_getfasta_command = "$BEDTOOLS getfasta -bed $WORK\/$dataset\/merged\/macs2_$replicate_i_ID\_$replicate_j_ID\/top_peaks.bed -fi $GENOME_FASTA > $WORK\/$dataset\/merged\/macs2_$replicate_i_ID\_$replicate_j_ID\/top_peaks.fa";
				print "\n### bedtools_getfasta_command:\n$bedtools_getfasta_command\n";
				`time $bedtools_getfasta_command`;
				
				# Motif discovery with meme-chip
				my $meme_chip_dir = 'meme_chip';
				#my $meme_chip_dir_num = 1;
				#while(-d "$WORK\/$dataset\/merged\/macs2_$replicate_i_ID\_$replicate_j_ID\/$meme_chip_dir") {
				#	$meme_chip_dir = 'meme_chip_' . $meme_chip_dir_num;
				#	$meme_chip_dir_num++;
				#}
				
				# Remove the directory if it already exists
				if(-d "$WORK\/$dataset\/merged\/macs2_$replicate_i_ID\_$replicate_j_ID\/$meme_chip_dir") {
					`rm -r $WORK\/$dataset\/merged\/macs2_$replicate_i_ID\_$replicate_j_ID\/$meme_chip_dir`;
				}
				
				my $meme_chip_command = "$MEME_CHIP -meme-nmotifs 5 -ccut $CCUT $WORK\/$dataset\/merged\/macs2_$replicate_i_ID\_$replicate_j_ID\/top_peaks.fa -o $WORK\/$dataset\/merged\/macs2_$replicate_i_ID\_$replicate_j_ID\/$meme_chip_dir 2>&1 | tee $WORK\/$dataset\/merged\/macs2_$replicate_i_ID\_$replicate_j_ID\/MEME_CHIP.out"; # needs to create a new, non-existent directory
				print "\n### meme_chip_command:\n$meme_chip_command\n";
				`time $meme_chip_command`;
				
				### DONE with this PAIR
				print "We reached the end of this pair.\n";
				
				
				##########################################################################
				##########################################################################
				# Do it for this unique GROUP, if it is larger than a pair (i/j/j+1/j+2/...)
				##########################################################################
				##########################################################################
				if(@bam_file_group > 2) {
					
					# Print files, obtain ID's, and determine ready types concurrently
					my @current_replicate_IDs;
					my $all_PE = 'TRUE';
					my $all_SE = 'TRUE';
					
					print "### BAM FILES TO EXAMINE:\n";
					for(my $index = 0; $index < @bam_file_group; $index++) {
						"bam$index\=$bam_file_group[$index]\n";
						my $this_replicate_ID = $bam_file_group[$index]; # e.g., $WORK\/$dataset\/ENC*\/unique_alignment_sorted_rd.bam
						$this_replicate_ID =~ s/^$WORK\/$dataset\///; # trim leading
						$this_replicate_ID =~ s/\/unique_alignment_sorted_rd.bam\s*$//; # trim trailing
						push(@current_replicate_IDs, $this_replicate_ID);
						
						if($replicate2read_type{$this_replicate_ID} eq 'PE') {
							$all_SE = 'FALSE';
						} elsif($replicate2read_type{$this_replicate_ID} eq 'SE') {
							$all_PE = 'FALSE';
						}
					}
					
					# For the replicate files string
					my $replicates_underscored = join('_', @current_replicate_IDs);
					
					######################################################################
					# S3: MACS2 for PEAK CALLING
					# DIFFERENT procedure for SE and PE
					######################################################################
					
					# Python 2.7
					
					# MACS2
					if($all_PE eq 'TRUE') { # $fastq_file_count == 2) { # PAIRED-END (PE)
						
						# REPLICATE PAIR (PE)
					    my $macs2_command = "$MACS2 callpeak -t @bam_file_group -c $CONTROL_STRING -f BAMPE --gsize $MACS2_gsize --mfold $MFOLD_MIN $MFOLD_MAX --outdir $WORK\/$dataset\/merged\/macs2_$replicates_underscored 2>&1 | tee $WORK\/$dataset\/merged\/macs2_$replicates_underscored\/MACS2.out";
						print "\n### macs2_command:\n$macs2_command\n";
						`time $macs2_command`;
						
					} elsif ($all_SE eq 'TRUE') { # $fastq_file_count == 1) { # SINGLE-END (SE)
						
						# REPLICATE PAIR (SE)
					    my $macs2_command = "$MACS2 callpeak -t @bam_file_group -c $CONTROL_STRING -f BAM --gsize $MACS2_gsize --mfold $MFOLD_MIN $MFOLD_MAX --outdir $WORK\/$dataset\/merged\/macs2_$replicates_underscored 2>&1 | tee $WORK\/$dataset\/merged\/macs2_$replicates_underscored\/MACS2.out";
						print "\n### macs2_command:\n$macs2_command\n";
						`time $macs2_command`;
						
					} else {
						die "There were 0 or >2 reads files or conflicting read types! TERMINATED.";
					}
					
					
					##########################################################################
					# S4: MOTIF DISCOVERY
					# Same procedure for SE and PE
					##########################################################################
					
					#my $awk_on_macs2_command = "awk 'BEGIN {WIDTH=$PEAK_RADIUS} {if(\$2<=WIDTH) print \$1 \"\\t1\\t\" \$2+WIDTH \"\\t\" \$4 \"\\t\" \$5; else print \$1 \"\\t\" \$2-WIDTH \"\\t\" \$2+WIDTH \"\\t\" \$4 \"\\t\" \$5}' $WORK\/$dataset\/merged\/macs2_$replicates_underscored\/NA_summits.bed > $WORK\/$dataset\/merged\/macs2_$replicates_underscored\/extended_peaks.bed";
					my $awk_on_macs2_command = "awk 'BEGIN {WIDTH=$PEAK_RADIUS} {if(\$2<=WIDTH) print \$1 \"\\t1\\t\" \$2+WIDTH \"\\t\" \$4 \"\\t\" \$5; else print \$1 \"\\t\" \$2-WIDTH \"\\t\" \$2+WIDTH \"\\t\" \$4 \"\\t\" \$5}' $WORK\/$dataset\/merged\/macs2_$replicates_underscored\/NA_summits.bed | awk '\$5 \>\= $MIN_PEAK_SCORE' > $WORK\/$dataset\/merged\/macs2_$replicates_underscored\/extended_peaks.bed";
					
					print "\n### awk_on_macs2_command:\n$awk_on_macs2_command\n";
					`time $awk_on_macs2_command`;
					# Outputs: $1=seq_name, $2=extended_start, $3=extended_end, $4=peak_ID, $5=peak_score
					
					# Remove peaks in blacklist with bedtools
					my $bedtools_subtract_command = "$BEDTOOLS subtract -a $WORK\/$dataset\/merged\/macs2_$replicates_underscored\/extended_peaks.bed -b $blacklist -A > $WORK\/$dataset\/merged\/macs2_$replicates_underscored\/extended_bk_removal_peaks.bed";
					print "\n### bedtools_subtract_command:\n$bedtools_subtract_command\n";
					`time $bedtools_subtract_command`;
					
					# Sort peaks by highest-to-lowest score and keep top 500 peaks only
					my $sort_peaks_command = "sort -r -k5 -n $WORK\/$dataset\/merged\/macs2_$replicates_underscored\/extended_bk_removal_peaks.bed | head -n $NUM_TOP_PEAKS > $WORK\/$dataset\/merged\/macs2_$replicates_underscored\/top_peaks.bed";
					print "\n### sort_peaks_command:\n$sort_peaks_command\n";
					`time $sort_peaks_command`;
					
					# Get FASTA sequences of 500 top peaks
					my $bedtools_getfasta_command = "$BEDTOOLS getfasta -bed $WORK\/$dataset\/merged\/macs2_$replicates_underscored\/top_peaks.bed -fi $GENOME_FASTA > $WORK\/$dataset\/merged\/macs2_$replicates_underscored\/top_peaks.fa";
					print "\n### bedtools_getfasta_command:\n$bedtools_getfasta_command\n";
					`time $bedtools_getfasta_command`;
					
					# Motif discovery with meme-chip
					my $meme_chip_dir = 'meme_chip';
					#my $meme_chip_dir_num = 1;
					#while(-d "$WORK\/$dataset\/merged\/macs2_$replicates_underscored\/$meme_chip_dir") {
					#	$meme_chip_dir = 'meme_chip_' . $meme_chip_dir_num;
					#	$meme_chip_dir_num++;
					#}
					
					# Remove the directory if it already exists
					if(-d "$WORK\/$dataset\/merged\/macs2_$replicates_underscored\/$meme_chip_dir") {
						`rm -r $WORK\/$dataset\/merged\/macs2_$replicates_underscored\/$meme_chip_dir`;
					}
					
					my $meme_chip_command = "$MEME_CHIP -meme-nmotifs 5 -ccut $CCUT $WORK\/$dataset\/merged\/macs2_$replicates_underscored\/top_peaks.fa -o $WORK\/$dataset\/merged\/macs2_$replicates_underscored\/$meme_chip_dir 2>&1 | tee $WORK\/$dataset\/merged\/macs2_$replicates_underscored\/MEME_CHIP.out"; # needs to create a new, non-existent directory
					print "\n### meme_chip_command:\n$meme_chip_command\n";
					`time $meme_chip_command`;
					
					### DONE with this GROUP
					print "We reached the end of this GROUP.\n";
					
				}
			}
		}
		
		print "We reached the end of EXPERIMENT $dataset.\n";
		
	} # end this dataset (passed test)
} # end all datasets (experiments)


exit;

