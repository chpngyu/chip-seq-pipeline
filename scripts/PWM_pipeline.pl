#! /nas3/cnelson/bin/anaconda2/bin/perl
# /usr/local/miniconda3/bin/perl

# DESCRIPTION: Perl script to automate the PWM pipeline

#########################################################################################
# EXAMPLE CALL. Must be called from a "Target" name directory, e.g., for GATA1, /home/cnelson/ChIP-seq/human/GATA1
#########################################################################################
# EXAMPLE: /home/cnelson/nas3/ChIP-seq/scripts/PWM_pipeline.pl --ncpus=${NCPUS} --working_directory=${WORK} 2>&1
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
### CHANGE THIS! MANUALLY SET GLOBAL VARIABLES WITH PATH NAMES TO SOFTWARE AND INPUT. ###
#########################################################################################
#########################################################################################
# Modules already loaded by job shell script that calls this one, run_chip_seq_pipeline_S*.sh
my $CORR = '/home/cnelson/nas3/ChIP-seq/scripts/correlation.py';
my $MOTIF_CORR = '/home/cnelson/nas3/ChIP-seq/scripts/motif_cluster.py';
my $TRIM = '/home/cnelson/nas3/ChIP-seq/scripts/pwm.py';
my $MEME2MEME = '/nas3/cnelson/bin/anaconda2/bin/meme2meme'; # meme 5.0.5
my $FIMO = '/nas3/cnelson/bin/anaconda2/bin/fimo'; # 5.0.5
my $PCC_threshold = 0; # 0.8 # isolates not having a PCC at least this value will be eliminated***
my $BEDTOOLS = '/nas3/cnelson/bin/anaconda2/bin/bedtools'; # bedtools v2.28.0
my $MACS2 = '/nas3/cnelson/bin/anaconda2/bin/macs2'; # macs2 2.1.2
my $MACS2_gsize = 'hs'; # 'dm'; # 'hs';
my $PEAK_RADIUS = 100; # default=100 400
my $MIN_PEAK_SCORE = 0; # default, nothing # 20
my $MEME_CHIP = '/nas3/cnelson/bin/anaconda2/bin/meme-chip'; # 5.0.5
my $CCUT = 100; # 400; # default=100
my $PYTHON3 = '/usr/local/software/Python/python-3.6.1/bin/python3.6';
my $NUM_TOP_PEAKS = 500;
my $MFOLD_MIN = 5; #2; # 5
my $MFOLD_MAX = 50; # 50
my $blacklist = '/home/cnelson/nas3/ChIP-seq/blacklist/ENCFF023CZC_sorted.bed';
my $GENOME_FASTA = '/home/cnelson/nas3/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa';
#########################################################################################
#########################################################################################
### END "CHANGE THIS" ###
#########################################################################################
#########################################################################################


#########################################################################################
# INITIALIZE
print "\n################################################################################" .
"\n##                                                                            ##" .
"\n##                        Motif PWM Analysis Initiated!                       ##" .
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
# /home/cnelson/nas3/ChIP-seq/human/ATF2
# /home/cnelson/nas3/ChIP-seq/human/ATF3
# /home/cnelson/nas3/ChIP-seq/human/GATA1

if($working_directory) { # $working_directory passed as argument
	if(-d $working_directory) { # $working_directory is a real directory
		print "\n### DIRECTORY TO BEGIN ANALYSIS:\n";
		chomp($working_directory);
		print "WORKING_DIRECTORY=$working_directory\n";
	} else { # $working_directory is NOT a real directory
		print "\n### WARNING: DIRECTORY $working_directory DOESN'T EXIST. USING WORKING DIRECTORY TO BEGIN ANALYSIS:\n";
		$working_directory = `pwd`; # /home/cnelson/nas3/ChIP-seq/human/ATF2 
		chomp($working_directory);
		print "WORKING_DIRECTORY=$working_directory\n";
	}
} else {
	print "\n### BEGIN ANALYSIS IN WORKING DIRECTORY:\n";
	$working_directory = `pwd`; # /home/cnelson/nas3/ChIP-seq/human/ATF2 
	chomp($working_directory);
	print "WORKING_DIRECTORY=$working_directory\n";
}

my $WORK = $working_directory; # to keep nomenclature consistent
chdir("$WORK"); # no longer need to change manually


# Extract name of TF
my $TF_name;
if($WORK =~ /.*\/(.+)$/) {
	$TF_name = $1;
} else {
	die "### WARNING: could not detect TF name. TERMINATED.\n\n";
}


#########################################################################################
# Loop through subdirectories and run jobs
my %replicate2read_type;
my @wd_contents = glob "*";


print "\n################################################################################\n";
print "EXPERIMENTS TO EXAMINE: @wd_contents\n";

print "\n################################################################################\n";
print "################################################################################\n";
print "SELECTING THE TOP PWM FOR EACH REPLICATE\n";

my %replicate_topPWM;

# LOOP ALL EXPERIMENTS (TFs and also controls)
EXPERIMENT: foreach my $experiment (@wd_contents) { # i.e., experiment ID
	
	if(-d "$WORK\/$experiment" && $experiment =~ /^ENC/) { # e.g., /home/cnelson/ChIP-seq/human/ENCSR000BKC
		
		print "\n################################################################################\n";
		print "################################################################################\n";
		print "Examining experiment $experiment\n";
		
		chdir("$WORK\/$experiment");
		
		
		##########################################################################
		# STEP 1: "SELECTION OF TOP 1 PWM" for each REPLICATE
		##########################################################################
		
		
		my @wd_contents = glob "*"; 
		print "ISOGENIC REPLICATES TO EXAMINE: @wd_contents\n";
		
		# LOOP ALL REPLICATES, each of which should contain a populated ~/meme_chip/ directory
		REPLICATE: foreach my $replicate (@wd_contents) { # i.e., isogenic replicate ID
	
			if(-f "$WORK\/$experiment\/$replicate\/meme_chip\/meme_out\/meme.txt" && $replicate =~ /^ENC/) { # e.g., /home/cnelson/ChIP-seq/human/ENCSR000BKC/ENCLB559JTU
				unless($replicate =~ /control/) { # NOT controls
					print "\n################################################################################\n";
					print "Examining replicate $replicate\, file $WORK\/$experiment\/$replicate\/meme_chip\/meme_out\/meme.txt\...\n";
					
					#chdir("$WORK\/$experiment\/$replicate"); # no need
					
					# USE MEME2MEME: "takes meme motifs in many forms and writes out a single database in minimal meme format to standard output."
					# -numbers: use numbers instead of strings for motif names; default: use existing ID
				    my $meme2meme_command = "$MEME2MEME -numbers $WORK\/$experiment\/$replicate\/meme_chip\/meme_out\/meme.txt > $WORK\/$experiment\/$replicate\/PWM.meme 2>&1 | tee $WORK\/$experiment\/$replicate\/MEME2MEME.out";
					print "\n### meme2meme_command:\n$meme2meme_command\n";
					`time $meme2meme_command`;
					
					
					# USE SED (Stream EDitor) and ECHO (to add a newline)
					# This will change the file, e.g., "4 MEME-1" will become "4_MEME-1"
					# -i: "Edit files in-place, saving backups with the specified extension"
					# -e: "Append the editing commands specified by the command argument to the list of commands.""
					#     # 0x10c    268    UNMAP,MUNMAP,SECONDARY
				    my $sed_command = "sed -i -e 's/ MEME/_MEME/g' $WORK\/$experiment\/$replicate\/PWM.meme 2> $WORK\/$experiment\/$replicate\/SED.out; echo >> $WORK\/$experiment\/$replicate\/PWM.meme";
					print "\n### sed_command:\n$sed_command\n";
					`time $sed_command`;
					
					
					# USE FIMO: count occurrences of each PWM
					my $fimo_command = "$FIMO --text $WORK\/$experiment\/$replicate\/PWM.meme $WORK\/$experiment\/$replicate\/meme_chip\/top_peaks.fa > $WORK\/$experiment\/$replicate\/fimo.txt 2> $WORK\/$experiment\/$replicate\/FIMO.out";
					print "\n### fimo_command:\n$fimo_command\n";
					`time $fimo_command`;
					
					
					# USE UNIX loop to count occurrences of each PWM
					my $unix_loop_command = "echo 'MOTIF Occur' > $WORK\/$experiment\/$replicate\/occurrences.txt; " .
						"for i in {1..5}; do\n" .
							"occur=\$(grep \${i}_MEME $WORK\/$experiment\/$replicate\/fimo.txt | awk '{print \$2}' | sort | uniq | wc -l); " .
							"motif=\$(grep \${i}_MEME $WORK\/$experiment\/$replicate\/fimo.txt | awk '{print \$1}' | sort | uniq); " .
							"echo \$motif \$occur >> $WORK\/$experiment\/$replicate\/occurrences.txt\n" .
						"done";
					print "\n### unix_loop_command:\n$unix_loop_command\n";
					`time $unix_loop_command`;
					
					
					# USE Python script to trim PWMs from the end; also an optional --take argument
					my $TRIM_command = "/usr/local/software/Python/python-3.6.1/bin/python3.6 $TRIM -i $WORK\/$experiment\/$replicate\/PWM.meme --trim 0.3 -o $WORK\/$experiment\/$replicate\/PWM_trimmed.meme 2> $WORK\/$experiment\/$replicate\/TRIM.out";
					print "\n### TRIM_command:\n$TRIM_command\n";
					`time $TRIM_command`;
					
					
					# USE Python script correlation.py to determine correlations between all PWMs; ***REQUIRES PYTHON 3 (not Python 2)
					#my $CORR_command = "/usr/local/software/Python/python-3.6.1/bin/python3.6 $CORR --method pcc -m1 $WORK\/$experiment\/$replicate\/PWM.meme -o $WORK\/$experiment\/$replicate\/motif_pcc.txt 2> $WORK\/$experiment\/$replicate\/CORR.out";
					my $CORR_command = "/usr/local/software/Python/python-3.6.1/bin/python3.6 $CORR --method pcc -m1 $WORK\/$experiment\/$replicate\/PWM_trimmed.meme -o $WORK\/$experiment\/$replicate\/motif_pcc.txt 2> $WORK\/$experiment\/$replicate\/CORR.out";
					print "\n### CORR_command:\n$CORR_command\n";
					`time $CORR_command`;
					
					
					# USE Python script motif_cluster.py to determine PWM clusters.
					my $MOTIF_CORR_command = "/usr/local/software/Python/python-3.6.1/bin/python3.6 $MOTIF_CORR $WORK\/$experiment\/$replicate\/motif_pcc.txt $WORK\/$experiment\/$replicate\/occurrences.txt > $WORK\/$experiment\/$replicate\/motif_cluster.txt 2> $WORK\/$experiment\/$replicate\/MOTIF_CORR.out";
					print "\n### MOTIF_CORR_command:\n$MOTIF_CORR_command\n";
					`time $MOTIF_CORR_command`;
					
					
					# SAVE this replicate's TOP PWM, identified in motif_cluster.txt
					my $top_PWM_ID;
					my $top_PWM_occurrences;
					
					# FOR NOW, we'll just take the first line, even if it's the same cluster
					# COMEBACK: LATER WE WILL KEEP ANYTHING OVER 100 OCCURRENCES
					open(MOTIF_CLUSTER, "$WORK\/$experiment\/$replicate\/motif_cluster.txt");
					while(<MOTIF_CLUSTER>) {
						if($_ =~ /^Cluster1\s+([\w\d\-]+)\s+(\d+)$/) {
							$top_PWM_ID = $1;
							$top_PWM_occurrences = $2;
							last;
						}
					}
					close MOTIF_CLUSTER;
					
					
					### SAVE TOP REPLICATES, NOT KEEPING THE HEADER OF EACH
					# Now find the actual PWM, located in PWM.meme
					#my $top_PWM = '';
					my $read_flag = 0; # starts at 0 because we do NOT need the top of the file
					
					open(PWM_MEME, "$WORK\/$experiment\/$replicate\/PWM.meme");
					while(<PWM_MEME>) {
						if($_ =~ /^MOTIF\s+([\w\d\-]+)/) {
							my $curr_motif_ID = $1;
							
							if($curr_motif_ID eq $top_PWM_ID) {
								$read_flag = 1;
							} else {
								$read_flag = 0;
							}
						}
						
						if($read_flag == 1) {
							$replicate_topPWM{"$experiment\_$replicate"} .= $_; # should include newline
						}
						
						
					}
					close PWM_MEME;
					
					
					### DONE
					print "We reached the end of replicate $replicate.\n";
					
				} # end this replicate (passed test)
			} # replicate wasn't a control
		} # end this experiment (all replicates)
		
		print "We reached the end of experiment $experiment.\n";
		
	} # end this dataset (passed test)
} # end all datasets (experiments)


###########################################################################
## STEP 2: print the top motif of each replicate to a file, find all PCCs
###########################################################################

# Print to file
open(TOP_PWMS_REPLICATES, ">$WORK\/top_PWMs_replicates.txt");

# Print a general header
print TOP_PWMS_REPLICATES "MEME version 5.0.5\n";
print TOP_PWMS_REPLICATES "\n";
print TOP_PWMS_REPLICATES "ALPHABET= ACGT\n";
print TOP_PWMS_REPLICATES "\n";
print TOP_PWMS_REPLICATES "strands: + -\n";
print TOP_PWMS_REPLICATES "\n";
print TOP_PWMS_REPLICATES "Background letter frequencies (from unknown source):\n";
print TOP_PWMS_REPLICATES " A 0.25 C 0.25 G 0.25 T 0.25\n";
print TOP_PWMS_REPLICATES "\n";

foreach my $replicate (sort keys %replicate_topPWM) {
	# The motif ID's are not unique so we need to add the replicate name to them
	my $this_PWM = $replicate_topPWM{$replicate};
	$this_PWM =~ s/MOTIF /MOTIF $TF_name\_$replicate\_/g;
	print TOP_PWMS_REPLICATES $this_PWM;
	#print TOP_PWMS_REPLICATES $replicate_topPWM{$replicate};
}
close TOP_PWMS_REPLICATES;


# USE Python script to trim PWMs from the end; also an optional --take argument
my $TRIM_command = "/usr/local/software/Python/python-3.6.1/bin/python3.6 $TRIM -i $WORK\/top_PWMs_replicates.txt --trim 0.3 -o $WORK\/top_PWMs_replicates_trimmed.txt 2> $WORK\/TRIM_replicates.out";
print "\n### TRIM_command:\n$TRIM_command\n";
`time $TRIM_command`;


# Find all PCCs between PWMs; USE Python script correlation.py to determine correlations between all PWMs; ***REQUIRES PYTHON 3 (not Python 2)
#my $TOP_CORR_command = "/usr/local/software/Python/python-3.6.1/bin/python3.6 $CORR --method pcc -m1 $WORK\/top_PWMs_replicates.txt -o $WORK\/top_PWMs_replicates_trimmed_pcc.txt 2> $WORK\/TOP_CORR.out";
my $TOP_CORR_command = "/usr/local/software/Python/python-3.6.1/bin/python3.6 $CORR --method pcc -m1 $WORK\/top_PWMs_replicates_trimmed.txt -o $WORK\/top_PWMs_replicates_trimmed_pcc.txt 2> $WORK\/TOP_CORR.out";
print "\n### TOP_CORR_command:\n$TOP_CORR_command\n";
`time $TOP_CORR_command`;


# Identify the motifs to keep, i.e., those with PCC at least $PCC_threshold
my %good_replicates;
open(PWM_MEME, "$WORK\/top_PWMs_replicates_trimmed_pcc.txt");
while(<PWM_MEME>) {
	if($_ =~ /^([\w\d\-]+)\s+([\w\d\-]+)\s+([\d\.]+)\s*$/) {
		#GATA1_ENCSR000EFT_ENCLB209ADA_5_MEME-1	GATA1_ENCSR000EFT_ENCLB209AJT_2_MEME-1	0.9403
		my $ID1 = $1;
		my $ID2 = $2;
		my $PCC = $3;
		
		if($ID1 ne $ID2 && $PCC >= $PCC_threshold) {
			if($ID1 ne $ID2) { # make sure they're not the same!
				$good_replicates{$ID1}++;
				$good_replicates{$ID2}++;
			}
		}
		
	}
}

my %good_replicates_short;
print "REPLICATES WITH AT LEAST ONE PCC >= $PCC_threshold (ID / COUNT / SHORT ID)\:\n";
foreach my $replicate (sort keys %good_replicates) {
	print "$replicate\t" . $good_replicates{$replicate} . "\t";
	
	if($replicate =~ /^[\w\-]+_[\w\-]+_([\w\-]+)_\d+\_MEME/) { # VERY specific
		$good_replicates_short{$1} = 1;
		print "$1\n";
	}
	
}


##########################################################################
# STEP 3: CALCULATE THE MERGED PWM for each EXPERIMENT, USING ONLY GOOD REPLICATES
##########################################################################

print "\n################################################################################\n";
print "################################################################################\n";
print "CALCULATING MERGED PWM FOR EACH EXPERIMENT USING GOOD REPLICATES\n";


# LOOP ALL EXPERIMENTS (TFs and also controls)
EXPERIMENT: foreach my $experiment (@wd_contents) { # i.e., experiment ID
	
	if(-d "$WORK\/$experiment" && $experiment =~ /^ENC/) { # e.g., /home/cnelson/ChIP-seq/human/ENCSR000BKC
		
		print "\n################################################################################\n";
		print "################################################################################\n";
		print "Examining experiment $experiment\n";
		
		chdir("$WORK\/$experiment");
		
		my @wd_contents = glob "*";
		
		print "ISOGENIC REPLICATES TO CONSIDER: @wd_contents\n";
		
		my @wd_contents_kept;
		
		foreach my $item (@wd_contents) {
			if($good_replicates_short{$item}) {
				push(@wd_contents_kept, $item);
			}
		}
		
		print "ISOGENIC REPLICATES TO USE: @wd_contents_kept\n";
		
		if(@wd_contents_kept > 0) { # at least one good replicate
			
			# ANALYZE BAM FILES: all good replicates
			
			# There should always be one BAM file, by the name 'unique_alignment_sorted_rd.bam'
			print "FINDING BAM FILES: unique_alignment_sorted_rd.bam\n";
			
			# Make output directory? What does MACS2 do if files already exist?
			unless(-d "$WORK\/$experiment\/macs2_merged") {
				mkdir("$WORK\/$experiment\/macs2_merged");
			}
			
	
			##########################################################################
			# S3a: MACS2 for PEAK CALLING
			# DIFFERENT procedure for SE and PE
			##########################################################################
			
			# This required a CIRCUS of Python 2.7 installation MADNESS
			
			# Gather names of control files, which are .bam files that have been placed in $WORK/$experiment/controls/
			my @control_fastq_file_list = glob "$WORK\/$experiment\/*_control\/unique_alignment_sorted_rd.bam";
			#my $CONTROL_STRING = '';
			#foreach my $control_fastq_file (@control_fastq_file_list) {
			#	if($control_fastq_file =~ /^.+\/(.+).bam$/) {
			#		#$control_fastq_file = "$1\.bam"; # we want the whole path this time
			#		$CONTROL_STRING .= "$control_fastq_file "; # space at the end!
			#	}
			#}
			
			my $CONTROL_STRING = "@control_fastq_file_list"; # this should do it
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
			
			# Identify FASTQ file names; there will be 
			
			
			# Find out if SE (1 fastq) or PE (2 fastq) and build pooled BAM argument
			my $BAM_argument = '';
			my $read_type = '';
			foreach my $replicate (@wd_contents_kept) {
				my @fastq_file_list = glob "$WORK\/$experiment\/$replicate\/*.fastq.gz";
				my $fastq_file_count = scalar(@fastq_file_list);
				
				if($fastq_file_count == 2) { # PAIRED-END ANALYSIS (if two inputs, or read_type =~ /^PE/)
				
					if($read_type eq '') {
						$read_type = 'PE';
					} elsif($read_type eq 'SE') {
						die "### WARNING: conflicting read types for $WORK\/$experiment\/$replicate\! TERMINATED.\n";
					}
				} elsif($fastq_file_count == 1) { # SINGLE-END ANALYSIS (if read_type =~ /^SE/)
					
					if($read_type eq '') {
						$read_type = 'SE';
					} elsif($read_type eq 'PE') {
						die "### WARNING: conflicting read types for $WORK\/$experiment\/$replicate\! TERMINATED.\n";
					}
					
				} else {
					die "There were 0 or >2 reads files! TERMINATED.";	
				}
				
				$BAM_argument .= "$WORK\/$experiment\/$replicate\/unique_alignment_sorted_rd.bam ";
			}
			chop($BAM_argument);
			
			if($read_type eq 'PE') { # $fastq_file_count == 2) { # PAIRED-END (PE)
				
				# Call peaks for this replicate ### WILL THE BELOW CREATE A macs2 DIRECTORY??
				
				# INDIVIDUAL REPLICATE
			    my $macs2_merged_command = "$MACS2 callpeak -t $BAM_argument -c $CONTROL_STRING -f BAMPE --gsize $MACS2_gsize --mfold $MFOLD_MIN $MFOLD_MAX --outdir $WORK\/$experiment\/macs2_merged 2>&1 | tee $WORK\/$experiment\/MACS2_MERGED.out";
				print "\n### macs2_merged_command:\n$macs2_merged_command\n";
				`time $macs2_merged_command`;
				
			} elsif ($read_type eq 'SE') { # $fastq_file_count == 1) { # SINGLE-END (SE)
				
				# Call peaks for this replicate
				
				# INDIVIDUAL REPLICATE
			    my $macs2_merged_command = "$MACS2 callpeak -t $BAM_argument -c $CONTROL_STRING -f BAM --gsize $MACS2_gsize --mfold $MFOLD_MIN $MFOLD_MAX --outdir $WORK\/$experiment\/macs2_merged 2>&1 | tee $WORK\/$experiment\/MACS2_MERGED.out";
				print "\n### macs2_merged_command:\n$macs2_merged_command\n";
				`time $macs2_merged_command`;
				
				
			} else {
				die "There were 0 or >2 reads files! TERMINATED.";
			}
			
			
			
			##########################################################################
			# S3b: MOTIF DISCOVERY
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
			# TOP OF: ATF3/ENCSR000BKC/ENCLB559JTU/macs2/NA_summits.bed
			#1	629938  629939  NA_peak_1       26.02448
			#1	634028  634029  NA_peak_2       35.47879
			#1	778706  778707  NA_peak_3       46.39153
			#1	1000913 1000914 NA_peak_4       7.42275
			my $awk_on_macs2_command = "awk 'BEGIN {WIDTH=$PEAK_RADIUS} {if(\$2<=WIDTH) print \$1 \"\\t1\\t\" \$2+WIDTH \"\\t\" \$4 \"\\t\" \$5; else print \$1 \"\\t\" \$2-WIDTH \"\\t\" \$2+WIDTH \"\\t\" \$4 \"\\t\" \$5}' $WORK\/$experiment\/macs2_merged/NA_summits.bed | awk '\$5 \>\= $MIN_PEAK_SCORE' > $WORK\/$experiment\/extended_peaks.bed"; # COMEBACK to figure out how to tee this, just the error
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
			
			# Remove peaks in blacklist with bedtools
			my $bedtools_subtract_command = "$BEDTOOLS subtract -a $WORK\/$experiment\/extended_peaks.bed -b $blacklist -A > $WORK\/$experiment\/extended_bk_removal_peaks.bed";
			print "\n### bedtools_subtract_command:\n$bedtools_subtract_command\n";
			`time $bedtools_subtract_command`;
			
			# Sort peaks by highest-to-lowest score and keep top 500 peaks only
			# -r reverse the output, sorting highest to lowest 
			# -k specifies which column(s) to sort by (here, 5, which is scores)
			my $sort_peaks_command = "sort -r -k5 -n $WORK\/$experiment\/extended_bk_removal_peaks.bed | head -n $NUM_TOP_PEAKS > $WORK\/$experiment\/top_peaks.bed";
			print "\n### sort_peaks_command:\n$sort_peaks_command\n";
			`time $sort_peaks_command`;
			
			# Get FASTA sequences of 500 top peaks
			my $bedtools_getfasta_command = "$BEDTOOLS getfasta -bed $WORK\/$experiment\/top_peaks.bed -fi $GENOME_FASTA > $WORK\/$experiment\/top_peaks.fa";
			print "\n### bedtools_getfasta_command:\n$bedtools_getfasta_command\n";
			`time $bedtools_getfasta_command`;
			
#			# Motif discovery with meme-chip
#			my $meme_chip_merged_dir = 'meme_chip_merged';
#			my $meme_chip_merged_dir_num = 1;
#			while(-d "$WORK\/$experiment\/$meme_chip_merged_dir") {
#				$meme_chip_merged_dir = 'meme_chip_merged_' . $meme_chip_merged_dir_num;
#				$meme_chip_merged_dir_num++;
#			}
			
#			my $meme_chip_command = "$MEME_CHIP -meme-nmotifs 5 $WORK\/$experiment\/top_peaks.fa -o $WORK\/$experiment\/$meme_chip_merged_dir 2>&1 | tee MEME_CHIP.out"; # needs to create a new, non-existent directory
			my $meme_chip_command = "$MEME_CHIP -meme-nmotifs 5 -ccut $CCUT $WORK\/$experiment\/top_peaks.fa -o $WORK\/$experiment\/meme_chip_merged 2>&1 | tee $WORK\/$experiment\/MEME_CHIP_MERGED.out"; # directory CAN'T ALREADY EXIST; needs to create a new, non-existent directory
			print "\n### meme_chip_command:\n$meme_chip_command\n";
			`time $meme_chip_command`;
			
			print "We reached the end of experiment $experiment.\n";
			
			
		} # this experiment has at least one good replicate
	} # end this experiment (passed test)
} # end all experiments


##########################################################################
# STEP 4: CALCULATE ALL PAIRWISE PCCs BETWEEN MERGED PWMs
##########################################################################

### EACH EXPERIMENT now has its own (MERGED) meme-chip results in /meme_chip_merged/

my %experiment_topPWM;

# LOOP ALL EXPERIMENTS (TFs and also controls)
EXPERIMENT: foreach my $experiment (@wd_contents) { # i.e., experiment ID
	
	if(-d "$WORK\/$experiment" && $experiment =~ /^ENC/) { # e.g., /home/cnelson/ChIP-seq/human/ENCSR000BKC
		
		print "\n################################################################################\n";
		print "################################################################################\n";
		print "Examining experiment $experiment\n";
		
		chdir("$WORK\/$experiment");
		
		
		##########################################################################
		# STEP 4-1: "SELECTION OF TOP 1 PWM" for each REPLICATE
		##########################################################################
		
		if(-f "$WORK\/$experiment\/meme_chip_merged\/meme_out\/meme.txt") { # e.g., /home/cnelson/ChIP-seq/human/ENCSR000BKC/
		
			print "\n################################################################################\n";
			print "Examining file $WORK\/$experiment\/meme_chip_merged\/meme_out\/meme.txt\...\n";
			
			#chdir("$WORK\/$experiment"); # no need
			
			# USE MEME2MEME: "takes meme motifs in many forms and writes out a single database in minimal meme format to standard output."
			# -numbers: use numbers instead of strings for motif names; default: use existing ID
		    my $meme2meme_command = "$MEME2MEME -numbers $WORK\/$experiment\/meme_chip_merged\/meme_out\/meme.txt > $WORK\/$experiment\/PWM.meme 2>&1 | tee $WORK\/$experiment\/MEME2MEME.out";
			print "\n### meme2meme_command:\n$meme2meme_command\n";
			`time $meme2meme_command`;
			
			
			# USE SED (Stream EDitor) and ECHO (to add a newline)
			# This will change the file, e.g., "4 MEME-1" will become "4_MEME-1"
			# -i: "Edit files in-place, saving backups with the specified extension"
			# -e: "Append the editing commands specified by the command argument to the list of commands.""
			#     # 0x10c    268    UNMAP,MUNMAP,SECONDARY
		    my $sed_command = "sed -i -e 's/ MEME/_MEME/g' $WORK\/$experiment\/PWM.meme 2> $WORK\/$experiment\/SED.out; echo >> $WORK\/$experiment\/PWM.meme";
			print "\n### sed_command:\n$sed_command\n";
			`time $sed_command`;
			
			
			# USE FIMO: count occurrences of each PWM
			my $fimo_command = "$FIMO --text $WORK\/$experiment\/PWM.meme $WORK\/$experiment\/meme_chip_merged\/top_peaks.fa > $WORK\/$experiment\/fimo.txt 2> $WORK\/$experiment\/FIMO.out";
			print "\n### fimo_command:\n$fimo_command\n";
			`time $fimo_command`;
			
			
			# USE UNIX loop to count occurrences of each PWM
			my $unix_loop_command = "echo 'MOTIF Occur' > $WORK\/$experiment\/occurrences.txt; " .
				"for i in {1..5}; do\n" .
					"occur=\$(grep \${i}_MEME $WORK\/$experiment\/fimo.txt | awk '{print \$2}' | sort | uniq | wc -l); " .
					"motif=\$(grep \${i}_MEME $WORK\/$experiment\/fimo.txt | awk '{print \$1}' | sort | uniq); " .
					"echo \$motif \$occur >> $WORK\/$experiment\/occurrences.txt\n" .
				"done";
			print "\n### unix_loop_command:\n$unix_loop_command\n";
			`time $unix_loop_command`;
			
			
			# USE Python script to trim PWMs from the end; also an optional --take argument
			my $TRIM_command = "/usr/local/software/Python/python-3.6.1/bin/python3.6 $TRIM -i $WORK\/$experiment\/PWM.meme --trim 0.3 -o $WORK\/$experiment\/PWM_trimmed.meme 2> $WORK\/$experiment\/TRIM.out";
			print "\n### TRIM_command:\n$TRIM_command\n";
			`time $TRIM_command`;

			
			# USE Python script correlation.py to determine correlations between all PWMs; ***REQUIRES PYTHON 3 (not Python 2)
			#my $CORR_command = "/usr/local/software/Python/python-3.6.1/bin/python3.6 $CORR --method pcc -m1 $WORK\/$experiment\/PWM.meme -o $WORK\/$experiment\/motif_pcc.txt 2> $WORK\/$experiment\/CORR.out";
			my $CORR_command = "/usr/local/software/Python/python-3.6.1/bin/python3.6 $CORR --method pcc -m1 $WORK\/$experiment\/PWM_trimmed.meme -o $WORK\/$experiment\/motif_pcc.txt 2> $WORK\/$experiment\/CORR.out";
			print "\n### CORR_command:\n$CORR_command\n";
			`time $CORR_command`;
			
			
			# USE Python script motif_cluster.py to determine PWM clusters.
			my $MOTIF_CORR_command = "/usr/local/software/Python/python-3.6.1/bin/python3.6 $MOTIF_CORR $WORK\/$experiment\/motif_pcc.txt $WORK\/$experiment\/occurrences.txt > $WORK\/$experiment\/motif_cluster.txt 2> $WORK\/$experiment\/MOTIF_CORR.out";
			print "\n### MOTIF_CORR_command:\n$MOTIF_CORR_command\n";
			`time $MOTIF_CORR_command`;
			
			
			# SAVE this replicate's TOP PWM, identified in motif_cluster.txt
			my $top_PWM_ID;
			my $top_PWM_occurrences;
			
			# FOR NOW, we'll just take the first line, even if it's the same cluster
			# COMEBACK: LATER WE WILL KEEP ANYTHING OVER 100 OCCURRENCES
			open(MOTIF_CLUSTER_EXP, "$WORK\/$experiment\/motif_cluster.txt");
			while(<MOTIF_CLUSTER_EXP>) {
				if($_ =~ /^Cluster1\s+([\w\d\-]+)\s+(\d+)$/) {
					$top_PWM_ID = $1;
					$top_PWM_occurrences = $2;
					last;
				}
			}
			close MOTIF_CLUSTER_EXP;
			
			### SAVE TOP REPLICATES, KEEPING THE HEADER OF EACH
			# Now find the actual PWM, located in PWM.meme
#					#my $top_PWM;
#					my $read_flag = 1; # starts at 1 because we need the top of the file
#					
#					open(PWM_MEME, "$WORK\/$experiment\/PWM.meme");
#					while(<PWM_MEME>) {
#						if($_ =~ /^MOTIF\s+([\w\d\-]+)/) {
#							my $curr_motif_ID = $1;
#							
#							if($curr_motif_ID eq $top_PWM_ID) {
#								$read_flag = 1;
#							} else {
#								$read_flag = 0;
#							}
#						}
#						
#						if($read_flag == 1) {
#							$replicate_topPWM{"$experiment\_$replicate"} .= $_; # should include newline
#						}
#						
#						
#					}
#					close PWM_MEME;
			
			
			### SAVE TOP REPLICATES, NOT KEEPING THE HEADER OF EACH
			# Now find the actual PWM, located in PWM.meme
			#my $top_PWM = '';
			my $read_flag = 0; # starts at 0 because we do NOT need the top of the file
			
			open(PWM_MEME, "$WORK\/$experiment\/PWM.meme");
			while(<PWM_MEME>) {
				if($_ =~ /^MOTIF\s+([\w\d\-]+)/) {
					my $curr_motif_ID = $1;
					
					if($curr_motif_ID eq $top_PWM_ID) {
						$read_flag = 1;
					} else {
						$read_flag = 0;
					}
				}
				
				if($read_flag == 1) {
					$experiment_topPWM{"$experiment"} .= $_; # should include newline
				}
				
			}
			close PWM_MEME;
			
		} else { # meme.txt existed
			print "### WARNING: no file $WORK\/$experiment\/meme_chip_merged\/meme_out\/meme.txt\! This experiment must not have enough good replicates.";
		}
		
		print "We reached the end of experiment $experiment.\n";
		
	} # end this dataset (passed test)
} # end all datasets (experiments)



###########################################################################
## STEP 5: print the top motif of each EXPERIMENT to a file, find all pairwise PCCs
###########################################################################

# Print to file
open(TOP_PWMS_EXPERIMENTS, ">$WORK\/top_PWMs_experiments.txt");

# Print a general header
print TOP_PWMS_EXPERIMENTS "MEME version 5.0.5\n";
print TOP_PWMS_EXPERIMENTS "\n";
print TOP_PWMS_EXPERIMENTS "ALPHABET= ACGT\n";
print TOP_PWMS_EXPERIMENTS "\n";
print TOP_PWMS_EXPERIMENTS "strands: + -\n";
print TOP_PWMS_EXPERIMENTS "\n";
print TOP_PWMS_EXPERIMENTS "Background letter frequencies (from unknown source):\n";
print TOP_PWMS_EXPERIMENTS " A 0.25 C 0.25 G 0.25 T 0.25\n";
print TOP_PWMS_EXPERIMENTS "\n";

foreach my $experiment (sort keys %experiment_topPWM) {
	# The motif ID's are not unique so we need to add the experiment name to them
	my $this_PWM = $experiment_topPWM{$experiment};
	$this_PWM =~ s/MOTIF /MOTIF $TF_name\_$experiment\_/g;
	print TOP_PWMS_EXPERIMENTS $this_PWM;
	#print TOP_PWMS_EXPERIMENTS $experiment_topPWM{$experiment};
}
close TOP_PWMS_EXPERIMENTS;


# USE Python script to trim PWMs from the end; also an optional --take argument
my $TRIM_command = "$PYTHON3 $TRIM -i $WORK\/top_PWMs_experiments.txt --trim 0.3 -o $WORK\/top_PWMs_experiments_trimmed.txt 2> $WORK\/TRIM_experiments.out";
print "\n### TRIM_command:\n$TRIM_command\n";
`time $TRIM_command`;


# Find all PCCs between PWMs; USE Python script correlation.py to determine correlations between all PWMs; ***REQUIRES PYTHON 3 (not Python 2)
#my $TOP_CORR_command = "/usr/local/software/Python/python-3.6.1/bin/python3.6 $CORR --method pcc -m1 $WORK\/top_PWMs_experiments.txt -o $WORK\/top_PWMs_experiments_trimmed_pcc.txt 2> $WORK\/TOP_CORR.out";
my $TOP_CORR_command = "$PYTHON3 $CORR --method pcc -m1 $WORK\/top_PWMs_experiments_trimmed.txt -o $WORK\/top_PWMs_experiments_trimmed_pcc.txt 2> $WORK\/TOP_CORR_experiments.out";
print "\n### TOP_CORR_command:\n$TOP_CORR_command\n";
`time $TOP_CORR_command`;

exit;


