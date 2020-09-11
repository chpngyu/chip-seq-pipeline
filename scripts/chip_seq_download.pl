#! /usr/bin/env perl

# DESCRIPTION: Perl script to automate the DOWNLOADING of ChIP-seq ENCODE data

# MUST ENTER MANUALLY FOR PBS:
#### ENTER the environment in INTERACTIVE mode, and go to your own node like so:
# qsub -I -N FASTQ_download -l ncpus=1,walltime=24:00:00

#########################################################################################
# EXAMPLE CALL. Must be called from the directory you wish to populate, e.g., /home/cnelson/ChIP-seq/species/
#########################################################################################
# EXAMPLE: chip_seq_download.pl --input_file=ENCODE_download_data.txt
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
#use warnings;
use Data::Dumper;
use Getopt::Long;
#use List::Util qw(max);
STDOUT->autoflush(1);


#########################################################################################
# INITIALIZE VARIABLE
my $input_file;

GetOptions( "input_file=s" => \$input_file)
or die "\n### WARNING: Error in command line arguments. Script terminated.\n\n";

# Get the time
my $time1 = time;
my $local_time1 = localtime;

print "\n##########################################################################################";
print "\nENCODE ChIP-seq data downloading initiated at local time $local_time1";
print "\n##########################################################################################\n\n";

unless(-f "$input_file") {
	die "\n### WARNING: An existing --input_file option must be provided\n".
		"### Script terminated.\n\n";
}


#########################################################################################
# Determine working directory
print "\n### BEGIN ANALYSIS IN WORKING DIRECTORY:\n";
my $working_directory = `pwd`; # /home/cnelson/nas3/ChIP-seq/dm6/
chomp($working_directory);
my $WORK = $working_directory; # to keep nomenclature consistent
print "WORKING_DIRECTORY=$WORK\n";


#########################################################################################
# LOOP FILE TO CREATE HASH TREE

# EXAMPLE FIRST RECORDS IN INPUT FILE:
#Target name	directory_name	AC	type	Library	url
#abd-A	abd_A	ENCSR609JDR	replicate	ENCLB367MZH	https://www.encodeproject.org/files/ENCFF036RZV/@@download/ENCFF036RZV.fastq.gz
#abd-A	abd_A	ENCSR609JDR	replicate	ENCLB506QBD	https://www.encodeproject.org/files/ENCFF999PWE/@@download/ENCFF999PWE.fastq.gz
#abd-A	abd_A	ENCSR609JDR	replicate	ENCLB589GMA	https://www.encodeproject.org/files/ENCFF646QDZ/@@download/ENCFF646QDZ.fastq.gz
#abd-A	abd_A	ENCSR609JDR	control	ENCLB428GIT_control	https://www.encodeproject.org/files/ENCFF120UZS/@@download/ENCFF120UZS.fastq.gz
#Abd-B	Abd_B	ENCSR465HPZ	replicate	ENCLB009GTL	https://www.encodeproject.org/files/ENCFF469WCT/@@download/ENCFF469WCT.fastq.gz
#Abd-B	Abd_B	ENCSR465HPZ	replicate	ENCLB205LSV	https://www.encodeproject.org/files/ENCFF969SEG/@@download/ENCFF969SEG.fastq.gz
#Abd-B	Abd_B	ENCSR465HPZ	control	ENCLB730TLW_control	https://www.encodeproject.org/files/ENCFF730IEL/@@download/ENCFF730IEL.fastq.gz
#Abd-B	Abd_B	ENCSR692UBK	replicate	ENCLB526KET	https://www.encodeproject.org/files/ENCFF820QJV/@@download/ENCFF820QJV.fastq.gz
#Abd-B	Abd_B	ENCSR692UBK	replicate	ENCLB588NTL	https://www.encodeproject.org/files/ENCFF252IVG/@@download/ENCFF252IVG.fastq.gz
#Abd-B	Abd_B	ENCSR692UBK	control	ENCLB728VIS_control	https://www.encodeproject.org/files/ENCFF354CHN/@@download/ENCFF354CHN.fastq.gz
#achi	achi	ENCSR959SWC	replicate	ENCLB061SFK	https://www.encodeproject.org/files/ENCFF979YUJ/@@download/ENCFF979YUJ.fastq.gz
#achi	achi	ENCSR959SWC	replicate	ENCLB064AEO	https://www.encodeproject.org/files/ENCFF881ITO/@@download/ENCFF881ITO.fastq.gz
#achi	achi	ENCSR959SWC	replicate	ENCLB722DLY	https://www.encodeproject.org/files/ENCFF721FQK/@@download/ENCFF721FQK.fastq.gz
#achi	achi	ENCSR959SWC	control	ENCLB240LXI_control	https://www.encodeproject.org/files/ENCFF548BRW/@@download/ENCFF548BRW.fastq.gz


#########################################################################################
# LOOP CCDS FILE A FIRST TIME TO DETERMINE WHICH CCDS_ID's WE WILL USE (ONLY THE LATEST)

# Determine header/column indices
my @header_names_arr;
my $line = 0;
open(INPUT_FILE, "$input_file") or die "Could not open data file $input_file\n";

while(<INPUT_FILE>) {
	chomp;

	if($line == 0) { # it's the header
		$_ =~ s/^#//g;
		@header_names_arr = split(/\t/, $_);
		$line++;
		last;
	}
}

close INPUT_FILE;

my %header_indices;

# Determine the index of each column
for (my $i = 0; $i < scalar(@header_names_arr); $i++) {
	my $curr_header = $header_names_arr[$i];
	$header_indices{$curr_header} = $i;
}

# Expected headers
my @required_headers = ('directory_name', 'AC', 'Library', 'url');

foreach(@required_headers) {
	unless(exists $header_indices{$_}) {
		die "### DIE: the header name $_ is not present in the input file.\n\n";
	}
}

# Store each directory_name -> AC -> Library -> url
my %url_data;

open(INPUT_FILE, "$input_file") or die "Could not open input file $input_file\n";

print "\n\nReading in records from input file...\n\n";

$line = 0;

while(<INPUT_FILE>) {
	chomp;
	
	if($line > 0) { # unless it's the header
		
		my @line_arr = split(/\t/, $_);
		
		# Store the data
		my $directory_name = $line_arr[$header_indices{directory_name}];
		my $AC = $line_arr[$header_indices{AC}];
		my $Library = $line_arr[$header_indices{Library}];
		my $url = $line_arr[$header_indices{url}];
		
		$url_data{$directory_name}->{$AC}->{$Library}->{$url} = 1;
		
	} # not header
	
	$line++;
} # end input file

close INPUT_FILE;

print "### WE FOUND THE FOLLOWING NUMBER OF UNIQUE DIRECTORIES (TFs):\n";
print "### " . scalar(keys %url_data) . "\n";


#########################################################################################
# DOWNLOAD ENCODE DATA
print "\n\nCreating directories (if not yet existent) and downloading the URLs (if not yet downloaded)...\n\n";

foreach my $directory_name (sort keys %url_data) {
	
	# Make TF directory if not existent
	unless(-d "$WORK\/$directory_name") {
		mkdir("$WORK\/$directory_name");
	}
	
	# Change to TF directory, which will contain experiment(s)
	print "### Entering directory $WORK\/$directory_name\n";
	chdir("$WORK\/$directory_name");
	
	foreach my $experiment (sort keys %{$url_data{$directory_name}}) {
	
		# Make experiment directory if not existent
		unless(-d "$WORK\/$directory_name\/$experiment") {
			mkdir("$WORK\/$directory_name\/$experiment");
		}
		
		# Change to experiment directory, which will contain replicate(s)
		print "### Entering directory $WORK\/$directory_name\/$experiment\n";
		chdir("$WORK\/$directory_name\/$experiment");
		
		foreach my $replicate (sort keys %{$url_data{$directory_name}->{$experiment}}) {
	
			# Make replicate directory if not existent
			unless(-d "$WORK\/$directory_name\/$experiment\/$replicate") {
				mkdir("$WORK\/$directory_name\/$experiment\/$replicate");
			}
			
			# Change to replicate directory, which will contain downloaded FASTQ(s)
			print "### Entering directory $WORK\/$directory_name\/$experiment\/$replicate\n";
			chdir("$WORK\/$directory_name\/$experiment\/$replicate");
			
			foreach my $url (sort keys %{$url_data{$directory_name}->{$experiment}->{$replicate}}) {
				
				my $url_file_name;
				if($url =~ /.+\/(.+)/) {
					$url_file_name = $1;
				}
				
				# Download FASTQ file if not existent
				if(-f "$WORK\/$directory_name\/$experiment\/$replicate\/$url_file_name" && -f "$WORK\/$directory_name\/$experiment\/$replicate\/wget.done") {
					print "\n### The following URL has already been downloaded: $url\n### Moving on!\n\n";
				} else {
					my $wget_command = 'wget ' . $url . ' > wget.log'; # WORKS. Guess you juts can't combine wget with a second (or preceding 'time') command.
					
					print "\n### wget_command:\n$wget_command\n\n";
					    
					# RUN WGET TO DOWNLOAD
					#`time $wget_command`; # doesn't work
					`$wget_command`;
					
					if($? == 0) { # produce done file if successful
						`touch wget.done`;
					}
					
					print "\n";
					
				} # finish downloading block
			}
		}
	}
}





