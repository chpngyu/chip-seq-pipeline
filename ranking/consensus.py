import os
import subprocess
class Consensus:
	def __init__(self, filename, combined_peaks, maxd, script_path):
		self.temp_group_bed = 'tmp_group.bed'
		self.temp_rmbk_bed = 'tmp_rmbk_group.bed'
		self.temp_score_bed = 'tmp_score_group.bed'
		self.temp_combined_bed = 'tmp_combined.bed'
		self.temp_combined_sorted_bed = 'tmp_combined_scorted.bed'
		self.temp_cluster = 'cluster.bed'
		self.top_combined_peaks_bed = combined_peaks
		self.max_distance = maxd
		self.path = script_path
		self.group_path = {}
		with open(filename) as infile:
			infile.readline()
			for line in infile:
				tokens = line.rstrip().split('\t')
				if tokens[0] in self.group_path:
					raise RuntimeError('Duplicate group found: %s' % tokens[1])
				else:
					self.group_path[tokens[0]] = tokens[1]
		print('Got %d groups in %s' % (len(self.group_path), filename))

	def add(self, group, path, blacklist, width=100):
		if width:
			value = '''BEGIN {WIDTH=%d} {if($2<=WIDTH) print $1 "\t1\t" $3+WIDTH-1 "\t" $4 "\t" $5; else print $1 "\t" $2-WIDTH "\t" $3+WIDTH-1 "\t" $4 "\t" $5}''' % width
		else:
			value = '''{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}'''
		ofile = open(self.temp_group_bed, 'w')
		command = ['awk', repr(value).strip("'''"), path]
		if not subprocess.run(command, stdout=ofile, check=True):
			print('Failed to run command: %s' % ' '.join(command))
			ofile.close()
			return None
		ofile.close()
		# remove peald in blacklist
		ofile = open(self.temp_rmbk_bed, 'w')
		command = ['bedtools', 'subtract', '-a', self.temp_group_bed, '-b',  blacklist, '-A']
		if not subprocess.run(command, stdout=ofile, check=True):
			print('Failed to run command: %s' % ' '.join(command))
			ofile.close()
			return None
		ofile.close()
		#scroing peaks
		SCORING = os.path.join(self.path, 'score_quantile.py')
		command = ['python', SCORING, '-s', self.temp_rmbk_bed, '-g', group, '-o', self.temp_score_bed]
		if not subprocess.run(command, check=True):
			print('Failed to run command: %s' % ' '.join(command))
			return None
		return open(self.temp_score_bed).readlines()

	def run_all(self, outfile, NoTopPeak, get_sequence, score_type, width, genome, blacklist):
		ofile = open(self.temp_combined_bed, 'w')
		for group in self.group_path:
			content = self.add(group, self.group_path[group], blacklist, width)
			if content:
				ofile.write('\n'.join(content))
		ofile.close()
		# sorting peaks by positions
		command = ['bedtools', 'sort', '-i', self.temp_combined_bed]
		ofile = open(self.temp_combined_sorted_bed, 'w')
		if not subprocess.run(command, stdout=ofile, check=True):
			print('Failed to run command: %s' % ' '.join(command))
			return None
		ofile.close()
		# clustering peaks
		if self.max_distance == 0:
			command = ['bedtools', 'cluster', '-i', self.temp_combined_sorted_bed]
		else:
			distance = '%d' % self.max_distance
			command = ['bedtools', 'cluster', '-d', distance, '-i', self.temp_combined_sorted_bed]
		ofile = open(self.temp_cluster, 'w')
		if not subprocess.run(command, stdout=ofile, check=True):
			print('Failed to run command: %s' % ' '.join(command))
			return None
		ofile.close()
		SEL_TOP = os.path.join(self.path, 'top_peak_sel.py')
		command = ['python', SEL_TOP, '-i', self.temp_cluster, '-o', self.top_combined_peaks_bed]
		if score_type == 1:
			command += ['--score', '1']
		else:
			command += ['--score', '2']

		if NoTopPeak:
			command += ['--top', str(NoTopPeak)]
		if not subprocess.run(command, check=True):
			print('Failed to run command: %s' % ' '.join(command))
			return None
		# get a sequence for each peak
		if get_sequence == False:
			return 1
		command = ['bedtools', 'getfasta', '-bed', self.top_combined_peaks_bed, '-fi', genome]
		ofile = open(outfile, 'w')
		if not subprocess.run(command, stdout=ofile, check=True):
			print('Failed to run command: %s' % ' '.join(command))
			return None
		ofile.close()
		return 1

	def __del__(self):
		tmp_files = [self.temp_group_bed, self.temp_rmbk_bed, self.temp_score_bed, self.temp_combined_bed, self.temp_combined_sorted_bed]
		for tmpf in tmp_files:
			if os.path.isfile(tmpf):
				os.remove(tmpf)

import argparse
parser = argparse.ArgumentParser(description='obtaining consenus peaks')
parser.add_argument('-i', '--input', type=str, required=True,
					help='give a table which includes groupped name and its paths of folder')
parser.add_argument('--genome', type=str, required=True,
					help ='genome assembly')
parser.add_argument('--blacklist', type=str, required=True,
					help ='black list for removing artifact and noise peaks')
parser.add_argument('--top', type=int,
					help='no. of top peaks be chose (default: all)')
parser.add_argument('--width', type=int, default=100,
					help='extension of peak (default: 100)')
parser.add_argument('--nosequence', action='store_false',
					help='whether to get sequence (default: Yes)')
parser.add_argument('--maxd', type=int, default=0,
					help='maximum distance between features to be clustered (default: 0)')
parser.add_argument('--score', type=int, default=1,
					help ='ouput score for individual score (1; default) or sumation score (2)')
parser.add_argument('--combined_peak', type=str, default='top_combined.bed',
					help ='ouput name for the combined peaks (default: top_combined.bed)')
parser.add_argument('--path', type=str, default='.',
					help ='path for external python scripts')
parser.add_argument('-o', '--output', type=str, default='top_combined_peaks.fa',
					help='output sequences for top peaks (default: top_combined_peaks.fa')
args = parser.parse_args()

consensus = Consensus(args.input, args.combined_peak, args.maxd, args.path)
consensus.run_all(args.output, args.top, args.nosequence, args.score, args.width, args.genome, args.blacklist)
