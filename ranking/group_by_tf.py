
def write_table(filename, records):
	with open(filename, 'w') as ofile:
		ofile.write('Group\tPath\n')
		for record in records:
			ofile.write('%s\t%s\n' % record)

import pandas as pd
import sys
import os

import argparse
parser = argparse.ArgumentParser(description='generate a script for infering a TF PWM that has mutiple biosamples')
parser.add_argument('-t', '--table', type=str, required=True,
					help='a table that includes TF and its called peak file')
parser.add_argument('-p', '--path', type=str, required=True,
					help='root path for the python scripts')
parser.add_argument('-g', '--genome', type=str, required=True,
					help='genome assebly file')
parser.add_argument('--blacklist', type=str, required=True,
					help='black list file to remove unwanted peaks')
parser.add_argument('-o', '--output', type=str, default='batch_run.sh',
					help='the ouput filename for the generated bash script (default: batch_run.sh)')
args = parser.parse_args()

consensus_script = os.path.join(args.path, 'consensus.py')
command1 = 'python %s -i table --width 100 --top 500 --path %s' % (consensus_script, args.path)
command1 += ' --genome %s --blacklist %s' % (args.genome, args.blacklist)
command2 = 'meme-chip -ccut 0 -meme-nmotifs 5 -oc meme_ranking_out top_combined_peaks.fa'
command3 = 'meme2meme meme_ranking_out/meme_out/meme.txt > top.meme'
command4 = '''BEGIN {WIDTH=100} {if($2<=WIDTH) print $1 "\t1\t" $3+WIDTH-1 "\t" $4 "\t" $5; else print $1 "\t" $2-WIDTH "\t" $3+WIDTH-1 "\t" $4 "\t" $5}'''

data = pd.read_csv(args.table, sep='\t', index_col=0)
batch_ofile = open(args.output, 'w')
pwd = os.path.abspath(os.getcwd())
batch_ofile.write('BASE=%s\n\n' % pwd)
for tf, subdata in data.groupby('TF'):
	if not os.path.isdir(tf):
		os.mkdir(tf)
	batch_ofile.write('cd $BASE/%s\n' % tf)
	if len(subdata) < 2:
		get_region = repr(command4).strip("'''")
		top_peak = "head -n 500 %s | awk '%s' > top_combined_peaks.bed" % (subdata['Path'].values[0], get_region)
		batch_ofile.write('%s\n' % top_peak)
		command = ['bedtools', 'getfasta', '-bed', 'top_combined_peaks.bed', '-fi', args.genome]
		batch_ofile.write('%s > top_combined_peaks.fa\n' % ' '.join(command))
		batch_ofile.write('%s\n%s\n\n' % (command2, command3))
	else:
		batch_ofile.write('%s\n%s\n%s\n\n' % (command1, command2, command3))
		biosample = subdata['Biosample'].values
		path = subdata['Path'].values
		records = []
		for x, y in zip(biosample, path):
			records.append( (x, y) )
		grp_file = os.path.join(tf, 'table')
		write_table(grp_file, records)
	print('%s done' % tf)
batch_ofile.close()
print('Please see %s' % sys.argv[2])
