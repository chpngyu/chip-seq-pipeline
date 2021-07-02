def write_table(filename, records):
	with open(filename, 'w') as ofile:
		ofile.write('Group\tPath\n')
		for record in records:
			ofile.write('%s\t%s\n' % record)

import pandas as pd
import sys
import os
import argparse
parser = argparse.ArgumentParser(description='generate a script for running represnetative peaks')
parser.add_argument('-t', '--table', type=str, required=True,
					help='a table that includes TF and its called peaks')
parser.add_argument('-p', '--path', type=str, required=True,
					help='root path for the python scripts')
parser.add_argument('-g', '--genome', type=str, required=True,
					help='genome assebly file')
parser.add_argument('--blacklist', type=str, required=True,
					help='black list file to remove unwanted peaks')
parser.add_argument('-b', '--biosample', type=str, default='biosample',
					help='generate a table for TFs with their biosamples')
parser.add_argument('-o', '--output', type=str, default='batch_run.sh',
					help='the ouput filename for the generated bash script (default: batch_run.sh)')
args = parser.parse_args()

consensus_script = os.path.join(args.path, 'consensus.py')
command = 'python %s -i table --width 0 --maxd 200 --nosequence --score 2 --path %s' % (consensus_script, args.path)
command += ' --genome %s --blacklist %s' % (args.genome, args.blacklist)

data = pd.read_csv(args.table, sep='\t')
batch_ofile = open(args.output, 'w')
pwd = os.path.abspath(os.getcwd())
batch_ofile.write('BASE=%s\n\n' % pwd)
TF = []
Biosample = []
Path = []
for tf, subdata in data.groupby('TF'):
	if len(subdata) < 2:
		continue
	for biosample, X in subdata.groupby('Biosample'):
		TF.append(tf)
		Biosample.append(biosample)
		if len(X) < 2:
			Path.append(X['Path'].values[0])
			continue
		records = []
		paths = X['Path'].values
		groups = X['Experiment'].values
		new_file = '%s.bed' % biosample
		batch_ofile.write('cd $BASE/%s\n%s\n' % (biosample, command))
		batch_ofile.write('echo "%s done"\n\n' % biosample)
		for s, p in zip(groups, paths):
			records.append( (s, p) )
		y = '%s/%s/top_combined.bed' % (pwd, biosample)
		Path.append(y)
		if not os.path.isdir(biosample):
			os.mkdir(biosample)
		grp_file = os.path.join(biosample, 'table')
		write_table(grp_file, records)
		print('%s done!' % biosample)
batch_ofile.close()
print('Please see %s' % args.output)
merged_data = {'TF': TF, 'Biosample': Biosample, 'Path': Path}
merged_data = pd.DataFrame(merged_data)
merged_data.to_csv(args.biosample, sep='\t')
print('Please see %s' % args.biosample)
print('Done')

