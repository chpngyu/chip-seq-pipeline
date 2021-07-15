# python score_quantile.py input.bed output.bed
#import numpy as np
from scipy.stats import rankdata
import argparse

parser = argparse.ArgumentParser(description='normalizing macs2 scores')
parser.add_argument('-s', '--score', type=str, required=True,
					help='summits.bed file by macs2')
parser.add_argument('-g', '--group', type=str, required=True,
					help='sample ID for identifying this peaks')
parser.add_argument('-e', '--extent', type=int, default=0,
					help='length for the extent peak region (default: 0)')
parser.add_argument('-o', '--output', type=str, default='score.bed',
					help='output filename (default: score.bed)')
# parse command line
args = parser.parse_args()

peaks = []
scores = []
with open(args.score) as infile:
	for no, line in enumerate(infile):
		tokens = line.rstrip().split('\t')
		start = int(tokens[1])
		start = max(1, start-args.extent)
		stop = int(tokens[2]) + args.extent
		peak = '%s_peak_%d' % (args.group, no+1)
		peaks.append( (tokens[0], start, stop, peak) )
		scores.append( float(tokens[4]) )
percentile = rankdata(scores)/len(scores)

ofile = open(args.output,'w')
for pos, pert in zip(peaks, percentile):
	ofile.write('%s\t%s\t%s\t%s\t' % pos)
	ofile.write('%g\n' % pert)
ofile.close()
