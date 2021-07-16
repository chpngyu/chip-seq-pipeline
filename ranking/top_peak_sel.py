import argparse

parser = argparse.ArgumentParser(description='selection of top peaks in cluster')
parser.add_argument('-i', '--input', type=str, required=True,
					help='a input file which includes peak regions be scoring and clutering')
parser.add_argument('--top', type=int,
					help='no. of top peaks to putput (default: all)')
parser.add_argument('--score', type=int, choices=[1,2], default=1,
					help='output score for individual (1) peak or summation (2) of peak cluster (default: individual score)')
parser.add_argument('-o', '--output', type=str, default='top_peaks.bed',
					help='output filename (default: top_peaks.bed)')
parser.add_argument('--log', type=str,
					help='show summation score for each cluster (default: no)')

# parse command line
args = parser.parse_args()

cluster_peak = {}
with open(args.input) as infile:
	for line in infile:
		tokens = line.rstrip().split('\t')
		pos = tuple(tokens[:4])
		score = float(tokens[4])
		cluster = tokens[5]
		# get group name
		group_pos = tokens[3].find('_peak')
		group = tokens[3][:group_pos]
		if cluster not in cluster_peak:
			cluster_peak[cluster] = {group: (pos, score)}
		elif group not in cluster_peak[cluster]:
			cluster_peak[cluster][group] = (pos, score)
		elif score > cluster_peak[cluster][group][1]:
			cluster_peak[cluster][group] = (pos, score)

cluster_score = []
for cls in cluster_peak:
	sum_score = 0
	for grp in cluster_peak[cls]:
		sum_score += cluster_peak[cls][grp][1]
	cluster_score.append( (cls, sum_score) )
cluster_score = sorted(cluster_score, key=lambda c: c[1], reverse=True)

ofile = open(args.output, 'w')
if args.top:
	no_peak = args.top
else:
	no_peak = len(cluster_score)
for cluster, sum_score in cluster_score[:no_peak]:
	best_pos = None
	max_score = -1
	for grp in cluster_peak[cluster]:
		if cluster_peak[cluster][grp][1] > max_score:
			best_pos, max_score = cluster_peak[cluster][grp]
	ofile.write('%s\t%s\t%s\t%s\t' % best_pos)
	if args.score == 1:
		ofile.write('%g\t' % max_score)
	else:
		ofile.write('%g\t' % sum_score)
	ofile.write('%s\n' % cluster)
ofile.close()

if args.log:
	logf = open(args.log, 'w')
	logf.write('Peak cluster\tSummation score\tNo. peak(s)\n')
	for cluster, sum_score in cluster_score:
		logf.write('%s\t%g\t%d\n' % (cluster, sum_score, len(cluster_peak[cluster])))
	logf.close()
