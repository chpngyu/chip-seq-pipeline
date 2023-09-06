# -*- coding: utf-8 -*-
"""
Created on Mon Oct 4 2021

@author: Chun-Ping Yu
"""

def read_motif_IC(filename):
	tf_motif2IC = {}
	with open(filename) as infile:
		infile.readline()
		for line in infile:
			tokens = line.rstrip().split('\t')
			motif, IC = tokens[0], float(tokens[1])
			tf = motif.split('_')[0]
			if tf not in tf_motif2IC:
				tf_motif2IC[tf] = {}
			tf_motif2IC[tf][motif] = IC
	return tf_motif2IC

import gzip
def read_correlation(filename):
	data = {}
	with gzip.open(filename, 'rt') as infile:
		for line in infile:
			tokens = line.rstrip('\n').split('\t')
			if len(tokens) < 3:
				continue
			data[(tokens[0], tokens[1])] = float(tokens[2])
	return data

import networkx as nx
def grouping(motifs, corr, minc):
	if len(motifs) == 1:
		return [(motifs[0], 'G1')]
	elif len(motifs) == 2:
		if corr[(motifs[0], motifs[1])] > minc:
			return [(motifs[0], 'G1'), (motifs[1], 'G1')]
		else:
			return [(motifs[0], 'G1'), (motifs[1], 'G2')]
	G = nx.Graph()
	G.add_nodes_from(motifs)
	for i, m1 in enumerate(motifs):
		for m2 in motifs[i+1:]:
			if corr[(m1, m2)] > minc:
				G.add_edge(m1, m2)
	res = []
	for no, subg in enumerate(nx.connected_components(G)):
		for m in subg:
			res.append((m, f'G{no+1}'))
	return res

import argparse
parser = argparse.ArgumentParser('grouping motifs')
# add options
parser.add_argument('-i', '--input', type=str, required=True,
					help='motif ICs')
parser.add_argument('-c', '--correlation', type=str, required=True,
					help='pairwise correlations between motifs')
parser.add_argument('--minc', type=float, default=0.8,
					help='minimum correlation to group a motif (default: 0.8)')
parser.add_argument('-o', '--output', type=str, default='motif_group.txt',
					help='filename for the output')
# parse arguments
args = parser.parse_args()

tf_motif2IC = read_motif_IC(args.input)
print(f'Got {len(tf_motif2IC)} TF(s) from {args.input}')

corr = read_correlation(args.correlation)
print(f'Got {len(corr)} pairs of corrlation from {args.correlation}')

ofile = open(args.output, 'w')
for tf in tf_motif2IC:
	motifs = [m[0] for m in tf_motif2IC[tf].items()]
	for motif, group in grouping(motifs, corr, args.minc):
		ofile.write(f'{tf}\t{motif}\t{group}\n')
ofile.close()
print(f'Done. Please see {args.output}')
