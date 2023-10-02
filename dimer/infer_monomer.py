"""
Created on Aug. 07, 2023

@author: Chun-Ping Yu
@mail: chpngyu@gmail.com
"""
import numpy as np
import motif
import networkx as nx
from pwm import MEME


def read_PPI(filename):
    biogrid = nx.Graph()
    with open(filename) as infile:
        infile.readline()
        for line in infile:
            tokens = line.split('\t')
            p1, p2 = tokens[:2]
            biogrid.add_edge(p1.upper(), p2.upper())
    print(f'There are {biogrid.number_of_nodes()} nodes and {biogrid.number_of_edges()} edges in {filename}')
    return biogrid


class Monomer:
    def __init__(self, _pwm, _init_kmer=3, _upper_kmer=6, overlap=0, _minIC=0.3):
        self.pwm = _pwm
        self.kmer = _init_kmer
        self.upper_kmer = _upper_kmer
        self.overlap = overlap
        self.minIC = _minIC

    def search(self, min_corr=0.8):
        pwm_width = self.pwm.shape[0]
        aligned_motif = {}
        for window in range(self.upper_kmer, self.kmer-1, -1):
            minIC = self.minIC * window
            for start in range(0, pwm_width - 2 * window + self.overlap + 1):
                end = start + window
                sub_pwm1 = self.pwm[start:end, :]
                IC1 = motif.information_content(sub_pwm1)
                if IC1 < minIC:
                    continue
                start2 = start + window - self.overlap
                sub_pwm2 = self.pwm[start2:, :]
                if sub_pwm2.shape[0] < window:
                    continue
                corr, aligned_strand = motif.piece_correlation(sub_pwm2, sub_pwm1)
                # print(start, start2, window)
                # print(corr, aligned_strand)
                if max(corr) >= min_corr:
                    idx = np.argmax(corr)
                    strand = aligned_strand[idx]
                    if strand == 1:
                        aligned_motif['Direction'] = 'repeat'
                    else:
                        aligned_motif['Direction'] = 'palindrome'
                    aligned_motif['Correlation'] = corr[idx]
                    aligned_motif['Width'] = window
                    aligned_motif['Aligned position (Left)'] = start+1
                    aligned_start2 = start2 + idx
                    aligned_end2 = aligned_start2 + window
                    aligned_motif['Aligned position (Right)'] = aligned_start2 +1
                    IC2 = motif.information_content(self.pwm[aligned_start2:aligned_end2])
                    aligned_motif['Left IC'] = IC1
                    aligned_motif['Right IC'] = IC2
                    if IC2 > IC1:
                        aligned_motif['sub-motif'] = self.pwm[aligned_start2:aligned_end2]
                        aligned_motif['Sub-motif used'] = 'Right'
                    else:
                        aligned_motif['sub-motif'] = sub_pwm1
                        aligned_motif['Sub-motif used'] = 'Left'
                if len(aligned_motif) > 0:
                    return aligned_motif
        return aligned_motif


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog=f'{__file__}', description='Infer a monomer in a motif')
    # add options
    parser.add_argument('-m', '--motif', type=str, required=True,
                        help='give a set of motif(s) to be inferred their monomer(s)')
    parser.add_argument('-i', '--interact', type=str, required=True,
                        help='external database to check whether a TF can be formed a homodimer')
    parser.add_argument('-c', '--correlation', type=float, default=0.8,
                        help='only consider a monomer that it is similar with its partner in a same motif'
                             ' (default: 0.8)')
    parser.add_argument('-k', '--kmer', type=int, nargs='+', default=[3],
                        help='set size(s) of k-mer form a to b (default: a=b=3)')
    parser.add_argument('--overlap', type=int, default=0,
                        help='allow two monomers has a maximum length of a overlap with the give value '
                             '(default: 0; no overlap)')
    parser.add_argument('--mIC', type=float, default=0.5,
                        help='only keep sub-motif which has a minimum IC per base (default: 0.5)')
    parser.add_argument('--monomer', type=str, default='monomer.meme',
                        help='output PWMs for the inferred monomers to a file (default: monomer.mem)')
    parser.add_argument('-o', '--output', type=str, default='monomer.txt',
                        help='output file that include name of monomer and its alignment (default: monomer.txt)')
    # parse arguments
    args = parser.parse_args()

    input_motif = MEME(args.motif)
    print(f'Got {len(input_motif)} motifs from {args.motif}')

    ppi = read_PPI(args.interact)

    if len(args.kmer) == 1:
        min_kmer = args.kmer[0]
        max_kmer = min_kmer
    else:
        min_kmer = args.kmer[0]
        max_kmer = args.kmer[1]

    output_motif = MEME()
    columns = ['Width', 'Aligned position (Left)', 'Aligned position (Right)', 'Correlation',
               'Left IC', 'Right IC', 'Sub-motif used', 'Direction']
    out_file = open(args.output, 'w')
    out_file.write('TF\tPWM ID of TF\tMonomer ID\t')
    out_file.write('\t'.join(columns))
    out_file.write('\tHomodimer\n')
    for PWM_ID in input_motif.data:
        pwm = input_motif.get(PWM_ID)
        monomer = Monomer(pwm, min_kmer, max_kmer, args.overlap, args.mIC)
        sub_motif = monomer.search(args.correlation)
        if len(sub_motif) == 0:
            continue
        TF_name = PWM_ID.split('_')[0]
        is_dimer = 'No'
        if TF_name in ppi:
            if TF_name in ppi.neighbors(TF_name):
                monomerID = f'monomer_{TF_name}'
                is_dimer = 'Yes'
            else:
                monomerID = f'sub_motif_{TF_name}'
        else:
            monomerID = f'sub_motif_{TF_name}'

        if is_dimer:
            output_motif.add(monomerID, sub_motif['sub-motif'], no_sites=2)
        out_file.write(f'{TF_name}\t{PWM_ID}\t{monomerID}')
        for column in columns:
            out_file.write(f'\t{sub_motif[column]}')
        out_file.write(f'\t{is_dimer}\n')
    out_file.close()
    output_motif.write(args.monomer)
    print(f'Done. Please see {args.output} and {args.monomer}')
