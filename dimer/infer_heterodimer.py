"""
Created on Aug. 08, 2023

@author: Chun-Ping Yu
@mail: chpngyu@gmail.com
"""

import numpy as np
import pandas as pd
import motif
from infer_monomer import read_PPI
from pwm import MEME


def non_overlap_regions(region, width, excluded_region=None, max_overlap=0):
    x1, x2 = region
    for start in range(x1, x2-width+1):
        if excluded_region is not None:
            r1 = set(range(start, start+width))
            r2 = set(range(*excluded_region))
            if len(r1 & r2) > max_overlap:
                continue
        yield start, start+width


class Heterodimer:
    def __init__(self, _protein_net, _tf_family, _monomeric_motifs, _core_motifs,
                 _init_kmer=3, _upper_kmer=6, overlap=0, _minIC=0.3):
        """
        initialized  function

        Parameters
        ----------
        _protein_net:
            a network for a protein-protein interaction
        _tf_family: dict
            TF and its family
        _monomeric_motifs: MEME
            monomeric motifs
        _core_motifs: MEME
            core motifs
        _init_kmer: int
            minimum k-mer (included)
        _upper_kmer: int
            maximum k-mer (included)
        overlap: int
            allow two sub-motifs has an overlap in a motif. default: no (0)
        _minIC: float
            minimum IC for a sub-motif

        Returns
        -------
        Heterodimer
        """
        self.binding_network = _protein_net
        self.TF_family = _tf_family
        self.monomeric_motifs = _monomeric_motifs
        # collect monomeric motifs in TF families
        self.monomeric_TF = {}
        for motif_name in self.monomeric_motifs:
            gene_name = motif_name.split('_')[-1]
            self.monomeric_TF[gene_name] = motif_name
        print(f'There are {len(self.monomeric_motifs)} monomeric motifs in {len(self.monomeric_TF)} TFs')

        self.core_motifs = _core_motifs
        # collect core motifs into TF families
        self.core_motifs_in_family = {}
        # print(self.core_motifs.description)
        for motif_name in self.core_motifs:
            _tf_family = self.core_motifs.get_description(motif_name)
            if _tf_family not in self.core_motifs_in_family:
                self.core_motifs_in_family[_tf_family] = set()
            self.core_motifs_in_family[_tf_family].add(motif_name)
        print(f'There are {len(self.core_motifs)} core motifs in {len(self.core_motifs_in_family)} TF families')

        self.kmer = _init_kmer
        self.upper_kmer = _upper_kmer
        self.overlap = overlap
        self.minIC = _minIC

    def search(self, pwm_target, min_corr=0.8):
        """
        Finding a motif in pwm_reference that is also occurred partially in one of motif of pwm_target

        Parameters
        ----------
        pwm_target : MEME
            a target set of motif
        min_corr : float
            minimum correlation value set to compare the similarity between two motifs

        Returns
        -------
        pandas.DataFrame
            store the results for two motifs
        """
        predictions = []
        total_no = len(pwm_target)
        for no, motif_name1 in enumerate(pwm_target):
            tf1 = motif_name1.split('_')[0]
            print(f'[{no+1}/{total_no}] {tf1}')
            motif1 = pwm_target.get(motif_name1)
            for motif_name2 in self.monomeric_motifs:
                tf2 = motif_name2.split('_')[-1]
                if tf1 == tf2:
                    continue
                if tf1 not in self.binding_network.nodes:
                    interact = 'No'
                elif tf2 in self.binding_network.neighbors(tf1):
                    # no evidence to support that two TFs can be formed a hetero-dimer
                    interact = 'Yes'
                else:
                    interact = 'No'
                motif2 = self.monomeric_motifs.get(motif_name2)
                width = motif2.shape[0]
                if motif1.shape[0] < width:
                    continue
                _correlation, _strands, _aligned_position = motif.auto_correlation(motif1, motif2, fill=False)
                # print(tf1, motif1.shape, tf2, motif2.shape, _correlation)
                if _correlation < min_corr:
                    continue
                # find the other sub-motif that may bind by itself DBD
                _excluded_region = _aligned_position, _aligned_position+width
                other_region = self.find_the_other_sub_motif(tf1, motif1, _excluded_region, min_corr)
                if other_region[1] - other_region[0] < self.kmer:
                    continue
                aligned_position = (_aligned_position+1, _aligned_position+width)
                other_region = (other_region[0]+1, other_region[1])
                record = (tf1, motif_name1, tf2, motif_name2, aligned_position, _correlation, _strands,
                          other_region, interact)
                predictions.append(record)
        columns = ['TF1', 'PWM1 ID', 'TF2', 'PWM2 ID', 'Aligned position of PWM1',
                   'Correlation (cos)', 'Aligned direction', 'self bound', 'BioGRID supported']
        predictions = pd.DataFrame(predictions, columns=columns)
        return predictions

    def find_the_other_sub_motif(self, target_name, target_motif, excluded_region, min_corr=0.8):
        """

        Parameters
        ----------
        target_name : str
            gene symbol
        target_motif : numpy.ndarray
            a motif (PWM) to find the other binding sub-motif that is similar in core motifs or monomeric motifs
        excluded_region : tuple
            do not find the sub-motif in the give region
        min_corr : float
            set minimum value of similarity for comparing two sub-motifs. default=0.8
        Returns
        -------
        tuple
            a possible region that may be bound by a DNA-binding domain. return (-1, -1) if not found
        """

        if target_name in self.monomeric_TF:
            # check the position of the monomeric motif in the target motif
            monomeric_motif = self.monomeric_motifs.get(self.monomeric_TF[target_name])
            _corr, _strands = motif.piece_correlation(target_motif, monomeric_motif)
            max_idx = np.argmax(_corr)
            if _corr[max_idx] < min_corr:
                return -1, -1
            region = max_idx, max_idx + monomeric_motif.shape[0]
            if len(set(range(*region)) & set(range(*excluded_region))) > self.overlap:
                return -1, -1
            return region

        target_family = self.TF_family[target_name]
        if target_family not in self.core_motifs_in_family:
            return -1, -1
        target_motif_width = target_motif.shape[0]
        for window in range(self.upper_kmer, self.kmer - 1, -1):
            for x1, x2 in non_overlap_regions((0, target_motif_width), window, excluded_region, self.overlap):
                sub_motif_target = target_motif[x1:x2]
                for core_motif_name in self.core_motifs_in_family[target_family]:
                    core_motif = self.core_motifs.get(core_motif_name)
                    if core_motif.shape[0] < window:
                        continue
                    _corr, _strand, _aligned_pos = motif.auto_correlation(core_motif, sub_motif_target)
                    if _corr >= min_corr:
                        return _aligned_pos, _aligned_pos+window
        return -1, -1

    def search_old(self, pwm_set, min_corr=0.8):
        predictions = []
        for motif_name1 in pwm_set:
            tf1 = motif_name1.split('_')[0]
            motif1 = pwm_set.get(motif_name1)
            for motif_name2 in pwm_set:
                tf2 = motif_name2.split('_')[0]
                if motif_name1 == motif_name2:
                    continue
                if tf2 not in self.binding_network.neighbors(tf1):
                    # no evidence to support that two TFs can be formed a hetero-dimer
                    continue
                motif2 = pwm_set.get(motif_name2)
                for k in range(self.kmer, self.upper_kmer+1):
                    for i in range(0, motif1.shape[0]-k):
                        sub_motif1 = motif1[i:i+k, :]
                        if motif.information_content(sub_motif1) < self.minIC:
                            continue
                        _correlation, _strands, _aligned_position = motif.auto_correlation(motif2, sub_motif1)
                        if _correlation >= min_corr:
                            record = (tf1, motif_name1, tf2, motif_name2, i+1, _aligned_position+1,
                                      k, _correlation, _strands)
                            predictions.append(record)
        columns = ['TF1', 'PWM1 ID', 'TF2', 'PWM2 ID', 'Aligned position of PWM1', 'Aligned position of PWM2',
                   'Aligned width (k-mer)', 'Correlation (cos)', 'Aligned direction']
        predictions = pd.DataFrame(predictions, columns=columns).sort_values(
                        by=['PWM1 ID', 'PWM2 ID', 'Aligned width (k-mer)', 'Correlation (cos)'],
                        ascending=[True, True, False, False]
                    )
        return predictions.drop_duplicates(['PWM1 ID', 'PWM2 ID'])


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog=f'{__file__}', description='Inferring heterodimer')
    # add options
    parser.add_argument('-m', '--motif', type=str, required=True,
                        help='give a set of motif(s) to be inferred their monomer(s)')
    parser.add_argument('--monomer', type=str, required=True,
                        help='give monomeric motifs to find them in a hetero-motif')
    parser.add_argument('--core', type=str, required=True,
                        help='give core motifs to make sure that one sub-motif is bound by it self in a hetero-motif')
    parser.add_argument('-i', '--interact', type=str, required=True,
                        help='external database to check whether a TF can be formed a homodimer')
    parser.add_argument('-f', '--family', type=str, required=True,
                        help='give a dataframe (csv) which includes genes and their families')
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
    parser.add_argument('-o', '--output', type=str, default='heterodimer.csv',
                        help='give a filename to write the results (default: heterodimer.csv)')
    # parse arguments
    args = parser.parse_args()

    motif_target = MEME(args.motif)
    print(f'Got {len(motif_target)} motifs from {args.motif}')
    monomeric_motifs = MEME(args.monomer)
    print(f'Got {len(monomeric_motifs)} motifs from {args.monomer}')
    core_motifs = MEME(args.core)
    print(f'Got {len(core_motifs)} motifs from {args.core}')

    tf_family = pd.read_csv(args.family, index_col=0)
    tf_family = tf_family.to_dict()[tf_family.columns[0]]
    print(f'Got {len(tf_family)} genes and their families from {args.family}')

    ppi = read_PPI(args.interact)

    if len(args.kmer) == 1:
        min_kmer = args.kmer[0]
        max_kmer = min_kmer
    else:
        min_kmer = args.kmer[0]
        max_kmer = args.kmer[1]

    tf_family = pd.read_csv(args.family, index_col=0)
    tf_family = tf_family.to_dict()[tf_family.columns[0]]
    hetero_dimer = Heterodimer(ppi, tf_family, monomeric_motifs, core_motifs,
                               min_kmer, max_kmer, args.overlap, args.mIC)
    inferred_dimers = hetero_dimer.search(motif_target, args.correlation)

    # inferred_dimers = pd.DataFrame(inferred_dimers)
    inferred_dimers.to_csv(args.output)
    print(f'Done. Please see {args.output}')
