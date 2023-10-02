# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 16:02:43 2021.

@author: Chun-Ping Yu
"""
import re
import numpy as np
import copy
import motif

def read(file_path):
    """
    Read PWMs from a file.

    Parameters
    ----------
    file_path : str
        (full path/)filename for reading motifs.

    Raises
    ------
    RuntimeError
        the format is not simple MEME format.

    Yields
    ------
    name : str
        name of a motif.
    _pwm : a list of a list
        position weight matrix  of the motif.

    """
    with open(file_path) as _infile:
        eof = False
        no_sites = 0
        e_value = 0
        pattern = '[^ ]+'
        while not eof:  # first, get motif name
            while not eof:
                _line = _infile.readline()
                if 'MOTIF' in _line:
                    tokens = re.findall(pattern, _line.rstrip())
                    if len(tokens) > 2:
                        name = tokens[1]
                        description = ' '.join(tokens[2:])
                    else:
                        name = tokens[1]
                        description = ''
                    break
                if len(_line) == 0:
                    eof = True
                    break
            if eof:
                break  # next, get width of the motif

            while not eof:
                _line = _infile.readline()
                if len(_line) == 0:
                    eof = True
                    break
                if 'letter-probability' in _line:
                    # read motif's width
                    seek = re.findall(r'w=\s*(\d+)', _line)
                    if len(seek) != 1:
                        raise RuntimeError('format error in search of the motif width')
                    width = int(seek[0])
                    # read no_sites
                    seek = re.findall(r'nsites=\s*(\d+)', _line)
                    if len(seek) > 0:
                        no_sites = int(seek[0])
                    # and e-value
                    seek = re.findall(r'E=\s*(.+)', _line)
                    if len(seek) > 0:
                        e_value = float(seek[0])
                    break
            # and read a PWM
            count = 0
            _pwm = np.zeros((width, 4))
            while not eof and count < width:
                _line = _infile.readline()
                if len(_line) == 0:
                    raise RuntimeError('Expect the width of motif %s is %d but got %d' % (name, width, count))
                _line = _line.rstrip()
                if len(_line) == 0:
                    continue
                _pwm[count] = [float(x) for x in _line.split()]
                count += 1
            yield name, _pwm, no_sites, e_value, description


def iupac_nucleotide(probabilities, min_prob=0.1):
    """
    Parameters
    ----------
    probabilities : list
        probability for occurring A, C, G, T, respectively
    min_prob : float, optional
        a nucleotide preference has to be greater than it. The default is 0.1.

    Returns
    -------
    str
        IUPAC code.

    """
    nucleotide = list('ACGT')
    pref_nc = []
    for idx, prob in enumerate(probabilities):
        if prob > min_prob:
            pref_nc.append(nucleotide[idx])
    if len(pref_nc) == 1:
        return pref_nc[0]
    temp = ''.join(sorted(pref_nc))
    if temp == 'AG':
        return 'R'
    elif temp == 'CT':
        return 'Y'
    elif temp == 'CG':
        return 'S'
    elif temp == 'AT':
        return 'W'
    elif temp == 'GT':
        return 'K'
    elif temp == 'AC':
        return 'M'
    elif temp == 'CGT':
        return 'B'
    elif temp == 'AGT':
        return 'D'
    elif temp == 'ACT':
        return 'H'
    elif temp == 'ACG':
        return 'V'
    return 'N'


class MEME:
    """
    A class for reading, writing, and handling MEME's type motifs.

    Members
    -------
        read(self, file_path) :
            read PWMs from a file
        add(self, name, pwm, no_sites=0, e_value=0) :
            add a new PWM
        open(self, file_path) :
            open a file for reading position weight matrices
        write(self, file_path, motif_names=None) :
            write PWMs to a file
    """

    def __init__(self, file_path=None):
        """
        Set a path that will read motifs from a file.

        Parameters
        ----------
        file_path : str, optional
            give a file that includes motifs. The default is None.

        Returns
        -------
        None.

        """
        self.data = {}
        self.description = {}
        self.no_sites = {}
        self.e_value = {}
        self.name_in_order = []
        if file_path:
            self.open(file_path)

    def open(self, file_path):
        """
        Read PWMs from a file.

        Parameters
        ----------
        file_path : str
            read a set of PWMs from an external file.
            Note that the available format is supported only for simple MEME format

        Returns
        -------
        None.

        """
        for name, pwm, no_sites, e_value, description in read(file_path):
            self.data[name] = pwm
            self.no_sites[name] = no_sites
            self.e_value[name] = e_value
            self.description[name] = description
            self.name_in_order.append(name)

    def write(self, file_path, motif_names=None):
        """
        Write PWMs to a file.

        Parameters
        ----------
        file_path : str
            write to an external file by giving a (full path/)filename.
        motif_names : list of str
            only write motifs for the giving motif IDs

        Returns
        -------
        None.

        """
        _no_motif = 0
        with open(file_path, 'w') as out_file:
            out_file.write('MEME version 4\n\nALPHABET= ACGT\n\n')
            out_file.write('strands: + -\n\nBackground letter frequencies\n')
            out_file.write('A 0.25 C 0.25 G 0.25 T 0.25\n\n')
            if motif_names is None:
                motif_names = self.data.keys()
            for name in motif_names:
                out_file.write('MOTIF %s\n\n' % name)
                pwm = self.data[name]
                out_file.write('letter-probability matrix: alength= {1} w= {0} nsites= {2} E= {3}\n'.format(
                    *pwm.shape, self.no_sites[name], self.e_value[name]))
                for row in pwm:
                    for col in row:
                        out_file.write('{:>10.6f}'.format(col))
                    out_file.write('\n')
                _no_motif += 1
                out_file.write('\n')
        return _no_motif

    def to_transfac(self, file_path):
        """
        Write PWMs to a file in TRNSFAC format.

        Parameters
        ----------
        file_path : str
            filename for the output data.

        Returns
        -------
        None.

        """
        with open(file_path, 'w') as out_file:
            for name in self.data:
                out_file.write(f'ID {name}\n')
                out_file.write('P0\tA\tC\tG\tT\n')
                pwm = self.data[name]
                for pos, row in enumerate(pwm):
                    out_file.write(f'{pos + 1:02d}')
                    consensus = iupac_nucleotide(row)
                    for prob in row:
                        prob = round(prob * 1e6)
                        out_file.write(f'\t{prob}')
                    out_file.write(f'\t{consensus}\n')
                out_file.write('XX\n//\n')

    def add(self, name, pwm, no_sites=0, e_value=0):
        """
        Add a new PWM into dataset.

        Parameters
        ----------
        name : str
            name of the motif to be included .
        pwm : list of list
            position weight matrix.
        no_sites : int
            number of sites that is used to build the motif
        e_value : float
            E value of the motif

        Returns
        -------
        None.

        """
        self.data[name] = copy.deepcopy(pwm)
        self.no_sites[name] = no_sites
        self.e_value[name] = e_value
        self.name_in_order.append(name)

    def get(self, name):
        """Get a PWM for giving a motif name."""
        return self.data[name]

    def get_consensus(self, name):
        """Get consensus sequence for giving a motif name."""
        max_indices = np.argmax(self.data[name], axis=1)
        nucleotides = list('ACGT')
        consensus = [nucleotides[i] for i in max_indices]
        return ''.join(consensus)

    def trim(self, min_prob=0.3, isIC=True):
        for name in self.data:
            PWM = self.data[name]
            self.data[name] = motif.trim(PWM, min_prob, isIC)

    def get_description(self, name):
        return self.description[name]

    def __len__(self):
        """Return number of PWM(s)."""
        return len(self.data)

    def __iter__(self):
        """Iterate names of PWMs."""
        for name in self.name_in_order:
            yield name

    def name(self):
        """Return all names in the dataset."""
        return self.data.keys()

    def __contains__(self, name):
        """Test a name existence in the current PWMs."""
        return name in self.data


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser('Read/Write motifs from a MEME file')
    parser.add_argument('-m', '--motif', type=str, required=True,
                        help='give the name of a file which includes one or more motifs in simple MEME format')
    parser.add_argument('-n', '--name', type=str,
                        help="if give a list of PWM's names in a file, the output will only show for the given PWMs")
    parser.add_argument('--trim', type=float, default=0,
                        help='trimming bases in the two ends of motifs which have lower ICs (default=0, no trimming)')
    parser.add_argument('--info', action='store_true',
                        help='show the information for each PWM instead of writing the PWM')
    parser.add_argument('-o', '--output', type=str, default='output.meme',
                        help='give a filename for the output (default: output.meme)')
    args = parser.parse_args()

    motifs = MEME(args.motif)
    print(f'Got {len(motifs)} motif(s) from {args.motif}')

    # names = None
    if args.name is not None:
        names = set()
        with open(args.name) as infile:
            for line in infile:
                names.add(line.rstrip())
        print(f'Got {len(names)} motif name(s) from {args.name}')
    else:
        names = None
    if 0 <= args.trim < 2:
        motifs.trim(args.trim)
        print(f'All motifs are trimmed by IC >= {args.trim}')
    else:
        err_msg = f'Invalid value of IC was set ({args.trim}). The available range is  0 <= IC < 2.'
        raise RuntimeError(err_msg)
    if args.info is not None:
        no_motif = 0
        with open(args.output, 'w') as out_file:
            out_file.write('Motif ID\tIC\tMotif length\tConsensus sequence\n')
            for name in motifs:
                if names is not None and name not in names:
                    continue
                pwm = motifs.get(name)
                length = len(pwm)
                ic = sum(motif.IC(pwm))
                consensus = motif.consensus(pwm)
                out_file.write(f'{name}\t{ic}\t{length}\t{consensus}\n')
                no_motif += 1
        print(f'A total of {no_motif} motifs has been processed, please see {args.output}')
    else:
        no_motif = motifs.write(args.output, names)
        print(f'Totally {no_motif} motifs have been written to {args.output}')
    print(f'Done.')
