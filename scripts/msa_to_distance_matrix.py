#!/usr/bin/env python3
"""
This script takes a FASTA MSA as input and outputs a PHYLIP distance matrix.

Copyright 2021 Ryan Wick (rrwick@gmail.com)

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version. This program is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should
have received a copy of the GNU General Public License along with this program. If not, see
<https://www.gnu.org/licenses/>.
"""

import gzip
import sys


def main():
    assert len(sys.argv) == 2
    seqs = load_fasta(sys.argv[1])
    print(len(seqs))
    names = []
    distances = {}
    for i, a in enumerate(seqs):
        name_a, seq_a = a
        names.append(name_a)
        for j in range(i, len(seqs)):
            name_b, seq_b = seqs[j]
            assert len(seq_a) == len(seq_b)  # it's an MSA so they must be equal length
            distance = 0
            for k, base_a in enumerate(seq_a):
                base_b = seq_b[k]
                if base_a != base_b:
                    distance += 1
            distances[(name_a, name_b)] = distance
            distances[(name_b, name_a)] = distance
    for name_a in names:
        print(name_a, end='\t')
        distance_to_others = [str(distances[(name_a, name_b)]) for name_b in names]
        print('\t'.join(distance_to_others))


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)
    unknown_file = open(filename, 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('Error: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('Error: cannot use zip format - use gzip instead')
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open


def load_fasta(filename):
    fasta_seqs = []
    with get_open_func(filename)(filename, 'rt') as fasta_file:
        name = ''
        sequence = ''
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    fasta_seqs.append((name.split()[0], sequence.upper()))
                    sequence = ''
                name = line[1:]
            else:
                sequence += line
        if name:
            fasta_seqs.append((name.split()[0], sequence.upper()))
    return fasta_seqs


if __name__ == '__main__':
    main()
