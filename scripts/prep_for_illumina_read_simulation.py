#!/usr/bin/env python3
"""
This script takes a FASTA as input and outputs the same sequences multiplied to higher copy
numbers. This is to facilitate:
* circular Illumina read simulation (i.e. reads over the start/end break).
* different read depths for different replicons

The number of times each sequence is copied is 10 times its depth, rounded to the nearest integer.
So a chromosome with 1x depth will get 10 copies, a plasmid with 3.11x depth will get 31 copies,
and so on.

Run it like this:
python3 prep_for_illumina_read_simulation.py input.fasta > multicopy.fasta


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
import re
import sys


def main():
    if len(sys.argv) != 2:
        sys.exit('Error: this script must be given one FASTA file as an argument')

    for name, seq in load_fasta(sys.argv[1]):
        assert 'circular=true' in name
        depth_strings = re.findall(r'depth=[\d.]+', name)
        assert len(depth_strings) == 1
        depth = float(re.findall(r'[\d.]+', depth_strings[0])[0])
        repeat_count = int(round(depth * 10.0))
        multi_seq = ''.join([seq] * repeat_count)
        print(f'>{name}')
        print(multi_seq)


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
        sequence = []
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    fasta_seqs.append((name, ''.join(sequence)))
                    sequence = []
                name = line[1:]
            else:
                sequence.append(line)
        if name:
            fasta_seqs.append((name, ''.join(sequence)))
    return fasta_seqs


if __name__ == '__main__':
    main()
