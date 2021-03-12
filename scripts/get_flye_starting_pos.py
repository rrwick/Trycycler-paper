#!/usr/bin/env python3
"""
This script takes two arguments:
* a full Flye assembly before I've extracted and rotated the chromosome
* a chromosome-only Flye assembly that I've rotated.

It outputs the position in the rotated assembly which corresponds to the original (unrotated)
starting position.


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
import subprocess
import sys
import tempfile

MARGIN_SIZE = 30000


def main():
    if len(sys.argv) != 3:
        sys.exit('Error: this script must be given two FASTA files as arguments')

    full_assembly_filename, chromosome_assembly_filename = sys.argv[1], sys.argv[2]
    full_assembly = load_fasta(full_assembly_filename)
    chromosome_assembly = load_fasta(chromosome_assembly_filename)

    assert len(full_assembly) >= 1 and len(chromosome_assembly) == 1

    rotated_chromosome_seq = chromosome_assembly[0]
    unrotated_chromosome_seq = sorted(full_assembly, key=lambda s: 1/len(s))[0]
    unrotated_chromosome_seq_start = unrotated_chromosome_seq[:MARGIN_SIZE]

    with tempfile.TemporaryDirectory() as temp_dir:
        seq_a_filename = temp_dir + '/a.fasta'
        seq_b_filename = temp_dir + '/b.fasta'
        write_to_fasta(rotated_chromosome_seq, 'a', seq_a_filename)
        write_to_fasta(unrotated_chromosome_seq_start, 'b', seq_b_filename)
        stdout, _ = run_command(f'minimap2 -c -x map-ont {seq_a_filename} {seq_b_filename}')
    print(stdout)

    alignments = [Alignment(line) for line in stdout.splitlines()]
    full_alignments = [a for a in alignments if a.query_start == 0 and a.query_end == MARGIN_SIZE]
    assert len(full_alignments) == 1
    alignment = full_alignments[0]
    assert alignment.ref_end > alignment.ref_start

    if alignment.strand == '+':
        print(alignment.ref_start)
    else:
        assert alignment.strand == '-'
        print(alignment.ref_end - 1)


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
    try:
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
                        fasta_seqs.append(sequence)
                        sequence = ''
                    name = line[1:]
                else:
                    sequence += line
            if name:
                fasta_seqs.append(sequence)
        return fasta_seqs
    except FileNotFoundError:
        return []


def write_to_fasta(seq, seq_name, filename):
    with open(filename, 'wt') as fasta:
        fasta.write(f'>{seq_name}\n')
        fasta.write(f'{seq}\n')


def run_command(command):
    result = subprocess.run(command, shell=True, check=True, capture_output=True)
    return result.stdout.decode(), result.stderr.decode()


class Alignment(object):

    def __init__(self, paf_line):
        line_parts = paf_line.strip().split('\t')
        if len(line_parts) < 11:
            sys.exit('\nError: alignment file does not seem to be in PAF format')

        self.query_name = line_parts[0]
        self.query_length = int(line_parts[1])
        self.query_start = int(line_parts[2])
        self.query_end = int(line_parts[3])
        self.strand = line_parts[4]

        self.ref_name = line_parts[5]
        self.ref_length = int(line_parts[6])
        self.ref_start = int(line_parts[7])
        self.ref_end = int(line_parts[8])

        self.matching_bases = int(line_parts[9])
        self.num_bases = int(line_parts[10])
        self.percent_identity = 100.0 * self.matching_bases / self.num_bases

        self.query_cov = 100.0 * (self.query_end - self.query_start) / self.query_length

        self.cigar, self.alignment_score = None, None
        for part in line_parts:
            if part.startswith('cg:Z:'):
                self.cigar = part[5:]
            if part.startswith('AS:i:'):
                self.alignment_score = int(part[5:])

    def __repr__(self):
        return self.query_name + ':' + str(self.query_start) + '-' + str(self.query_end) + \
               '(' + self.strand + '), ' + \
               self.ref_name + ':' + str(self.ref_start) + '-' + str(self.ref_end) + \
               ' (' + ('%.3f' % self.percent_identity) + '%)'


if __name__ == '__main__':
    main()
