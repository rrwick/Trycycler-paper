#!/usr/bin/env python3
"""
This script takes two or more FASTA as input, aligns all pairwise comparisons between the
sequences and outputs a summary of the types of differences found.

Copyright 2021 Ryan Wick (rrwick@gmail.com)

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version. This program is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should
have received a copy of the GNU General Public License along with this program. If not, see
<https://www.gnu.org/licenses/>.
"""

import collections
import edlib
import itertools
import gzip
import re
import sys


MAX_INDEL = 16


def main():
    assert len(sys.argv) >= 3
    fasta_filenames = sys.argv[1:]
    names, seqs = load_all_seqs(fasta_filenames)
    print_header()
    for name_a, name_b in itertools.combinations(names, 2):
        process_seq_pair(name_a, name_b, seqs[name_a], seqs[name_b])


def load_all_seqs(fasta_filenames):
    seqs = []
    for filename in fasta_filenames:
        s = load_fasta(filename)
        assert len(s) == 1
        seqs.append(s[0])
    names = [s[0] for s in seqs]
    assert len(names) == len(set(names))  # ensure no duplicate names
    return names, dict(seqs)


def process_seq_pair(name_a, name_b, seq_a, seq_b):
    counts = collections.defaultdict(int)
    assert 0.9 < len(seq_a) / len(seq_b) < 1.1  # sequences should be similar length
    seq_a, seq_b = align_sequences(seq_a, seq_b)
    assert len(seq_a) == len(seq_b)  # now aligned so must be equal length

    sub_positions = []
    indel_positions_by_size = collections.defaultdict(list)

    i = 0
    while True:
        a, b = seq_a[i], seq_b[i]
        if a != b:
            if a != '-' and b != '-':
                sub_positions.append(i)
            else:
                assert (a == '-') != (b == '-')
                indel_start, indel_end = i, i+1
                if a == '-':
                    while seq_a[indel_end] == '-':
                        indel_end += 1
                elif b == '-':
                    while seq_b[indel_end] == '-':
                        indel_end += 1
                indel_len = indel_end - indel_start
                indel_positions_by_size[indel_len].append(indel_start)
                i += (indel_len - 1)
        i += 1
        if i >= len(seq_a):
            break

    counts['sub'] += len(sub_positions)
    for indel_size, indel_positions in indel_positions_by_size.items():
        for indel_position in indel_positions:
            if is_indel_homopolymer(indel_position, indel_size, seq_a, seq_b):
                counts[f'{indel_size}_homo'] += 1
            else:
                counts[str(indel_size)] += 1
                if is_methyl_motif(indel_position, seq_a, seq_b):
                    counts['methyl'] += 1
                else:
                    counts['non-methyl'] += 1

    for sub_position in sub_positions:
        if is_methyl_motif(sub_position, seq_a, seq_b):
            counts['methyl'] += 1
        else:
            counts['non-methyl'] += 1

    print_counts(name_a, name_b, counts)


def align_sequences(seq_a, seq_b):
    result = edlib.align(seq_a, seq_b, mode="NW", task="path")
    cigar = result['cigar']
    expanded_cigar = get_expanded_cigar(cigar)
    i, j = 0, 0
    aligned_a, aligned_b = [], []
    for c in expanded_cigar:
        if c == '=':
            aligned_a.append(seq_a[i])
            aligned_b.append(seq_b[j])
            assert seq_a[i] == seq_b[j]
            i += 1
            j += 1
        if c == 'X':
            aligned_a.append(seq_a[i])
            aligned_b.append(seq_b[j])
            assert seq_a[i] != seq_b[j]
            i += 1
            j += 1
        if c == 'I':
            aligned_a.append(seq_a[i])
            aligned_b.append('-')
            i += 1
        if c == 'D':
            aligned_a.append('-')
            aligned_b.append(seq_b[j])
            j += 1
    assert i == len(seq_a) and j == len(seq_b)
    return aligned_a, aligned_b


def is_indel_homopolymer(i, indel_size, seq_a, seq_b):
    """
    Checks to see if the indel in question is a homopolymer. In this function, the deletion can be
    on either of the two sequences.
    """
    if seq_a[i] == '-':
        for j in range(i, i+indel_size):
            assert seq_a[i] == '-'
            assert seq_b[i] != '-'
            return is_indel_homopolymer_2(i, indel_size, seq_a, seq_b)
    elif seq_b[i] == '-':
        for j in range(i, i+indel_size):
            assert seq_a[i] != '-'
            assert seq_b[i] == '-'
            return is_indel_homopolymer_2(i, indel_size, seq_b, seq_a)
    else:
        assert False


def is_indel_homopolymer_2(i, indel_size, seq_a, seq_b):
    """
    Checks to see if the indel in question is a homopolymer. In this function, the deletion is on
    sequence A.
    """
    bases = []
    for j in range(i, i+indel_size):
        assert seq_a[j] == '-'
        assert seq_b[j] != '-'
        bases.append(seq_b[j])

    # If there are multiple different bases in this indel, it can't be a homopolymer.
    if len(set(bases)) > 1:
        return False

    base = bases[0]
    assert all(b == base for b in bases)

    if seq_b[i-1] == base and seq_b[i-2] == base and seq_b[i-3] == base:
        return True
    j = i + indel_size - 1
    if seq_b[j+1] == base and seq_b[j+2] == base and seq_b[j+3] == base:
        return True

    return False


def is_methyl_motif(i, seq_a, seq_b):
    """
    Checks to see if the error in question (can be a substitution or an indel) is in a methylation
    motif. This function assumes that the second sequence is the polished one, so it looks for the
    motif there.
    """
    middle = seq_b[i]
    
    before = ''
    j = i - 1
    while len(before) < 4:
        before = seq_b[j] + before
        before = before.replace('-', '')
        j -= 1

    after = ''
    j = i + 1
    while len(after) < 4:
        after = after + seq_b[j]
        after = after.replace('-', '')
        j += 1

    ref_seq = before + middle + after
    ref_seq = ref_seq.replace('-', '')

    return 'CCAGG' in ref_seq or 'CCTGG' in ref_seq or 'GATC' in ref_seq


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


def get_expanded_cigar(cigar):
    expanded_cigar = []
    cigar_parts = re.findall(r'\d+[IDX=]', cigar)
    for p in cigar_parts:
        size = int(p[:-1])
        letter = p[-1]
        expanded_cigar.append(letter * size)
    return ''.join(expanded_cigar)


def print_header():
    print(f'assembly_1\tassembly_2\tsubstitutions', end='')
    for i in range(1, MAX_INDEL+1):
        print(f'\t{i}_bp_homopolymer_indels', end='')
        print(f'\t{i}_bp_non-homopolymer_indels', end='')
    print(f'\tmethylation_motif_errors', end='')
    print(f'\tnon-methylation_motif_errors', end='')
    print('', flush=True)


def print_counts(name_a, name_b, counts):
    indel_sizes = []
    for difference_type in counts:
        if difference_type.isdigit():
            indel_sizes.append(int(difference_type))
        if '_homo' in difference_type:
            indel_sizes.append(int(difference_type.replace('_homo', '')))
    max_indel = sorted(indel_sizes)[-1]

    if max_indel > MAX_INDEL:
        sys.exit(f'Error: maximum indel too large: {max_indel}')

    print(f'{name_a}\t{name_b}', end='')

    sub_count = counts['sub']
    print(f'\t{sub_count}', end='')
    for i in range(1, MAX_INDEL+1):
        indel_count = counts[f'{i}_homo']
        print(f'\t{indel_count}', end='')
        indel_count = counts[str(i)]
        print(f'\t{indel_count}', end='')

    methyl_count = counts['methyl']
    print(f'\t{methyl_count}', end='')
    non_methyl_count = counts['non-methyl']
    print(f'\t{non_methyl_count}', end='')

    print('', flush=True)


if __name__ == '__main__':
    main()
