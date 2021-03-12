#!/usr/bin/env python3
"""
This script aligns two sequences and reports lots of info about them:
* Where repeats are
* The position of errors (places where the sequences differ)
* The number of mismatches, insertions and deletions
* Identity and q score in: the entire sequence, repeat regions, non-repeat regions

All positions given are positions in the alignment, not positions in the sequence (i.e. shifted
by insertions and deletions). Optionally, up to two extra values can be given (positions in the
first and second sequence) and these will be reported in terms of the alignment. This was used for
starting positions - when a contig was rotated from its original starting position, I needed to
know that position in the alignment.

The two input FASTA files are assumed to have one sequence each and to have already been normalised
for strand and starting position.


Copyright 2021 Ryan Wick (rrwick@gmail.com)

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version. This program is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should
have received a copy of the GNU General Public License along with this program. If not, see
<https://www.gnu.org/licenses/>.
"""

import edlib
import gzip
import math
import re
import subprocess
import sys


def main():
    if len(sys.argv) < 3 or len(sys.argv) > 5:
        sys.exit('Error: this script must be given two FASTA files as arguments')

    assembly_1_filename, assembly_2_filename = sys.argv[1], sys.argv[2]
    assembly_1, assembly_2 = load_fasta(assembly_1_filename), load_fasta(assembly_2_filename)

    try:
        starting_pos_1 = int(sys.argv[3])
    except IndexError:
        starting_pos_1 = None
    try:
        starting_pos_2 = int(sys.argv[4])
    except IndexError:
        starting_pos_2 = None

    assembly_1_repeat_ranges = get_repeat_ranges(assembly_1_filename)
    assembly_2_repeat_ranges = get_repeat_ranges(assembly_2_filename)

    if len(assembly_1) == 0 or len(assembly_2) == 0:
        sys.exit(f'Error: empty assembly')

    if len(assembly_1) > 1:
        sys.exit(f'Error: {assembly_1_filename} has more than one sequence')
    if len(assembly_2) > 1:
        sys.exit(f'Error: {assembly_2_filename} has more than one sequence')

    seq_1, seq_2 = assembly_1[0], assembly_2[0]
    assert seq_1[:10] == seq_2[:10]
    assert seq_1[-10:] == seq_2[-10:]

    result = edlib.align(seq_1, seq_2, mode="NW", task="path")
    cigar = result['cigar']
    expanded_cigar = get_expanded_cigar(cigar)

    length = len(expanded_cigar)
    print(f'length\t{length}')

    repeats, unaligned_to_aligned_1, unaligned_to_aligned_2 = \
        get_alignment_repeat_ranges(expanded_cigar, seq_1, seq_2,
                                    assembly_1_repeat_ranges, assembly_2_repeat_ranges)
    for start, end in repeats.ranges:
        print(f'repeat\t{start}\t{end}')

    errors, mismatches, insertions, deletions = get_error_positions(expanded_cigar)
    for e in errors:
        print(f'error\t{e}')
    print(f'mismatches\t{mismatches}')
    print(f'insertions\t{insertions}')
    print(f'deletions\t{deletions}')

    print_accuracy_stats(repeats, errors, length)

    if starting_pos_1 is not None:
        aligned_starting_pos_1 = unaligned_to_aligned_1[starting_pos_1]
        print(f'starting_pos_1\t{aligned_starting_pos_1}')
    if starting_pos_2 is not None:
        aligned_starting_pos_2 = unaligned_to_aligned_2[starting_pos_2]
        print(f'starting_pos_2\t{aligned_starting_pos_2}')


def get_alignment_repeat_ranges(expanded_cigar, seq_1, seq_2,
                                assembly_1_repeat_ranges, assembly_2_repeat_ranges):
    i, j, k = 0, 0, 0
    unaligned_to_aligned_1, unaligned_to_aligned_2 = {}, {}
    for c in expanded_cigar:
        if c == '=':
            unaligned_to_aligned_1[i] = k
            unaligned_to_aligned_2[j] = k
            assert seq_1[i] == seq_2[j]
            i += 1
            j += 1
        elif c == 'X':
            unaligned_to_aligned_1[i] = k
            unaligned_to_aligned_2[j] = k
            assert seq_1[i] != seq_2[j]
            i += 1
            j += 1
        elif c == 'I':
            unaligned_to_aligned_1[i] = k
            i += 1
        elif c == 'D':
            unaligned_to_aligned_2[j] = k
            j += 1
        k += 1

    aligned_repeats = IntRange()
    for start, end in assembly_1_repeat_ranges.ranges:
        aligned_start, aligned_end = unaligned_to_aligned_1[start], unaligned_to_aligned_1[end]
        aligned_repeats.add_range(aligned_start, aligned_end)
    for start, end in assembly_2_repeat_ranges.ranges:
        aligned_start, aligned_end = unaligned_to_aligned_2[start], unaligned_to_aligned_2[end]
        aligned_repeats.add_range(aligned_start, aligned_end)
    return aligned_repeats, unaligned_to_aligned_1, unaligned_to_aligned_2


def get_error_positions(expanded_cigar):
    errors = []
    mismatches, insertions, deletions = 0, 0, 0
    for i, c in enumerate(expanded_cigar):
        if c == '=':
            pass
        elif c == 'X':
            mismatches += 1
            errors.append(i)
        elif c == 'I':
            insertions += 1
            errors.append(i)
        elif c == 'D':
            deletions += 1
            errors.append(i)
    return errors, mismatches, insertions, deletions


def get_expanded_cigar(cigar):
    expanded_cigar = []
    cigar_parts = re.findall(r'\d+[IDX=]', cigar)
    for p in cigar_parts:
        size = int(p[:-1])
        letter = p[-1]
        expanded_cigar.append(letter * size)
    return ''.join(expanded_cigar)


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


def get_repeat_ranges(filename):
    _, _ = run_command(f'nucmer --maxmatch --nosimplify --prefix=seq_seq {filename} {filename}')
    out, _ = run_command('show-coords -H -T -r seq_seq.delta | awk \'{if ($2-$1 < 1000000 && '
                         '$2-$1 > 1000 && $1 != $3 && $2 != $4 && $7 > 95) print $1"\t"$2;}\'')
    repeat_ranges = IntRange()
    for line in out.splitlines():
        parts = line.split('\t')
        assert len(parts) == 2
        start, end = int(parts[0]), int(parts[1])
        repeat_ranges.add_range(start, end)
    _, _ = run_command('rm seq_seq.delta')
    return repeat_ranges


def print_accuracy_stats(repeats, errors, length):
    error_count = len(errors)
    print(f'overall_errors\t{error_count}/{length}')

    overall_identity = 1.0 - (error_count / length)
    print(f'overall_identity\t{100.0 * overall_identity}%')
    print(f'overall_q_score\t{get_q_score(overall_identity)}')

    repeat_error_count, not_repeat_error_count = 0, 0
    for e in errors:
        if repeats.contains(e):
            repeat_error_count += 1
        else:
            not_repeat_error_count += 1
    repeat_length = repeats.total_length()
    not_repeat_length = length - repeat_length

    print(f'not_repeat_errors\t{not_repeat_error_count}/{not_repeat_length}')
    not_repeat_identity = 1.0 - (not_repeat_error_count / not_repeat_length)
    print(f'not_repeat_identity\t{100.0 * not_repeat_identity}%')
    print(f'not_repeat_q_score\t{get_q_score(not_repeat_identity)}')

    print(f'repeat_errors\t{repeat_error_count}/{repeat_length}')
    repeat_identity = 1.0 - (repeat_error_count / repeat_length)
    print(f'repeat_identity\t{100.0 * repeat_identity}%')
    print(f'repeat_q_score\t{get_q_score(repeat_identity)}')


def get_q_score(identity):
    if identity == 1.0:
        return 'Qinf'
    else:
        q_score = -10.0 * math.log10(1 - identity)
        return f'Q{q_score:.1f}'


def run_command(command):
    result = subprocess.run(command, shell=True, check=True, capture_output=True)
    return result.stdout.decode(), result.stderr.decode()


class IntRange(object):
    """
    This class contains one or more integer ranges. Overlapping ranges will be merged together.
    It stores its ranges in a Python-like fashion where the last value in each range is
    exclusive.
    """
    def __init__(self, ranges=None):
        if not ranges:
            ranges = []
        self.ranges = []
        self.add_ranges(ranges)
        self.simplify()

    def __repr__(self):
        return str(self.ranges)

    def add_range(self, start, end):
        """Adds a single range."""
        self.add_ranges([(start, end)])

    def add_ranges(self, ranges):
        """Adds multiple ranges (list of tuples)."""
        self.ranges += ranges
        self.simplify()

    def total_length(self):
        """Returns the number of integers in the ranges."""
        return sum([x[1] - x[0] for x in self.ranges])

    def simplify(self):
        """Collapses overlapping ranges together."""
        fixed_ranges = []
        for int_range in self.ranges:
            if int_range[0] > int_range[1]:
                fixed_ranges.append((int_range[1], int_range[0]))
            elif int_range[0] < int_range[1]:
                fixed_ranges.append(int_range)
        starts_ends = [(x[0], 1) for x in fixed_ranges]
        starts_ends += [(x[1], -1) for x in fixed_ranges]
        starts_ends.sort(key=lambda z: z[0])
        current_sum = 0
        cumulative_sum = []
        for start_end in starts_ends:
            current_sum += start_end[1]
            cumulative_sum.append((start_end[0], current_sum))
        prev_depth = 0
        start = 0
        combined = []
        for pos, depth in cumulative_sum:
            if prev_depth == 0:
                start = pos
            elif depth == 0:
                combined.append((start, pos))
            prev_depth = depth
        self.ranges = combined

    def overlaps(self, other):
        """Returns True if the other IntRange overlaps with this IntRange."""
        for other_range in other.ranges:
            other_start, other_end = other_range
            for this_start, this_end in self.ranges:
                if (this_start <= other_start < this_end) or (this_start < other_end <= this_end):
                    return True
        return False

    def contains(self, i):
        """Returns True if i is contained in any of the ranges."""
        for start, end in self.ranges:
            if start <= i < end:
                return True
        return False


if __name__ == '__main__':
    main()
