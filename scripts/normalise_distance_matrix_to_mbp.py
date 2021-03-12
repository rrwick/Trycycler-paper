#!/usr/bin/env python3
"""
This script takes a PHYLIP distane matrix and a genome size as input and outputs another PHYLIP
distance matrix, where the values are per-Mbp.

Copyright 2021 Ryan Wick (rrwick@gmail.com)

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version. This program is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should
have received a copy of the GNU General Public License along with this program. If not, see
<https://www.gnu.org/licenses/>.
"""

import sys


def main():
    assert len(sys.argv) == 3
    matrix_filename = sys.argv[1]
    genome_size = int(sys.argv[2])
    ratio = 1000000 / genome_size
    with open(matrix_filename, 'rt') as matrix:
        for line in matrix:
            parts = line.strip().split('\t')
            if len(parts) == 1:  # first line:
                print(parts[0])
            else:
                print(parts[0], end='')
                for part in parts[1:]:
                    new_val = int(part) * ratio
                    print(f'\t{new_val:.1f}', end='')
                print()


if __name__ == '__main__':
    main()
