#!/usr/bin/env python3
"""
This script takes in a bridged Unicycler GFA (004_bridges_applied.gfa) and outputs the length of
the longest bridge in the chromosome. It assumes the largest connected component of the graph is
the chromosome.


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
import sys


def main():
    if len(sys.argv) != 2:
        sys.exit('Error: this script must be given one GFA file as an argument')

    gfa_filename = sys.argv[1]
    segments, links, bridges, lengths = load_gfa(gfa_filename)
    connected_components = get_connected_components(segments, links)
    chromosome_segs = get_largest_component(connected_components, lengths)
    chromosome_bridges = [s for s in chromosome_segs if s in bridges]
    longest_bridge = max(lengths[b] for b in chromosome_bridges)
    print(longest_bridge)


def load_gfa(gfa_filename):
    segments, links, bridges, lengths = set(), collections.defaultdict(set), set(), dict()
    with open(gfa_filename, 'rt') as gfa:
        for line in gfa:
            parts = line.strip().split('\t')
            if line.startswith('S\t'):
                num = int(parts[1])
                segments.add(num)
                lengths[num] = len(parts[2])
                if 'bridge' in line:
                    bridges.add(num)
            elif line.startswith('L\t'):
                num_1 = int(parts[1])
                num_2 = int(parts[3])
                links[num_1].add(num_2)
                links[num_2].add(num_1)
    return segments, links, bridges, lengths


def get_connected_components(segments, links):
    visited = set()
    connected_components = []
    for s in sorted(segments):
        if s not in visited:
            component = dfs(links, s)
            connected_components.append(component)
            for c in component:
                assert c not in visited
                visited.add(c)
    return connected_components


def get_largest_component(connected_components, lengths):
    largest_component, largest_component_length = None, 0
    for component in connected_components:
        component_length = sum(lengths[s] for s in component)
        if component_length > largest_component_length:
            largest_component = component
            largest_component_length = component_length
    return largest_component


def dfs(graph, start):
    visited, stack = set(), [start]
    while stack:
        vertex = stack.pop()
        if vertex not in visited:
            visited.add(vertex)
            stack.extend(graph[vertex] - visited)
    return visited


if __name__ == '__main__':
    main()
