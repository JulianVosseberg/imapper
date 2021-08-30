#!/usr/bin/env python

import argparse
import sys

# Parse arguments
parser = argparse.ArgumentParser(description = "This script combines multiple intron tables to create one large intron table, which can be used as input for Malin.")
parser.add_argument("table", nargs = '+', help = 'multiple intron tables')
parser.add_argument("-o", "--output", metavar = 'filename', help = 'name of output file (default: combined_introns.table)', default = 'combined_introns.table')
args = parser.parse_args()

family_files = args.table

species_introns = {}
intron_info = []
for i, family_file in enumerate(family_files):
    sys.stderr.write(f'\r{i+1:4} / {len(family_files)} families processing...')
    with open(family_file) as family_introns:
        species_covered = []
        for line in family_introns:
            line = line.rstrip()
            if line[:4] != '#MAP':
                species, introns = line.split('\t')
                species_covered.append(species)

                # Add to previous intron presences or fill up with ? if previous genes were missing
                species_introns[species] = species_introns.get(species, '?'*len(intron_info)) + introns
            else:
                fields = line.split('\t')
                og = fields[2]
                position = fields[3]
                intron_info.append((og, position))

        # Check if all species covered in this family, if not fill up with ?
        for species in species_introns:
            if species not in species_covered:
                species_introns[species] += '?' * len(introns)

sys.stderr.write(f'\nCreating combined intron table of {len(intron_info)} intron positions...\n')
with open(args.output, 'w') as comb_table:
    for species, introns in species_introns.items():
        if len(introns) != len(intron_info):
            sys.exit(f'Error: unequal number of intron positions detected: {species} ({len(introns)})')
        comb_table.write(f'{species}\t{introns}\n')
    for i, intron in enumerate(intron_info):
        og, position = intron
        comb_table.write(f'#MAP\t{i}\t{og}\t{position}\n')
sys.stderr.write('Done!\n')
