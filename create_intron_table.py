#!/usr/bin/env python

import argparse
import sys
import os

# Parse arguments
parser = argparse.ArgumentParser(description = "This script creates an intron table from a file of mapped introns, which can be used as input for Malin.")
parser.add_argument("intronmap", help = 'file of mapped introns created by Intron mapper (*_introns.tsv)')
parser.add_argument("features", help = "Jalview features file as created by Intron mapper (*_jalview.sff)")
parser.add_argument("sequences", help = "file containing the sequence IDs for each OG in FASTA format")
parser.add_argument("-i", "--ignore", metavar = 'species', help = 'species to ignore (comma-separated)')
parser.add_argument("-o", "--output", metavar = 'filename', help = 'name of output file (default: introns.table)', default = 'introns.table')
args = parser.parse_args()

family_file = args.intronmap
species_ignore = []
if args.ignore:
    species_ignore = args.ignore.split(',')
    for species in species_ignore:
        if len(species) != 4:
            sys.exit(f'Error: {species} not recognised!')

all_species = set()
og_sequences = {}
with open(args.sequences) as sequences_file:
    for line in sequences_file:
        if line[0] == '>':
            og_seqid = line[1:].rstrip()
            i_sep = og_seqid.find('_')
            og = og_seqid[:i_sep]
            seqid = og_seqid[i_sep + 1:]
            all_species.add(seqid[:4])
            og_sequences[og] = og_sequences.get(og, []) + [seqid]
    if len(og_sequences) == 0:
        sys.exit(f'Error: {args.sequences} not in FASTA format!')

failed_seqs = []
with open(args.features) as features_file:
    for i in range(4):
        features_file.readline()
    for line in features_file:
        fields = line.rstrip().split('\t')
        if fields[5] == "NA":
            failed_seqs.append(fields[1].split('_')[1])

sp_introns = {sp: [] for sp in all_species if sp not in species_ignore}
intron_info = []

family_base = os.path.basename(family_file)
if family_base.endswith('_merged_introns.tsv'):
    family = family_base[:-19] # Removed _merged_introns.tsv
elif family_base.endswith('_introns.tsv'):
    family = family_base[:-12]
else:
    family = family_base[:family_base.rfind('.')]

with open(family_file) as intron_file:
    pr_og = ""
    intron_file.readline()
    for line in intron_file:
        fields = line.rstrip().split('\t')
        if len(fields) != 4:
            sys.exit(f'Error: format of mapped introns file {family_file} not correct!')
        og = fields[2]
        if og != pr_og:
            species_present = set([seqid[:4] for seqid in set(og_sequences[og]).difference(set(failed_seqs))])
            pr_og = og
        seqids = fields[3].split(',')
        intron_species = set([seqid[:4] for seqid in seqids if seqid[:4] not in species_ignore])
        for sp in sp_introns:
            if sp in intron_species:
                sp_introns[sp].append("1")
            elif sp not in species_present:
                sp_introns[sp].append("?")
            else:
                sp_introns[sp].append("0")
        position = 3 * (int(fields[0]) - 1) + int(fields[1])
        intron_info.append([f'{family}/{og}', position])

sys.stderr.write(f'Creating intron table of {len(intron_info)} intron positions...')
with open(args.output, 'w') as out_file:
    for sp, introns in sp_introns.items():
        if len(introns) != len(intron_info):
            sys.exit(f'\nError: unequal number of intron positions detected: {species} ({len(introns)})')
        print(sp, ''.join(introns), sep = '\t', file = out_file)
    for i, intron in enumerate(intron_info):
        out_file.write(f'#MAP\t{i}\t{intron[0]}\t{intron[1]}\n')
#MAP<tab>#intron<tab>alignment name<tab>intron site position

sys.stderr.write(' Done!\n')
