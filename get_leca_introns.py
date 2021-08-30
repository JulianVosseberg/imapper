#!/usr/bin/env python

import argparse
parser = argparse.ArgumentParser(description = "This script extracts the LECA introns from an intron table and a corresponding site histories file.")
parser.add_argument("table", help = 'intron table')
parser.add_argument("histories", help = 'site histories file')
parser.add_argument("-t", metavar = 'threshold', help = 'LECA posterior probability (default 0.5)', default = 0.5)
parser.add_argument("-p", help = 'include Pfam in OG name', action = 'store_true')
args = parser.parse_args()

introns = {}
with open(args.table) as intron_table:
    for line in intron_table:
        if line[:4] != '#MAP':
            continue
        fields = line.rstrip().split('\t')
        if args.p:
            og = fields[2]
        else:
            og = fields[2].split('/')[1]
        position = int(fields[3])
        phase = position % 3
        if phase == 0:
            location = f"{position // 3}'3"
        else:
            location = f"{position // 3 + 1}'{phase}"
        introns[fields[1]] = {'Family': og, 'Location': location}

with open(args.histories) as site_histories:
    for i in range(3):
        site_histories.readline()
    colnames = site_histories.readline().split('\t')
    i = colnames.index('LECA-present')
    #check = site_histories.readline().split('\t')[317]
    #if check != 'LECA-present':
    #    print(check)
    #    sys.exit('Check if file is correct.')
    for line in site_histories:
        fields = line.rstrip().split('\t')
        site = fields[0]
        leca_pr = float(fields[i])
        if leca_pr >= args.t:
            info = introns[site]
            print(f"{site}\t{info['Family']}\t{info['Location']}\t{round(leca_pr, 1)}")
