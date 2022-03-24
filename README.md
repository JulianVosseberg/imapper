# Spread of introns in proto-eukaryotic paralogs
  
This repository contains scripts to map introns onto protein alignments and create intron tables that can be used for ancestral intron reconstructions in Malin (Csűrös, 2008, http://www.iro.umontreal.ca/~csuros/introns/malin/). We used this for "The spread of the first introns in proto-eukaryotic paralogs" (Vosseberg _et al._, https://doi.org/10.1101/2021.09.28.462139).

## Dependencies
- Python 3.6 or higher

## Map intron positions onto alignment

### Input
- Protein alignment in fasta format
  - Fasta headers are in the format {OG}_{sequence ID}, e.g. >KOG1735_HSAP050599
  - The used alignments are available in figshare: https://doi.org/10.6084/m9.figshare.16601744.
- A database of eukaryotic genomes
  - Genome features files in directory `data_set/gff_files/`
  - Proteome metadata files in directory `data_set/proteomes_metadata` and a concatenated one (`data_set/eukarya_proteomes_metadata.tsv`)
  - These files are available upon request (mailto:b.snel@uu.nl)
- In case of Pfam domains: a file containing the coordinates of the domains in the full-length proteins

### Usage
```
usage: map_introns_alignment.py [-h] [-o outdir] [-e eukarya_path] [-i]
                                [-s shifts] [-p gene%] [-t species%]
                                [-n species] [-d domains]
                                alignment

This script maps intron positions onto a protein alignment.

positional arguments:
  alignment        protein alignment in fasta format

optional arguments:
  -h, --help       show this help message and exit
  -o outdir        directory for output files (default: current)
  -e eukarya_path  directory containing the Eukarya database (default:
                   ~julian/julian2/snel-clan-genomes/eukarya_new)
  -i               infer LECA introns and introns predating duplications
                   (default: off)
  -s shifts        number of nucleotide shifts allowed for an intron (default:
                   0)
  -p gene%         percentage of genes in which an intron at least has to
                   occur to call it a LECA intron (default: 7.5)
  -t species%      percentage of species in which an intron at least has to
                   occur to call it a LECA intron (default: 15)
  -n species       minimum number of Opimoda and Diphoda species that should
                   have an intron to call it a LECA intron (default: 2)
  -d domains       coordinates file
                   (<seqid>\t<start>..<end>(,<start>..<end>,..) for domains
                   instead of full-length genes
```

### Output files
- Table with per intron position and phase in the alignment the corresponding sequence IDs (*_introns.tsv)
- Jalview features file (*_jalview.sff) to show the intron positions on the protein alignment in Jalview (https://www.jalview.org).
- Log file of analysis

## Create intron table

### Input
- File with mapped introns (see above)
- Jalview features file (see above)
- Protein alignment

### Usage
```
usage: create_intron_table.py [-h] [-i species] [-o filename]
                              intronmap features sequences

This script creates an intron table from a file of mapped introns, which can
be used as input for Malin.

positional arguments:
  intronmap             file of mapped introns created by Intron mapper
                        (*_introns.tsv)
  features              Jalview features file as created by Intron mapper
                        (*_jalview.sff)
  sequences             file containing the sequence IDs for each OG in FASTA
                        format

optional arguments:
  -h, --help            show this help message and exit
  -i species, --ignore species
                        species to ignore (comma-separated)
  -o filename, --output filename
                        name of output file (default: introns.table)
```

## Concatenate tables

### Usage
```
usage: concatenate_intron_tables.py [-h] [-o filename] table [table ...]

This script combines multiple intron tables to create one large intron table,
which can be used as input for Malin.

positional arguments:
  table                 multiple intron tables

optional arguments:
  -h, --help            show this help message and exit
  -o filename, --output filename
                        name of output file (default: combined_introns.table)
```

The combined introns tables that were used in the manuscript are available in figshare.
