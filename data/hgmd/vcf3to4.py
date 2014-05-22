#!/usr/bin/env python

"""

"""

# Author: Orion Buske
# Date:   246 Sept 2013
import os
import sys
import pickle

from collections import defaultdict


def load_genome(filename):
    """Return dict: sequence -> str from multiFASTA file"""
    genome = {}
    name = None
    sequences = []
    print("Loading genome from FASTA file: {}".format(filename), 
          file=sys.stderr)
    with open(filename) as ifp:
        for line in ifp:
            line = line.strip()
            if line.startswith('>'):
                if name:
                    genome[name] = ''.join(sequences)

                name = line.lstrip('>').strip()
                print("- {}".format(name), file=sys.stderr)
                sequences = []
            elif line:
                sequences.append(line)

        if name:
            genome[name] = ''.join(sequences)

    print("Found {} sequences.".format(len(genome)), file=sys.stderr)
    return genome

def convert3to4(ifp, ofp, genome):
    n_dropped = 0
    for line in ifp:
        line = line.strip()
        if line.startswith('#'):
            if line.startswith('##fileformat=VCFv3'):
                print('##fileformat=VCFv4.0', file=ofp)
            else:
                print(line, file=ofp)
            continue

        tokens = line.split('\t')
        if tokens[6] == '0':
            tokens[6] = 'PASS'

        chrom = tokens[0]
        chrom = chrom if chrom.startswith('chr') else ('chr' + chrom)
        pos = int(tokens[1])
        ref = tokens[3]
        alt = tokens[4]
        if ',' in alt:
            if '2' in tokens[9].split(':', 1)[0]:
                # Drop sites with two alternative haploytypes
                n_dropped += 1
                continue
            alt = alt.split(',')[0]


        if alt.startswith('I'):
            alt = ref + alt.lstrip('I')
        elif alt.startswith('D'):
            n = int(alt.lstrip('D'))
            pos -= 1
            chromosome = genome.get(chrom)
            if not chromosome:
                n_dropped += 1
                continue

            ref = chromosome[pos - 1: pos + n]
            try:
                alt = ref[0]
            except IndexError:
                n_dropped += 1
                pass

        tokens[0] = chrom
        tokens[1] = str(pos)
        tokens[3] = ref
        tokens[4] = alt
        print('\t'.join(tokens), file=ofp)

    print('Dropped {} entries'.format(n_dropped), file=sys.stderr)

def script(genome_fasta, outdir, vcfs, pickle_genome=None):
    if not pickle_genome or not os.path.isfile(pickle_genome):
        genome = load_genome(genome_fasta)
        if pickle_genome:
            with open(pickle_genome, 'wb') as ofp:
                print('Saving genome: {}'.format(pickle_genome), 
                      file=sys.stderr)
                pickle.dump(genome, ofp, -1)
    elif pickle_genome and os.path.isfile(pickle_genome):
        with open(pickle_genome, 'rb') as ifp:
            print('Loading saved genome: {}'.format(pickle_genome), 
                  file=sys.stderr)
            genome = pickle.load(ifp)

    if not os.path.isdir:
        os.makedirs(outdir)
        print("Created directory: {}".format(outdir), file=sys.stderr)

    for vcf in vcfs:
        outfile = os.path.join(outdir, os.path.basename(vcf))
        with open(vcf) as ifp:
            if os.path.isfile(outfile):
                print("Already exists: {}".format(outfile), file=sys.stderr)
                continue

            with open(outfile, 'w') as ofp:
                print("{} -> {}".format(vcf, outfile), file=sys.stderr)
                convert3to4(ifp, ofp, genome)

def parse_args(args):
    from argparse import ArgumentParser
    description = __doc__.strip()
    
    parser = ArgumentParser(description=description)
    parser.add_argument('genome_fasta', metavar='GENOME', help="FASTA file")
    parser.add_argument('outdir', metavar='OUTDIR')
    parser.add_argument("vcfs", nargs='+', metavar='VCF3')
    parser.add_argument("--pickle-genome", metavar='PKL')

    return parser.parse_args(args)

def main(args=sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
