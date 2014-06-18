#!/usr/bin/env python

#Given a directory, get the average number of genes listed in ezr ranking files

__author__ = 'Tal Friedman'

import os
import sys
import logging

def script(ezr_paths):
    for ezr_path in ezr_paths:
        logging.basicConfig(filename = os.path.join(ezr_path, 'genes_info.log'), level = logging.INFO)
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        formatter = logging.Formatter('%(levelname)s - %(message)s')
        ch.setFormatter(formatter)
        logging.getLogger().addHandler(ch)

        num_genes = 0
        num_variants = 0
        max_genes = 0
        min_genes = 100000000
        min_variants = 1000000
        max_variants = 0
        contents = os.listdir(ezr_path)
        ezr_files = filter(lambda f: f.endswith('.ezr'),contents)

        for ezr in ezr_files:
            s = set()
            variants = 0
            with open(os.path.join(ezr_path, ezr)) as file:
                for line in file:
                    if line.startswith('#'):continue
                    variants += 1
                    s.add(line.split('\t')[7].split(';')[0].split('=')[1])
            num_genes += len(s)
            num_variants += variants
            max_genes = max(len(s), max_genes)
            min_genes = min(len(s), min_genes)
            max_variants = max(variants, max_variants)
            min_variants = min(variants, min_variants)

        logging.info('Total Patients: ' + str(len(ezr_files)))
        logging.info('Average # of genes: ' + str(float(num_genes)/len(ezr_files)))
        logging.info('Average # of variants: ' + str(float(num_variants)/len(ezr_files)))
        logging.info('Max # of genes: ' + str(max_genes))
        logging.info('Max # of variants: ' + str(max_variants))
        logging.info('Min # of genes: ' + str(min_genes))
        logging.info('Min # of variants: ' + str(min_variants))

def parse_args(args):
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Annotate an ezr filled directory with the average number of genes and variants per file')
    parser.add_argument('ezr_paths',metavar='DIR',nargs='+',help='ezr directories')
    return parser.parse_args(args)

def main(args = sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
