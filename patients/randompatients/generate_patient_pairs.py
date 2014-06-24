#!/usr/bin/env python

#Generate pairs of patients which have inserted variants and associated HPO lists from the same OMIM diseases

__author__ = 'Tal Friedman'

import os
import sys

from collections import defaultdict
from argparse import ArgumentParser


def parse_args(args):
    parser = ArgumentParser(description='Generate randomly sampled pairs of sick patients')
    parser.add_argument('pheno_file', metavar='PHENO', help='phenotype annotation .tab file')
    parser.add_argument('hgmd_file', metavar='HGMD', help='Annotated HGMD file in .vcf format')
    parser.add_arguemnt('vcf_path', metavar='IN', help='Directory from which to take .vcf and .vcf.gz')
    parser.add_argument('orphanet_lookup',metavar='ORPHLOOK', help='Orphanet XML file which crossreferences to genotypic OMIMs')
    parser.add_argument('orphanet_inher',metavar='ORPHINHER', help='Orphanet XML file giving inheritance patterns.')
    parser.add_argument('orphanet_geno_pheno',metavar='ORPHGENOPHENO', help='Orphanet XML file relation genotypes to phenotypes')
    parser.add_argument('-I', '--Inheritance',nargs='+',choices=['AD','AR'], help='Which inheritance pattern sampled diseases should have')
    return parser.parse_args(args)

def main(args = sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
