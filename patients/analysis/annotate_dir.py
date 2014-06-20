#!/usr/bin/env python

#Given a directory, will add a corresponding text file for each patient containing info on Genotypic OMIM, Phenotypic OMIM, Variant, and rank as given by exomizer 

__author__ = 'Tal Friedman'

import os
import sys

from hgmd import HGMD
from argparse import ArgumentParser

def script(path):
    hgmd = HGMD('/dupa-filer/talf/matchingsim/patients/hgmd_correct.jv.vcf')

def main(args = sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name == '__main__':
    sys.exit(main())
