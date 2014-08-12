#!/usr/bin/env python

"""
The purpose of this script is to generate patients for evaluating
the preformance of exomiser on a disease by disease basis. It will
generate a certain amount of patients for each disease
"""


import os
import sys
import logging

import generate_patient_pairs as gp

from argparse import ArgumentParser
from orpha import Orphanet
from hgmd import HGMD
from omim import MIM

def script(data_path, vcf_path, out_path, num_per, drop_intronic, inheritance=None, **kwargs):
    try:
        hgmd, omim_dict, orph, hp = gp.load_data(data_path)
    except IOError, e:
        logging.error(e)
        sys.exit(1)

    # If we are dropping intronic variants from hgmd, do it now
    if drop_intronic:
        gp.drop_intronic_variants(hgmd)

    # Set up our corrected lookup
    rev_hgmd = hgmd.get_by_omim()
    orph_diseases = orph.filter_lookup(orph.lookup, omim_dict, rev_hgmd, inheritance)
    
    # If vcf dir given, need to check there are at least 2 vcf files
    if vcf_path:
        contents = os.listdir(vcf_path)
        vcf_files = filter(lambda x: x.endswith('.vcf'), contents)
        assert len(vcf_files) > 2, "Need at least 2 vcf files"
    

    for num, dis in orph_diseases.iteritems():
        for i in range(num_per):
            if vcf_path:
                new_patient = gp.copy_vcf(vcf_files, vcf_path, out_path, num, i, 1)[0]
            else:
                new_patient = orphanum + '_' + str(i) + '.vcf'

            if vcf_path:
                gp.infect_geno(os.path.join(out_path, new_patient), dis, rev_hgmd)
            gp.infect_pheno(os.path.join(out_path, new_patient), dis, omim_dict,
                    hp, False, False, 1.0)

def parse_args(args):
    parser = ArgumentParser(description=__doc__.strip())

    parser.add_argument('data_path', metavar='DATA',
            help='Directory from which to grab data files')
    parser.add_argument('--vcf_path', 
            help='Directory from which to take .vcf and .vcf.gz')
    parser.add_argument('out_path', metavar='OUT', 
            help='Directory where to put the generated patient files')
    parser.add_argument('num_per', type=int, metavar='NUM',
            help='Number of patients to generate for each disease')
    parser.add_argument('-I', '--inheritance', nargs='+',
            choices=['AD','AR'],
            help='Which inheritance pattern sampled diseases should have')
    parser.add_argument('--drop_intronic', action='store_true',
            help='Drop intronic variants from HGMD')
    return parser.parse_args(args)

def main(args = sys.argv[1:]):
    args = parse_args(args)
    logging.basicConfig(level=logging.INFO)
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
