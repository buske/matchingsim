#!/usr/bin/env python

#Generate pairs of patients which have inserted variants and associated HPO lists from the same OMIM diseases

__author__ = 'Tal Friedman'

import os
import sys
import shutil

from random import Random
from collections import defaultdict
from argparse import ArgumentParser
from orpha import Orphanet
from hgmd import HGMD
from omim import MIM

#Assume orphanet disease has a phenotype entry
def sample_phenotypes(omim, orph_disease):
    phenotypes = []
    omimd = next(x for x in omim if x.id == orph_disease.pheno[0])
    for pheno, freq in omimd.phenotype_freqs.iteritems():
        if not freq:
            phenotypes.append(pheno)
        else:
            if random.random() < freq:
                phenotypes.append(pheno)
    if phenotypes:
        return phenotypes
    else:
        #we don't want an empty phenotype list
        return sample_phenotypes(omim, orph_disease)


def infect_patient(patient,disease,omim,rev_hgmd):
    phenotypes = sample_phenotypes(omim, disease)


#Choices and weights should be parallel corresponding lists
def weighted_choice(choices, weights):
    total = sum(weights)
    threshold = random.uniform(0,total)
    for k, weight in enumerate(weights):
        total -= weight
        if total < threshold:
            return choices[k]

def load_data(data_path):
    #load hgmd
    hgmd = HGMD(os.path.join(data_path, 'hgmd_correct.jv.vcf'))
  
    #load and filter OMIM
    mim = MIM(os.path.join(data_path, 'phenotype_annotation.tab'))
    omim = filter(lambda d:d.db == 'OMIM', mim.diseases)

    #grab orphanet file names and then load orphanet data
    orphanet_lookup = os.path.join(data_path, 'orphanet_lookup.xml')
    orphanet_inher = os.path.join(data_path, 'orphanet_inher.xml')
    orphanet_geno_pheno = os.path.join(data_path, 'orphanet_geno_pheno.xml')
    orph = Orphanet(orphanet_lookup, orphanet_inher, orphanet_geno_pheno)

    return hgmd, omim, orph

def script(data_path, vcf_path, out_path, N, I=None):
    try:
        hgmd, omim, orph = load_data(data_path)
    except IOError, e:
        print >> sys.stderr, e
        sys.exit(1)
   
    #Set up our corrected lookup
    rev_hgmd = hgmd.get_by_omim()
    lookup = orph.correct_lookup(orph.lookup, omim, rev_hgmd, Inheritance)

    #First, need to check there are at least 2 vcffiles
    contents = os.listdir(vcf_path)
    vcf_files = filter(lambda x: x.endswith('.vcf'), contents)
    assert len(vcf_files) > 2, "Need at least 2 vcf files"
 
    for i in range(N):
        #first, copy pair
        old_pair = Random.sample(vcf_files, 2)
        new_pair = map(lambda x: x[:-4] + '_' + str(i) + '.vcf', old_pair)
        for old, new in zip(old_pair, new_pair):
            shutil.copy(os.path.join(vcf_path, old), os.path.join(out_path, new))         

        #next, get a disease
        #We do a weighted sample based on the number of associated harmful variants
        disease = lookup[weighted_choice(lookup.keys(), [len(rev_hgmd[x.geno[0]]) for x in lookup.itervalues()])]

def parse_args(args):
    parser = ArgumentParser(description='Generate randomly sampled pairs of sick patients')
    praser.add_argument('data_path', metavar='DATA', help='Directory from which to grab data files')
    parser.add_arguemnt('vcf_path', metavar='IN', help='Directory from which to take .vcf and .vcf.gz')
    parser.add_argument('out_path', metavar='OUT', help='Directory where to put the generated patient files')
    parser.add_argument('-N', type=int, help='Number of pairs of patients to generate', required=True)
    parser.add_argument('-I', '--Inheritance',nargs='+',choices=['AD','AR'], help='Which inheritance pattern sampled diseases should have')
    return parser.parse_args(args)

def main(args = sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
