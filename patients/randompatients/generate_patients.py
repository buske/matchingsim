#!/usr/bin/env python


"""
Generate random patients (phenotype and genotype) with OMIM disorders.
Phenotypes are sampled randomly, using frequency information where available.
Genotypes are generated using control VCF files, with harmful variants added
from HGMD.
"""


# OMIM, Orphanet  parsing code taken/adapted from Orion Buske, original can be found at
# https://github.com/buske/udp-dating/blob/master/


import os
import sys
import random
import gzip
import logging

from collections import defaultdict

from orpha import Orphanet
from hgmd import HGMD
from omim import MIM


__author__ = 'Tal Friedman'

def weighted_choice(choices, weights):
    """Return a random choice, given corresponding weights.

    Args:
      choices: a list of options to choose from
      weights: a list of numeric weights, one per choice
      choices and weights must be parallel corresponding lists

    Returns:
      a random element from choices
    """
    assert len(choices) == len(weights)
    total = sum(weights)
    threshold = random.uniform(0, total)
    for k, weight in enumerate(weights):
        threshold -= weight
        if threshold <= 0:
            return choices[k]

def sample_phenotypes(omim, orph_disease):
    """Sample phenotypes randomly for an orphanet disease.

    Args:
      omim: a MIM object
      orph_disease: the orphanet disease to sample (type?)

    Each orphanet disease must have at least one phenotype entry
    """
    phenotypes = []

    # Lookup OMIM entry for orphanet disease
    # Is the OMIM ID always the first?
    omim_id = orph_disease.pheno[0]

    # Ideally, omim.py would be fixed so that the MIM object
    # can be indexed with an ID.
    omimd = next(x for x in omim if x.id == omim_id)
    assert omimd, "Could not find OMIM entry for: %s" % omim_id
    
    # Sample phenotypes randomly (weighted by frequency info if present)
    phenotype_freqs = omimd.phenotype_freqs
    assert phenotype_freqs, "Missing phenotypes for: %s" % omim_id
    for pheno, freq in phenotype_freqs.iteritems():
        if not freq:
            phenotypes.append(pheno)
        else:
            if random.random() <= freq:
                phenotypes.append(pheno)

    if phenotypes:
        return phenotypes
    else:
        logging.warning("Random phenotype sampling for %s resulted in"
            " empty set" % omim_id)
        return sample_phenotypes(omim, orph_disease)

def sample_variants(hgmd_variants, inheritance):
    """Return list of (variant, homozygous?) tuples, sampled according to inheritance."""
    # What about when multiple inheritance types are specified?
    variants = []
    if inheritance == 'Autosomal recessive':
        if len(hgmd_variants) == 1 or random.random() < 0.5:
            # One hom variant
            variants.append((hgmd_variants[0], True))
        else:
            # Two het variants
            vars = random.sample(hgmd_variants, 2)
            for var in vars:
                variants.append((var, False))
    elif inheritance == 'Autosomal dominant':
        # If AD then one heterozygous mutation added
        variants.append((random.choice(hgmd_variants), False))
    else:
        raise NotImplementedError('Unexpected inheritance: %s' % inheritance)

    return variants
 
def generate_vcf_line(var, hom=False):
    """Return newline-terminated VCF line for given variant"""
    gt = '1/1' if hom else '0/1'
    return '%s\n' % '\t'.join([var.chrom, var.loc, '.', var.ref, 
        var.alt, '255', 'PASS', '.', 'GT', gt])

def annotate_patient(patient, rev_hgmd, omim, lookup, by_variant, produce_omim=False):
    """Generate a random disease patient.

    Modify the patient VCF file with a variant causal of a random disease,
    and output sampled phenotypes to a file.

    Args:
      patient:
      rev_hgmd:
      omim:
      lookup:
      produce_omim:
    """
    assert patient.endswith('.vcf') or patient.endswith('.vcf.gz'), \
        "Incorrect file format, use .vcf.gz or .vcf"
    
    # Open file
    if patient.endswith('.vcf'):
        file = open(patient, 'a')
        name = patient[:-4]
    else:  # .vcf.gz
        file = gzip.open(patient, 'ab')
        name = patient[:-7]


    # Sample random variant and corresponding OMIM disease
    if by_variant:
        diseases = list(lookup)
        variants_per_disease = [len(rev_hgmd[lookup[disease].geno[0]]) for disease in diseases]
        orph_disease = lookup[weighted_choice(diseases, variants_per_disease)]
    else:
        diseases = list(lookup)
        orph_disease = lookup[random.choice(diseases)]
    phenotypes = sample_phenotypes(omim, orph_disease)
    
    # If autosomal recessive, if we only have one variant available use it as 
    # homozygous otherwise randomly (50/50) pick two and use as heterozygous 
    # or pick one as homozygous
    inheritance = orph_disease.inheritance[0]
    hgmd_variants = rev_hgmd[orph_disease.geno[0]]

    variants = sample_variants(hgmd_variants, inheritance)

    # Save all the writing logic for one place:
    for variant, hom in variants:
        file.write(generate_vcf_line(variant, hom=hom))

    if produce_omim:
        with open(name + '_omim.txt', 'w') as hpo:
            hpo.write(orph_disease.pheno[0])
    else:
        with open(name + '_hpo.txt', 'w') as hpo:
            hpo.write(','.join(phenotypes))

    file.close()

def annotate_patient_dir(patient_dir, rev_hgmd, omim, lookup, produce_omim, by_variant):
    """Annotate all patient VCF files in directory patient_dir."""
    for filename in os.listdir(patient_dir):
        filepath = os.path.join(patient_dir, filename)
        if os.path.isfile(filepath) and (filepath.endswith('.vcf') or filepath.endswith('.vcf.gz')):
            annotate_patient(filepath, rev_hgmd, omim, lookup, produce_omim,by_variant) 

def has_pattern(patterns, o):
    return any(x in patterns for x in o.inheritance)

def has_pheno(omim, o):
    return any(x.id == o.pheno[0] for x in omim)

def filter_lookup(lookup, omim, rev_hgmd, inheritance=None):
    """Return a new, filtered lookup based on inheritance and mapping ability.

    Args:
      lookup: dict ??? -> ???
      omim: ?
      rev_hgmd: ?
      inheritance: ?

    Returns:
      a new lookup dict ??? -> ???
    """
    # Get diseases with straightforward orphanet mappings
    # Better variable name?
    newlook = {k:v for k,v in lookup.iteritems() 
               if len(v.pheno) == 1 and len(v.inheritance) == 1 and len(v.geno) == 1}
    
    # Filter diseases to any specified inheritance models
    if inheritance:
        patterns = []
        if 'AD' in inheritance:
             patterns.append('Autosomal dominant')
        if 'AR' in inheritance:
            patterns.append('Autosomal recessive')

        newlook = {k:v for k,v in newlook.iteritems() if has_pattern(patterns, v)}
    
    # Ensure all orphanet cases have phenotypic annotations
    lookup = {k:v for k,v in newlook.iteritems() if has_pheno(omim, v)}

    # Ensure all orphanet cases have at least one associated variant
    newlook = {}
    for k, o in lookup.iteritems():
        try:
            a = rev_hgmd[o.geno[0]]
            newlook[k] = o
        except KeyError:
            pass

    return newlook

def script(pheno_file, hgmd_file, patient_path, orphanet_lookup, 
           orphanet_inher, orphanet_geno_pheno,  
           produce_omim=False, by_variant=False, inheritance=None, **kwargs):
    # I wouldn't actually bother with these try/excepts, since they make the 
    # code a lot harder to read, and don't make things that much clearer for
    # the user in the case of an error. This is debatable though.
    try:
        mim = MIM(pheno_file)
    except IOError:
        logging.error("OMIM file not found or invalid")
        raise
 
    try:    
        hgmd = HGMD(hgmd_file)
    except IOError:
        logging.error("HGMD file not found or invalid")
        raise
   
    try:
        orph = Orphanet(orphanet_lookup, orphanet_inher, orphanet_geno_pheno)
    except IOError:
        logging.error("Orphanet files not found or invalid")
        raise


    omim = filter(lambda d: d.db == 'OMIM', mim.diseases)

    # Get hgmd variants by omim
    rev_hgmd = hgmd.get_by_omim()

    logging.debug(orph.lookup)
    lookup = filter_lookup(orph.lookup, omim, rev_hgmd, inheritance)
    
    if os.path.isdir(patient_path):
        # If we are given a directory, annotate each vcf.gz or vcf file in the 
        # directory assuming it is a patient
        logging.info("Processing directory of patient VCF files...")
        annotate_patient_dir(patient_path, rev_hgmd, omim, lookup, produce_omim, by_variant)
    elif os.path.isfile(patient_path):
        # If we are given a single file just annotate it normally
        logging.info("Processing single patient VCF file...")
        annotate_patient(patient_path, rev_hgmd, omim, lookup, produce_omim, by_variant)
    else:
        logging.error("Patient file/folder not found or invalid")

def parse_args(args):
    from argparse import ArgumentParser
    parser = ArgumentParser(description=__doc__.strip())

    parser.add_argument('pheno_file', metavar='PHENO', 
                help='OMIM->HPO phenotype annotation tab file')
    parser.add_argument('hgmd_file', metavar='HGMD', 
                help='Annotated HGMD file in .vcf format')
    parser.add_argument('patient_path', metavar='PATH', 
                help='Path to a .vcf, .vcf.gz or a directory with'
                ' multiple of these in it')
    parser.add_argument('orphanet_lookup', metavar='ORPH_LOOK', 
                help='Orphanet XML file which crossreferences to OMIM')
    parser.add_argument('orphanet_inher', metavar='ORPH_INHER', 
                help='Orphanet XML file giving inheritance patterns')
    parser.add_argument('orphanet_geno_pheno', metavar='ORPH_GENO_PHENO', 
                help='Orphanet XML file relating genotypic inheritance'
                ' to phenotypic outcome')
    
    parser.add_argument('-V', dest='by_variant', action='store_true',
                help='Sample diseases weighted by variant, default is uniform') 
    parser.add_argument('-I', '--inheritance', nargs='+', choices=['AD','AR'], 
                help='Which inheritance pattern sampled diseases'
                ' should have, default is any, including unknown.'
                ' May be specified multiple times.')
    parser.add_argument('-O', dest="produce_omim", default=False, 
                action='store_true', help='Produce OMIM IDs rather'
                ' than a list of phenotypes')
    parser.add_argument('--logging', default='WARNING',
                help='Logging level (e.g. DEBUG, ERROR)')

    return parser.parse_args(args)

def main(args = sys.argv[1:]):
    args = parse_args(args)
    logging.basicConfig(level=args.logging)
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
