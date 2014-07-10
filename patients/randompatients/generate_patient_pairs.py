#!/usr/bin/env python


"""
Generate random pairs of patients (with genotype and phenotype).
Phenotypes are sampled randomly, using frequency information when 
possible. Genotypes are generated by using control VCFs and adding
harmful variants from HGMD. Bot patients in a pair are infected by
the same genotypic/phenotypic OMIM pair, but do not necessarily
have the same variants added or the exact same phenotypes.
"""


import os
import sys
import shutil
import logging
import random

import hpo

from collections import defaultdict
from argparse import ArgumentParser
from orpha import Orphanet
from hgmd import HGMD
from omim import MIM


__author__ = 'Tal Friedman (talf301@gmail.com)'

def add_imprecision(hp, phenotypes):
    """Randomly add some imprecision by moving terms up the
    tree at random to simulate different language usage

    Args:
        hp: an hpo.HPO instance
        phenotypes: list of phenotypes to add imprecision to

    Returns:
        A list of the phenotypes with imprecision added. 
    Note that the list may be shorter if intersection happens
    """
    # We're going to use a set, since scaling phenotypes up could
    # potentially make them intersect
    new_pheno = set()
    
    for pheno in phenotypes:
        try:
            # Key may not be found since there are things like
            # Inheritance patterns in phenotypic annotations
            ancestors = list(hpo.get_ancestors(hp[pheno]))
        except KeyError:
            continue

        # Get as string ids and remove root
        ancestors = map(lambda x: x.id, ancestors)
        ancestors.remove('HP:0000118')

        # Now randomly add an ancestor
        new_pheno.add(random.choice(ancestors))

    # Return phenotypes as a list
    return list(new_pheno)

def sample_phenotypes(omim_dict, orph_disease, hp, imprecision, default_freq=1.0):
    """Sample phenotypes randomly from an orphanet disease

    Args:
        omim: a dict of OMIM number -> omim.Disease
        orph_disease: an orpha.Disease
        hp: an hpo.HPO instance 
        imprecision: whether or not to add imprecision
        default_freq: default frequency to use if not specified

    Each disease should have at least one phenotype entry

    Returns:
        A list of sampled phenotypes
    """
    phenotypes = []

    # Lookup phenotypic OMIM for orphanet disease
    omim_id = orph_disease.pheno[0]

    # Grab the Disease entry for this OMIM
    try:
        omim_dis = omim_dict[omim_id]
    except KeyError:
        logging.warning('Could not find OMIM entry for %s' % omim_id)

    # If frequency available, we will sample, otherwise use default freq 
    phenotype_freqs = omim_dis.phenotype_freqs
    assert phenotype_freqs, "Missing phenotypes for: %s" % omim_id
    for pheno, freq in phenotype_freqs.iteritems():
        if not freq and random.random() < default_freq:
            phenotypes.append(pheno)
        else:
            if random.random() < freq:
                phenotypes.append(pheno)
    if phenotypes:
        # Add imprecision if necessary
        if imprecision:
           phenotypes = add_imprecision(hp, phenotypes) 
        return phenotypes
    else:
        logging.warning("Random phenotype sampling for %s resulted in"
                " empty set" % omim_id)
        return sample_phenotypes(omim_dict, orph_disease, default_freq)

def sample_variants(rev_hgmd, orph_disease):
    """Sample variants randomly from an orphanet disease

    Args:
        rev_hgmd: A dict of OMIM number -> list(hgmd.Entry)
        orph_disease: An orpha.Disease instance to be sampled off of

    Returns:
        A list of (variant, homozygous?) tuples
    """
    variants = []

    # Grab inheritance pattern and list of associated variants
    inheritance = orph_disease.inheritance[0]
    hgmd_variants = rev_hgmd[orph_disease.geno[0]]

    if inheritance == 'Autosomal recessive':
        if len(hgmd_variants) == 1 or random.random() < 0.5:
            # One hom variant
            variants.append((random.choice(hgmd_variants), True))
        else:
            # Two het variants
            vars = random.sample(hgmd_variants, 2)
            for var in vars:
                variants.append((var, False))
    elif inheritance == 'Autosomal dominant':
        # If AD then one het variant added
        variants.append((random.choice(hgmd_variants), False))
    else:
        raise NotImplementedError('Unexpected inheritance: %s' % inheritance)

    return variants

def generate_vcf_line(var, hom=False):
    """Return newline-terminated VCF line for given variant"""
    gt = '1/1' if hom else '0/1'
    return '%s\n' % '\t'.join([var.chrom, var.loc, '.', var.ref, 
        var.alt, '255', 'PASS', var.info_line, 'GT', gt])

def infect_pheno(patient, orph_disease, omim_dict, hp, imprecision, default_freq):
    """Do phenotypic infection by producing a file with a
    list of phenotypes associated with disease

    Args:
        patient: a string path to the patient vcf
        orph_disease: an orpha.Disease instant to infect patient
        omim_dict: a dict of OMIM number -> omim.Disease
        hp: an hpo.HPO instance
        imprecision: whether or not to add imprecision to pheno sampling
        default_freq: default frequency to use if info not found
    """
    # Sample phenotypes
    phenotypes = sample_phenotypes(omim_dict, orph_disease, hp,
            imprecision, default_freq)

    assert patient.endswith('.vcf')

    # If patient is HG01.vcf, phenotypes in HG01_hpo.txt
    pheno_file = patient[:-4] + '_hpo.txt'
    with open(pheno_file, 'w') as hpo:
        hpo.write(','.join(phenotypes))

def infect_geno(patient, orph_disease, rev_hgmd):
    """Do genotypic infection by inserting harmful variants 
    associated with disease.

    Args:
        patient: a string path to the patient vcf
        orph_disease: an orpha.Disease instance to infect patient with
        rev_hgmd: a dict of OMIM number -> list(hgmd.Entry)
    """
    # Sample variants
    variants = sample_variants(rev_hgmd, orph_disease)

    assert patient.endswith('.vcf')

    with open(patient, 'a') as vcf:
        for variant, hom in variants:
            vcf.write(generate_vcf_line(variant, hom=hom))

def weighted_choice(choices, weights):
    """Return a random choice, given corresponding weights
    
    Args:
        choices: a list of options to choose from
        weights: a list of numeric weights, one per choice
        choices and weights must be parallel corresponding lists

    Returns:
        a random element from choices
    """
    assert len(choices) == len(weights)
    total = sum(weights)
    threshold = random.uniform(0,total)
    for k, weight in enumerate(weights):
        total -= weight
        if total < threshold:
            return choices[k]

def load_data(data_path):
    """Load all the required data files the program needs
    
    Args:
        data_path: String file path to the directory files are in

    Returns
        (HGMD, OMIM number -> omim.Disease, Orphanet, HPO)
    """
    # Load hgmd
    hgmd = HGMD(os.path.join(data_path, 'hgmd_correct.jv.vcf'))
  
    # Load and filter OMIM into an OMIM number: Disease dict
    mim = MIM(os.path.join(data_path, 'phenotype_annotation.tab'))
    omim = filter(lambda d:d.db == 'OMIM', mim.diseases)
    omim_dict = {dis.id:dis for dis in omim}

    # Grab orphanet file names and then load orphanet data
    orphanet_lookup = os.path.join(data_path, 'orphanet_lookup.xml')
    orphanet_inher = os.path.join(data_path, 'orphanet_inher.xml')
    orphanet_geno_pheno = os.path.join(data_path, 'orphanet_geno_pheno.xml')
    orph = Orphanet(orphanet_lookup, orphanet_inher, orphanet_geno_pheno)

    # Load hpo
    hp = hpo.HPO(os.path.join(data_path, 'hp.obo'))
    # Filter to phenotypic abnormailities
    hp.filter_to_descendants('HP:0000118')

    return hgmd, omim_dict, orph, hp

def copy_vcf(vcf_files, vcf_path, out_path, orphanum, i, num_vcf):
    """Copy in a new vcf group, sampling at random from vcf_files

    Args:
        vcf_files: list of vcf files in source directory
        vcf_path: path to vcf_files
        out_path: path where new files should be put
        orphanum: orphanet disease number do sign new file with
        i: iteration to sign new file with 
        num_vcf: number of vcf's to sample and copy
    Returns:
        A list of the new patient locations
    """

    old_pair = random.sample(vcf_files, num_vcf)
    new_pair = map(lambda x: x[:-4] + '_' + orphanum + '_' + str(i) + '.vcf', old_pair)
    for old, new in zip(old_pair, new_pair):
        shutil.copy(os.path.join(vcf_path, old), os.path.join(out_path, new))         
    return new_pair

def drop_intronic_variants(hgmd):
    """Drop all intronic variants in the given hgmd instance
    
    Args:
        hgmd: an hgmd instance
    """
    dropped_counter = 0
    # We make a list of accepted effects
    accepted = ['FS_DELETION', 'FS_SUBSTITUTION', 'NON_FS_DELETION', 
            'NON_FS_SUBSTITUTION', 'NONSYNONYMOUS', 'SPLICING', 
            'STOPGAIN', 'STOPLOSS']

    # Also make a list of rejected effects, so we know when something is totally wrong
    rejected = ['ERROR', 'INTERGENIC', 'INTRONIC', 'ncRNA_EXONIC', 
            'ncRNA_SPLICING', 'SYNONYMOUS', 'UTR3', 'UTR5']

    new_entries = []
    for entry in hgmd.entries:
        if entry.effect in accepted:
            new_entries.append(entry)
        elif entry.effect in rejected:
            dropped_counter += 1
        else:
            logging.error("Entry did not have matched effect %s\n" % entry)
    logging.warning("%d HGMD variants dropped as intronic\n" % dropped_counter)

    #Finally, reassign
    hgmd.entries = new_entries

def script(data_path, vcf_path, out_path, generate, num_samples, by_variant, default_freq, 
        drop_intronic, imprecision, inheritance=None, **kwargs):
    try:
        hgmd, omim_dict, orph, hp = load_data(data_path)
    except IOError, e:
        logging.error(e)
        sys.exit(1)
  
    # If we are dropping intronic variants from hgmd, do it now
    if drop_intronic:
        drop_intronic_variants(hgmd)

    # Set up our corrected lookup
    rev_hgmd = hgmd.get_by_omim()
    orph_diseases = orph.filter_lookup(orph.lookup, omim_dict, rev_hgmd, inheritance)

    # If vcf dir given, need to check there are at least 2 vcf files
    if vcf_path:
        contents = os.listdir(vcf_path)
        vcf_files = filter(lambda x: x.endswith('.vcf'), contents)
        assert len(vcf_files) > 2, "Need at least 2 vcf files"
    
    # Dealing with pairs
    if generate == 'PAIRS':
        for i in range(num_samples):
            # First, get a disease
            if by_variant:
                # We do a weighted sample based on the number of associated harmful variants
                orphanum = weighted_choice(orph_diseases.keys(), 
                    [len(rev_hgmd[x.geno[0]]) for x in orph_diseases.itervalues()])
                disease = orph_diseases[orphanum]
            else:
                # A uniform sample over all diseases
                orphanum = random.choice(orph_diseases.keys())
                disease = orph_diseases[orphanum]

            # Next, if we have a vcf dir copy pair
            if vcf_path:
                new_pair = copy_vcf(vcf_files, vcf_path, out_path, orphanum, i, 2)
            else:
                # Otherwise, name patients based on just disease and iteration (with fake vcf)
                new_pair = ['First_' + orphanum + '_' + str(i) + '.vcf',
                        'Second_' + orphanum + '_' + str(i) + '.vcf']

            # Finally, infect both patients with disease
            for patient in new_pair:
                if vcf_path:
                    infect_geno(os.path.join(out_path, patient), disease, rev_hgmd)
                infect_pheno(os.path.join(out_path, patient), disease, omim_dict, 
                        hp, imprecision,  default_freq)

    # Dealing with individual patients
    if generate == 'PATIENTS':
        for i in range(num_samples):
            # First, get a disease
            if by_variant:
                # We do a weighted sample based on the number of associated harmful variants
                orphanum = weighted_choice(orph_diseases.keys(), 
                    [len(rev_hgmd[x.geno[0]]) for x in orph_diseases.itervalues()])
                disease = orph_diseases[orphanum]
            else:
                # A uniform sample over all diseases
                orphanum = random.choice(orph_diseases.keys())
                disease = orph_diseases[orphanum]
            
            # Next, if we have a vcf dir copy one over 
            if vcf_path:
                new_patient = copy_vcf(vcf_files, vcf_path, out_path, orphanum, i, 1)[0]
            else:
                # Otherwise, name patient based on just disease and iteration (with fake vcf)
                new_patient = orphanum + '_' + str(i) + '.vcf'

            # Finally, infect patient with geno and pheno
            if vcf_path:
                infect_geno(os.path.join(out_path, new_patient), disease, rev_hgmd)
            infect_pheno(os.path.join(out_path, new_patient), disease, omim_dict, 
                    hp, imprecision,default_freq)

def parse_args(args):
    parser = ArgumentParser(description=__doc__.strip())

    parser.add_argument('data_path', metavar='DATA', 
            help='Directory from which to grab data files')
    parser.add_argument('--vcf_path', 
            help='Directory from which to take .vcf and .vcf.gz')
    parser.add_argument('out_path', metavar='OUT', 
            help='Directory where to put the generated patient files')
    parser.add_argument('--generate', dest='generate',
            choices=['PATIENTS', 'PAIRS'], default='PAIRS',
            help='Generate pairs or individuals (defult is pairs)')
    parser.add_argument('-N', type=int, dest='num_samples',
            help='Number of samples to generate (patients or pairs)')
    parser.add_argument('-I', '--inheritance',nargs='+',
            choices=['AD','AR'], 
            help='Which inheritance pattern sampled diseases should have')
    parser.add_argument('-D', '--default_freq', type=float,
            help='Default frequency for phenotype if info not found (default'
            'is 1.0)', default=1.0)
    parser.add_argument('--drop_intronic', action='store_true',
            help='Drop intronic variants from HGMD')
    parser.add_argument('--imprecision', action='store_true',
            help='Add imprecision of hpo term selection (randomly send'
            'terms to ancestors')
    parser.add_argument('-V', dest='by_variant', action='store_true',
            help='Sample diseases weighted by variant, default is uniform')
    parser.add_argument('--logging', default='WARNING',
            choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
            help='Logging level')
    return parser.parse_args(args)

def main(args = sys.argv[1:]):
    args = parse_args(args)
    logging.basicConfig(level=args.logging)
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
