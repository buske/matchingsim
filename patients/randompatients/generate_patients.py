#!/usr/bin/env python

# OMIM, Orphanet  parsing code taken/adapted from Orion Buske, original can be found at
# https://github.com/buske/udp-dating/blob/master/

__author__ = 'Tal Friedman'

import os
import sys
import random
import gzip

from collections import defaultdict
from orpha import Orphanet
from hgmd import HGMD
from omim import MIM

#Assume orphanet disease has a phenotype entry
def sample_phenotypes(omim, orph_disease):
    phenotypes = []
    omimd = next(x for x in omim if x.id == orph_disease[0][0])
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

def annotate_patient(patient,rev_hgmd,omim,lookup):
    try:
        if not patient[-4:] == '.vcf' and not patient[-7:] == '.vcf.gz':
            print >> sys.stderr, "Incorrect file format, use .vcf.gz or .vcf"
            return
    except IndexError:
            print >> sys.stderr, "Incorrect file format, use .vcf.gz or .vcf"
    
    #find if .vcf or .vcf.gz
    if patient.endswith('.vcf'):
        file = open(patient, 'a')
        name = patient[:-4]
    elif patient.endswith('.vcf.gz'):
        file = gzip.open(patient, 'ab')
        name = patient[:-7]        

    orph_disease = lookup[random.choice(lookup.keys())]
    phenotypes = sample_phenotypes(omim, orph_disease)
    
    #if autosomal recessive, if we only have one variant available use it as homozygous otherwise randomly (50/50) pick two and use as heterozygous or pick one as homozygous
    if orph_disease[1][0] == 'Autosomal recessive':
        if len(rev_hgmd[orph_disease[2][0]]) == 1 | random.random() < 0.5:
            var = rev_hgmd[orph_disease[2][0]][0]
            file.write('\t'.join([var.chrom,var.loc,'.',var.ref,var.alt,'100','PASS','.','GT','1|1'])+'\n')
        else:
            vars = random.sample(rev_hgmd[orph_disease[2][0]], 2)
            for var in vars:
                file.write('\t'.join([var.chrom,var.loc,'.',var.ref,var.alt,'100','PASS','.','GT','0|1'])+'\n')
                    
    #If AD (or something else, but those should never happen) then one heterozygous mutation added
    else orph_disease[1][0] == 'Autosomal dominant':
        var = random.choice(rev_hgmd[orph_disease[2][0]])
        file.write('\t'.join([var.chrom,var.loc,'.',var.ref,var.alt,'100','PASS','.','GT','0|1'])+'\n')
    

    hpo = open(name + '_hpo.txt','w')
    hpo.write(phenotypes[0])
    for p in phenotypes[1:]:
        hpo.write(','+p) 
    file.close() 

def annotate_patient_dir(pdir,rev_hgmd,omim,lookup):
    for f in os.listdir(pdir):
        if os.path.isfile(os.path.join(pdir,f)) and (f.endswith('.vcf') or f.endswith('.vcf.gz')):
            annotate_patient(os.path.join(pdir,f),rev_hgmd,omim,lookup) 

def has_pattern(patterns, o):
    return any(x in patterns for x in o[1])

def has_pheno(omim, o):
    return any(x.id == o[0][0] for x in omim)

#ensure that all elements of lookup are entirely useable
def correct_lookup(lookup, omim,rev_hgmd, Inheritance=None):
    #get ideal orphanet cases
    newlook = {}
    for k,v in lookup.iteritems():
        if len(v[0]) == 1 and len(v[1]) == 1 and len(v[2]) == 1:
            newlook[k] = v
    #get the right disease set based on inheritance
    if Inheritance:
        patterns = []
        if 'AD' in Inheritance:
             patterns.append('Autosomal dominant')
        if 'AR' in Inheritance:
            patterns.append('Autosomal recessive')
        newlook = {k:v for k,v in newlook.iteritems() if has_pattern(patterns, v)}
    
    #ensure all orphanet cases have phenotypic annotations
    lookup = {k:v for k,v in newlook.iteritems() if has_pheno(omim, v)}
    #ensure all orphanet cases have at least one associated variant
    newlook = {}
    for k, o in lookup.iteritems():
        try:
            a = rev_hgmd[o[2][0]]
            newlook[k] = o
        except KeyError:
            pass
    return newlook

def script(pheno_file, hgmd_file, patient_path, orphanet_lookup, orphanet_inher, orphanet_geno_pheno,  Inheritance=None):
    try:
        mim = MIM(pheno_file)
    except IOError:
        print >> sys.stderr, "OMIM file not found or invalid"
 
    omim = filter(lambda d:d.db == 'OMIM',mim.diseases)

    try:    
        hgmd = HGMD(hgmd_file)
    except IOError:
        print >> sys.stderr, "HGMD file not found or invalid"
   
    try:
        orph = Orphanet(orphanet_lookup,orphanet_inher, orphanet_geno_pheno)
    except IOError:
        print >> sys.stderr, "Orphanet files not found or invalid"

    #get hgmd variants by omim
    rev_hgmd = hgmd.get_by_omim()

    lookup = correct_lookup(orph.lookup,omim,rev_hgmd,Inheritance)
    #if we are given a directory, annotate each vcf.gz or vcf file in the directory assuming it is a patient 
    if os.path.isdir(patient_path):
        print "Given a directory full of patients!"
        annotate_patient_dir(patient_path, rev_hgmd, omim, lookup)
    #if we are given a single file just annotate it normally
    elif os.path.isfile(patient_path):
        print "Given a single patient file!"
        annotate_patient(patient_path, rev_hgmd, omim, lookup)
    else:
        print >> sys.stderr, "Patient file/folder not found or invalid"   
  
def parse_args(args):
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Generate randomly sampled sick patients')
    parser.add_argument('pheno_file', metavar='PHENO', help='phenotype annotation tab file')
    parser.add_argument('hgmd_file',metavar='HGMD', help='Annotated HGMD file in .vcf format')
    parser.add_argument('patient_path',metavar='PATH', help='Path to a .vcf, .vcf.gz or a directory with multiple of these in it')
    parser.add_argument('orphanet_lookup',metavar='ORPHLOOK', help='Orphanet XML file which crossreferences to OMIM')
    parser.add_argument('orphanet_inher',metavar='ORPHINHER', help='Orphanet XML file giving inheritance patterns')
    parser.add_argument('orphanet_geno_pheno',metavar='ORPHGENOPHENO', help = 'Orphanet XML file relating genotypic inheritance to phenotypic outcome')
    parser.add_argument('-I', '--Inheritance',nargs='+', choices=['AD','AR'], help='Which inheritance pattern sampled diseases should have, default is any, including unknown')
    return parser.parse_args(args)

def main(args = sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
