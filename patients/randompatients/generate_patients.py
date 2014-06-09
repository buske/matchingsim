#!/usr/bin/env python

# OMIM, Orphanet  parsing code taken/adapted from Orion Buske, original can be found at
# https://github.com/buske/udp-dating/blob/master/

__author__ = 'Tal Friedman'

import os
import sys
import re
import logging
import random
import gzip

from collections import defaultdict
import xml.etree.ElementTree as ET
from orpha import Orphanet
from hgmd import HGMD

FREQUENCIES = {'very rare':  0.01, 
               'rare':       0.05, 
               'occasional': 0.075, 
               'frequent':   0.33, 
               'typical':    0.5, 
               'variable':   0.5, 
               'common':     0.75, 
               'hallmark':   0.9, 
               'obligate':   1.0}
fraction_frequency_re = re.compile(r'of|/')

class Disease:
    def __init__(self, db, id, name, phenotype_freqs):
        self.db = db
        self.id = id
        self.name = name
        self.phenotype_freqs = phenotype_freqs

class MIM:
    def __init__(self, filename):
        self.diseases = list(self.iter_diseases(filename))

    def __iter__(self):
        return iter(self.diseases)

    @classmethod
    def iter_disease_lines(cls, filename):
        with open(filename) as ifp:
            cur_disease = None
            cur_lines = []
            for line in ifp:
                line = line.rstrip()
                tokens = line.split('\t')
                if len(tokens) == 1: continue
                disease = (tokens[0].strip(), tokens[1].strip())

                if disease == cur_disease:
                    cur_lines.append(tokens)
                else:
                    if cur_disease:
                        yield cur_disease, cur_lines

                    cur_lines = [tokens]
                    cur_disease = disease
            if cur_disease:
                yield cur_disease, cur_lines

    @classmethod
    def parse_frequency(cls, s, default=None):
        """Return float parsed frequency or default if problem or absent"""
        s = s.lower()
        if not s:
            freq = default
        elif s in FREQUENCIES:
            freq = FREQUENCIES[s]
        elif s.endswith('%'):
            s = s.replace('%', '')
            if '-' in s:
                # Average any frequency ranges
                low, high = s.split('-')
                freq = (float(low) + float(high)) / 2 / 100
            else:
                freq = float(s) / 100
        else:
            try:
                num, denom = fraction_frequency_re.split(s)
            except:
                logging.error("Error parsing frequency: {!r}".format(s))
                freq = default
            else:
                freq = float(num) / float(denom)

        return freq

    @classmethod
    def iter_diseases(cls, filename, default_freq=None):
        for disease, tokens_list in cls.iter_disease_lines(filename):
            db, id = disease
            raw_phenotypes = defaultdict(list)
            name = None
            for tokens in tokens_list:
                freq = cls.parse_frequency(tokens[8])
                hp_term = tokens[4].strip()
                raw_phenotypes[hp_term].append(freq)
                if not name:
                    name = tokens[2].strip()

            phenotype_freqs = {}
            for hp_term, freqs in raw_phenotypes.items():
                non_null = [x for x in freqs if x is not None]
                if non_null:
                    freq = sum(non_null) / len(non_null)
                else:
                    freq = default_freq

                phenotype_freqs[hp_term] = freq
                
            disease = Disease(db, id, name, phenotype_freqs)
            yield disease

def annotate_patient(patient,hgmd,omim):
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

    dis, phenotypes = hgmd.sample_disease(omim)
    file.write('\t'.join([dis.chrom,dis.loc,'.',dis.ref,dis.alt,'100','PASS','.','GT','1|1'])+'\n')
    hpo = open(name + '_hpo.txt','w')
    hpo.write(phenotypes[0])
    for p in phenotypes[1:]:
        hpo.write(','+p) 
    file.close() 

def annotate_patient_dir(pdir,hgmd,omim):
    for f in os.listdir(pdir):
        if os.path.isfile(os.path.join(pdir,f)) and (f.endswith('.vcf') or f.endswith('.vcf.gz')):
            annotate_patient(os.path.join(pdir,f),hgmd,omim) 

def has_pattern(patterns, orphanet, e):
    try:
        return any(x in patterns for x in orphanet.inheritance[e.omimid])
    except KeyError:
        return false


def script(pheno_file, hgmd_file, patient_path, orphanet_lookup, orphanet_inher, Inheritance=None):
    try:
        mim = MIM(pheno_file)
    except FileNotFoundError:
        print >> sys.stderr, "OMIM file not found or invalid"
 
    omim = filter(lambda d:d.db == 'OMIM',mim.diseases)

    try:    
        hgmd = HGMD(hgmd_file)
    except FileNotFoundError:
        print >> sys.stderr, "HGMD file not found or invalid"
   
    try:
        orph = Orphanet(orphanet_lookup,orphanet_inher)
    except FileNotFoundError:
        print >> sys.stderr, "Orphanet files not found or invalid"

    hgmd.entries = filter(lambda d:d.omimid,hgmd.entries)         
     
    #get the right disease set based on inheritance
    if inheritance:
        patterns = []
        if 'AD' in inheritance:
             patterns.append('Autosomal dominant')
        if 'AR' in inheritance:
            patterns.append('Autosomal recessive')
        hgmd = [x for x in hgmd if has_pattern(patterns, orph, x)]                    
  
  #if we are given a directory, annotate each vcf.gz or vcf file in the directory assuming it is a patient 
    if os.path.isdir(patient_path):
        print "Given a directory full of patients!"
        annotate_patient_dir(patient_path, hgmd, omim)
    #if we are given a single file just annotate it normally
    elif os.path.isfile(patient_path):
        print "Given a single patient file!"
        annotate_patient(patient_path, hgmd, omim)
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
    parser.add_argument('-I', '--Inheritance',nargs='+', choices=['AD','AR'], help='Which inheritance pattern sampled diseases should have, default is any, including unknown')
    return parser.parse_args(args)

def main(args = sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    print "begun"
    #sys.exit(main())
