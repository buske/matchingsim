#!/usr/bin/env python

# OMIM parsing code taken from Orion Buske, original can be found at
# https://github.com/buske/udp-dating/blob/master/mim.py

__author__ = 'Tal Friedman'

import os
import sys
import re
import logging
import random
import gzip

from collections import defaultdict

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
        for i in range(5):
            logging.debug(self.diseases[i].__dict__)

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

class Entry:
    def __init__(self, chrom, loc, ref, alt, effect, pmid, omimid):
        self.chrom = chrom
        self.loc = loc
        self.ref = ref
        self.alt = alt
        self.effect = effect
        self.omimid = omimid
        self.pmid = pmid
    
    def __str__(self):
        return [self.chrom, self.loc, self.ref, self.alt, self.effect, self.omimid, self.pmid].__str__()

    def __repr__(self):
        return self.__str__()
    
    def get_phenotypes(self, omim):
        if self.omimid: 
            om = next(x for x in omim if x.id == self.omimid)
        return list(om.phenotype_freqs)
     
class HGMD:
    def __init__(self, filename):
        self.entries = list(self.iter_lines(filename))
    
    def __iter__(self):
        return iter(self.entries)
   
    def iter_lines(self, filename):
         with open(filename) as hgmd:
            for line in hgmd:
                if line == '\n': continue
                if line[0] == '#': continue
                info = line.split('\t')
                effect = info[7].split(';')[0].split('=')[1]
                omimid = info[7].split(';')[2].split(':')[1]
                pmid = info[7].split(';')[3].split(':')[1]      
                yield Entry(info[0],info[1],info[3],info[4],effect,pmid,omimid)
 
    def __str__(self):
        return self.entries.__str__()
    
    def sample_disease(self, omim, effects=None):
        if effects:
            diseases = self.get_entries_effects(effects)
        else:
            diseases = self.entries
        phenotypes = []
        #run until we get a disease with an omim entry
        while True:
            dis = random.choice(diseases)
            omimd = next((x for x in omim if x.id == dis.omimid), None)
            if not omimd: continue
            for pheno in list(omimd.phenotype_freqs):
                if not omimd.phenotype_freqs[pheno]:
                    phenotypes.append(pheno)
                else:
                    if random.random < omimid.phenotype_freqs[pheno]:
                        phenotypes.append(pheno) 
            return dis, phenotypes
 
    def get_entries_effects(self, effects):
        returned = []
        for effect in effects:
            returned += filter(lambda e:e.effect==effect,self.entries)

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
    file.write('\t'.join([dis.chrom,dis.loc,'.',dis.ref,dis.alt,'50','PASS','.','GT','./.'])+'\n')
    hpo = open(name + '_hpo.txt','w')
    hpo.write(phenotypes[0])
    for p in phenotypes[1:]:
        hpo.write(','+p) 
    file.close() 

def annotate_patient_dir(pdir,hgmd,omim):
    for f in os.listdir(pdir):
        if os.path.isfile(os.path.join(pdir,f)) and (f.endswith('.vcf') or f.endswith('.vcf.gz')):
            annotate_patient(os.path.join(pdir,f),hgmd,omim) 

def script(pheno_file, hgmd_file, patient_path):
    try:
        mim = MIM(pheno_file)
    except FileNotFoundError:
        print >> sys.stderr, "OMIM file not found or invalid"
 
    omim = filter(lambda d:d.db == 'OMIM',mim.diseases)

    try:    
        hgmd = HGMD(hgmd_file)
    except FileNotFoundError:
        print >> sys.stderr, "HGMD file not found or invalid"
   
    hgmd = filter(lambda d:d.omimid,hgmd.entries)         
    
    #if we are given a directory, annotate each vcf.gz or vcf file in the directory assuming it is a patient 
    if os.path.isdir(patient_path):
        print "Given a directory full of patients!"
        annotate_patient_dir(patient_path, hgmd, omim)
    elif os.path.isfile(patient_path):
        print "Given a single patient file!"
        annotate_patient(patient_path, hgmd, omim)
    else:
        print >> sys.stderr, "Patient file/folder not found or invalid"   
   
if __name__ == '__main__':
    print "begun"
    mim = MIM('/dupa-filer/talf/matchingsim/patients/phenotype_annotation.tab')
    omim = filter(lambda d: d.db == 'OMIM',mim.diseases)
    hgmd = HGMD('/dupa-filer/talf/matchingsim/data/hgmd/hgmd_correct.jv.vcf')
    hgmd.entries = filter(lambda d:d.omimid,hgmd.entries)
