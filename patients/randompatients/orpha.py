#!/usr/bin/env python

import os
import sys
from collections import defaultdict
import xml.etree.ElementTree as ET

class Orphanet:
    def __init__(self, lookup_filename, inher_filename, geno_pheno_filename):
        self.lookup = self.parse_lookup(lookup_filename)
        self.inheritance = self.parse_inheritance(inher_filename, self.lookup)
        self.counter = self.parse_geno_pheno(geno_pheno_filename,self.lookup)
        
    @classmethod
    def parse_lookup(cls, filename):
        tree = ET.parse(filename)
        root = tree.getroot()
        lookup = defaultdict(lambda: [[],[],[]]) # orphanet -> omim
        for disorder in root.findall('.//Disorder'):
            orphanum = disorder.find('OrphaNumber').text
            for ref in disorder.findall('./ExternalReferenceList/ExternalReference'):
                if ref.find('Source').text == 'OMIM':
                    omim = ref.find('Reference').text
                    lookup[orphanum][0].append(omim)
        lookup = {k:v for k,v in lookup.items() if v != [[],[],[]]} 
        return lookup

    @classmethod
    def parse_inheritance(cls, filename, lookup):
        tree = ET.parse(filename)
        root = tree.getroot()
        inheritance = defaultdict(list) # omim -> inheritance pattern list
        for disorder in root.findall('.//Disorder'):
            orphanum = disorder.find('OrphaNumber').text
            #ensure that this disorder has an omim number
            try:
                ids = lookup[orphanum][0]
            except KeyError:
                continue
            for inher in  disorder.findall('./TypeOfInheritanceList/TypeOfInheritance'):
                pattern = inher.find('Name').text
                lookup[orphanum][1].append(pattern)
                for id in ids:
                    inheritance[id].append(pattern)
        return inheritance    
        
    @classmethod
    def parse_geno_pheno(cls, filename, lookup):
        tree = ET.parse(filename)
        root = tree.getroot()
        counter = 0
        for disorder in root.findall('.//Disorder'):
            restart = False
            orphanum = disorder.find('OrphaNumber').text
            for ref in disorder.findall('.//ExternalReferenceList/ExternalReference'):
                if ref.find('Source').text == 'OMIM':
                    omim = ref.find('Reference').text
                    try:
                        lookup[orphanum][2].append(omim)           
                    except KeyError:
                        restart = True
                        counter += 1
                        break
            if restart:
                continue
        return counter

    def write_file(self, filename):
        with open(filename, 'w') as out:
            for k in self.lookup.keys():
                out.write(str(self.lookup[k]) + '\n')

    def write_stats(self, filename):
        with open(filename, 'w') as out:
            n_ideal = 0
            n_miss_inher = 0
            n_miss_pheno = self.counter
            n_miss_geno = 0
            n_one_pheno = 0
            n_one_geno = 0
            n_many_geno = 0
            n_many_pheno = 0
            for k in self.lookup.keys():
                o = self.lookup[k]
                if len(o[0]) == 1 and len(o[1]) == 1 and len(o[2]) == 1: n_ideal += 1
                if len(o[1]) == 0: n_miss_inher += 1
                if len(o[2]) == 0: n_miss_geno += 1
                if len(o[0]) == 1: n_one_pheno += 1
                if len(o[2]) == 1: n_one_geno += 1
                if len(o[0]) > 1: n_many_pheno += 1
                if len(o[2]) > 1: n_many_geno += 1
            
            out.write("PHENO OMIM:" + '\n')
            out.write(str(n_miss_pheno) + " entries missing OMIM Pheno entry\n")
            out.write(str(n_one_pheno) + " entries with one OMIM Pheno entry\n")
            out.write(str(n_many_pheno) + " entries with many OMIM pheno entries\n")
            out.write("INERITANCE:" + '\n')
            out.write(str(n_miss_inher) + " entries missing inheritance pattern\n")
            out.write("GENO OMIM:" + '\n')
            out.write(str(n_miss_geno) + " entries missing OMIM Geno entry\n")
            out.write(str(n_one_geno) + " entries with one OMIM Geno entry\n")
            out.write(str(n_many_geno) + " entries with many OMIM Geno entries\n")
            out.write(str(n_ideal) + " ideal entries (1 of each)")

if __name__ == '__main__':
    orph = Orphanet('/dupa-filer/talf/matchingsim/patients/orphanet_lookup.xml', '/dupa-filer/talf/matchingsim/patients/orphanet_inher.xml', '/dupa-filer/talf/matchingsim/patients/orphanet_geno_pheno.xml')
    orph.write_file('orphanet_parsed.txt')
    orph.write_stats('orphanet_parsed.log')
