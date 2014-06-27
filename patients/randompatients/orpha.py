#!/usr/bin/env python

"""
Parser a collection of orphanet files together to get a useable lookup table of 
diseases with associated genotypic and phenotypic OMIMs as well as inheritance
patterns.
"""


import os
import sys
import logging

from collections import defaultdict
import xml.etree.ElementTree as ET


__author__ = 'Tal Friedman (talf301@gmail.com)'

class Disease:
    def __init__(self):
        self.pheno = []
        self.inheritance = []
        self.geno = []

class Orphanet:
    def __init__(self, lookup_filename, inher_filename, geno_pheno_filename):
        self.lookup = self.parse_lookup(lookup_filename)
        self.inheritance = self.parse_inheritance(inher_filename, self.lookup)
        # Counter for when we write stats
        self.counter = self.parse_geno_pheno(geno_pheno_filename,self.lookup)
    
    @classmethod
    def parse_lookup(cls, filename):
        """
        Parse the external reference file for orphanet

        Args:
            filename: string path of file

        Return:
            Orphanet # -> Disease dict, with each entry having nonempty
            pheno and empty inheritance and geno
        """
        tree = ET.parse(filename)
        root = tree.getroot()
        lookup = defaultdict(Disease) # orphanet -> omim
        counter = 0
        for disorder in root.findall('.//Disorder'):
            orphanum = disorder.find('OrphaNumber').text
            for ref in disorder.findall('./ExternalReferenceList/ExternalReference'):
                if ref.find('Source').text == 'OMIM':
                    omim = ref.find('Reference').text
                    # Ensure Orphanum and OMIM are numeric
                    try:
                        int(omim)
                        int(orphanum)
                        lookup[orphanum].pheno.append(omim)
                        counter += 1
                    except ValueError:
                        logging.error("Malformed OMIM or Orphanum at %s" % orphanum)
        logging.info("%d orphanet diseases initially found" % counter)
        return lookup

    @classmethod
    def parse_inheritance(cls, filename, lookup):
        """
        Parse the inheritance reference file for orphanet

        Args:
            filename: string path of file
            lookup: Orphanet # -> Disease dict, expect each pheno entry
            to be nonempty
        """
        tree = ET.parse(filename)
        root = tree.getroot()
        counter = 0
        for disorder in root.findall('.//Disorder'):
            orphanum = disorder.find('OrphaNumber').text

            # Ensure that this disorder has an omim number
            if orphanum in lookup:
                for inher in  disorder.findall('./TypeOfInheritanceList/TypeOfInheritance'):
                    pattern = inher.find('Name').text
                    assert pattern in ['X-linked dominant', 'Mitochondrial inheritance', \
                            'Unknown', 'Autosomal recessive', 'Multigenic/multifactorial', \
                            'X-linked recessive', 'Sporadic', 'Autosomal dominant', \
                            'No data available'], "Unrecognized inheritance pattern %s" % pattern
                    lookup[orphanum].inheritance.append(pattern)
            else:
                counter += 1
        logging.warning("%d unmatched orphanums from inheritance" % counter)

    @classmethod
    def parse_geno_pheno(cls, filename, lookup):
        """
        Parse the genotypic omim reference file for orphanet.

        Args: 
            filename: string path of file
            lookup: Orphanet # -> Disease dict, we expect each entry
            will have non-empty pheno and inheritance entries and empty
            geno entry

        Return:
            A count of the number of orphanet diseases found in this file without a 
            genotypic mapping.
        """
        tree = ET.parse(filename)
        root = tree.getroot()
        counter = 0
        for disorder in root.findall('.//Disorder'):
            restart = False
            orphanum = disorder.find('OrphaNumber').text
            for ref in disorder.findall('.//ExternalReferenceList/ExternalReference'):
                if ref.find('Source').text == 'OMIM':
                    omim = ref.find('Reference').text
                    # Ensure OMIM is numeric
                    try:
                        int(omim)
                    except ValueError:
                        logging.error("Malformed OMIM %s" % omim)

                    try:
                        lookup[orphanum].geno.append(omim)           
                    except KeyError:
                        restart = True
                        counter += 1
                        break
            if restart:
                continue
        logging.warning("%d Disorders were unmatched to a phenotypic omim" % counter)
        return counter

    def write_file(self, filename):
        """Serialize this Orphanet object to a file"""
        with open(filename, 'w') as out:
            for k in self.lookup.keys():
                out.write(str(self.lookup[k]) + '\n')

    def write_stats(self):
        """Logs some stats about this Orphanet object"""
        n_ideal = 0
        n_miss_inher = 0
        n_miss_pheno = self.counter
        n_miss_geno = 0
        n_one_pheno = 0
        n_one_geno = 0
        n_many_geno = 0
        n_many_pheno = 0
        for o in self.lookup.itervalues():
            if len(o.pheno) == 1 and len(o.inheritance) == 1 and len(o.geno) == 1:
                n_ideal += 1
            if len(o.inheritance) == 0: n_miss_inher += 1
            if len(o.geno) == 0: n_miss_geno += 1
            if len(o.pheno) == 1: n_one_pheno += 1
            if len(o.geno) == 1: n_one_geno += 1
            if len(o.pheno) > 1: n_many_pheno += 1
            if len(o.geno) > 1: n_many_geno += 1
        
        logging.info("PHENO OMIM:" + '\n')
        logging.info(str(n_miss_pheno) + " entries missing OMIM Pheno entry\n")
        logging.info(str(n_one_pheno) + " entries with one OMIM Pheno entry\n")
        logging.info(str(n_many_pheno) + " entries with many OMIM pheno entries\n")
        logging.info("INERITANCE:" + '\n')
        logging.info(str(n_miss_inher) + " entries missing inheritance pattern\n")
        logging.info("GENO OMIM:" + '\n')
        logging.info(str(n_miss_geno) + " entries missing OMIM Geno entry\n")
        logging.info(str(n_one_geno) + " entries with one OMIM Geno entry\n")
        logging.info(str(n_many_geno) + " entries with many OMIM Geno entries\n")
        logging.info(str(n_ideal) + " ideal entries (1 of each)")

    @classmethod
    def has_pattern(cls, patterns, dis):
        """Return if disease follows any of the given inheritance patterns"""
        return any(x in patterns for x in dis.inheritance)

    @classmethod
    def has_pheno(cls, omim_dict, dis):
        """Return if disease has HPO terms associated with it""" 
        return dis.pheno[0] in omim_dict

    @classmethod
    def filter_lookup(cls, lookup, omim_dict, rev_hgmd, Inheritance=None):
        """Return a new lookup table filtered based on inheritance and mapping ability

        Args:
            lookup: dict of Orphanet Number -> orph.Disease instances
            omim_dict: OMIM number -> omim.Disease
            rev_hgmd: OMIM number -> list(hgmd.Entry)
            Inheritance: list(string) containing accepted inheritance pattern
        
        Return:
            dict of Orphanet Number -> orph.Disease each with exactly one geno,pheno, 
            and inheritance entry, and the matching inheritance pattern, as well as
            having at least one corresponding variant in hgmd
        """
        # Get ideal orphanet cases
        newlook = {orph:dis for orph, dis in lookup.iteritems() 
                if len(dis.pheno) == 1 and len(dis.inheritance) == 1 and len(dis.geno) == 1}

        # Get the right disease set based on inheritance
        if Inheritance:
            patterns = []
            if 'AD' in Inheritance:
                 patterns.append('Autosomal dominant')
            if 'AR' in Inheritance:
                patterns.append('Autosomal recessive')
            newlook = {orph:dis for orph, dis in newlook.iteritems() 
                    if cls.has_pattern(patterns, dis)}
        
        # Ensure all orphanet cases have phenotypic annotations
        lookup = {orph:dis for orph, dis in newlook.iteritems() 
                if cls.has_pheno(omim_dict, dis)}
        
        # Ensure all orphanet cases have at least one associated variant
        newlook = {}
        for orph, dis in lookup.iteritems():
            if dis.geno[0] in rev_hgmd:
                newlook[orph] = dis

        return newlook

if __name__ == '__main__':
    orph = Orphanet('/dupa-filer/talf/matchingsim/patients/orphanet_lookup.xml', '/dupa-filer/talf/matchingsim/patients/orphanet_inher.xml', '/dupa-filer/talf/matchingsim/patients/orphanet_geno_pheno.xml')
    orph.write_file('orphanet_parsed.txt')
    orph.write_stats()
