#!/usr/bin/env python


"""
Parse an hgmd file in vcf format annotated (by Jannovar), returning an easy to use format.
"""


import os
import sys
import logging

from collections import defaultdict


__author__ = 'Tal Friedman (talf301@gmail.com)'

class Entry:
    def __init__(self, chrom, loc, ref, alt, effect, pmid, omimid, info_line):
        self.chrom = chrom
        self.loc = loc
        self.ref = ref
        self.alt = alt
        self.effect = effect
        self.omimid = omimid
        self.pmid = pmid
        self.info_line = info_line

        # Casts will raise an error if things that should be ints aren't ints
        int(loc)
        int(omimid)
    
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
        """Iterate through the lines in given hgmd file, with each
        line yielding a corresponding Entry object.
        """
        with open(filename) as hgmd:
            for line in hgmd:
                if line == '\n': continue
                if line[0] == '#': continue

                info = line.rstrip().split('\t')
                assert len(info) == 10, "Malformed line %s" % line
                
                try:
                    info_line = info[7]
                    effect = info_line.split(';')[0].split('=')[1]
                    omimid = info_line.split(';')[2].split(':')[1]
                    pmid = info_line.split(';')[3].split(':')[1]      
                except IndexError:
                    logging.error("Malformed line %s" % line)
                try:
                    yield Entry(info[0],info[1],info[3],info[4],effect,pmid,omimid,info[7].strip())
                except ValueError:
                    logging.error("Malformed line %s" % line)
 
    def __str__(self):
        return self.line

    def get_entries_effects(self, effects):
        """Get all hgmd variants with a certain effect
        
        Args:
            effects: A list of strings representing the desired effects

        Returns:
            A list a of hgmd variants with the given effect
        """
        returned = []
        for e in self.entries:
            if e.effect in effects:
                returned += e
        return returned

    def get_by_omim(self):
        oentries = filter(lambda x: x.omimid, self.entries)
        ret = {}
        for o in oentries:
            ret.setdefault(o.omimid, []).append(o)
        return ret

if __name__ == '__main__':
    try:
        hgmd = HGMD('/dupa-filer/talf/matchingsim/patients/hgmd_correct.jv.vcf')
    except IOError, e:
        logging.error(e)
    print "Success!"
