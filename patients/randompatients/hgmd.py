import os
import sys

from collections import defaultdict

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
    
    def get_entries_effects(self, effects):
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
    hgmd = HGMD('/dupa-filer/talf/matchingsim/patients/hgmd_correct.jv.vcf')
    print "Success!"
