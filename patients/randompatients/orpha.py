import os
import sys
from collections import defaultdict
import xml.etree.ElementTree as ET

class Orphanet:
    def __init__(self, lookup_filename, inher_filename, geno_pheno_filename):
        self.lookup = self.parse_lookup(lookup_filename)
        self.inheritance = self.parse_inheritance(inher_filename, self.lookup)
        self.parse_geno_pheno(geno_pheno_filename,self.lookup)

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
                        break
            if restart:
                continue

if __name__ == '__main__':
    orph = Orphanet('/dupa-filer/talf/matchingsim/patients/orphanet_lookup.xml', '/dupa-filer/talf/matchingsim/patients/orphanet_inher.xml', '/dupa-filer/talf/matchingsim/patients/orphanet_geno_pheno.xml')
