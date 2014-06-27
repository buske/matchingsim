import os
import sys
import re
import logging

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
        counter = 0
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
            
            if all(not freq for freq in phenotype_freqs.itervalues()) and db == 'OMIM':
                counter += 1

            disease = Disease(db, id, name, phenotype_freqs)
            yield disease
        
        logging.warning("%d OMIM diseases had no freq info at all" % counter)
if __name__ == '__main__':
    omim = MIM('/dupa-filer/talf/matchingsim/patients/phenotype_annotation.tab')
