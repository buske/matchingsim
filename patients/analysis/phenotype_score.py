#!/usr/bin/env python


import os
import sys
import logging

from argparse import ArgumentParser
from collections import defaultdict
__author__ = 'Tal Friedman (talf301@gmail.com)'


def script(res_file, annotate, **kwargs):
    logging.basicConfig(filename = os.path.join(os.path.dirname(res_file), 'pheno_score.log'), level = logging.INFO, filemode = 'w')
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logging.getLogger().addHandler(ch)

    # If we need to annotate, do it first
    if annotate:
        with open(res_file) as res:
            with open(res_file + '.annotated', 'w') as anno:
                for line in res:
                    if line.startswith('#'): continue
                    line = line.strip()
                    tokens = line.split('\t')
                    assert len(tokens) > 3, "%s" % line
                    # Check by name signature if they come from same disease
                    if tokens[0].split('_')[-2] == tokens[1].split('_')[-2]:
                        anno.write('\t'.join(['1',line]) + '\n')
                    else:
                        anno.write('\t'.join(['0',line]) + '\n')

        # Next we find out how many patients have one with the same disease as top hit
    with open(res_file) as res: 
        top_counter = 0
        # dict with Patient ID -> list of tuples (patient, score)
        scores = defaultdict(set)
        # Parse pairwise scores
        for line in res:
            if line.startswith('#'): continue
            line = line.strip()
            tokens = line.split('\t')
            assert len(tokens) >= 3
            first = tokens[0]
            second = tokens[1]
            score = float(tokens[2])
            scores[first].add((second, score))
            scores[second].add((first, score))
        
        for name, scoreset in scores.iteritems():
            scoreset = list(scoreset)
            scoreset.sort(key=lambda x: x[1], reverse=True)
            if scoreset[0][0].split('_')[-2] == name.split('_')[-2]:
                top_counter += 1

        logging.info("Total patients: %d\n" % len(scores))
        logging.info("Patients where top hit was the same disease: %d\n" % top_counter)
        logging.info("Total accuracy of top hit: %f\n" % (float(top_counter)/len(scores)))

def parse_args(args):
    parser = ArgumentParser()
    parser.add_argument('res_file', metavar='RESULTS')
    parser.add_argument('-A', dest='annotate', action='store_true')
    return parser.parse_args(args)

def main(args = sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
