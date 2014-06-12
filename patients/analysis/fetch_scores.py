#!/usr/bin/env python

#Given a directory, this script will check for which vcf/ezr files have the added variant (assumed to be the last) as the top ranked variant in the ezr file.

__author__ = 'Tal Friedman'

import os
import sys
import logging

def get_last_line(path):
    with open(path) as file:
        for line in file:
            s = line
    return s

def get_actual_lines(path):
    with open(path) as file:    
        return filter(lambda x: not x.startswith('#'), list(file))

def is_match(linev, linee):
    infov = linev.split('\t')
    infoe = linee.split('\t')
    return infov[0] == infoe[0][3:] and infov[1] == infoe[1]

def script(vcf_ezr_path, A, N=None):
    logging.basicConfig(filename = os.path.join(vcf_ezr_path, 'score.log'), level = logging.INFO)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logging.getLogger().addHandler(ch)

    if N:
        N=int(N[0])
    
    counter = 0
    Acounter = 0
    Ncounter = 0
    contents = os.listdir(vcf_ezr_path)
    vcf_files = filter(lambda f: f.endswith('.vcf'),contents)
    ezr_files = filter(lambda f: f.endswith('.ezr'),contents)
    #give warning if we're missing ezr files, and take only matched vcf's; likely job didn't complete 
    if len(vcf_files) > len(ezr_files):
        logging.warning("vcf/ezr numbers do not align, your job probably didn't finish")
        vcf_files = filter(lambda f: ''.join([f[:-4],'.ezr']) in ezr_files, vcf_files) 
   #sort file lists so that when we loop through them in parallel they are in same order
    vcf_files.sort()
    ezr_files.sort()
    
    for vcf, ezr in zip(vcf_files, ezr_files):
        v = get_last_line(os.path.join(vcf_ezr_path,vcf))
        elines = get_actual_lines(os.path.join(vcf_ezr_path,ezr))
        efirst = elines[0] 
        if is_match(v,efirst):
            counter += 1
        if A:
            if any(is_match(v,x) for x in elines):
                Acounter += 1
        if N:
            if any(is_match(v,x) for x in elines[:N]):
                Ncounter += 1

    logging.info('Total Patients: ' + str(len(vcf_files)))
    logging.info('Total Patients exomizer ranked inserted variant #1: ' + str(counter))
    logging.info('Accuracy of top hits: ' + str(float(counter)/len(vcf_files)))
    if N:
        logging.info('Total patients exomizer ranked inserted variant top #' + str(N) + ': ' + str(Ncounter))
        logging.info('Accuracy of top ' + str(N) + ' hits: ' + str(float(Ncounter)/len(vcf_files)))
    if A:
        logging.info('Total patients exomizer ranked inserted variant at all: ' + str(Acounter))
        logging.info('Accuracy of entire file: ' + str(float(Acounter)/len(vcf_files)))


def parse_args(args):
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Annotate a vcf/ezr file with the top hit success rate')
    parser.add_argument('-A',help='check entire ezr ranking for a hit',action='store_true')
    parser.add_argument('-N',help='check if hit is in the top N entries',nargs=1)
    parser.add_argument('vcf_ezr_path',metavar='DIR', help='the directory were vcf/ezr files are located')
    return parser.parse_args(args)

def main(args = sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
