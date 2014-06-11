#!/usr/bin/env python

#Given a directory, this script will check for which vcf/ezr files have the added variant (assumed to be the last) as the top ranked variant in the ezr file.

__author__ = 'Tal Friedman'

import os
import sys

def get_last_line(path):
    with open(path) as file:
        for line in file:
            s = line
    return s

def get_first_line(path):
    with open(path) as file:
        for line in file:
            if line.startswith('#'): continue
            return line

def script(vcf_ezr_path):
    counter = 0
    contents = os.listdir(vcf_ezr_path)
    vcf_files = filter(lambda f: f.endswith('.vcf'),contents)
    ezr_files = filter(lambda f: f.endswith('.ezr'),contents)
    #give warning if we're missing ezr files, and take only matched vcf's; likely job didn't complete 
    if len(vcf_files) > len(ezr_files):
        print >> sys.stderr, "WARNING: vcf/ezr numbers do not align, your job probably didn't finish"
        vcf_files = filter(lambda f: ''.join([f[:-4],'.ezr']) in ezr_files, vcf_files) 
   #sort file lists so that when we loop through them in parallel they are in same order
    vcf_files.sort()
    ezr_files.sort()
    
    for vcf, ezr in zip(vcf_files, ezr_files):
        last = get_last_line(os.path.join(vcf_ezr_path,vcf))
        first = get_first_line(os.path.join(vcf_ezr_path,ezr))
        linfo = last.split('\t')
        finfo = first.split('\t')
        if linfo[0] == finfo[0][3:] and linfo[1] == finfo[1]:
            counter += 1

    with open(os.path.join(vcf_ezr_path, 'score.txt'), 'w') as file:
        file.write('Total Patients: ' + str(len(vcf_files)) + '\n')
        file.write('Total Patients exomizer ranked inserted variant #1: ' + str(counter) + '\n')
        file.write('Accuracy: ' + str(float(counter)/len(vcf_files)) + '\n')

def parse_args(args):
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Annotate a vcf/ezr file with the top hit success rate')
    parser.add_argument('vcf_ezr_path',metavar='DIR', help='the directory were vcf/ezr files are located')
    return parser.parse_args(args)

def main(args = sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
