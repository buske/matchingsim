#!/usr/bin/env python

#Given a directory, this script will check for which vcf/ezr files have the added variant (assumed to be the last) as the top ranked variant in the ezr file.

__author__ = 'Tal Friedman'

import os
import sys
import logging
import subprocess

from collections import defaultdict
from hgmd import HGMD

def get_last_line(path):
    with open(path) as file:
        for line in file:
            s = line
    return [s]

def get_last_recessive(path):
    with open(path) as file:
        cont = list(file)
        info = cont[-1].split('\t')
        if info[-1].strip() == '1|1':
            return cont[-1:]
        else:
            return cont[-2:]

def get_actual_lines(path):
    with open(path) as file:    
        return filter(lambda x: not x.startswith('#'), list(file))

def is_match(linevs, linee):
    for linev in linevs:
        infov = linev.split('\t')
        infoe = linee.split('\t')
        if infov[0] == infoe[0][3:] and infov[1] == infoe[1]:
            return True
    return False

def script(vcf_ezr_paths, A, R, E, D, RD, V, N=None):
    if E or D:
        hgmd = HGMD('/dupa-filer/talf/matchingsim/patients/hgmd_correct.jv.vcf')
    if D:
        rev_hgmd = hgmd.get_by_omim()
    if N:
        N=int(N[0])
    for vcf_ezr_path in vcf_ezr_paths:
        logging.basicConfig(filename = os.path.join(vcf_ezr_path, 'score.log'), level = logging.INFO, filemode = 'w')
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        formatter = logging.Formatter('%(levelname)s - %(message)s')
        ch.setFormatter(formatter)
        logging.getLogger().addHandler(ch)

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
        
        if D:
            Dcounter = defaultdict(lambda:[0,0]) 

        if E:
            Ecounter = defaultdict(lambda:[0,0])
        
        if RD:
            singlecounter = 0
            doublecounter = 0
            singlecorr = 0
            doublecorr = 0

        for vcf, ezr in zip(vcf_files, ezr_files):
            #if it is a recessive disease we may be pulling 2 lines
            if R:
                v = get_last_recessive(os.path.join(vcf_ezr_path,vcf))
            else:
                v = get_last_line(os.path.join(vcf_ezr_path,vcf))

            elines = get_actual_lines(os.path.join(vcf_ezr_path,ezr))
            efirst = elines[0]
            if RD:
                singlecounter += len(v) == 1
                doublecounter += len(v) == 2
                singlecorr += len(v) == 1 and is_match(v,efirst)
                doublecorr += len(v) == 2 and is_match(v,efirst)

            if is_match(v,efirst):
                counter += 1

            if D:
                id = next(x for x in hgmd.entries if x.chrom == v[0].split('\t')[0] and x.loc == v[0].split('\t')[1]).omimid
                Dcounter[id][1] += 1
                if is_match(v,efirst):
                    Dcounter[id][0] += 1

            if E:
                #find the effect of the inserted variant
                effect = next(x for x in hgmd.entries if x.chrom == v[0].split('\t')[0] and x.loc == v[0].split('\t')[1]).effect
                Ecounter[effect][1] += 1
                if is_match(v,efirst):
                    Ecounter[effect][0] += 1
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
            logging.info('Total patients exomizer ranked inserted variant top ' + str(N) + ': ' + str(Ncounter))
            logging.info('Accuracy of top ' + str(N) + ' hits: ' + str(float(Ncounter)/len(vcf_files)))
        if A:
            logging.info('Total patients exomizer ranked inserted variant at all: ' + str(Acounter))
            logging.info('Accuracy of entire file: ' + str(float(Acounter)/len(vcf_files)))
        if E:
            for k,v in Ecounter.iteritems():
                logging.info('Total patients for effect ' + k + ': ' + str(v[1]))
                logging.info('Accuracy for effect ' + k + ': ' + str(float(v[0])/v[1]))
        if D: 
            dis = Dcounter.items()
            if V:
                dis.sort(key=lambda d: len(rev_hgmd[d[0]]), reverse=True)
            else:
                #we are sorting first by the number of cases, and within that by the success rate
                dis.sort(key=lambda d: d[1][1] + max(float(d[1][0])/d[1][1]-0.001, 0), reverse=True)
            highcounter = 0
            lowcounter = 0
            zerocounter = 0
            totalcounter = 0
            for d in dis:
                logging.info('Total patients for disease ' + d[0] + ': ' + str(d[1][1]))
                logging.info('Accuracy for disease ' + d[0] + ': ' + str(float(d[1][0])/d[1][1]))
                command = "grep " + d[0] + " /dupa-filer/talf/matchingsim/patients/orphanet_geno_pheno.xml | wc -l"
                hits = subprocess.check_output(command,shell=True)
                logging.info('Orphanet hits for genotypic omim ' + d[0] + ': ' + hits.strip())
                logging.info('#of variants for genotypic omim ' + d[0] + ': ' + str(len(rev_hgmd[d[0]]))+ '\n')
                #check what percentage of total is either below 30 or above 60%
                if d[1][1] > 2:
                    totalcounter += d[1][1]
                    if float(d[1][0])/d[1][1] < 0.3: 
                        lowcounter += d[1][1]
                    elif float(d[1][0])/d[1][1] > 0.6: 
                        highcounter += d[1][1]
                    if d[1][0] == 0:
                        zerocounter += d[1][1]

            logging.info(str(lowcounter) + ' patients have a disease with accuracy < 30% (from diseases with > 2 patients)')
            logging.info(str(highcounter) + ' patients have a disease with accuracy > 60% (from diseases with > 2 patients)')
            logging.info(str(zerocounter) + ' patients a have a disease with accuracy 0% (from diseases with > 2 patients)')
            logging.info(str(float(lowcounter + highcounter) / totalcounter) + 'of patients have a disease in one of these ranges of diseases with > 2 patients')
        
        if RD:
            logging.info('Total patients with a single inserted mutation: ' + str(singlecounter))
            logging.info('Accuracy of single inserted mutation: ' + str(float(singlecorr)/ singlecounter))
            logging.info('Total patients with two inserted mutations: ' + str(doublecounter))
            logging.info('Accuracy of two inserted mutations: ' + str(float(doublecorr)/doublecounter))

def parse_args(args):
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Annotate a vcf/ezr filled directory with the top hit success rate')
    parser.add_argument('-R',help='files to analyze were infected with autosomal recessive diseases (default is AD)', action='store_true')
    parser.add_argument('-A',help='check entire ezr ranking for a hit',action='store_true')
    parser.add_argument('-N',help='check if hit is in the top N entries',nargs=1)
    parser.add_argument('-E',help='give info about accuracy per variant effect type',action='store_true')
    parser.add_argument('-RD',help='for AR, give info about accuracy for single gene vs. 2 gene', action='store_true')
    parser.add_argument('-D',help='give info about accuracy per disease',action='store_true')
    parser.add_argument('-V',help='when per disease flag is given, sort diseases by number of associated variants',action='store_true')
    parser.add_argument('vcf_ezr_paths',metavar='DIR',nargs='+',help='the directory were vcf/ezr files are located')
    return parser.parse_args(args)

def main(args = sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
