#!/usr/bin/env python

#Given a directory, this script will check for which vcf/ezr files have the added variant (assumed to be the last) as the top ranked variant in the ezr file.

__author__ = 'Tal Friedman'

import os
import sys
import logging
import subprocess

import annotate_dir
from collections import defaultdict
from hgmd import HGMD
from orpha import Orphanet

def get_last_line(path):
    with open(path) as file:
        for line in file:
            s = line
    return [s]

def get_last_recessive(path):
    with open(path) as file:
        cont = list(file)
        info = cont[-1].split('\t')
        if info[-1].strip() == '1|1' or info[-1].strip() == '1/1':
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

def script(vcf_ezr_paths, A, R, D, RD, V, N=None):
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
        ezr_files = filter(lambda f: f.endswith('.ezr'), contents)

        # If we're missing some annotation files, just run annotater
        if any(not any(s == f[:-4] + '.txt' for s in contents) for f in ezr_files):
            if R:
                annotate_dir.main([vcf_ezr_path, '-R'])
            else:
                annotate_dir.main([vcf_ezr_path])

        contents = os.listdir(vcf_ezr_path)
        txt_files = filter(lambda f: f.endswith('.txt') and any(s == f[:-4] + '.ezr' for s in contents),contents)
        
        if D:
            hgmd = HGMD('/dupa-filer/talf/matchingsim/patients/hgmd_correct.jv.vcf')
            # Grab orphanet file names and then load orphanet data
            orphanet_lookup = '/dupa-filer/talf/matchingsim/patients/orphanet_lookup.xml'
            orphanet_inher = '/dupa-filer/talf/matchingsim/patients/orphanet_inher.xml'
            orphanet_geno_pheno = '/dupa-filer/talf/matchingsim/patients/orphanet_geno_pheno.xml'
            orph = Orphanet(orphanet_lookup, orphanet_inher, orphanet_geno_pheno)
            lookup = orph.lookup
     
            rev_hgmd = hgmd.get_by_omim()
            Dcounter = defaultdict(lambda:[0,0]) 

        if RD:
            singlecounter = 0
            doublecounter = 0
            singlecorr = 0
            doublecorr = 0
        
        use_orph = False
        for txt in txt_files:
            info = list(open(os.path.join(vcf_ezr_path,txt)))
            rank = info[0].split(' ')[-1].strip() 
            has_orph = info[-1].startswith('Orphanum')
            if has_orph:
                use_orph = True

            if RD:
                single = len(info) == 4 or len(info) == 5 and has_orph
                double = len(info) == 6 or len(info) == 5 and not has_orph
                singlecounter += single
                doublecounter += double
                singlecorr += single and rank == '1'
                doublecorr += double and rank == '1'

            counter += rank == '1'

            if D:
                id = info[-1].split(' ')[-1].strip()
                Dcounter[id][1] += 1
                Dcounter[id][0] += rank == '1'
                if len(Dcounter[id]) == 2:
                    Dcounter[id].append(info[-1].split(' ')[-1].strip())

            if A:
                Acounter += rank.isdigit()
            if N:
                Ncounter += rank.isdigit() and int(rank) <= N

        logging.info('Total Patients: ' + str(len(txt_files)))
        logging.info('Total Patients exomizer ranked inserted variant #1: ' + str(counter))
        if txt_files:
            logging.info('Accuracy of top hits: ' + str(float(counter)/len(txt_files)) + '\n')
        if N:
            logging.info('Total patients exomizer ranked inserted variant top ' + str(N) + ': ' + str(Ncounter))
            if txt_files:
                logging.info('Accuracy of top ' + str(N) + ' hits: ' + str(float(Ncounter)/len(txt_files))+'\n')
        if A:
            logging.info('Total patients exomizer ranked inserted variant at all: ' + str(Acounter))
            if txt_files:
                logging.info('Accuracy of entire file: ' + str(float(Acounter)/len(txt_files))+'\n')
        if D: 
            dis = Dcounter.items()
            if use_orph:
                for i in range(len(dis)):
                    dis[i] = (lookup[dis[i][0]].geno[0], dis[i][1])
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
                logging.info('Phenotypic omim used in association with this genotypic omim: ' + d[1][2] + '\n')
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
            if totalcounter:
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
    parser.add_argument('-RD',help='for AR, give info about accuracy for single gene vs. 2 gene', action='store_true')
    parser.add_argument('-D',help='give info about accuracy per disease',action='store_true')
    parser.add_argument('-V',help='when per disease flag is given, sort diseases by number of associated variants',action='store_true')
    parser.add_argument('vcf_ezr_paths',metavar='DIR',nargs='+',help='the directory where vcf/ezr files are located')
    return parser.parse_args(args)

def main(args = sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
