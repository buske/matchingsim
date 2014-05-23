import os
import sys
import gzip

from collections import defaultdict


def load_genome(filename):
    """Return dict: sequence -> str from multiFASTA file"""
    genome = {}
    name = None
    sequences = []
    with open(filename) as ifp:
        for line in ifp:
            line = line.strip()
            if line.startswith('>'):
                if name:
                    genome[name] = ''.join(sequences)

                name = line.lstrip('>').strip()
                sequences = []
            elif line:
                sequences.append(line)

        if name:
            genome[name] = ''.join(sequences)
    return genome



#load hgmd and dummy values we'll need
def load_hgmd(filename):
    hgmd = []
    with open(filename) as file:
        for line in file:
            if line[0] == '#':continue
            info = line.split('\t')
            hgmd.append([info[0]] + [info[1]] + ['.'] + info[2:4] + ['50','PASS','.','GT','./.']+info[5:7])
    return hgmd 
def load_vcf(filename):
    vcf = {}
    with gzip.open(filename) as file:
        for line in file:
            if line[0] == '#':continue
            info = line.split('\t')
            chr = 'chr' + info[0]
            vcf[chr].append([chr] + [info[1]] + info[3:4])
    return vcf

def create_fa(filename,hgmd):
    long = filter(lambda x:len(x[3]) >= 20, hgmd)
    out = open(filename,'w')
    for l in long:
        out.write('>'+l[1]+'\n')
        out.write(l[3] + '\n')
    out.close()

def get_correct(hgmd, genome):
    return filter(lambda x: x[3] and genome[x[0]][int(x[1])-1:int(x[1])+len(x[3])-1].upper() == x[3].upper(), hgmd)

def test_accuracy(hgmd, genome):
    counter = 0
    for h in hgmd:
        if not h[3] or h[3] == '-' or  genome[h[0]][int(h[1])-1:int(h[1])+len(h[3])-1].upper() == h[3].upper():
            counter = counter + 1
    return float(counter)/float(len(hgmd))
    
if __name__ == '__main__':
    #genome = load_genome('/filer/hg18/hg18.fa')
    hgmd = load_hgmd('hgmd_pro_allmut_2013.4')
    #genome = load_genome('/filer/hg19/hg19.fa') 
    #dbsnp = load_vcf('dbSnp.vcf')   
    #vcf = load_vcf('/dupa-filer/talf/matchingsim/data/1000gp/samples/complete/HG00096.vcf.gz')
    #for v in vcf:
    #    if genome[v[0]][int(v[1])-1:int(v[1])+len(v[2])-1].upper() != v[2]:
    #        counter = counter + 1

   # print counter
   # print float(counter)/float(len(vcf))
   
