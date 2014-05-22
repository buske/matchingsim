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
            hgmd.append([info[0]] + [info[1]] + ['.'] + info[2:4] + ['50','PASS','.','GT','./.'])
    return hgmd 
def load_vcf(filename):
    vcf = []
    with gzip.open(filename) as file:
        for line in file:
            if line[0] == '#':continue
            info = line.split('\t')
            vcf.append(['chr'+info[0]] + [info[1]] + [info[3]])
    return vcf

if __name__ == '__main__':
    genome = load_genome('/filer/hg18/hg18.fa')
    hgmd = load_hgmd('hgmd_pro_allmut_2013.4')
    counter = 0
    for h in hgmd:
        if h[3] and h[3] != '-' and genome[h[0]][int(h[1])-1:int(h[1])+len(h[3])-1].upper() != h[3]:
           #print h[3]
           #print genome[h[0]][int(h[1])-1:int(h[1])+len(h[3])-1] 
           counter = counter + 1
       # if i > 20: break
    print 'hg18:'
    print counter
    print float(counter)/float(len(hgmd))
    #counter = 0 
    #genome2 = load_genome('/filer/hg19/hg19.fa') 
    #for h in hgmd:
    #    if genome2[h[0]][int(h[1])-1:int(h[1])+len(h[3])-1] != h[3]:
    #        counter = counter + 1
       
    #print 'hg19:'
    #print counter
    #print float(counter)/float(len(hgmd))
    #vcf = load_vcf('/dupa-filer/talf/matchingsim/data/1000gp/samples/complete/HG00096.vcf.gz')
    #for v in vcf:
    #    if genome[v[0]][int(v[1])-1:int(v[1])+len(v[2])-1].upper() != v[2]:
    #        counter = counter + 1

   # print counter
   # print float(counter)/float(len(vcf))
