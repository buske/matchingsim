import os
import sys
import gzip

from collections import defaultdict
from sets import Set

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

class Entry:
    def __init__(self, chrom, loc, ref, alt, pmid,omimid):
        self.chrom = chrom
        self.loc = loc
        self.ref = ref
        self.alt = alt
        self.pmid = pmid
        self.omimid = omimid
    
    def __str__(self):
        return [self.chrom, self.loc, self.ref, self.alt, self.pmid, self.omimid].__str__()

    def __repr__(self):
        return self.__str__()
 
class HGMD:
    def __init__(self, filename):
        self.entries = list(self.iter_lines(filename))
    
    def __iter__(self):
        return iter(self.entries)
    
    #NOTE: Jannovar messes up on chr17 41243453
    @classmethod
    def write_vcf(cls, filename, entries):
        with open(filename, 'w') as out:
            out.write('##fileformat=VCFv4.1\n')
            out.write('\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER', 'INFO', 'FORMAT', 'FIRST']) + '\n') 
            for e in entries:
                if e.loc:
                    out.write('\t'.join([e.chrom[3:],e.loc,'.',e.ref,e.alt,'50','PASS','OMIM:' + e.omimid + ';PM:' + e.pmid,'GT','./.'])+'\n')
            out.close()
    @classmethod
    def write_vcf_17(cls, filename, entries):
        final = filter(lambda x: x.chrom == 'chr17',entries)
        HGMD.write_vcf(filename, final)

    def iter_lines(self, filename):
         with open(filename) as hgmd:
            for line in hgmd:
                if line[0] == '#':continue
                info = line.split('\t')
                if len(info) < 8 or not info[7]:continue
                if info[7][-1] == '\n': info[7] = info[7][:-1]
                yield Entry(info[0],info[1],info[2],info[3],info[6], info[7])

    def __str__(self):
        return self.entries.__str__()

def load_hgmd(filename):
    return HGMD(filename).entries 

def load_vcf_gz(filename):
    vcf = {}
    for i in range(1,23):
        vcf['chr'+str(i)] = []
    vcf['chrX'] = []
    vcf['chrY'] = []   
    with gzip.open(filename) as file:
        for line in file:
            if line[0] == '#':continue
            info = line.split('\t')
            chr = 'chr' + info[0]
            if chr in list(vcf):
                vcf[chr].append([chr] + [info[1]] + info[3:4])
    return vcf

def load_hgmd_vcf(filename):
    entries = []
    with open(filename) as file:
        for line in file:
            info = line.split('\t')
            if len(info) < 10: continue 
            print line
            omimid = info[7].split(';')[0].split(':')[1]
            pmid = info[7].split(';')[1].split(':')[1]
            entries.append(Entry('chr'+info[0],info[1],info[3],info[4],pmid,omimid))        
    return entries    

def load_vcf(filename):
    vcf = {}
    for i in range(1,23):
        vcf['chr'+str(i)] = []
    vcf['chrX'] = []
    vcf['chrY'] = []   
    with open(filename) as file:
        for line in file:
            if line[0] == '#':continue
            info = line.split('\t')
            chr = 'chr' + info[0]
            if chr in list(vcf):
                vcf[chr].append([chr] + [info[1]] + info[3:5])
    return vcf

def create_fa(filename,hgmd):
    long = filter(lambda x:len(x.ref) >= 20, hgmd)
    out = open(filename,'w')
    for l in long:
        out.write('>'+out.loc+'\n')
        out.write(out.ref + '\n')
    out.close()

def remove_dash(hgmd, genome):
    new = []
    for h in hgmd:
        if h.ref == '-' or h.alt == '-':
            nloc = str(int(h.loc) - 1)
            nref = h.ref.strip('-')
            nalt = h.alt.strip('-')
            add = genome[h.chrom][int(nloc) - 1]
            nref = add + nref
            nalt = add + nalt
            new.append(Entry(h.chrom,nloc,nref,nalt,h.pmid,h.omimid))
        else:
            new.append(Entry(h.chrom,h.loc,h.ref,h.alt,h.pmid,h.omimid))
    return new
   #[start, end)
def intersect_neg(e, refbed):
    for r in refbed[e.chrom]:
        if max(int(e.loc)-1, r[1]) < min(int(e.loc)+max(len(e.ref), len(e.alt)) - 1, r[2]):
            return True 
    return False            

def load_refbed(filename):
    ref = {}
    for i in range(1,23):
        ref['chr'+str(i)] = []
    ref['chrX'] = []
    ref['chrY'] = []
    with open(filename) as file:
        for line in file:
            info = line.split('\t')
            if not info[0] in list(ref):
                return ref
            if info[5].strip('\n') == '-':
                ref[info[0]].append([info[0],int(info[1]),int(info[2])])


def reverse_complement(seq):
    if seq == '-': return seq
    seq = seq[::-1]
    sw = {'A':'T', 'T':'A','G':'C','C':'G'}
    seq = map(lambda x: sw[x],seq)
    return ''.join(seq) 
 
def get_correct(hgmd, genome, refbed):
    ret = []
    counter = 0
    for e in hgmd:
        counter += 1
        if counter % 1000 == 0: print counter
        #if it fits naively
        if genome[e.chrom][int(e.loc)-1:int(e.loc)+len(e.ref)-1].upper() == e.ref.upper():
            ret.append(e)
            continue    
        #if it is snp, also try reversing ref and alt (for rare ref cases)
        if len(e.ref) == 1 and len(e.alt) == 1 and e.ref != '-' and e.alt != '-':
            if genome[e.chrom][int(e.loc)-1].upper() == e.alt.upper():
                ret.append(Entry(e.chrom, e.loc, e.alt, e.ref, e.pmid, e.omimid))
                continue
        #anything involving negative strand; only want to check strand once
        if intersect_neg(e,refbed):
            if e.ref != '-':
                if genome[e.chrom][int(e.loc)-1:int(e.loc)+len(e.ref)-1].upper() == reverse_complement(e.ref.upper()):
                    ret.append(Entry(e.chrom, e.loc, reverse_complement(e.ref), reverse_complement(e.alt), e.pmid, e.omimid))
                    continue
    return ret 

def get_long(hgmd):
    return filter(lambda x: len(x.ref) >= 5, hgmd) 

def get_snp(hgmd):
    return filter(lambda x: len(x.ref) == 1 and len(x.alt) == 1 and x.ref != '-' and x.alt != '-',hgmd) 

def test_accuracy(hgmd, genome, refbed):
    return float(len(get_correct(hgmd,genome,refbed)))/float(len(hgmd))    

def found_vcf(entry,vcf):
    for v in vcf[entry.chrom]:
        if v[1] == entry.loc and v[2] == entry.ref and v[3] == entry.alt:
            return True
    return False

def get_found_vcf(hgmd, vcf):
    return filter(lambda x: found_vcf(x,vcf), hgmd)    
 
def get_unique_diseases(hgmd):
    s = Set()
    for h in hgmd:
        s.add(h.pmid)
    return len(s) 

def get_unique_omim(hgmd):
    s = Set()
    for h in hgmd:
        s.add(h.omimid)
    return len(s) 

#what is currently being used for hgmd
def current_hgmd(hgmd,genome,refbed):
    hgmd = filter(lambda x: x.ref and x.alt and x.ref != '-', hgmd)
    return get_correct(hgmd,genome,refbed)

if __name__ == '__main__':
    hgmd = load_hgmd('hgmd_pro_allmut_2013.4')
    #genome = load_genome('/filer/hg19/hg19.fa') 
    #refbed = load_refbed('/dupa-filer/talf/matchingsim/data/refgene.bed')
    #dbsnp = load_vcf('dbSnp.vcf')   
    #vcf = load_vcf('/dupa-filer/talf/matchingsim/data/1000gp/samples/complete/HG00096.vcf.gz')
    #for v in vcf:
    #    if genome[v[0]][int(v[1])-1:int(v[1])+len(v[2])-1].upper() != v[2]:
    #        counter = counter + 1

   # print counter
   # print float(counter)/float(len(vcf))
   
