import sys


if __name__ == '__main__':
    output = open('final_hgmd.vcf', 'w')
    with open('hgmd_pro_allmut_2013.4') as hgmd:
        with open('out.jv.vcf') as vcf:
            iter = hgmd.__iter__()
            iter.next() 
            for vline in vcf:
                if vline[0] == '#': continue
                hline = iter.next()
                hinfo = hline.split('\t')
                vinfo = vline.split('\t')
                output.write('\t'.join(vinfo[:-1]+hinfo[5:8])+'\n')
