#!/usr/bin/python
import gzip
import sys

__author__ = 'Tal Friedman'

MAN_COL = 9

def get_names(file):
	for line in file:
		if line[0:3] == '#CH':
			info = line.split('\t')
			return info[MAN_COL:]

def write_line(info, out, index):
	out.write('\t'.join(info[0:9] + [info[index]])+'\n')

if __name__ == '__main__':
	#with gzip.open(sys.argv[1],'rb') as file:
	file = gzip.open(sys.argv[1],'rb')
	names = get_names(file)		
	outs = []
	for i,n in enumerate(names):
		f = gzip.open(sys.argv[2] + "/" + n + '.vcf.gz', 'wb')
		outs.append(f) 
	for line in file:
		if line[0] == '#': continue
		info = line.split('\t')
		for i,n in enumerate(info[MAN_COL:]):
			if n[0] == '1' or n[2] == '1':
				write_line(info,outs[i],i+MAN_COL)
			
