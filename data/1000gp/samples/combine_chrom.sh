#!/usr/bin/env bash
# Flags:
#$ -V
#$ -N "combine"
#$ -pe parallel "1"
#$ -l h_vmem="4G"
#$ -e /home/friedman/sge_logs/samples/
#$ -o /home/friedman/sge_logs/samples/
set -eu
set -o pipefail

temp=$TMPDIR
mkdir -p $temp
for s in `cat /dupa-filer/talf/matchingsim/data/1000gp/samples/names.txt`
do
	cat /dupa-filer/talf/matchingsim/data/1000gp/samples/vcf_header.txt >> $temp/$s.vcf
	for c in {1..22} X
	do
		zcat "/dupa-filer/talf/matchingsim/data/1000gp/samples/chr$c/$s.vcf.gz" >> $temp/$s.vcf
	done
	gzip $temp/$s.vcf
done
mv -v $temp/* /dupa-filer/talf/matchingsim/data/1000gp/samples/complete/	
