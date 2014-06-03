#!/usr/bin/
#$ -V
#$ -N "getlast"
#$ -pe parallel "1"
#$ -l h_vmem="4G"
#$ -e /home/friedman/sge_logs/samples/
#$ -o /home/friedman/sge_logs/samples/

set -eu
set -o pipefail

for chr in {1..22} X
do
    vcf-subset -e -c NA20828 /dupa-filer/talf/matchingsim/data/1000gp/exome/$chr.vcf.gz | gzip > /dupa-filer/talf/matchingsim/data/1000gp/samples/chr$chr/NA20828.vcf.gz
    touch /dupa-filer/talf/matchingsim/data/1000gp/samples/chr$chr/finished
done
     
