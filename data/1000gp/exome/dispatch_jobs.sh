#!/usr/bin/env bash

# A few bash settings that make the script much safer
# (they cause it to die when things go wrong, rather than
# powering on)
set -eu
set -o pipefail

# Memory to reserve on the node, e.g. 100M, 4.2G, etc.
memory=4G
# Number of processors to reserve on the node
processors=1
logdir=~/sge_logs/

mkdir -pv $logdir

for chrom in {1..22} X; do
    # Create a bash script in the current directory
    # for the specific job we want to run
    script="dispatch_$chrom.sh"
    cat > "$script" <<EOF
#!/usr/bin/env bash
# Flags for qsub:
#$ -V
#$ -N "chr$chrom"
#$ -pe parallel "$processors"
#$ -l h_vmem="$memory"
#$ -e $logdir
#$ -o $logdir

set -x
set -eu
set -o pipefail

temp=\$TMPDIR/$chrom.temp
out=/dupa-filer/talf/matchingsim/data/1000gp/exome/$chrom.vcf.gz
zcat /data/1000genomes/release/20110521/ALL.chr$chrom.phase1_release_v2.20101123.snps_indels_svs.vcf.gz \
  | bedtools intersect -u -header -a stdin -b /dupa-filer/talf/matchingsim/data/nrefgene.bed \
  | gzip \
  >> \$temp

# Two-step move for safety.
# The first is a copy across the network (slow)
# The second is a rename on the same file system (very fast)
# Thus, the output file will only exist if everything was successful
mv -v \$temp \$out.temp
mv -v \$out.temp \$out
EOF

    # Submit the script to the cluster
    qsub -S /bin/sh "$script"
done
