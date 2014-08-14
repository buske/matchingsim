#!/usr/bin/env bash

set -eu
set -o pipefail

function usage {
    cat <<EOF
Usage: $0 dir

Dispatch Exomiser on vcf files in dir.
Script is idempotent, so can be run again to finish failed jobs.
EOF
    exit 1
}


if [ $# -ne 1 ]; then
    usage
fi

out=$1
max_jobs=20
#sge
memory=10g
processors=1



function sge_wait {
    sleep_time=1  # seconds
    # Check SGE for number of jobs
    local sge_jobs="$(qstat | grep "EZR2_" | wc -l)" || true
    while [[ $sge_jobs -ge $max_jobs ]]; do
	sleep $sleep_time
	sleep_time=$(expr $sleep_time "*" 2)
	sge_jobs="$(qstat | grep "EZR2_" | wc -l)" || true
    done
}


logdir=~/sge_logs/gen_exomise/"$out"
mkdir -pv ~/sge_logs/gen_exomise/"$out"
mkdir -pv "$out/scripts"
for file in `ls $out | grep .vcf | grep -v .vcf.results.vcf`; do
    f=`basename $file .vcf`
    outfile="$out/$f.vcf.results.vcf"
    if [[ -s $outfile ]]; then
	echo "Output file already exists: $outfile" >&2
	continue
    fi

    script="$out/scripts/dispatch_$f.sh"
    cat > "$script" <<EOF
#!/usr/bin/env bash
#$ -V
#$ -N "EZR2_$f"
#$ -pe parallel "$processors"
#$ -l h_vmem="$memory"
#$ -e $logdir
#$ -o $logdir
#$ -l hostname="supa*"

set -eu
set -o pipefail
temp=\$TMPDIR/$f.ezr

echo "Current directory: \$(pwd)" >&2
echo "Temp directory: \$temp" >&2
echo "Input file: $out/$f.vcf" >&2

test -s $out/$f.vcf

#only unzip if the unziped file doesn't already exist (i.e. only unzip on first run)

java -Xms5g -Xmx5g -jar /filer/tools/exomiser/2.0.0/exomiser-2.0.0.jar -D /filer/tools/exomiser/2.0.0/data/ucsc_hg19.ser -I AD -F 1 -W /filer/tools/exomiser/2.0.0/data/rw_string_9_05.gz -X /filer/tools/exomiser/2.0.0/data/rw_string_9_05_id2index.gz --hpo_ids `cat $out/"$f"_hpo.txt` -v $out/$f.vcf --vcf_output -o \$temp -P

mv -v \$temp $out/$f.ezr.temp
mv -v $out/$f.ezr.temp $outfile
EOF
    # Wait for space on cluster
    sge_wait
    # Submit job
    qsub -S /bin/sh "$script"
done

