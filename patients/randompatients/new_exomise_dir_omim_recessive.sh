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
    local sge_jobs="$(qstat | grep "EZR_" | wc -l)" || true
    while [[ $sge_jobs -ge $max_jobs ]]; do
	sleep $sleep_time
	sleep_time=$(expr $sleep_time "*" 2)
	sge_jobs="$(qstat | grep "EZR_" | wc -l)" || true
    done
}


logdir=~/sge_logs/gen_exomise/"$out"
mkdir -pv ~/sge_logs/gen_exomise/"$out"
mkdir -pv "$out/scripts"
for file in $out/*.vcf; do
    f=`basename $file .vcf`
    omim=`echo $f | cut -f 3 -d '-'`
    outfile="$out/$f.ezr"
    if [[ -s $outfile ]]; then
	echo "Output file already exists: $outfile" >&2
	continue
    fi

    script="$out/scripts/dispatch_$f.sh"
    cat > "$script" <<EOF
#!/usr/bin/env bash
#$ -V
#$ -N "EZR_$f"
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
java -Xmx5000m -Xms1000m -jar /data/Exomiser/Exomizer.jar --db_url jdbc:postgresql://supa01.biolab.sandbox/nsfpalizer -D /data/Exomiser/ucsc.ser -I AR -F 1 -A $omim -v $out/$f.vcf --vcf_output -o \$temp -P

mv -v \$temp $out/$f.ezr.temp
mv -v $out/$f.ezr.temp $outfile
EOF
    # Wait for space on cluster
    sge_wait
    # Submit job
    qsub -S /bin/sh "$script"
done

