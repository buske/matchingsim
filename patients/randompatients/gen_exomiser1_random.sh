#!/usr/bin/env bash

set -eu
set -o pipefail
#Use date and time as a signature for the generated data
sig=`date +%F-%H-%M-%S`
#sge stuff
memory=4G
processors=1
logdir=~/sge_logs/gen_exomise/$sig/
mkdir -pv $logdir
out=/dupa-filer/talf/matchingsim/patients/$sig/

#location of files given first, number of files to generate is given as second argument
loc=$1
num=$2
#step 1: make copies of some files (randomly choose the number specified)
mkdir -pv $out
for s in `ls $loc/*.vcf.gz | sort -R | head -n $num`
do
    cp $s $out 
done
#step 2: run our patient generator to make "infect" all of these patients
python /dupa-filer/talf/matchingsim/patients/randompatients/generate_patients.py /dupa-filer/talf/matchingsim/patients/phenotype_annotation.tab /dupa-filer/talf/matchingsim/data/hgmd/hgmd_correct.jv.vcf $out
#step 3: create a script and dispatch exomizer job
mkdir -pv $out/scripts
for file in $out/*.vcf.gz; do
    #create a bash script
    #get only ending to name script
    f=`echo $file | rev | cut -d '/' -f1 | rev | cut -d '.' -f1`
    script="$out/scripts/dispatch_$f.sh"   
    cat > "$script" <<EOF
#!/usr/bin/env bash
#$ -V
#$ -N "$f"
#$ -pe parallel "$processors"
#$ -l h_vmem="$memory"
#$ -e $logdir
#$ -o $logdir

set -eu
set -o pipefail
temp=\$TMPDIR/

#Run exomizer, only if the required file doesn't already exist
gunzip $out/$f.vcf.gz | java -Xms5g -Xmx5g -jar /filer/tools/exomiser/2.0.0/exomiser-2.0.0.jar -D /filer/tools/exomiser/2.0.0/data/ucsc_hg19.ser -I AD -F 1 -W /filer/tools/exomiser/2.0.0/data/rw_string_9_05.gz -X /filer/tools/exomiser/2.0.0/data/rw_string_9_05_id2index.gz --hpo_ids `cat $out/"$f"_hpo.txt` -v $out/$f.vcf --vcf_output 

touch $out/$f.complete
EOF
    #Submit
    #qsub -S /bin/sh "$script"
    gunzip $out/$f.vcf.gz 
    java -Xms5g -Xmx5g -jar /filer/tools/exomiser/2.0.0/exomiser-2.0.0.jar -D /filer/tools/exomiser/2.0.0/data/ucsc_hg19.ser -I AD -F 1 -W /filer/tools/exomiser/2.0.0/data/rw_string_9_05.gz -X /filer/tools/exomiser/2.0.0/data/rw_string_9_05_id2index.gz --hpo_ids `cat $out/"$f"_hpo.txt` -v $out/$f.vcf --vcf_output 


done

