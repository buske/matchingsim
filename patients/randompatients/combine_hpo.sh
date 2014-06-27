#!/usr/bin/env bash
usage="usage: $0 dir"

if [ $# -lt 1 ]; then
    echo "$usage"
    exit 1
fi

loc=$1

rm -v -f $loc/combined_hpo.txt
for f in $loc/*_hpo.txt; do
    s=`basename $f _hpo.txt`
    echo -e $s$'\t'`cat $f` >> $loc/combined_hpo.txt
done


