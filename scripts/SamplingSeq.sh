#!/bin/bash

PathIn=$1
PathOut=$2

for file in $PathIn/*.fa; do
	if [ -f "$file" ]; then
		echo $file
		echo "------"
		seqkit seq -m 1300 -M 1750 $file | seqkit sample -n 500 -s 991 > $PathOut/`echo $file | cut -d/ -f3`
	fi
done
