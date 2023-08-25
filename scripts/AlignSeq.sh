#!/bin/bash

PathIn=$1
PathOut=$2

for file in $PathIn/*.fa; do
	if [ -f "$file" ]; then
		echo $file
		echo "------"
		~/Software/clustalo-1.2.4-Ubuntu-x86_64 -i $file --threads 20 -o $PathOut/`echo $file | cut -d/ -f3`
	fi
done
