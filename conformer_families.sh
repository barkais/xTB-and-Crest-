#!/bin/bash

declare -i a=2
declare -i num_of_atoms=$(head $1 -n 1 | grep -o '[[:digit:]]*')
declare -i num_of_lines=$(($a+$num_of_atoms))

split -l "$num_of_lines" $1 --numeric-suffixes structure
sed -i '2d' structure* 
sed -i '1 a \\n' structure*
sed -i '2d' structure*
find . -type f -not -name "*.*" -exec mv "{}" "{}".xyz \;

for file in structure*
do
	clean_name=$(basename $file .xyz)
	mkdir "$clean_name"
	mv "$clean_name".xyz "$clean_name"/basic.xyz

done
