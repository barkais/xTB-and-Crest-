#!/bin/bash

mkdir crest_files
mkdir xyz_originals

read -p 'Constrain atoms: ' atoms 

cp /gpfs0/gaus/projects/crest ..

for file in *.xyz
do
        clean_name=$(basename $file)
	../crest $file --xnam /gpfs0/gaus/users/barkais/home/xtb-6.4.1/bin/xtb --constrain $atoms
	../crest $file --xnam /gpfs0/gaus/users/barkais/home/xtb-6.4.1/bin/xtb --cinp .xcontrol.sample
	mv $file xyz_originals
	cat crest_best.xyz | sed '2d' | sed '1 a \\n' | sed '2d' > "$clean_name"
	mv "$clean_name" crest_files
	rm cre* gfnff_* coord* MRMSD/* struc.xyz wbo

done

