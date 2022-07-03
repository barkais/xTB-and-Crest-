#!/bin/bash

mkdir crest_files
mkdir xyz_originals
mkdir maximum_dipole
mkdir minimum_dipole

read -p 'Constrain atoms: ' atoms

cp ~/../../../projects/crest .
cp -R -u -p //gpfs0/gaus/xtb-6.4.1/bin/xtb .

for file in *.xyz
do
        clean_name=$(basename $file .xyz) 
        crest $file --xnam /gpfs0/gaus/xtb-6.4.1/bin/xtb --constrain $atoms
        crest $file --xnam /gpfs0/gaus/xtb-6.4.1/bin/xtb --cinp .xcontrol.sample --cluster
        mv $file xyz_originals
	rm anmr_* coord* cregen_* cre_* crest_best.xyz crest_conformers.xyz crest_rotamers.xyz gfnff_* struc.xyz wbo  crest.energies
	rm -r MRMSD/
	mkdir crest_files/"$clean_name"
	mkdir Errors
	mv crest_clustered.xyz crest_files/"$clean_name"
	cd crest_files/"$clean_name"

	# split the structure cluster file to separate xyz files

	declare -i a=2
	declare -i num_of_atoms=$(head crest_clustered.xyz -n 1 | grep -o '[[:digit:]]*')
	declare -i num_of_lines=$(($a+$num_of_atoms))
	split -l "$num_of_lines" crest_clustered.xyz --numeric-suffixes structure
	sed -i '2d' structure*
	sed -i '1 a \\n' structure*
	sed -i '2d' structure*
	find . -type f -not -name "*.*" -exec mv "{}" "{}".xyz \;
	rm crest_clustered.xyz
	cluster_check=`ls -l structure* | wc -l`
        if [ "$cluster_check" -eq 1 ]
	then
		name_mol=`basename "$PWD"`
		mv structure00.xyz "$name_mol".xyz
		cp "$name_mol".xyz ../../maximum_dipole
		cp "$name_mol".xyz ../../minimum_dipole
	else
		rm ../../cluster.order
		mkdir dipole_moments
		for file in *.xyz
		do	
			clean_name_2=$(basename $file .xyz)
			../../xtb $file --dipole > xtb.log
			grep -A 3 dipole xtb.log | tail -1 > dipole
			mv dipole dipole_moments/dipole_"$clean_name_2".csv
			rm wbo xtbrestart xtbtopo.mol charges xtb.log
		done
		cat dipole_moments/* | awk '{ print $5 }' > tot_dip
		cat tot_dip | sort -nk1,1 | tail -1 > max
		cat tot_dip | sort -nk1,1 | head -1 > min
		a=`awk 'FNR==NR{l[$0]=NR; next}; $0 in l{print $0, l[$0], FNR}' max tot_dip | awk '{ print $3 }'`
		b=`awk 'FNR==NR{l[$0]=NR; next}; $0 in l{print $0, l[$0], FNR}' min tot_dip | awk '{ print $3 }'`
		ls -l structure* | awk -v max="$a" -v min="$b" 'BEGIN{min=min;max=max};{ if (NR==min) print $9 ;{ if(NR==max) print $9}}' > keep_list
		ls | egrep -vf keep_list | xargs -n 1 echo rm | grep -v dipole_moments | grep -v keep_list | grep -v tot_dip > rmscript
		sh -x rmscript
		name_min=`head -n 1 keep_list`
		name_max=`tail -n 1 keep_list`
		name_mol=`basename "$PWD"`
		mv "$name_min" "$name_mol"_min.xyz
		mv "$name_max" "$name_mol"_max.xyz
		cp "$name_mol"_min.xyz ../../minimum_dipole
		cp "$name_mol"_max.xyz ../../maximum_dipole
		rm keep_list
		length_dipole=`cat tot_dip | wc -l`
		if ! [ "$cluster_check" -eq "$length_dipole" ]
		then
			cd ..
			cp "$clean_name" ../Errors
			cd "$clean_name"
		fi
	fi
	cd ../..
done
rm xtb crest
