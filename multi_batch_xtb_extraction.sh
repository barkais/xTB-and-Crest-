#!/bin/bash


find . -name '*_xtb' -exec rm -r {} +

for f in *
do
    	if [ -d "$f" ]; then
		cd $f
		cp -R -u -p //gpfs0/gaus/xtb-6.4.1/bin/xtb .
		cp -R -u -p //gpfs0/gaus/projects/obabel/build/bin/obabel .

		mkdir ../"$f"_xtb
		mkdir ../"$f"_xtb/extracted_info
		mkdir ../"$f"_xtb/xtb_optimized_xyz
		for file in *.xyz
			do 
				clean_name=$(basename $file .xyz)
				./obabel $file -Omol2 --gen3d		
				xtb $file --ohess --dipole > xtb.log
				
				mkdir ../"$f"_xtb/extracted_info/"$clean_name"
				grep -A 3 dipole xtb.log | tail -1 > dipole
				cat mol2 | grep @\<TRIPOS\>ATOM -A 100 | grep @\<TRIPOS\>BOND -B 100 | head -n -1 | tail -n +2 | awk '{ print $6 }' > atypes_"$clean_name".csv
				mv atypes_"$clean_name".csv ../"$f"_xtb/extracted_info/"$clean_name"
				rm mol2
				mv dipole ../"$f"_xtb/extracted_info/"$clean_name"/dipole_"$clean_name".csv
				mv g98.out ../"$f"_xtb/extracted_info/"$clean_name"/movec_"$clean_name".csv
				mv hessian ../"$f"_xtb/extracted_info/"$clean_name"/hessian_"$clean_name".csv
				mv wbo ../"$f"_xtb/extracted_info/"$clean_name"/wbo_"$clean_name".csv
				sed -e '2d' xtbopt.xyz | sed '2i\\' > temp.xyz
				sed 's/ \{1,\}/,/g' temp.xyz > temp2.xyz
				cp temp2.xyz ../"$f"_xtb/extracted_info/"$clean_name"/xyz_"$clean_name".csv
				cp temp2.xyz ../"$f"_xtb/xtb_optimized_xyz/"$clean_name".xyz
				mv vibspectrum ../"$f"_xtb/extracted_info/"$clean_name"/vibspectrum_"$clean_name".csv
				mv xtb.log ../"$f"_xtb/extracted_info/"$clean_name"
			done
			

		rm temp.xyz
		rm temp2.xyz
		rm obabel
		rm xtb
		rm xtb*
		rm charges
			
		cd ..
	fi 
done

mkdir Finished_multi_batch_xtb
mv *_xtb Finished_multi_batch_xtb
	
	
