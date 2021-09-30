#!/bin/bash -l

# Correct the xy phase of U/V images

for obsid in 1*
do
	cat xycorr.tmpl.sh | sed "s;OBSID;${obsid};g" > ${obsid}.xycorr.sh
	mv ${obsid}.xycorr.sh ${obsid}/
	cd ${obsid}/
	# rename "U.fits" "U_b.fits" *U.fits
	# rename "V.fits" "V_b.fits" *V.fits
	sbatch ${obsid}.xycorr.sh
	cd ../
	
done
