#!/bin/bash -l

for obsid in 1*
do
	cat pbcorr.tmpl.sh | sed "s;OBSID;${obsid};g" > ${obsid}.pbcorr.sh
	mv ${obsid}.pbcorr.sh ${obsid}/
	cd ${obsid}/
	# rename "${obsid}-" "${obsid}_" *.fits
	sbatch ${obsid}.pbcorr.sh
	cd ../
	
done
