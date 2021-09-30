#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --cluster=garrawarla
#SBATCH --account=mwasci
#SBATCH -p workq

obsid=OBSID
MFS_freq=184.955

module load singularity
export SINGULARITY_BINDPATH=/astro/mwasci/zxiang

rename "${obsid}-" "${obsid}_" *.fits

# intervals=$(singularity exec /astro/mwasci/tgalvin/gleamx_testing_small.img ../../ms_stats.py ${obsid}.ms)

# for (( ts=0; ts<intervals; ts++ ))
# do

for chan in {0..11}
do
	modu=$((3%16))
	finechan=$(printf '%04d' $chan)
	freq=$(python -c "print $chan*2.56+170.875")	# For GLEAM phase I band 169
	# timestep=$(printf 't%04d' $ts)
	
	if (($modu==0)) || (($modu==15))
	then
		# rm ${obsid}_${finechan}*
		continue
	else
		# 
		singularity exec /astro/mwasci/tgalvin/gleamx_testing_small.img lookup_jones.py ${obsid} _${finechan}-XX-image.fits ${obsid}-${finechan}-beam- .fits -f $freq -vv --beam_path /astro/mwasci/zxiang/database/gleam_jones.hdf5
		singularity exec /astro/mwasci/tgalvin/gleamx_testing_small.img pbcorrect ${obsid}_${finechan} image.fits ${obsid}-${finechan}-beam ${obsid}-${finechan}-pbcorr-image
		rm ${obsid}-${finechan}-beam*.fits
		# rm ${obsid}_${finechan}*image.fits
	fi

	# if [[ $(jobs -r -p | wc -l) -gt 38 ]]; then
	# 	# wait only for first job
	# 	wait -n
	# fi
	
done
	# wait

# done

# Adding pbcorrect for MFS images

singularity exec /astro/mwasci/tgalvin/gleamx_testing_small.img lookup_jones.py ${obsid} _MFS-XX-image.fits ${obsid}-MFS-beam- .fits -f $MFS_freq -vv --beam_path /astro/mwasci/zxiang/database/gleam_jones.hdf5

singularity exec /astro/mwasci/tgalvin/gleamx_testing_small.img pbcorrect ${obsid}_MFS image.fits ${obsid}-MFS-beam ${obsid}-MFS-pbcorr-image

# rm *beam*.fits

# rm *image.fits
