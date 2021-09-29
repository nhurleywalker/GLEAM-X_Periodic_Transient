#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --cluster=garrawarla
#SBATCH --account=mwasci
#SBATCH -p workq

rename "OBSID_" "OBSID-" *.fits
rename "XY-image.fits" "XY-image_b.fits" *XY-image.fits
rename "XYi-image.fits" "XYi-image_b.fits" *XYi-image.fits

module load singularity
export SINGULARITY_BINDPATH=/astro/mwasci/zxiang

singularity exec /astro/mwasci/tgalvin/gleamx_testing_small.img ../xyphase_corr.py OBSID

# rm *_b.fits

