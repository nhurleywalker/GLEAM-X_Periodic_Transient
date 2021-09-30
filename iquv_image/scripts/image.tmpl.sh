#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=38
#SBATCH --time=12:00:00
#SBATCH --cluster=garrawarla
#SBATCH --account=mwasci
#SBATCH -p workq
#SBATCH --mem=340G

module load singularity
export SINGULARITY_BINDPATH=/astro/mwasci/zxiang

obsid=OBSID

basescale=0.06

chan=`singularity exec /astro/mwasci/tgalvin/gleamx_testing_small.img pyhead.py -p CENTCHAN ${obsid}.metafits | awk '{print $3}'`

scale=`echo "$basescale / $chan" | bc -l` # At least 4 pix per synth beam for each channel

# intervals=$(singularity exec /astro/mwasci/tgalvin/gleamx_testing_small.img ../../ms_stats.py ${obsid}.ms)

singularity exec /astro/mwasci/tgalvin/gleamx_testing_small.img wsclean \
		-name ${obsid} \
		-abs-mem 300 \
		-size 1600 1600 \
		-mgain 1 \
		-nmiter 0 \
		-niter 50000 \
		-weight briggs 0.5 \
		-pol xx,yy,xy,yx \
		-join-polarizations \
		-join-channels \
		-channels-out 96 \
		-auto-mask 3 \
		-auto-threshold 1 \
		-scale 12asec \
		-no-dirty \
		${obsid}.ms | tee wsclean.log

# singularity exec /astro/mwasci/tgalvin/gleamx_testing_small.img wsclean \
# 		-name ${obsid}-MFS \
# 		-abs-mem 320 \
# 		-size 250 250 \
# 		-mgain 1 \
# 		-nmiter 0 \
# 		-weight natural \
# 		-pol xx,yy,xy,yx \
# 		-join-polarizations \
# 		-auto-mask 3 \
# 		-auto-threshold 1 \
# 		-scale ${scale:0:8} \
# 		-no-dirty \
# 		${obsid}.ms | tee wsclean.log
