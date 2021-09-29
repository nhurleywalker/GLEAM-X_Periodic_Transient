#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --cluster=garrawarla
#SBATCH --account=mwasci
#SBATCH -p workq

module load singularity
export SINGULARITY_BINDPATH=/astro/mwasci/zxiang

# WORKDIR='/scratch1/zha292/GLEAMX/galactic_centre/gleam_0.0_-27_169/corrected_chan/'

# if [ ! -d "./Cubes" ]; then
#     echo "* INFO: Making directory: ./Cubes"
#     mkdir ./Cubes/
# fi

# cd $LOCALDIR
# mkdir ./Cubes

# python -u ${WORKDIR}/makeFitsCubeMosaic.py -m -s -i  ${WORKDIR}/Mosaics/ -o ./Cubes #./Cubes_mosaic

# XZ: core dump issue. Try not using localdir?

# python -u makeFitsCubeMosaic.py -m -s -i  ./Mosaics/ -o ./Cubes

singularity exec /astro/mwasci/tgalvin/gleamx_testing_small.img python -u makeGiantCubeMosaic.py -i  ./chan_imgs/ -o ./Cubes

# echo "Copying cubes..."
# cp -r Cubes/* ${WORKDIR}/Cubes/.
echo "** DONE **"
