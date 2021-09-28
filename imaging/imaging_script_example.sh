# This script shows how WSClean (Offringa et al. 2014) was run to produce deep images, then individual timestep images at 0.5s time resolution, for each observation
# If you have access to the MWA primary beam, the line
#        -apply-primary-beam \
# can be added to the WSClean commands
# Alternatively, this can be applied afterwards via "beam" and "applybeam", found in https://github.com/ICRAR/mwa-reduce

msigma=3
tsigma=1
imsize=8000
robust=0.5
basescale=0.6
# Pixel scale depends on the central frequency of the observation, so set accordingly
# E.g. a 154 MHz observation is at central channel 121
# This observation can be read from ASVO databse
chan=121
 # At least 4 pix per synth beam for each channel
scale=$(echo "$basescale / $chan" | bc -l)

# Create a template image for the mask
wsclean \
        -abs-mem ${GXMEMORY} \
        -mgain 1.0 \
        -nmiter 1 \
        -niter 0 \
        -name ${obsnum}_template \
        -size ${imsize} ${imsize} \
        -scale ${scale:0:8} \
        -pol XX \
        -channel-range 4 5 \
        -interval 4 5 \
        -no-update-model-required \
        "${obsnum}.ms"

# Make a mask around the transient
python ../make_mask.py ${obsnum}_template.fits

# Create deep images
wsclean \
        -abs-mem ${GXMEMORY} \
        -fits-mask ${obsnum}_mask.fits \
        -multiscale -mgain 0.85 -multiscale-gain 0.15
        -nmiter 5 \
        -niter 10000000 \
        -auto-mask ${msigma} \
        -auto-threshold ${tsigma} \
        -name ${obsnum}_deep \
        -size ${imsize} ${imsize} \
        -scale ${scale:0:8} \
        -weight briggs ${robust} \
        -pol I \
        -join-channels \
        -channels-out 4 \
        -save-source-list \
        -fit-spectral-pol 2 \
        "${obsnum}.ms"

# Create timestep images with deep model subtracted
# $nscans should be set depending on the length of the observation
wsclean \
        -subtract-model \
        -intervals-out ${nscans} \
        -mgain 1 \
        -nmiter 0 \
        -channels-out 96 \
        -auto-mask ${msigma} \
        -auto-threshold ${tsigma} \
        -name ${obsnum}_time \
        -size ${imsize} ${imsize} \
        -scale ${scale:0:8} \
        -weight natural \
        -pol I,Q,U,V \
        -no-update-model-required \
        ${obsnum}.ms
