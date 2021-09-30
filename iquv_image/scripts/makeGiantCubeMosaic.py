#!/usr/bin/env python
"""
makeFitsCube.py

This script is intended for merging channel maps into a FITS cube.

Written by Sarrvesh S. Sridhar

Modified by Xiang Zhang to work on TB scale cubes

To do list:
* Change all pyFits dependencies to AstroPy
"""
import optparse
import glob
import os
try:
    import numpy as np
except ImportError:
    raise Exception('Unable to import Numpy')
try:
    from astropy.io import fits
except ImportError:
    raise Exception('Unable to import pyFits.')

version_string = 'v1.0, 8 June 2015\nWritten by Sarrvesh S. Sridhar'
print('makeFitsCube.py', version_string)
print('')

def getValidFitsList(fileList):
    """
    Extracts fits files from a list of files using its magic number.
    """
    validFitsList = []
    for name in fileList:
        #print "* INFO: Validating {0}".format(name)
        if 'FITS' in os.popen('file {0}'.format(name)).read():
            validFitsList.append(name)
    return validFitsList

def checkFitsShape(fitsList):
    """
    Checks if the list of fits files all have the same shape.
    If True, return the shape and memory in bytes of a single fits file
    If False, raises an exception causing the execution to terminate
    """
    for i, name in enumerate(fitsList):
        if i == 0:
            templateShape = fits.open(name, readonly=True)[0].data.shape
        elif templateShape != fits.open(name, readonly=True)[0].data.shape:
            raise Exception('Fits file {0} has an incompatible shape'.format(name))
    return templateShape

# def concatenateWithMemMap(validFitsList, shape, memMapName, FLAG):
#     """
#     Concatenate a given list of fits files into a single cube using memory map.
#     Return the concatenated array and frequency list
#     """
#     concatCube = np.memmap(memMapName, dtype='float32', mode='w+',\
#                            shape=(len(validFitsList), shape[-2], shape[-1]))
#     freqList = []
#     for i, name in enumerate(validFitsList):
#         print "* INFO: Getting data for {0}".format(name)
#         if len(shape) == 2:
#             tempData = np.squeeze(fits.open(name, readonly=True)[0].data)
#         else:
#             tempData = np.squeeze(fits.open(name, readonly=True)[0].data[0])
#         tempHead = fits.open(name, readonly=True)[0].header
#         freqList.append(tempHead[FLAG])
#         concatCube[i] = np.copy(tempData)
#     return concatCube, np.array(freqList)

def makefreqlist(validFitsList, shape, FLAG):
    """
    Return frequency list
    """
    freqList = []
    for i, name in enumerate(validFitsList):
        tempHead = fits.open(name, readonly=True)[0].header
        freqList.append(tempHead[FLAG])
    return np.array(freqList)

def main(options):
    """
    Main function
    """
    # Check user input
    if options.inp == '':
        raise Exception('An input glob string must be specified.')
    if options.out == '':
        raise Exception('An output filename must be specified.')

    # Iterate over Q and U
    for pol in ['U']:
        memMapName = 'memMap_cube_{0}'.format(pol)
        # Get the list of FITS files; XZ: modify line 86 according to the name of corrected fine channel images
        # fileList = sorted(glob.glob(options.inp+'/eor*{0}.fits'.format(pol))) # For eor data
        fileList = sorted(glob.glob(options.inp+'/1*{0}.fits'.format(pol)))

        # print fileList
        validFitsList = getValidFitsList(fileList)
        print('INFO: Identified {0} fits files from {1} files selected by input string'.\
              format(len(validFitsList), len(fileList)))

        # Proceed with the execution if we have non-zero FITS files
        if len(validFitsList) == 0:
            raise Exception('No valid fits files were selected by the glob string')

        # Check if the list of supplied fits files have the same shape
        shape = checkFitsShape(validFitsList)
        print('INFO: All fits files have shape {0}'.format(shape))
        if len(shape) not in [2, 3, 4]:
            raise Exception('Fits files have unknown shape')

        # Merge the cubes
        if options.mwafiles:
            FLAG = 'FREQ'
        elif options.restfrq:
            FLAG = 'RESTFREQ'
        else:
            FLAG = 'CRVAL3'
        # finalCube, freqList = concatenateWithMemMap(validFitsList, shape, memMapName, FLAG)
        freqList = makefreqlist(validFitsList, shape, FLAG)


        if options.mwafiles: freqList *= 1.e6

        # Create outdir if it doesn't already exist
        if not os.path.exists( options.out ):
            os.system("mkdir {0}".format(options.out))

        # Write the frequency list to disk
        f = open(options.out+'/freq-{0}-mosaic.txt'.format(pol), "w")
        for line in freqList:
            f.write(str(line)+'\n')
        f.close()
        print("INFO: Frequency range is {0} - {1} MHz".format(freqList[0]/1.e6, freqList[-1]/1.e6))

        # Get a template header from a fits file
        header = fits.open(validFitsList[0], readonly=True)[0].header
        # print 'INFO: Writing the concatenated fits file to {0}'.format(options.out)

        # create dummy axis 3 header entries if there are only 2 axes

        # XZ: change this part so we create rotated cubes by default

        newHeader = fits.open(validFitsList[0], readonly=True)[0].header

        # while len(newHeader) < (60 * 4 - 1):
        #     newHeader.append()  # Adds a blank card to the end
        
        newHeader['NAXIS']=3
        newHeader['NAXIS3'] = len(validFitsList)
        newHeader['CTYPE3']='FREQ'
        newHeader['CRPIX3']=1
        newHeader['CRVAL3']=freqList[0]
        newHeader['CDELT3']=freqList[1]-freqList[0]

        # newHeader['NAXIS'] = 3

        # newHeader['NAXIS1'] = len(validFitsList)
        # newHeader['NAXIS3'] = header['NAXIS1']

        # newHeader['CDELT3'] = header['CD1_1']
        # newHeader['CRVAL3'] = header['CRVAL1']
        # newHeader['CRPIX3'] = header['CRPIX1']
        # newHeader['CTYPE3'] = header['CTYPE1']

        # newHeader['CTYPE1']='FREQ'
        # newHeader['CRPIX1']=1
        # newHeader['CRVAL1']=freqList[0]
        # newHeader['CDELT1']=freqList[1]-freqList[0]

        # newHeader['NAXIS1'] = 400
        # newHeader['NAXIS3'] = 646

        # newHeader.verify('fix')

        # print newHeader

        newHeader.tofile(options.out+'/cube-{0}-mosaic.fits'.format(pol))

        # print 'newheader to file'

        giantfitscube = options.out+'/cube-{0}-mosaic.fits'.format(pol)

        # XZ: fix the header

        hdu = fits.open(giantfitscube, mode='update')
        hdu[0].verify('fix')
        # # test data writing speed
        # hdu[0].header['NAXIS1'] = 400
        # hdu[0].header['NAXIS3'] = 646

        hdu.close()

        # print header

        with open(giantfitscube, 'rb+') as fobj:
            fobj.seek(len(newHeader.tostring()) + (newHeader['NAXIS1'] * newHeader['NAXIS2'] 
                * newHeader['NAXIS3'] * np.abs(header['BITPIX']//8)) - 1)
            fobj.write(b'\0')

        # print 'header made fits'

        # from tempfile import mkdtemp
        # import os.path as path
        # memmapname = path.join(mkdtemp(), 'memmap.dat')

        # temp_array = np.memmap(memMapName, dtype='float32', mode='w+', shape=(newHeader['NAXIS3'], newHeader['NAXIS2'], newHeader['NAXIS1']))

        # XZ: now write the giant cube
        for i in range(0,newHeader['NAXIS1']):
            print("* INFO: Getting data for {0}".format(validFitsList[i]))
            hdu = fits.open(giantfitscube, mode='update')
            hdu_chan = fits.open(validFitsList[i])
            # temp_array[:,:,i] = hdu_chan[0].data.T
            # XZ: this is a rotated cube. NAXIS1 is freq, NAXIS2 is DEC, NAXIS3 is RA.
            hdu[0].data[i,:,:] = hdu_chan[0].data
            # hdu[0].data[:,:,i] = hdu_chan[0].data
            hdu.close()
            hdu_chan.close()

        # rotate the cube?
        # temp_rotated = np.swapaxes(temp_array, 0, 2)

        # hdu = fits.open(giantfitscube, mode='update')
        # hdu[0].data = temp_array
        # hdu.close()


        # Rotate the cube if -s is used
        # XZ: this part is no longer used, since we are not able to rotate a TB size cube.
        # if options.swapaxis:
        #    rotCube = np.memmap(memMapName+"_rot", dtype='float32', mode='w+',\
        #                        shape=(shape[-2], shape[-1], len(validFitsList)))
        #    rotCube = np.swapaxes(finalCube, 0, 2)
        #    print "INFO: Rotated cube has shape ", rotCube.shape
        #    newHeader = fits.open(validFitsList[0], readonly=True)[0].header
        #    newHeader['CDELT1'] = header['CDELT3']
        #    newHeader['CRVAL1'] = header['CRVAL3']
        #    newHeader['CRPIX1'] = header['CRPIX3']
        #    newHeader['CTYPE1'] = header['CTYPE3']

        #    newHeader['CDELT3'] = header['CD1_1']
        #    newHeader['CRVAL3'] = header['CRVAL1']
        #    newHeader['CRPIX3'] = header['CRPIX1']
        #    newHeader['CTYPE3'] = header['CTYPE1']

        #    hdu = fits.PrimaryHDU(data=rotCube, header=newHeader)
        # else:
        #    hdu = fits.PrimaryHDU(data=finalCube, header=header)

        # hdu = fits.PrimaryHDU(data=finalCube, header=newHeader)
        # hdu.writeto(options.out+'/cube-{0}-mosaic.fits'.format(pol), overwrite=True)
        # os.remove(memMapName)
        # os.remove(memMapName+"_rot")
        #os.system('cp -r {0}/Cubes {1}/.'.format(options.out, options.inp))
        #os.remove(options.out)

if __name__ == '__main__':
    opt = optparse.OptionParser()
    opt.add_option('-i', '--inp', help='Glob selection string for input files '+
                   '[no default]', default='')
    opt.add_option('-o', '--out', help='Destination for output cube(s) '+
                   '[no default]', default='')
    opt.add_option('-f', '--freq',
                   help='Filename to write the frequency list [default: frequency.txt]',
                   default='frequency.txt')
    opt.add_option('-r', '--restfrq', help='Frequency is stored in RESTFRQ '+
                   'instead of CRVAL3 [default: False]', default=False,
                   action='store_true')
    opt.add_option('-m', '--mwafiles', help='FITS files are MWA based (they '+
                   'have FREQ header keyword in MHz); this overrides --restfrq '+
                   '[default False]', default=False, action='store_true')
    # opt.add_option('-s', '--swapaxis', help='Make frequency axis as the first '+
    #                'axis [default: False]', default=False, action='store_true')
    inOpts, arguments = opt.parse_args()
    main(inOpts)

