#! /usr/bin/env python

import pyraf
import pyfits
import iraf
import glob
import re

iraf.onedspec(Stdout=1)

# global header

idlfitss = glob.glob('spec1d.*.fits')
idlfitss.sort()

#idlfitss = ['']

#this extracts spectra from the binary fits which the deimos pipeline produces and 
# turns it into the sort of fits files iraf likes. It the resulting spectra are also linearly dispersed

#After running this, run velcor.py for perform heliocentric velocity corrections


# writes a pair of arrays as a fits file using info from the header of the binary fits file
# does this by writing the data to an .asc file, using onedspec.rspectext to create a fits file then using pyfits 
# to update the header
# if I was more clever I would use pyfits to write the fits dirrectly but this way the wcs is sort out for me
def write(name, x, y, header):
  #
#  try:
    l = min(len(x), len(y))
    f = open(name+'.asc', 'w')
    for i in range(l):
    	f.write(repr(x[i] + 0.) + '\t' + repr(y[i] + 0.) + '\n')
    f.close()
    iraf.rspectext(input=name+'.asc', output=name+'.fits', dtype='interp', crval1='INDEF', cdelt1='INDEF')
    #
    fits = pyfits.open(name+'.fits', mode='update')
    fitshead = fits[0].header
    fitshead.set('DATE-OBS', header['DATE-OBS'])
    fitshead.set('TELESCOP', 'Keck II')
    fitshead.set('INSTRUME', 'DEIMOS')
    fitshead.set('OBSERVER', header['OBSERVER'])
    fitshead.set('OBJECT', header['OBJNO'])
    fitshead.set('EQUINOX', header['EQUINOX'])	
    fitshead.set('RA', header['RA_OBJ'])
    fitshead.set('DEC', header['DEC_OBJ'])
    fitshead.set('UTC', header['UTC'])
    fitshead.set('MJD-OBS', header['MJD-OBS'])
    fitshead.set('AIRMASS', header['AIRMASS'])
    fitshead.set('EXPTIME', header['EXPTIME'])
    fitshead.set('SYNOPSIS', header['SYNOPSIS'])
    fitshead.set('CHIPNO', header['CHIPNO'])
    fits.close()
    return name + ' ' + header['RA_OBJ'] + ' ' + header['DEC_OBJ']
    #
	# there is always bad data  
  # except Exception:
  #   iraf.delete(files=name+"*",verify='no')
  #   if verbose: print 'Something wrong writing ' + name
		
def runFitsfitter(verbose=False):
  output = ''
  #
  # go through all the 1d spectra and create fits files
  for idlfits in idlfitss:
  #	
    name = idlfits.replace('spec1d.','').replace('.fits','')
    if verbose: print idlfits + ' -> '+ name
    #
   # if True:
    try:
      f = pyfits.open(idlfits)	
      bdata = f[3].data
      blambda = bdata.field('LAMBDA')[0]
      bspec = bdata.field('SPEC')[0]
      bskyspec = bdata.field('SKYSPEC')[0]
      bivar = bdata.field('IVAR')[0]
      #
      if 'SKYIVAR' in bdata.dtype.names:
        bskyivar = bdata.field('SKYIVAR')[0]
      else:
        bskyivar = None
    #
      rdata = f[4].data
      rlambda = rdata.field('LAMBDA')[0]
      rspec = rdata.field('SPEC')[0]
      rskyspec = rdata.field('SKYSPEC')[0]
      rivar = rdata.field('IVAR')[0]
      if 'SKYIVAR' in rdata.dtype.names:
        rskyivar = rdata.field('SKYIVAR')[0]
      else:
        rskyivar = None
    #
      header = f[3].header
      f.close()
    #
      write(name + '.b', blambda, bspec, header)
      write(name + '.b.ivar', blambda, bivar, header)
      write(name + '.b.sky', blambda, bskyspec, header)
      if bskyivar != None:
        write(name + '.b.skyivar', blambda, bskyivar, header)
      #
      info = write(name + '.r', rlambda, rspec, header)
      output += info + '\n'
      write(name + '.r.ivar', rlambda, rivar, header)
      write(name + '.r.sky', rlambda, rskyspec, header)
      if rskyivar != None:
        write(name + '.r.skyivar', rlambda, rskyivar, header)
    #
  	# there is always bad slits
    except Exception:
      if verbose: print 'Something wrong with ' + idlfits + ' Skipping this file'
  #
  open('index', 'w').write(output)	
