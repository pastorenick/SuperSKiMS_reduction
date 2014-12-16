#!/usr/bin/env python
#Filename: plotFLUX.py

from Nicola import *

#Retrieve dictionary
lib_path = os.path.abspath('/Users/npastorello/Desktop/Galaxies/General_studies/')
sys.path.append(lib_path)
from galaxyParametersDictionary_v5 import *



#Main
def plotFlux(Deimos_mask, inputFile='DATA2.DAT'):
  namegal = Deimos_mask.galaxy
  ra0, dec0, pa0, ba = Deimos_mask.gal_RA, Deimos_mask.gal_Dec, Deimos_mask.gal_PA0, Deimos_mask.gal_ba
  # ra0, dec0 = CentreCoordinates[namegal]
  # dummy = convAngCoord(ra0)
  # ra0 = dummy[4]*15. #gal centre coordinates in deg
  # dummy = convAngCoord(dec0)
  # dec0 = dummy[4] #gal centre coordinates in deg
  # pa0,ba = PA0[namegal], b_a[namegal]
  #
  inputTab = numpy.loadtxt('DATA2.DAT', 
            dtype={'names': ('filename', 'action', 'RA', 'Dec', 'cont1','band','cont2','eqw'),
                                 'formats': ('S50', 'S3', 'f10', 'f10', 'f10', 'f10', 'f10', 'f10')})
  #
  #Measuring distance
  RArel, Decrel = (inputTab['RA']-ra0)*3600., (inputTab['Dec']-dec0)*3600. #in arcsec
  #
  Xrot = RArel*numpy.cos(numpy.pi*(pa0)/180.) - Decrel*numpy.sin(numpy.pi*(pa0)/180.) 
  Yrot = RArel*numpy.sin(numpy.pi*(pa0)/180.) + Decrel*numpy.cos(numpy.pi*(pa0)/180.)
  #
  distance = numpy.sqrt(((Xrot/numpy.sqrt(ba))**2)+((Yrot*numpy.sqrt(ba))**2))
  #
  plt.ion()
  fig = figure()
  clf()
  ax = subplot(111)
  p = plot(distance, inputTab['eqw'], 'r.')
  ylim([0,max(inputTab['eqw'])+0.2])
  xlim([0,max(distance)+10])
  ylabel('Sky index')
  xlabel('Radius (arcsec)')
  grid()
  #
  for ii in numpy.arange(len(distance)):
    annotate(inputTab['filename'][ii], xy=(distance[ii], inputTab['eqw'][ii]))
  #
  show()
  raw_input('Select at least 3 good sky spectra, change their label in DATA2.DAT from "TAR" to "SKY" and press a key to continue')
  