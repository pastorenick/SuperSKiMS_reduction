#!/usr/bin/env python
# Filename: SKiMS_red.py
# Run the SKiMS procedure on the fits files obtained from the deep2 pipeline
# 
# Written by Nicola Pastorello
# 5/12/14

# v 1.1 (15/12/14) - create the objects 'mask' and 'slit'
# v 1.2 (????????) - save output in a pickle dictionary
# 
#
#
# The code requests the following files:
#   - Nicola.py
#   - flux_level.py
#   - flux_level2.py
#   - plotFLUX.py
#   - fitsfitter_MOD.py

# It needs SKiMS_red__def__.py, fitsfitter_MOD.py, flux_level.py, flux_level2.py, plotFLUX.py

from SKiMS_red__def__ import *
__builtins__.c = 299792.458

SKiMSreduction = True
check_spectra_by_eye = True
drawingDistro = True
pPXF_fitting = True
saveKinPlot = True
__builtins__.verbose = True

masktype = 'Major'
galname = 'NGC1023'

# RETRIEVE MASK ALIGNMENT
listPlanFiles = glob.glob('*plan*')
if len(listPlanFiles) == 1:
  fileIn = open(listPlanFiles[0], 'r')
  flag = True
  while flag:
    tmp = fileIn.readline()
    if tmp[:11] == 'SCIENCENAME': 
      nameFile = tmp[13:-1]
      flag = False
    if tmp == '':
      flag = False
  fileIn.close()
else:
  print "ERROR"
  raw_input()

#### RUNNING SKiMS REDUCTION

Deimos_mask = mask()
Deimos_mask.assignGalaxy(galname)
Deimos_mask.assignMaskPars(nameFile)

if SKiMSreduction:
  runSKiMSreduction(Deimos_mask)
 

#### PLOTTING SUBTRACTED SPECTRA

#if check_spectra_by_eye:
runCheckSpectra(Deimos_mask, galname=galname)
  
#if drawingDistro:
drawDistro(Deimos_mask, galname=galname, outputname='slitPosition_Major.pdf')
 

### CHECK FROM FITS HEADER THE LENGTH AND THE INCLINATION OF THE SLITS
Deimos_mask.assignExtraParametersSlits(nameFile)

# dicOut = {'name':numpy.array(listNAME), 'RA':numpy.array(listRA), 'DEC':numpy.array(listDEC), 'sel':numpy.array(listCheck)}
# fileOut = open(outputname, 'wb')
# pickle.dump(dicOut, fileOut)
# fileOut.close()


fileOut = open('Major_axis_mask_tmp.dat', 'wb')
pickle.dump(Deimos_mask, fileOut)
fileOut.close()

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fileIn = open('Major_axis_mask_tmp.dat', 'rb')
Deimos_mask = pickle.load(fileIn)
fileIn.close()
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#### RUN pPXF ON SUBTRACTED SPECTRA
if pPXF_fitting:
  listFail = run_pPXF(Deimos_mask)
  if len(listFail) > 0: 
    print "Failed with the following slits:"
    print '\n'.join(listFail)


if saveKinPlot:
  run_saveKinPlot(Deimos_mask, outputFile='kinplot_'+masktype+'.pdf')


#### SAVING OUTPUT Mask Object ####
fileOut = open('Major_axis_mask.dat', 'wb')
pickle.dump(Deimos_mask, fileOut)
fileOut.close()








