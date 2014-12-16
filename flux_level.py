#!/usr/bin/env python

#Filename: flux_level.py

from Nicola import *

#Copy of the flux_level.for wrote in Fortran by Caroline Foster


#FUNCTIONS
def countF(flux, wl, ll, hl):
  #finding position
  for ii in range(0,len(flux)-1):
    if ((wl[ii] <= ll) & (wl[ii+1] >= ll)):
      w1 = ii
    if ((wl[ii] <= hl) & (wl[ii+1] >= ll)):
      w2 = ii
  edge1 = flux[w1]*(wl[w1+1]-ll)/(wl[w1+1]-wl[w1])
  edge2 = flux[w2]*(hl-wl[w2])/(wl[w2+1]-wl[w2])    
  #
  mid = 0
  for ii in range(w1+1,w2):
    mid = mid+flux[ii]
  count = (mid+edge1+edge2)
  return count

#Main

def fluxLevel(Deimos_mask, inputFile='DATA.DAT', outputFile='DATA2.DAT', setLines='standard'):
  if setLines == 'standard':
  #Definition interval limits
    WL1, WL2 = 8478., 8489.
    WL3, WL4 = 8605., 8712.
    WL5, WL6 = 8813., 8822.
#
  inputTab= numpy.loadtxt(inputFile, dtype={'names': ('filename', 'action', 'RA', 'Dec'),
                                 'formats': ('S50', 'S3', 'f10', 'f10')})
#
  listDname = inputTab['filename']
  action = inputTab['action']
  listRA = inputTab['RA']
  listDec = inputTab['Dec']
#
#Definizione matrice per salvataggio risultati
  resMat = numpy.empty((len(listDname),4))
#
  listFails = []
  for ii in numpy.arange(len(listDname)):
    try:
      #Opening spectra
      tmp = numpy.loadtxt(listDname[ii][:-4]+'asc')
      WV, flux = tmp[:,0], tmp[:,1]
      #
      cont1 = countF(flux,WV,WL1,WL2)
      cont2 = countF(flux,WV,WL5,WL6)
      #
      cont1 = cont1/(WL2-WL1)
      cont2 = cont2/(WL6-WL5)
      #
      m = (cont2-cont1)/((WL5+WL6)/2-(WL1+WL2)/2)
      c = cont1-m*(WL1+WL2)/2
      #
      band = countF(flux,WV,WL3,WL4)
      band = band-(((m/2)*(WL4**2)+c*WL4)-((m/2)*(WL3**2)+c*WL3))
      band = band/(WL4-WL3)
      #
      eqw = 2*band/(cont1+cont2)
      #
      resMat[ii,:] = (cont1,band,cont2,eqw)
      #
      # Assigning values to objects 'slit'
      #
      for jj in Deimos_mask.listSlits:
        if listDname[ii] == jj.name:
          jj.skyindex = eqw
    except:
    #
      listFails.append(listDname[ii][:-4]+'asc')
      resMat[ii,:] = (nan,nan,nan,nan)
#
  outputMat = transpose(numpy.vstack((listDname,action,listRA,listDec,transpose(resMat))))
#
  numpy.savetxt(outputFile,outputMat,fmt='%s',
      header='file\taction\tRAdeg\tDecdeg\tcont1\tband\tcont2\teqw',newline='\n',delimiter='\t')
  #
  if len(listFails) > 0:
    print "Failed with the following:"
    print '\n'.join(listFails)
  #
  return "DONE"
