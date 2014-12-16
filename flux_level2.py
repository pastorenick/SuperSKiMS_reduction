#!/usr/bin/env python

#Filename: flux_level2.py

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

def run_flux_level2():
  #Definizione limiti intervalli
  WL1, WL2 = 8526., 8534.
  WL3, WL4 = 8605., 8695.5
  WL5, WL6 = 8795., 8820.
  
  sfile = asciidata.open('./temp.lis')
  outputMat = []
  for i in numpy.arange(len(sfile[0])):
    try:
      nfile=sfile[0][i]
      tmp=numpy.loadtxt('./'+nfile[:-4]+'asc')
      WV=tmp[:,0] #wavelenght
      flux=tmp[:,1]
      #
      cont1=countF(flux,WV,WL1,WL2)
      cont2=countF(flux,WV,WL5,WL6)
      #
      cont1=cont1/(WL2-WL1)
      cont2=cont2/(WL6-WL5)
      #
      m=(cont2-cont1)/((WL5+WL6)/2-(WL1+WL2)/2)
      c=cont1-m*(WL1+WL2)/2
      #
      band=countF(flux,WV,WL3,WL4)
      band=band-(((m/2)*(WL4**2)+c*WL4)-((m/2)*(WL3**2)+c*WL3))
      band=band/(WL4-WL3)
      #
      eqw=2*band/(cont1+cont2)
      #
      #Salvataggio
      outputMat.append(nfile+'  '+str(cont1)+'  '+str(band)+'   '+str(cont2)+'  '+str(eqw))
    except:
      print nfile
      outputMat.append(nfile+'  '+str(nan)+'  '+str(nan)+'   '+str(nan)+'  '+str(nan))
  
  
  numpy.savetxt('./temp_levels.lis',outputMat,fmt='%s')
  
