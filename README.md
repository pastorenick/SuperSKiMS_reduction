SuperSKiMS_reduction
====================

Python code to reduce SuperSKiMS data and get the stellar kinematics. 


This code runs the original SKiMS reduction on DEIMOS data, as used by Norris+08, Proctor+09, Foster+09. 
#
It is designed to retrieve the correct slit position from SuperSKiMS masks and to run pPXF on the stellar extracted spectra in 
order to obtain the stellar kinematics at large radii. 
#
The final output is a python object containing all the parameters and spectra of both the DEIMOS mask and the single 
slits. 
#

v1.0 16/12/2014 - code works. 


Notes:
  - the velocities in the output object are already heliocentric corrected
