import numpy as np

### PARAMETERS #################################################################################################

C_Ha = 15.539801701451697

K = 7.9e-42   #SFR factor in [solar mass * s / yr / erg]
d = 9.5   #distance of the galaxy in Mpc
dist = d * 3e24   #convert distance from Mpc to cm

sigma_sky = 3.4889957167289254e-13   #sky flux from our computataion

size = 10.2 * 60   #size of the target in arcsec
px = 0.44   #arcsec per pixel
radius = size / px / 2   #whole galaxy radius in pixels
area = np.pi * radius**2   #area of the whole galaxy in pixels

radius_sky = 168.1818181818182   #galaxy small radius in pixels
area_sky = np.pi * radius_sky**2   #small area of the galaxy in pixels

################################################################################################################

SFR = 0.5   #from literature (aka Fossati)

flux_ph = SFR / (4 * np.pi * dist**2) / K   #flux in physycal units
f = flux_ph * 10**C_Ha   #flux in electrons/s

f_sky = (sigma_sky * 10**C_Ha * 1200)**2 / 1200 / area_sky * area  #sky flux in electrons/s (normalized on #pix)

Texp = (5 * np.sqrt(f + f_sky) / f)**2   #exposure time in s

print(Texp / 60, 'min')