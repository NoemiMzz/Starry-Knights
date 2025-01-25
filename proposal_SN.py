import numpy as np

#%%
### FUNCTIONS ##################################################################################################

def m_to_flux(mag, lambda_c):
    f0 = F0 * 10**(-23)   #reference flux in erg/s/cm2/Hz
    f = f0 * 10**(-mag/2.5)   #flux in erg/s/cm2/Hz
    return f * 3*10**(18) / lambda_c**2


#%%
### PARAMETERS Ha ##############################################################################################

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


### Ha EXPOSURE TIME ###########################################################################################

SFR = 0.5   #from literature (aka Fossati)

flux_ph = SFR / (4 * np.pi * dist**2) / K   #flux in physycal units
f = flux_ph * 10**C_Ha   #flux in electrons/s

f_sky = (sigma_sky * 10**C_Ha * 1200)**2 / 1200 / area_sky * area  #sky flux in electrons/s (normalized on #pix)

Texp_Ha = (5 * np.sqrt(f + f_sky) / f)**2   #exposure time in s

print("\nHalpha\n", Texp_Ha / 60, 'min')


#%%
### PARAMETERS R ###############################################################################################

F0 = 3631   #reference flux in Junsky

lambda_r = 6165.0   #central wavelenght of the filter in A of SDSS
mag_r = 10.66   #magnetude from NED (SDSS)

C_r = 17.125059673553615

areaSDSS_as = 346 * 55   #SDSS aperture in arcsec
areaSDSS = areaSDSS_as / px**2   #SDSS aperture in pixels

### R EXPOSURE TIME ############################################################################################

flux_r = m_to_flux(mag_r, lambda_r)   #flux in physycal units
fr = flux_r * 10**C_r   #flux in electrons/s

f_sky_r = (sigma_sky * 10**C_r * 300)**2 / 300 / area_sky * areaSDSS  #sky flux in electrons/s (normalized on #pix)

Texp_r = (5 * np.sqrt(fr + f_sky_r) / fr)**2   #exposure time in s

print("\nR band\n", Texp_r / 60, 'min')



