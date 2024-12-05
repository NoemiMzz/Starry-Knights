import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import photutils.aperture as pha

#%%
### DATA #######################################################################################################

path='/Volumes/NOEMI USB/Lab data acquisition/'

im = fits.open(path+'M74_Rband.fits')[0].data

F0 = 3631   #reference flux in Junsky
lambda_r = 6580   #central wavelenght of the filter in A


#%%
### FUNCTIONS ##################################################################################################

def plotimage(data, minclim, maxclim, title):
    plt.figure()
    plt.imshow(data, cmap='viridis', clim=[minclim, maxclim])
    plt.title(title)
    plt.colorbar()
    plt.show()

def color_transformation(SDSS_r, SDSS_g):
    return SDSS_r - 0.0061 * (SDSS_g - SDSS_r) + 0.0001

def m_to_flux(mag, lambda_c):
    f0 = F0 * 10**(-23)   #reference flux in erg/s/cm2/Hz
    f = f0 * 10**(-mag/2.5)   #flux in erg/s/cm2/Hz
    return f * 3*10**(18) / lambda_c**2

def err_median(x):
    return 0.7413 * (np.percentile(x, 75) - np.percentile(x, 25))


#%%
### OPTIMAL APERTURE ###########################################################################################

f = open(path+'stars_coord.txt')
stars = np.genfromtxt(f, delimiter=';', skip_header=1)

radii = np.arange(1, 30)

int_flux = []

for n in range(len(stars)):
    flux = []
    for r in radii:
        aperture = pha.CircularAperture((stars[n,0], stars[n,1]), r)   #define aperture with different radii
        ff = pha.aperture_photometry(im, aperture)   #collect flux varying the aperture
        flux.append(ff['aperture_sum'][0])
    int_flux.append(flux)

#plot grow curves
plt.figure()
[plt.plot(radii, int_flux[n]/int_flux[n][20]) for n in range(len(stars))]
plt.title('Growth curves')
plt.xlabel('aperture radius')
plt.ylabel('integrated flux')
plt.show()


#%%
### ZERO POINT #################################################################################################

flux20 = []
apertures = []

for n in range(len(stars)):
    a = pha.CircularAperture((stars[n,0], stars[n,1]), 20)   #define aperture with choosen radius
    apertures.append(a)
    ff = pha.aperture_photometry(im, a)   #collect flux for each star
    flux20.append(ff['aperture_sum'][0])

#plot apertures
plt.figure(dpi=150)
plt.imshow(im, cmap='viridis', clim=[0, 1])
[apertures[n].plot(color='red', lw=7) for n in range(len(apertures))]
plt.title('Calibration stars')
plt.colorbar()
plt.show()

#extrapolate r band magnetude from SDSS r and g band
m_R = [color_transformation(stars[i,5], stars[i,4]) for i in range(len(stars))]

#compute the zero point for each star
zp_single = [m_R[i] + 2.5 * np.log10(flux20[i]) for i in range(len(stars))]

ZP = np.median(zp_single)
ZP_err = err_median(zp_single)

print('\nZero point r band:')
print(ZP)
print(ZP_err)


#%%
### C FACTOR ###################################################################################################

#convert Ha magnetudes in cgs fluxes
flux_SDSS = [m_to_flux(mag, lambda_r) for mag in m_R]

#compute the C factor for each star
c_single = np.log10(np.array(flux20) / np.array(flux_SDSS))

C = np.median(c_single)
C_err = err_median(c_single)

print('\nCalibration constant r band:')
print(C)
print(C_err)




