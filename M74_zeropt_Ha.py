import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import photutils.aperture as pha

#%%
### DATA #######################################################################################################

path='/Volumes/NOEMI USB/Lab data acquisition/'

im = fits.open(path+'M74_wcs_Ha.fits')[0].data


#%%
### FUNCTIONS ##################################################################################################

def plotimage(data, minclim, maxclim, title):
    plt.figure()
    plt.imshow(data, cmap='viridis', clim=[minclim, maxclim])
    plt.title(title)
    plt.colorbar()
    plt.show()

def color_transformation(SDSS_r, SDSS_g):
    return SDSS_r - 0.2568 * (SDSS_g - SDSS_r) + 0.1470


#%%
### OPTIMAL APERTURE ###########################################################################################

f = open(path+'stars_coord.txt')
stars = np.genfromtxt(f, delimiter=';', skip_header=1)

radii = np.arange(1, 25)

int_flux = []

for n in range(len(stars)):
    flux = []
    for r in radii:
        aperture = pha.CircularAperture((stars[n,0], stars[n,1]), r)
        ff = pha.aperture_photometry(im, aperture)
        flux.append(ff['aperture_sum'][0])
    int_flux.append(flux)

plt.figure()
[plt.plot(radii, int_flux[n]/int_flux[n][15]) for n in range(len(stars))]
plt.title('Growth curves')
plt.xlabel('aperture radius')
plt.ylabel('integrated flux')
plt.show()


#%%
### OPTIMAL APERTURE ###########################################################################################

flux15 = []
apertures = []

for n in range(len(stars)):
    a = pha.CircularAperture((stars[n,0], stars[n,1]), 15)
    apertures.append(a)
    ff = pha.aperture_photometry(im, a)
    flux15.append(ff['aperture_sum'][0])

plt.figure(dpi=150)
plt.imshow(im, cmap='viridis', clim=[0, 50])
[apertures[n].plot(color='red', lw=7) for n in range(len(apertures))]
plt.title('Calibration stars')
plt.colorbar()
plt.show()

m_Ha = [color_transformation(stars[i,5], stars[i,4]) for i in range(len(stars))]

zp_single = [m_Ha[i] + 2.5 * np.log10(flux15[i]) for i in range(len(stars))]

ZP = np.median(zp_single)
ZP_err = 0.7413 * (np.percentile(zp_single, 75) - np.percentile(zp_single, 25))



