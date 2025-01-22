import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import photutils.aperture as pha
from tqdm import tqdm

#%%
### DATA #######################################################################################################

path='/Volumes/NOEMI USB/Lab data acquisition/'

printfileinfo = False
approx_rms = False

imR = fits.open(path+'M74_Rband.fits')[0].data
imHa = fits.open(path+'M74_Ha.fits')[0].data
imnetHa = fits.open(path+'M74_netHa.fits')[0].data

with fits.open(path+'bpt/m74datacube.fits') as hdul:
    
    if printfileinfo:
        hdul.info()   #examine the infos about the extensions in the file

    data_cube = hdul[0].data   #extract the data cube from the first extension (usually 0 is the primary header)
    header = hdul[0].header

    imCubeHa = hdul['HA6562_FLUX'].data   #import fluxes
    imCubeHa_err = hdul['HA6562_FLUX_ERR'].data   #import errors

K = 7.9e-42   #SFR factor in [solar mass * s / yr / erg]
d = 9.5   #distance of the galaxy in Mpc
dist = d * 3e24   #convert distance from Mpc to cm

arcsec = 0.44   #arcsec per pixel for TOBi
arcsec_cube = 0.2   #arcsec per pixel for MUSE


#%%
### FUNCTIONS ##################################################################################################

def plotimage(data, minclim, maxclim, title):
    plt.figure(dpi=150)
    plt.imshow(data, cmap='viridis', clim=[minclim, maxclim])
    plt.title(title)
    plt.colorbar()
    plt.show()

def SFR(flux):
    return ((4 * np.pi * dist**2) * flux) * K

def rms(x):
    return np.sqrt(np.mean((x - np.mean(x))**2))


#%%
### DATACUBE APERTURE ##########################################################################################

gal_cx_cube = 700
gal_cy_cube = 700
radius_cube = 370

gal_aperture_cube = pha.CircularAperture((gal_cx_cube, gal_cy_cube), radius_cube)

plt.figure(dpi=150)
plt.imshow(imCubeHa, cmap='viridis', clim=[0, 1000])
plt.scatter(gal_cx_cube, gal_cy_cube, marker='x', s=70, color='r')
gal_aperture_cube.plot(color='red')
plt.title('Datacube net $H_{\\alpha}$')
plt.colorbar()
plt.show()


#%%
### RMS COMPUTATION MUSE #######################################################################################

print("\n--- MUSE data ---")

py, px = np.ogrid[:1426, :1412]
circle = (px - gal_cx_cube)**2 + (py - gal_cy_cube)**2 <= radius_cube**2

rms_pix_cube = imCubeHa_err[circle] * 10**(-20)

rms_sky_cube = np.sqrt(np.sum(np.square(rms_pix_cube)))

print('\nRMS taken from the sky:', rms_sky_cube)


#%%
### SFR ESTIMATION MUSE ########################################################################################

flc = pha.aperture_photometry(imCubeHa, gal_aperture_cube)
cube_flux = flc['aperture_sum'][0] * 10**(-20)

cube_SN = cube_flux / rms_sky_cube

print('\nS/N:', cube_SN)

cube_SFR_S = SFR(cube_flux)
cube_SFR_C = cube_SFR_S / 1.57
err_cube_SFR_S = rms_sky_cube * (cube_SFR_S / cube_flux)
err_cube_SFR_C = rms_sky_cube * (cube_SFR_C / cube_flux)
print('\nSFR (Salpeter):', np.round(cube_SFR_S, 5), '[solar mass / yr]')
print('with error:', np.round(err_cube_SFR_S, 5))
print('SFR (Chabrier):', np.round(cube_SFR_C, 6), '[solar mass / yr]')
print('with error:', np.round(err_cube_SFR_C, 6))


#%%
### GALAXY APERTURE ############################################################################################

gal_cx = 2287
gal_cy = 1941
radius = radius_cube * arcsec_cube / arcsec

gal_aperture = pha.CircularAperture((gal_cx, gal_cy), radius)

plt.figure(dpi=150)
plt.imshow(imR, cmap='viridis', clim=[0, 0.5])
plt.scatter(gal_cx, gal_cy, marker='x', s=70, color='r')
gal_aperture.plot(color='red')
plt.title('Zoom r band')
plt.xlim(1500, 3100)
plt.ylim(2800, 1200)
plt.colorbar()
plt.show()

plt.figure(dpi=150)
plt.imshow(imHa, cmap='viridis', clim=[0, 0.03])
plt.scatter(gal_cx, gal_cy, marker='x', s=70, color='r')
gal_aperture.plot(color='red')
plt.title('Zoom $H_{\\alpha}$')
plt.xlim(1500, 3100)
plt.ylim(2800, 1200)
plt.colorbar()
plt.show()

plt.figure(dpi=150)
plt.imshow(imnetHa, cmap='viridis', clim=[0, 1.7e-16])
plt.scatter(gal_cx, gal_cy, marker='x', s=70, color='r')
gal_aperture.plot(color='red')
plt.title('Zoom net $H_{\\alpha}$')
plt.xlim(1500, 3100)
plt.ylim(2800, 1200)
plt.colorbar()
plt.show()


#%%
### RMS COMPUTATION APPROX #####################################################################################

if approx_rms:
    print("\n--- TOBi data ---")
    
    rms_cx = 1000
    rms_cy = 2700
    rms_aperture = pha.CircularAperture((rms_cx, rms_cy), radius)
    
    plt.figure(dpi=150)
    plt.imshow(imnetHa, cmap='viridis', clim=[0, 1.7e-16])
    plt.scatter(gal_cx, gal_cy, marker='x', s=70, color='r')
    gal_aperture.plot(color='red')
    rms_aperture.plot(color='cyan')
    plt.title('Net $H_{\\alpha}$')
    plt.colorbar()
    plt.show()
    
    py, px = np.ogrid[:3599, :4499]
    circle = (px - rms_cx)**2 + (py - rms_cy)**2 <= radius**2
    
    rms_pix = np.std(imnetHa[circle])
    
    n_pixels = len(imnetHa[circle])
    rms_sky = rms_pix * np.sqrt(n_pixels)
    
    print('\nRMS per pixel:', rms_pix)
    print('RMS taken from the sky:', rms_sky)


#%%
### RMS COMPUTATION ############################################################################################

else:
    print("\n--- TOBi data ---")
    
    rms_cx = np.arange(radius+20, 4499-radius, int(radius))
    rms_cy = np.arange(radius+20, 3599-radius, int(radius))
    
    rms_aperture = []
    for rx in rms_cx:
        for ry in rms_cy:
            if np.sqrt((rx-gal_cx)**2 + (ry-gal_cy)**2) >= 430*2:
                rms_aperture.append(pha.CircularAperture((rx, ry), radius))
    
    plt.figure(dpi=150)
    plt.imshow(imnetHa, cmap='viridis', clim=[0, 1.7e-16])
    plt.scatter(gal_cx, gal_cy, marker='x', s=70, color='r')
    gal_aperture.plot(color='red')
    [rap.plot(color='cyan') for rap in rms_aperture]
    plt.title('Net $H_{\\alpha}$')
    plt.colorbar()
    plt.show()
    
    rms_fluxes = [pha.aperture_photometry(imnetHa, rap)['aperture_sum'][0] for rap in rms_aperture]
    
    sky_mean = np.mean(rms_fluxes)
    rms_sky = np.std(rms_fluxes)
    
    print('\nRMS taken from the sky:', rms_sky)


#%%
### SFR ESTIMATION #############################################################################################

fl = pha.aperture_photometry(imnetHa, gal_aperture)
gal_flux = fl['aperture_sum'][0]

sigma_dist = (gal_flux - sky_mean) / rms_sky

gal_SN = gal_flux / rms_sky

print('\nS/N:', gal_SN)

sky_vals = '$\\langle f_{sky} \\rangle$='+str(np.round(sky_mean, 14))+'  $\\sigma_{sky}$='+str(np.round(rms_sky, 14))
gal_vals = '$f_{gal}$='+str(np.round(gal_flux, 13))+'  dist='+str(np.round(sigma_dist, 1))+'$\\sigma$'

if not approx_rms:
    plt.hist(rms_fluxes, bins=int(np.sqrt(len(rms_fluxes))), density=True, color='c', label=sky_vals, alpha=0.6)
    plt.axvline(gal_flux, color='red', label=gal_vals, lw=2)
    plt.legend()
    plt.title("Sky flux distribution")
    plt.xlabel("flux (erg/s$^2$/cm$^2$)")
    plt.show()

gal_SFR_S = SFR(gal_flux)
gal_SFR_C = gal_SFR_S / 1.57
err_SFR_S = rms_sky * (gal_SFR_S / gal_flux)
err_SFR_C = rms_sky * (gal_SFR_C / gal_flux)
print('\nSFR (Salpeter):', np.round(gal_SFR_S, 3), '[solar mass / yr]')
print('with error:', np.round(err_SFR_S, 3))
print('SFR (Chabrier):', np.round(gal_SFR_C, 3), '[solar mass / yr]')
print('with error:', np.round(err_SFR_C, 3))