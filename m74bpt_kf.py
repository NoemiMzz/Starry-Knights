#!pip install photutils
#import os
#from google.colab import drive
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import photutils.aperture as pha
from tqdm import tqdm

#%%
### DATA #######################################################################################################

path = '/Volumes/NOEMI USB/Lab data acquisition/bpt/'

printfileinfo = False

with fits.open(path+'m74datacube.fits') as hdul:
    
    if printfileinfo:
        hdul.info()   #examine the infos about the extensions in the file

    data_cube = hdul[0].data   #extract the data cube from the first extension (usually 0 is the primary header)
    header = hdul[0].header

    imHa = hdul['HA6562_FLUX'].data   #import fluxes
    #imSii= hdul['SII6716_FLUX'].data
    imOiii= hdul['OIII5006_FLUX'].data
    imHb= hdul['HB4861_FLUX'].data
    imNii = hdul['NII6548_FLUX'].data
    im = [imHa, imNii, imHb, imOiii]
    im = np.array(im)
    im_names = ["$H\\alpha$", "NII", "$H\\beta$", "OIII"]

    imHa_err = hdul['HA6562_FLUX_ERR'].data   #import errors
    #imSii_err= hdul['SII6716_FLUX_ERR'].data
    imOiii_err= hdul['OIII5006_FLUX_ERR'].data
    imHb_err= hdul['HB4861_FLUX_ERR'].data
    imNii_err = hdul['NII6548_FLUX_ERR'].data
    im_err = [imHa_err, imNii_err, imHb_err, imOiii_err]
    im_err = np.array(im_err)

ARCSEC = 0.2   #arcesc per pixel


#%%
### FUNCTIONS ##################################################################################################

def plotimage(data, minclim, maxclim, title):
    plt.figure(dpi=150)
    plt.imshow(data, cmap='viridis', clim=[minclim, maxclim])
    plt.title(title)
    plt.colorbar()
    plt.show()
    
def plot4x4(datavec, minclim, maxclim, title, c, intpol):
    fig, axs = plt.subplots(2, 2, figsize=(9, 8), dpi=150, sharey=True, sharex=True)
    ax = axs.flatten()
    for i in range(4):
        color = ax[i].imshow(datavec[i], cmap=c, clim=[minclim, maxclim], interpolation=intpol)
        ax[i].set_title(im_names[i])
    fig.colorbar(color, ax=axs, orientation='vertical', fraction=.1)
    fig.suptitle(title, fontsize='16', y=0.96)
    plt.show()
    
def kewley(f):
  return 1.19 + 0.61 / (f - 0.47)   #central lambdas are slightly different

def kauffmann(f):
    return 1.3 + 0.61 / (f - 0.05)

def SB_pixel(f):
    return f / N_PIX / ARCSEC**2
    
    
#%%
### SHOW DATA ##################################################################################################

if printfileinfo:
    print("\n\n\n--- FILE HEADER ---\n")
    print(f"Header: {header}")

#plot all the filters
plot4x4(im, 0, 1000, "Datacube images", 'viridis', 'auto')


#%%
### SELECTING DATA #############################################################################################

SN = im/im_err
choosenSN = 5

#plot all the filters
plot4x4(SN, 0, choosenSN, "Signal to noise", 'coolwarm', 'None')

mask = SN >= choosenSN
im_SN = np.where(mask, im, np.nan)

#plot all the filters
plot4x4(im_SN, 0, 1000, "Points with SN>"+str(choosenSN), 'viridis', 'None')


#%%
### SELECTING VALID APERTURES ##################################################################################

# defining the grid with all the apertures
aperture_array = []
radius = 10
x = np.arange(0, 1401, 20)
y = np.arange(0, 1401, 20)
min_valid_pixels = 274   #~90% of 305, which is the n. of pixels in each aperture
aperture_array = [pha.CircularAperture((xi, yi), r=radius) for xi in x for yi in y]

N_PIX = 305   #numebr of pixels in an apertures defined in this way

valid_aperture_idx = []

#we chose to select the valid apertures using the OIII image since it's the one with less points with SN > 5

for j, ap in tqdm(enumerate(aperture_array)):
    f = pha.aperture_photometry(im[3], ap)['aperture_sum'][0]
    #we select apertures with surface brightness per pixel bigger than a certain treshold
    if SB_pixel(f) > 2200:
        valid_aperture_idx.append(j)

apertures = []
[apertures.append(aperture_array[i]) for i in valid_aperture_idx]   #selected apertures

#show the selected apertures
plt.figure(dpi=150)
plt.imshow(im_SN[3], cmap='viridis', clim=[0, 1000])
for ap in apertures:
    ap.plot(color='red', lw=0.8)
plt.title("Selected apertures")
plt.colorbar()
plt.show()

print("\nNumber of valid apertures:")
print(len(apertures))

#show apertures in each filter
fig, axs = plt.subplots(2, 2, figsize=(9, 8), dpi=150, sharey=True, sharex=True)
ax = axs.flatten()
for i in range(4):
    color = ax[i].imshow(im[i], cmap='viridis', clim=[0, 1000], interpolation='None')
    ax[i].set_title(im_names[i])
    for ap in apertures:
        ap.plot(color='red', lw=0.8, ax=ax[i])
fig.colorbar(color, ax=axs, orientation='vertical', fraction=.1)
fig.suptitle("Selected apertures", fontsize='16', y=0.96)
plt.show()


#%%
### FLUX COMPUTATION ###########################################################################################

flux = []

for i in range(4):
    f = []
    for ap in apertures:
        f.append(pha.aperture_photometry(im[i], ap)['aperture_sum'][0])
    f = np.array(f)
    flux.append(f)
    

#%%
### BPT DIAGRAM ################################################################################################
  
x = np.log10(flux[1] / flux[0])
y = np.log10(flux[3] / flux[2])

x_axis_kw = np.linspace(-1.6, 0.2, 100)
x_axis_kf = np.linspace(-1.6, -0.1, 100)

cc = y - kewley(x)

plt.figure(dpi=150)
#plt.scatter(x, y, s=15, color='royalblue')
plt.scatter(x, y, c=cc, s=15, cmap='RdYlBu_r', clim=[-1.5, 1.5])
plt.plot(x_axis_kw, kewley(x_axis_kw), color='gold', label='Kewley')
plt.plot(x_axis_kf, kauffmann(x_axis_kf), color='gold', ls='dashed', label='Kauffmann')
plt.legend(loc=4)
plt.colorbar()
plt.xlim(-1.7, 0.4)
plt.ylim(-1.2, 1)
plt.title("BPT diagram")
plt.xlabel('Log (NII/$H\\alpha$)')
plt.ylabel('Log (OIII/$H\\beta$)')
plt.text(-1.6, 0.45, "stellar ionization", color='cornflowerblue')
plt.text(-0.5, 0.7, "AGN ionization", color='indianred')
plt.show()

#convert from apertures to points
circles_x = []
circles_y = []
[circles_x.append(ap.positions[0]) for ap in apertures]
[circles_y.append(ap.positions[1]) for ap in apertures]

#show the ionization map on the galaxy
plt.figure(dpi=150)
plt.imshow(im_SN[0], cmap='Greys', clim=[0, 1000], alpha=0.4)
plt.scatter(circles_x, circles_y, c=cc, cmap='RdYlBu_r', clim=[-1.5, 1.5])
plt.title("Ionization map")
plt.colorbar(label='distance from Kewley law')
plt.show()


