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
    imSii= hdul['SII6716_FLUX'].data
    imOiii= hdul['OIII5006_FLUX'].data
    imHb= hdul['HB4861_FLUX'].data
    #imNii = hdul['NII6548_FLUX'].data
    im = [imHa, imSii, imHb, imOiii]
    im = np.array(im)
    im_names = ["$H_{\\alpha}$", "SII", "$H_{\\beta}$", "OIII"]

    imHa_err = hdul['HA6562_FLUX_ERR'].data   #import errors
    imSii_err= hdul['SII6716_FLUX_ERR'].data
    imOiii_err= hdul['OIII5006_FLUX_ERR'].data
    imHb_err= hdul['HB4861_FLUX_ERR'].data
    #imNii_err = hdul['NII6548_FLUX_ERR'].data
    im_err = [imHa_err, imSii_err, imHb_err, imOiii_err]
    im_err = np.array(im_err)


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
  return 1.30 + 0.72 / (f - 0.32)   #central lambdas are slightly different
    
    
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


valid_aperture_idx = []

#we chose to select the valid apertures using the OIII image since it's the one with less points with SN > 5

for j, ap in tqdm(enumerate(aperture_array)):
    mask = ap.to_mask(method='center')   #masking the pixels inside each aperture
    #adapting the mask to image dimensions (1=inside, 0=out, an intermediate value otherwise)
    cutout = mask.to_image(im_SN[3].shape)   #coordinates of the pixels inside each aperture
    covered_pixels = im_SN[3][cutout.astype(bool)]   
    if np.count_nonzero(~np.isnan(covered_pixels)) > min_valid_pixels:   #number of pixels inside each aperture
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
    color = ax[i].imshow(im_SN[i], cmap='viridis', clim=[0, 1000], interpolation='None')
    ax[i].set_title(im_names[i])
    for ap in apertures:
        ap.plot(color='red', lw=0.8, ax=ax[i])
fig.colorbar(color, ax=axs, orientation='vertical', fraction=.1)
fig.suptitle("Selected apertures", fontsize='16', y=0.96)
plt.show()


#%%
### FLUX COMPUTATION ###########################################################################################

flux = []

im_SN_nonan = np.where(np.isnan(im_SN), 0, im_SN)

for i in range(4):
    print("-------FILTER ", i)
    f = []
    for ap in apertures:
        #counting the number of pixels in each aperture
        mask = ap.to_mask(method='center')
        cutout = mask.to_image(im_SN[i].shape)
        covered_pixels = im_SN[i][cutout.astype(bool)]
        n_pix = np.count_nonzero(~np.isnan(covered_pixels))
        print('\n', n_pix)
        if n_pix < 305:   #if the number of valid pixels inside the aperture is <305
            ff = pha.aperture_photometry(im_SN_nonan[i], ap)['aperture_sum'][0]
            f.append( ff / n_pix * 305 )   #I normalize the flux on 305 pixels
            print("norm")
        else:
            f.append(pha.aperture_photometry(im_SN_nonan[i], ap)['aperture_sum'][0])
    f = np.array(f)
    flux.append(f)
    

#%%
### BPT DIAGRAM ################################################################################################
  
x = np.log10(flux[1] / flux[0])
y = np.log10(flux[3] / flux[2])

x_axis = np.linspace(-1.2, 0, 100)

plt.figure(dpi=150)
plt.scatter(x, y, s=15, color='royalblue')
plt.plot(x_axis, kewley(x_axis), color='crimson')
plt.title("BPT diagram")
plt.xlabel('Log (SII/$H\\alpha$)')
plt.ylabel('Log (OIII/$H\\beta$)')
plt.text(-1.1, 0.35, "stellar ionization", color='tomato')
plt.text(-0.6, 0.7, "AGN ionization", color='tomato')
plt.show()

