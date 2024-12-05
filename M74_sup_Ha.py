import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import photutils.centroids as cent

#%%
### DATA #######################################################################################################

path='/Volumes/NOEMI USB/Lab data acquisition/'

#collecting all the r unbiased and flattened images
print('...')
images = []
images.append(fits.open(path+'Ha01.fit')[0].data)
images.append(fits.open(path+'Ha02.fit')[0].data)
images.append(fits.open(path+'Ha03.fit')[0].data)
images = np.array(images)
print('Ha images imported \n')

n_im = len(images)


#%%
### FUNCTIONS ##################################################################################################

def plotimage(data, minclim, maxclim, title):
    plt.figure()
    plt.imshow(data, cmap='viridis', clim=[minclim, maxclim])
    plt.title(title)
    plt.colorbar()
    plt.show()
    
def plotimagebw(data, minclim, maxclim, title):
    plt.figure()
    plt.imshow(data, cmap='Greys_r', clim=[minclim, maxclim])
    plt.title(title)
    plt.colorbar()
    plt.show()
    
    
#%%
### FINDING CENTROIDS ##########################################################################################

xcent = []
ycent = []

for i in range(n_im):
    xycent = cent.centroid_quadratic(images[i, 1730:1830, 1670:1770])   #finding centroid of the chosen star
    xcent.append(xycent[0]+1670)
    ycent.append(xycent[1]+1730)

xcent = np.array(xcent)
ycent = np.array(ycent)

xoff = xcent[0]-xcent   #compute the displacements
yoff = ycent[0]-ycent

print('Displacements:')
print(xoff)
print(yoff)

#plot the star centroids
fig, axs = plt.subplots(1, n_im, figsize=(12, 12*n_im))
axs = axs.flatten()
for i in range(n_im):
    axs[i].imshow(images[i, 1730:1830, 1670:1770], cmap='viridis', clim=[0, 50])
    axs[i].scatter(xcent[i]-1670, ycent[i]-1730, marker='x', s=70, color='r')
    fig.suptitle('Selected star', y=0.55)


#%%
### TRANSLATING CENTROIDS ######################################################################################

xoff_int = np.rint(xoff).astype(int)   #round displacements to int
yoff_int = np.rint(yoff).astype(int)   #I need to translate for an integer number of pixels


for i in range(n_im):
    images[i,...] = np.roll(images[i,...], (xoff_int[i], yoff_int[i]), axis = (1,0))   #translation
    
fig, axs = plt.subplots(1, n_im, figsize=(12, 12*n_im))
axs = axs.flatten()
for i in range(n_im):
    axs[i].imshow(images[i, 1730:1830, 1670:1770], cmap='viridis', clim=[0, 50])
    axs[i].scatter(xcent[0]-1670, ycent[0]-1730, marker='x', s=70, color='r')
    fig.suptitle('Translation', y=0.55)


#%%
### SUPERIMPOSING IMAGES #######################################################################################

final = np.median(images/1200, axis=0)   #median of all the images

for i in range(n_im):
    plotimage(images[i], 0, 50, 'final image -'+str(i+1))

plotimage(final, 0, 0.03, 'M74 in $H\\alpha$')

out = fits.PrimaryHDU(final)
out.writeto(path+'M74_Ha.fits', overwrite=True)




