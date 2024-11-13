import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astroML.linear_model import PolynomialRegression
from tqdm import tqdm

#%%
### DATA #######################################################################################################

path='/Volumes/Noemi USB/Lab data acquisition/'

#collecting all the r unbiased and flattened images
print('...')
images = []
images.append(fits.open(path+'20241101/20241101_int_Ha01.fit')[0].data)
images.append(fits.open(path+'20241101/20241101_int_Ha02.fit')[0].data)
images.append(fits.open(path+'20241101/20241101_int_Ha03.fit')[0].data)
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
### ANALYSING IMAGES ###########################################################################################

#plot electron counting
plt.figure()
[plt.hist(images[i].flatten(), bins=150, histtype='step', lw=2, alpha=0.8, label='image '+str(i+1))
 for i in range(n_im)]
plt.yscale('log')
plt.xlim(-1000, 15000)
plt.ylabel('# pixels')
plt.xlabel('electron counts')
plt.title('Electron counts per pixel - with sky')
plt.legend()
plt.show()


#%%
### STAR MASK ##################################################################################################

#define the circle to mask the galaxy
cy = 2000
cx = 2300
radius = 430

mask = np.array(images)

py, px = np.ogrid[:3599, :4499]
circle = (px - cx)**2 + (py - cy)**2 <= radius**2

for i in range(n_im):
    mask[i] = np.where(circle, np.nan, images[i])   #mask the galaxy with a circle
    mask[i] = np.where(images[i] > 900, np.nan, mask[i])   #mask the stars above a certain intensity
    
#plot images
plotimage(mask[0], 340, 400, 'Mask - 1')
plotimage(mask[1], 330, 390, 'Mask - 2')
plotimage(mask[2], 330, 390, 'Mask - 3')


#%%
### SKY FIT ####################################################################################################

apx, apy = np.meshgrid(np.arange(4499), np.arange(3599))
all_pix = np.array([apy.ravel(), apx.ravel()])
all_pix = all_pix.T
model = PolynomialRegression(3)

sky_pix = []
sky_fit = []

print('Fit of the sky:')
for i in tqdm(range(n_im)):
    m = mask[i]
    sp = np.array(np.where(~np.isnan(m)))
    sky_pix.append(sp)
    sp = sp.T
    mask_not_nan = m[~np.isnan(m)]
    model.fit(sp, mask_not_nan)
    pred = model.predict(all_pix)
    sky_fit.append(pred.reshape((3599, 4499)))

#plot images
plotimage(sky_fit[0], 340, 400, 'Sky fit - 1')
plotimage(sky_fit[1], 330, 390, 'Sky fit - 2')
plotimage(sky_fit[2], 330, 390, 'Sky fit - 3')


#%%
### SKY SUBTRACTION ############################################################################################

im_nosky = [images[i] - sky_fit[i] for i in range(n_im)]

#plot images
plotimage(im_nosky[0], 0, 50, 'Subtracted sky - 1')
plotimage(im_nosky[1], 0, 50, 'Subtracted sky - 2')
plotimage(im_nosky[2], 0, 50, 'Subtracted sky - 3')
#and save
for i in range(n_im):
    out = fits.PrimaryHDU(im_nosky[i])
    out.writeto(path+'20241101/20241101_Ha0'+str(i+1)+'.fit', overwrite=True)

#plot electron counts
plt.figure()
[plt.hist(im_nosky[i].flatten(), bins=150, histtype='step', lw=2, alpha=0.8, label='image '+str(i+1))
 for i in range(n_im)]
plt.yscale('log')
plt.xlim(-11000, 14000)
plt.ylabel('# pixels')
plt.xlabel('electron counts')
plt.title('Electron counts per pixel - with sky')
plt.legend()
plt.show()










