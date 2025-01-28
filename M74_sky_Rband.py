import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astroML.linear_model import NadarayaWatson
from astroML.linear_model import PolynomialRegression
from tqdm import tqdm


#%%
### DATA #######################################################################################################

path='/Volumes/NOEMI USB/Lab data acquisition/'

#collecting all the r unbiased and flattened images
print('...')
images = []
images.append(fits.open(path+'int_R01.fit')[0].data)
images.append(fits.open(path+'int_R02.fit')[0].data)
images.append(fits.open(path+'int_R03.fit')[0].data)
images.append(fits.open(path+'int_R04.fit')[0].data)
images = np.array(images)
print('r images imported \n')

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

#choosing column indices
cin = 200
cfin = 800

#plot images with columns
for i in range(n_im):
    plt.figure()
    plt.imshow(images[i], cmap='viridis', clim=[7600, 8400])
    plt.vlines(cin, 0, 3590, colors='white', linestyles='dashed')
    plt.vlines(cfin, 0, 3590, colors='white', linestyles='dashed')
    plt.vlines(4499-cin, 0, 3590, colors='white', linestyles='dashed')
    plt.vlines(4499-cfin, 0, 3590, colors='white', linestyles='dashed')
    plt.title('Raw image - '+str(i+1))
    plt.colorbar()
    plt.show()

#plot electron counting
plt.figure()
[plt.hist(images[i].flatten(), bins=150, histtype='step', lw=2, alpha=0.8, label='image '+str(i+1))
 for i in range(n_im)]
plt.yscale('log')
plt.xlim(2000, 70000)
plt.ylabel('# pixels')
plt.xlabel('electron counts')
plt.title('Electron counts per pixel - with gradient')
plt.legend()
plt.show()


#%%
### FITTING GRADIENT ###########################################################################################

### find the gradient ###
med = []
gradient = []
pixels = np.arange(3599)

for i in range(n_im):
    #I take two slits in the sky from left and right
    mask = np.concatenate((images[i][:, cin:cfin], images[i][:, -cfin:-cin]), axis=1)
    med.append( np.median(mask, axis=1) )   #I take the median per row
    gradient.append( med[i] / np.mean(med[i]) )   #normalize different images
gradient = np.array(gradient)

#plot gradient
plt.figure()
[plt.scatter(pixels, med[i], s=1) for i in range(n_im)]
plt.title('Gradient')
plt.xlabel('pixels in a row')
plt.ylabel('electron counts')
plt.show()

#plot normalized gradient
plt.figure()
[plt.scatter(pixels, gradient[i], s=1) for i in range(n_im)]
plt.title('Normalized gradient')
plt.xlabel('pixels in a row')
plt.ylabel('electron ratio')
plt.show()

#collect togheter the data from different images
pixels_all = np.concatenate((pixels, pixels, pixels, pixels))
pixels_skl = pixels_all[:, np.newaxis]
gradient_all = gradient.flatten()


### fit the gradient ###
model = NadarayaWatson(kernel='gaussian', h=30)
model.fit(pixels_skl, gradient_all)
correction = model.predict(pixels[:, np.newaxis])

#plot the fit
plt.figure()
plt.scatter(pixels_all, gradient_all, color='gainsboro', s=1)
plt.plot(pixels, correction, color='r')
plt.title('KDE fit of the gradient')
plt.xlabel('pixels in a row')
plt.ylabel('electron ratio')
plt.show()


#%%
### REMOVE GRADIENT ############################################################################################

im_nograd = [x / correction[:, np.newaxis] for x in images]

#plot images
plotimage(im_nograd[0], 8050, 8300, 'Subtracted gradient - 1')
plotimage(im_nograd[1], 7900, 8100, 'Subtracted gradient - 2')
plotimage(im_nograd[2], 7800, 8000, 'Subtracted gradient - 3')
plotimage(im_nograd[3], 7650, 7850, 'Subtracted gradient - 4')

#plot zoom
plt.figure()
plt.imshow(im_nograd[0], cmap='viridis', clim=[8100, 8300])
plt.title('Subtracted gradient - 1 zoom')
plt.xlim(1500, 3000)
plt.ylim(2500, 1400)
plt.colorbar()
plt.show()

#plot electron counts
plt.figure()
[plt.hist(im_nograd[i].flatten(), bins=150, histtype='step', lw=2, alpha=0.8, label='image '+str(i+1))
 for i in range(n_im)]
plt.yscale('log')
plt.xlim(2000, 70000)
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
radius = 350

mask = np.array(im_nograd)

py, px = np.ogrid[:3599, :4499]
circle = (px - cx)**2 + (py - cy)**2 <= radius**2

for i in range(n_im):
    mask[i] = np.where(circle, np.nan, im_nograd[i])   #mask the galaxy with a circle
    mask[i] = np.where(im_nograd[i] > 8800, np.nan, mask[i])   #mask the stars above a certain intensity
    
#plot images
plotimage(mask[0], 8050, 8300, 'Mask - 1')
plotimage(mask[1], 7900, 8100, 'Mask - 2')
plotimage(mask[2], 7800, 8000, 'Mask - 3')
plotimage(mask[3], 7650, 7850, 'Mask - 4')


#%%
### SIGMA SKY (PROPOSAL) #######################################################################################

sigma_sky_pix = [np.nanstd(mask[i]) for i in range(n_im)]

sigma_sky_pix_mean = np.mean(sigma_sky_pix)


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
plotimage(sky_fit[0], 8050, 8300, 'Sky fit - 1')
plotimage(sky_fit[1], 7900, 8100, 'Sky fit - 2')
plotimage(sky_fit[2], 7800, 8000, 'Sky fit - 3')
plotimage(sky_fit[3], 7650, 7850, 'Sky fit - 4')


#%%
### SKY SUBTRACTION ############################################################################################

im_nosky = [im_nograd[i] - sky_fit[i] for i in range(n_im)]

#plot images
plotimage(im_nosky[0], 0, 300, 'Subtracted sky - 1')
plotimage(im_nosky[1], 0, 300, 'Subtracted sky - 2')
plotimage(im_nosky[2], 0, 300, 'Subtracted sky - 3')
plotimage(im_nosky[3], 0, 300, 'Subtracted sky - 4')
#and save
for i in range(n_im):
    out = fits.PrimaryHDU(im_nosky[i])
    out.writeto(path+'R0'+str(i+1)+'.fit', overwrite=True)

#plot electron counts
plt.figure()
[plt.hist(im_nosky[i].flatten(), bins=150, histtype='step', lw=2, alpha=0.8, label='image '+str(i+1))
 for i in range(n_im)]
plt.yscale('log')
plt.xlim(-10000, 60000)
plt.ylabel('# pixels')
plt.xlabel('electron counts')
plt.title('Electron counts per pixel - with sky')
plt.legend()
plt.show()








