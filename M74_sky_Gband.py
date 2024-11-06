import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astroML.linear_model import NadarayaWatson

#%%
### DATA #######################################################################################################

path='/Volumes/Noemi USB/Lab data acquisition/'

#collecting all the r unbiased and flattened images
print('...')
images = []
images.append(fits.open(path+'20241101/20241101_G01.fit')[0].data)
images = np.array(images)
print('Ha images imported \n')

n_im = len(images)


#%%
### FUNCTIONS ##################################################################################################

def error(X, y, fit_method):
    Y = fit_method.predict(X)
    return np.sqrt( np.sum( (y-Y)**2 ) / len(X) )


#%%
### ANALYSING IMAGES ###########################################################################################

#choosing column indices
cin = 200
cfin = 800

for i in range(n_im):
    plt.figure()
    plt.imshow(images[i], cmap='viridis', clim=[4950, 5100])
    plt.vlines(cin, 0, 3590, colors='white', linestyles='dashed')
    plt.vlines(cfin, 0, 3590, colors='white', linestyles='dashed')
    plt.vlines(4499-cin, 0, 3590, colors='white', linestyles='dashed')
    plt.vlines(4499-cfin, 0, 3590, colors='white', linestyles='dashed')
    plt.title('M74 g band - '+str(i+1))
    plt.colorbar()
    plt.show()
    
plt.figure()
[plt.hist(images[i].flatten(), bins=150, histtype='step', lw=2, alpha=0.8, label='image '+str(i+1))
 for i in range(n_im)]
plt.yscale('log')
plt.xlim(0, 70000)
plt.ylabel('# pixels')
plt.xlabel('electron counts')
plt.title('Electron counts per pixel')
plt.legend()
plt.show()


#%%
### SKY GRADIENT ###############################################################################################

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
plt.title('M74 g band - gradient')
plt.xlabel('pixels in a row')
plt.ylabel('electron counts')
plt.show()

#plot normalized gradient
plt.figure()
[plt.scatter(pixels, gradient[i], s=1) for i in range(n_im)]
plt.title('M74 g band - normalized gradient')
plt.xlabel('pixels in a row')
plt.ylabel('electron ratio')
plt.show()

#collect togheter the data from different images
pixels_all = pixels
pixels_skl = pixels_all[:, np.newaxis]
gradient_all = gradient.flatten()


### fit the gradient ###
model = NadarayaWatson(kernel='gaussian', h=30)
model.fit(pixels_skl, gradient_all)
correction = model.predict(pixels[:, np.newaxis])

#plot
plt.figure()
plt.scatter(pixels_all, gradient_all, color='gainsboro', s=1)
plt.plot(pixels, correction, color='r')
plt.title('Kernel estimation')
plt.xlabel('pixels in a row')
plt.ylabel('electron ratio')
plt.show()


#%%
### SKY GRADIENT ###############################################################################################

im_nograd = [x / correction[:, np.newaxis] for x in images]

#save
out = fits.PrimaryHDU(im_nograd[0])
out.writeto(path+'20241029/20241029_nograd_G01.fit', overwrite=True)
#plot
plt.figure()
plt.imshow(im_nograd[0], cmap='viridis', clim=[4950, 5100])
plt.title('M74 g band - 1')
plt.colorbar()
plt.show()

#plot zoom
plt.figure()
plt.imshow(im_nograd[0], cmap='viridis', clim=[4950, 5100])
plt.title('M74 g band - 1 zoom')
plt.xlim(1500, 3000)
plt.ylim(2500, 1400)
plt.colorbar()
plt.show()





