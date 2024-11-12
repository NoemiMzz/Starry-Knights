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
images.append(fits.open(path+'20241101/20241101_Ha01.fit')[0].data)
images.append(fits.open(path+'20241101/20241101_Ha02.fit')[0].data)
images.append(fits.open(path+'20241101/20241101_Ha03.fit')[0].data)
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

for i in range(n_im):
    plt.figure()
    plt.imshow(images[i], cmap='viridis', clim=[330, 390])
    plt.title('M74 H$\\alpha$ - '+str(i+1))
    plt.colorbar()
    plt.show()
    
plt.figure()
[plt.hist(images[i].flatten(), bins=150, histtype='step', lw=2, alpha=0.8, label='image '+str(i+1))
 for i in range(n_im)]
plt.yscale('log')
plt.xlim(-2500, 15000)
plt.ylabel('# pixels')
plt.xlabel('electron counts')
plt.title('Electron counts per pixel')
plt.legend()
plt.show()




