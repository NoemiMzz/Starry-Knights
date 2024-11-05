import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astroML.linear_model import NadarayaWatson


#%%
### DATA #######################################################################################################

path='/Volumes/Noemi USB/Lab data acquisition/'

#collecting all the r unbiased and flattened images
print('...')
imagesR = []
imagesR.append(fits.open(path+'20241029/20241029_R01.fit')[0].data)
imagesR.append(fits.open(path+'20241029/20241029_R02.fit')[0].data)
imagesR.append(fits.open(path+'20241029/20241029_R03.fit')[0].data)
imagesR.append(fits.open(path+'20241029/20241029_R04.fit')[0].data)
imagesR = np.array(imagesR)
print('r images imported \n')


#%%
### FUNCTIONS ##################################################################################################

def error(X, y, fit_method):
    Y = fit_method.predict(X)
    return np.sqrt( np.sum( (y-Y)**2 ) / len(X) )


#%%
### ANALYSING IMAGES ###########################################################################################

for i in range(4):
    plt.figure()
    plt.imshow(imagesR[i], cmap='viridis', clim=[7600, 8400])
    plt.title('M74 r - '+str(i+1))
    plt.colorbar()
    plt.show()
    
    plt.figure()
    plt.hist(imagesR[i].flatten(), bins=150)
    plt.yscale('log')
    plt.xlim(2000, 70000)
    plt.ylabel('# pixels')
    plt.xlabel('electron counts')
    plt.title('M74 r - '+str(i+1))
    plt.show()


#%%
### SKY GRADIENT ###############################################################################################

### find the gradient ###
med = []
gradient = []
pixels = np.arange(3599)

for i in range(4):
    #I take two slits in the sky from left and right
    mask = np.concatenate((imagesR[i][:, 100:400], imagesR[i][:, -400:-100]), axis=1)
    med.append( np.median(mask, axis=1) )   #I take the median per row
    gradient.append( med[i] / np.mean(med[i]) )   #normalize different images
gradient = np.array(gradient)

#plot gradient
plt.figure()
[plt.scatter(pixels, med[i], s=1) for i in range(4)]
plt.title('M74 r - gradient')
plt.xlabel('pixels in a row')
plt.ylabel('electron counts')
plt.show()

#plot normalized gradient
plt.figure()
[plt.scatter(pixels, gradient[i], s=1) for i in range(4)]
plt.title('M74 r - normalized gradient')
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

R_nograd = [r / correction[:, np.newaxis] for r in imagesR]

#save
out = fits.PrimaryHDU(R_nograd[0])
out.writeto(path+'20241029/20241029_nograd_R01.fit', overwrite=True)
#plot
plt.figure()
plt.imshow(R_nograd[0], cmap='viridis', clim=[8050, 8300])
plt.title('M74 r band - 1')
plt.colorbar()
plt.show()

#plot zoom
plt.figure()
plt.imshow(R_nograd[0], cmap='viridis', clim=[8070, 8300])
plt.title('M74 r band - 1 zoom')
plt.xlim(1500, 3000)
plt.ylim(2500, 1400)
plt.colorbar()
plt.show()

#save
out = fits.PrimaryHDU(R_nograd[1])
out.writeto(path+'20241029/20241029_nograd_R02.fit', overwrite=True)
#plot
plt.figure()
plt.imshow(R_nograd[1], cmap='viridis', clim=[7850, 8100])
plt.title('M74 r band - 2')
plt.colorbar()
plt.show()

#save
out = fits.PrimaryHDU(R_nograd[2])
out.writeto(path+'20241029/20241029_nograd_R03.fit', overwrite=True)
#plot
plt.figure()
plt.imshow(R_nograd[2], cmap='viridis', clim=[7750, 8000])
plt.title('M74 r band - 3')
plt.colorbar()
plt.show()

#save
out = fits.PrimaryHDU(R_nograd[3])
out.writeto(path+'20241029/20241029_nograd_R04.fit', overwrite=True)
#plot
plt.figure()
plt.imshow(R_nograd[3], cmap='viridis', clim=[7600, 7850])
plt.title('M74 r band - 4')
plt.colorbar()
plt.show()





