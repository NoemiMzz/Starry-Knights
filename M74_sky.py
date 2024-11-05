import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

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
### SKY REMOVAL ################################################################################################

med = []
gradient = []

for i in range(4):
    mask = np.concatenate((imagesR[i][:, 100:400], imagesR[i][:, -400:-100]), axis=1)
    med.append( np.median(mask, axis=1) )
    gradient.append( med[i] / np.mean(med[i]) )

plt.figure()
[plt.scatter(np.arange(3599), med[i], s=1) for i in range(4)]
plt.title('M74 r - gradient')
plt.xlabel('pixels in a row')
plt.ylabel('electron counts')
plt.show()

plt.figure()
[plt.scatter(np.arange(3599), gradient[i], s=1) for i in range(4)]
plt.title('M74 r - normalized gradient')
plt.xlabel('pixels in a row')
plt.ylabel('electron ratio')
plt.show()










