import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

#%%
### DATA #######################################################################################################

path='/Volumes/Noemi USB/Lab data acquisition/'

# collecting all the bias images
print('...')
bias = []
bias.append(fits.open(path+'20241101/Raw_data/calib_037_bias.fit')[0].data)
bias.append(fits.open(path+'20241101/Raw_data/calib_038_bias.fit')[0].data)
bias.append(fits.open(path+'20241101/Raw_data/calib_039_bias.fit')[0].data)
bias.append(fits.open(path+'20241101/Raw_data/calib_040_bias.fit')[0].data)
bias.append(fits.open(path+'20241101/Raw_data/calib_041_bias.fit')[0].data)
bias.append(fits.open(path+'20241101/Raw_data/calib_042_bias.fit')[0].data)
bias = np.array(bias)
print('biases imported \n')


# collecting all the dark images
print('...')
dark_300 = []
dark_300.append(fits.open(path+'20241101/Raw_data/calib_037_dark300.fit')[0].data)
dark_300.append(fits.open(path+'20241101/Raw_data/calib_038_dark300.fit')[0].data)
dark_300.append(fits.open(path+'20241101/Raw_data/calib_039_dark300.fit')[0].data)
dark_300.append(fits.open(path+'20241101/Raw_data/calib_040_dark300.fit')[0].data)
dark_300.append(fits.open(path+'20241101/Raw_data/calib_041_dark300.fit')[0].data)
dark_300.append(fits.open(path+'20241101/Raw_data/calib_042_dark300.fit')[0].data)
dark_300 = np.array(dark_300)
print('Darks imported \n')


# collecting all the flat for r band images
print('...')
flatG = []
flatG.append(fits.open(path+'20241101/Raw_data/calib_019_flat_g.fit')[0].data)
flatG.append(fits.open(path+'20241101/Raw_data/calib_020_flat_g.fit')[0].data)
flatG.append(fits.open(path+'20241101/Raw_data/calib_021_flat_g.fit')[0].data)
flatG.append(fits.open(path+'20241101/Raw_data/calib_022_flat_g.fit')[0].data)
flatG.append(fits.open(path+'20241101/Raw_data/calib_023_flat_g.fit')[0].data)
flatG.append(fits.open(path+'20241101/Raw_data/calib_024_flat_g.fit')[0].data)
flatG = np.array(flatG)
print('G flats imported \n')


#collecting the G exposure time
G_head = fits.open(path+'20241101/Raw_data/m74_g_014.fit')[0].header
texp_G = G_head['EXPTIME']
#collecting all the G images
print('...')
G = []
G.append(fits.open(path+'20241101/Raw_data/m74_g_014.fit')[0].data)
G = np.array(G)
print('exposure time: ' + str(texp_G))
print('G images imported \n')


#%%
### STACKING CALIBRATION FRAMES ################################################################################

### master bias ###
master_bias = np.median(bias, axis=0)
#save
out = fits.PrimaryHDU(master_bias)
out.writeto(path+'20241101/20241101_bias.fit', overwrite=True)
#plot
plt.figure()
plt.imshow(master_bias, cmap='viridis', clim=[300, 500])
plt.title('Master bias')
plt.colorbar()
plt.show()


### master dark ###
#the t_exp is just right, rescaling is not needed
master_dark300 = np.median(dark_300, axis=0)
#save
out = fits.PrimaryHDU(master_dark300)
out.writeto(path+'20241101/20241101_dark300.fit', overwrite=True)
#plot
plt.figure()
plt.imshow(master_dark300, cmap='viridis', clim=[300, 500])
plt.title('Master dark $t_{exp}=300$')
plt.colorbar()
plt.show()


### master flat for G ###
unbiased_flatG = [(x - master_bias) for x in flatG]   #I subtract the master bias for each flat image
norm_flatG = [(y / np.mean(y)) for y in unbiased_flatG]   #I normalize such as mean=1
master_flatG = np.median(norm_flatG, axis=0)
#save
out = fits.PrimaryHDU(master_flatG)
out.writeto(path+'20241101/20241101_flatG.fit', overwrite=True)
#plot
plt.figure()
plt.imshow(master_flatG, cmap='viridis', clim=[0.9, 1.1])
plt.title('Master flat for the g band$')
plt.colorbar()
plt.show()


#%%
### PRODUCING IMAGES ###########################################################################################

### R band ###
imagesG = [((h - master_dark300) / master_flatG) for h in G]

#save
out = fits.PrimaryHDU(imagesG[0])
out.writeto(path+'20241101/20241101_G01.fit', overwrite=True)
#plot
plt.figure()
plt.imshow(imagesG[0], cmap='viridis', clim=[4950, 5100])
plt.title('M74 g band - 1')
plt.colorbar()
plt.show()























