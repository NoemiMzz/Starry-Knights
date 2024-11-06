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
dark_1200 = []
dark_1200.append(fits.open(path+'20241101/Raw_data/calib_037_dark1200.fit')[0].data)
dark_1200.append(fits.open(path+'20241101/Raw_data/calib_038_dark1200.fit')[0].data)
dark_1200.append(fits.open(path+'20241101/Raw_data/calib_039_dark1200.fit')[0].data)
dark_1200.append(fits.open(path+'20241101/Raw_data/calib_040_dark1200.fit')[0].data)
dark_1200.append(fits.open(path+'20241101/Raw_data/calib_041_dark1200.fit')[0].data)
dark_1200.append(fits.open(path+'20241101/Raw_data/calib_042_dark1200.fit')[0].data)
dark_1200 = np.array(dark_1200)
print('Darks imported \n')


# collecting all the flat for r band images
print('...')
flatHa = []
flatHa.append(fits.open(path+'20241101/Raw_data/calib_031_flat_ha.fit')[0].data)
flatHa.append(fits.open(path+'20241101/Raw_data/calib_032_flat_ha.fit')[0].data)
flatHa.append(fits.open(path+'20241101/Raw_data/calib_033_flat_ha.fit')[0].data)
flatHa.append(fits.open(path+'20241101/Raw_data/calib_034_flat_ha.fit')[0].data)
flatHa.append(fits.open(path+'20241101/Raw_data/calib_035_flat_ha.fit')[0].data)
flatHa.append(fits.open(path+'20241101/Raw_data/calib_036_flat_ha.fit')[0].data)
flatHa = np.array(flatHa)
print('Ha flats imported \n')


#collecting the Ha exposure time
Ha_head = fits.open(path+'20241101/Raw_data/m74_Ha_011.fit')[0].header
texp_Ha = Ha_head['EXPTIME']
#collecting all the Ha images
print('...')
Ha = []
Ha.append(fits.open(path+'20241101/Raw_data/m74_Ha_011.fit')[0].data)
Ha.append(fits.open(path+'20241101/Raw_data/m74_Ha_012.fit')[0].data)
Ha.append(fits.open(path+'20241101/Raw_data/m74_Ha_013.fit')[0].data)
Ha = np.array(Ha)
print('exposure time: ' + str(texp_Ha))
print('Ha images imported \n')


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
master_dark1200 = np.median(dark_1200, axis=0)
#save
out = fits.PrimaryHDU(master_dark1200)
out.writeto(path+'20241101/20241101_dark1200.fit', overwrite=True)
#plot
plt.figure()
plt.imshow(master_dark1200, cmap='viridis', clim=[300, 500])
plt.title('Master dark $t_{exp}=1200$')
plt.colorbar()
plt.show()


### master flat for Ha ###
unbiased_flatHa = [(x - master_bias) for x in flatHa]   #I subtract the master bias for each flat image
norm_flatHa = [(y / np.mean(y)) for y in unbiased_flatHa]   #I normalize such as mean=1
master_flatHa = np.median(norm_flatHa, axis=0)
#save
out = fits.PrimaryHDU(master_flatHa)
out.writeto(path+'20241101/20241101_flatHa.fit', overwrite=True)
#plot
plt.figure()
plt.imshow(master_flatHa, cmap='viridis', clim=[0.9, 1.1])
plt.title('Master flat for H$\\alpha$')
plt.colorbar()
plt.show()


#%%
### PRODUCING IMAGES ###########################################################################################

### R band ###
imagesHa = [((h - master_dark1200) / master_flatHa) for h in Ha]

#save
out = fits.PrimaryHDU(imagesHa[0])
out.writeto(path+'20241101/20241101_Ha01.fit', overwrite=True)
#plot
plt.figure()
plt.imshow(imagesHa[0], cmap='viridis', clim=[340, 400])
plt.title('M74 H$\\alpha$ - 1')
plt.colorbar()
plt.show()

#save
out = fits.PrimaryHDU(imagesHa[1])
out.writeto(path+'20241101/20241101_Ha02.fit', overwrite=True)
#plot
plt.figure()
plt.imshow(imagesHa[1], cmap='viridis', clim=[330, 390])
plt.title('M74 H$\\alpha$ - 2')
plt.colorbar()
plt.show()

#save
out = fits.PrimaryHDU(imagesHa[2])
out.writeto(path+'20241101/20241101_Ha03.fit', overwrite=True)
#plot
plt.figure()
plt.imshow(imagesHa[2], cmap='viridis', clim=[330, 390])
plt.title('M74 H$\\alpha$ - 3')
plt.colorbar()
plt.show()























