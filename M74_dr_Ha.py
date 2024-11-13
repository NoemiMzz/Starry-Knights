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

n_im = len(Ha)


#%%
### FUNCTIONS ##################################################################################################

def plotimage(data, minclim, maxclim, title):
    plt.figure()
    plt.imshow(data, cmap='viridis', clim=[minclim, maxclim])
    plt.title(title)
    plt.colorbar()
    plt.show()


#%%
### STACKING CALIBRATION FRAMES ################################################################################

### master bias ###
master_bias = np.median(bias, axis=0)
#save
out = fits.PrimaryHDU(master_bias)
out.writeto(path+'20241101/20241101_bias.fit', overwrite=True)
#plot
plotimage(master_bias, 300, 500, 'Master bias')


### master dark ###
#the t_exp is just right, rescaling is not needed
master_dark1200 = np.median(dark_1200, axis=0)
#save
out = fits.PrimaryHDU(master_dark1200)
out.writeto(path+'20241101/20241101_dark1200.fit', overwrite=True)
#plot
plotimage(master_dark1200, 300, 500, 'Master dark $t_{exp}=1200$')


### master flat for R ###
unbiased_flatHa = [(x - master_bias) for x in flatHa]   #I subtract the master bias for each flat image
norm_flatHa = [(y / np.mean(y)) for y in unbiased_flatHa]   #I normalize such as mean=1
master_flatHa = np.median(norm_flatHa, axis=0)
#save
out = fits.PrimaryHDU(master_flatHa)
out.writeto(path+'20241101/20241101_flatHa.fit', overwrite=True)
#plot
plotimage(master_flatHa, 0.9, 1.1, 'Master flat for H$\\alpha$')


#%%
### PRODUCING IMAGES ###########################################################################################

### R band ###
imagesHa = [((r - master_dark1200) / master_flatHa) for r in Ha]

#plot images
plotimage(imagesHa[0], 340, 400, 'Raw image - 1')
plotimage(imagesHa[1], 330, 390, 'Raw image - 2')
plotimage(imagesHa[2], 330, 390, 'Raw image - 3')
#and save
for i in range(n_im):
    out = fits.PrimaryHDU(imagesHa[i])
    out.writeto(path+'20241101/20241101_int_Ha0'+str(i+1)+'.fit', overwrite=True)









