import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

#%%
### DATA #######################################################################################################

path='/Volumes/NOEMI USB/Lab data acquisition/'

#the calibration images are from 11nov since the ones from 29oct had some issues

# collecting all the bias images
print('...')
bias = []
bias.append(fits.open(path+'Raw_data/calib_037_bias.fit')[0].data)
bias.append(fits.open(path+'Raw_data/calib_038_bias.fit')[0].data)
bias.append(fits.open(path+'Raw_data/calib_039_bias.fit')[0].data)
bias.append(fits.open(path+'Raw_data/calib_040_bias.fit')[0].data)
bias.append(fits.open(path+'Raw_data/calib_041_bias.fit')[0].data)
bias.append(fits.open(path+'Raw_data/calib_042_bias.fit')[0].data)
bias = np.array(bias)
print('Biases imported \n')


# collecting all the dark images
print('...')
dark_300 = []
dark_300.append(fits.open(path+'Raw_data/calib_037_dark300.fit')[0].data)
dark_300.append(fits.open(path+'Raw_data/calib_038_dark300.fit')[0].data)
dark_300.append(fits.open(path+'Raw_data/calib_039_dark300.fit')[0].data)
dark_300.append(fits.open(path+'Raw_data/calib_040_dark300.fit')[0].data)
dark_300.append(fits.open(path+'Raw_data/calib_041_dark300.fit')[0].data)
dark_300.append(fits.open(path+'Raw_data/calib_042_dark300.fit')[0].data)
dark_300 = np.array(dark_300)
print('Darks imported \n')


# collecting all the flat for r band images
print('...')
flatR = []
flatR.append(fits.open(path+'Raw_data/calib_025_flat_r.fit')[0].data)
flatR.append(fits.open(path+'Raw_data/calib_026_flat_r.fit')[0].data)
flatR.append(fits.open(path+'Raw_data/calib_027_flat_r.fit')[0].data)
flatR.append(fits.open(path+'Raw_data/calib_028_flat_r.fit')[0].data)
flatR.append(fits.open(path+'Raw_data/calib_029_flat_r.fit')[0].data)
flatR = np.array(flatR)
print('r band flats imported \n')


#collecting the r band exposure time
R_head = fits.open(path+'Raw_data/M74_r_042.fit')[0].header
texp_R = R_head['EXPTIME']
#collecting all the r band images
print('...')
R = []
R.append(fits.open(path+'Raw_data/M74_r_042.fit')[0].data)
R.append(fits.open(path+'Raw_data/M74_r_043.fit')[0].data)
R.append(fits.open(path+'Raw_data/M74_r_044.fit')[0].data)
R.append(fits.open(path+'Raw_data/M74_r_045.fit')[0].data)
R = np.array(R)
print('exposure time: ' + str(texp_R))
print('r band images imported \n')

n_im = len(R)


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
out.writeto(path+'bias.fit', overwrite=True)
#plot
plotimage(master_bias, 300, 500, 'Master bias')


### master dark ###
#the t_exp is just right, rescaling is not needed
master_dark300 = np.median(dark_300, axis=0)
#save
out = fits.PrimaryHDU(master_dark300)
out.writeto(path+'dark300.fit', overwrite=True)
#plot
plotimage(master_dark300, 300, 500, 'Master dark $t_{exp}=300$')


### master flat for R ###
unbiased_flatR = [(x - master_bias) for x in flatR]   #I subtract the master bias for each flat image
norm_flatR = [(y / np.mean(y)) for y in unbiased_flatR]   #I normalize such as mean=1
master_flatR = np.median(norm_flatR, axis=0)
#save
out = fits.PrimaryHDU(master_flatR)
out.writeto(path+'flatR.fit', overwrite=True)
#plot
plotimage(master_flatR, 0.9, 1.1, 'Master flat for the r band')


#%%
### PRODUCING IMAGES ###########################################################################################

### R band ###
imagesR = [((r - master_dark300) / master_flatR) for r in R]

#plot images
plotimage(imagesR[0], 8000, 8400, 'Raw image - 1')
plotimage(imagesR[1], 7800, 8200, 'Raw image - 2')
plotimage(imagesR[2], 7800, 8200, 'Raw image - 3')
plotimage(imagesR[3], 7600, 8100, 'Raw image - 4')
#and save
for i in range(n_im):
    out = fits.PrimaryHDU(imagesR[i])
    out.writeto(path+'int_R0'+str(i+1)+'.fit', overwrite=True)



