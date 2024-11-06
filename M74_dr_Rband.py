import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

### DATA #######################################################################################################

path='/Volumes/Noemi USB/Lab data acquisition/'

# collecting all the bias images
print('...')
bias = []
bias.append(fits.open(path+'20241029/Raw_data/Calibration_052_Bias.fit')[0].data)
bias.append(fits.open(path+'20241029/Raw_data/Calibration_053_Bias.fit')[0].data)
bias.append(fits.open(path+'20241029/Raw_data/Calibration_054_Bias.fit')[0].data)
bias.append(fits.open(path+'20241029/Raw_data/Calibration_055_Bias.fit')[0].data)
bias.append(fits.open(path+'20241029/Raw_data/Calibration_056_Bias.fit')[0].data)
bias.append(fits.open(path+'20241029/Raw_data/Calibration_057_Bias.fit')[0].data)
bias.append(fits.open(path+'20241029/Raw_data/Calibration_058_Bias.fit')[0].data)
bias.append(fits.open(path+'20241029/Raw_data/Calibration_059_Bias.fit')[0].data)
bias.append(fits.open(path+'20241029/Raw_data/Calibration_060_Bias.fit')[0].data)
bias.append(fits.open(path+'20241029/Raw_data/Calibration_061_Bias.fit')[0].data)
bias = np.array(bias)
print('Biases imported \n')


# collecting all the dark images
print('...')
dark_300 = []
dark_300.append(fits.open(path+'20241029/Raw_data/Calibration_052_Dark300.fit')[0].data)
dark_300.append(fits.open(path+'20241029/Raw_data/Calibration_053_Dark300.fit')[0].data)
dark_300.append(fits.open(path+'20241029/Raw_data/Calibration_054_Dark300.fit')[0].data)
dark_300.append(fits.open(path+'20241029/Raw_data/Calibration_055_Dark300.fit')[0].data)
dark_300.append(fits.open(path+'20241029/Raw_data/Calibration_056_Dark300.fit')[0].data)
dark_300.append(fits.open(path+'20241029/Raw_data/Calibration_057_Dark300.fit')[0].data)
dark_300.append(fits.open(path+'20241029/Raw_data/Calibration_058_Dark300.fit')[0].data)
dark_300.append(fits.open(path+'20241029/Raw_data/Calibration_059_Dark300.fit')[0].data)
dark_300.append(fits.open(path+'20241029/Raw_data/Calibration_060_Dark300.fit')[0].data)
dark_300.append(fits.open(path+'20241029/Raw_data/Calibration_061_Dark300.fit')[0].data)
dark_300 = np.array(dark_300)
print('Darks imported \n')


# collecting all the flat for r band images
print('...')
flatR = []
flatR.append(fits.open(path+'20241029/Raw_data/Calib-01-flat_r_017.fit')[0].data)
flatR.append(fits.open(path+'20241029/Raw_data/Calib-01-flat_r_018.fit')[0].data)
flatR.append(fits.open(path+'20241029/Raw_data/Calib-01-flat_r_019.fit')[0].data)
flatR.append(fits.open(path+'20241029/Raw_data/Calib-01-flat_r_020.fit')[0].data)
flatR.append(fits.open(path+'20241029/Raw_data/Calib-01-flat_r_021.fit')[0].data)
flatR = np.array(flatR)
print('r band flats imported \n')


#collecting the r band exposure time
R_head = fits.open(path+'20241029/Raw_data/M74_r_042.fit')[0].header
texp_R = R_head['EXPTIME']
#collecting all the r band images
print('...')
R = []
R.append(fits.open(path+'20241029/Raw_data/M74_r_042.fit')[0].data)
R.append(fits.open(path+'20241029/Raw_data/M74_r_043.fit')[0].data)
R.append(fits.open(path+'20241029/Raw_data/M74_r_044.fit')[0].data)
R.append(fits.open(path+'20241029/Raw_data/M74_r_045.fit')[0].data)
R = np.array(R)
print('exposure time: ' + str(texp_R))
print('r band images imported \n')


#%%
### STACKING CALIBRATION FRAMES ################################################################################

### master bias ###
master_bias = np.median(bias, axis=0)
#save
out = fits.PrimaryHDU(master_bias)
out.writeto(path+'20241029/20241029_bias.fit', overwrite=True)
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
out.writeto(path+'20241029/20241029_dark300.fit', overwrite=True)
#plot
plt.figure()
plt.imshow(master_dark300, cmap='viridis', clim=[300, 500])
plt.title('Master dark $t_{exp}=300$')
plt.colorbar()
plt.show()


### master flat for R ###
unbiased_flatR = [(x - master_bias) for x in flatR]   #I subtract the master bias for each flat image
norm_flatR = [(y / np.mean(y)) for y in unbiased_flatR]   #I normalize such as mean=1
master_flatR = np.median(norm_flatR, axis=0)
#save
out = fits.PrimaryHDU(master_flatR)
out.writeto(path+'20241029/20241029_flatR.fit', overwrite=True)
#plot
plt.figure()
plt.imshow(master_flatR, cmap='viridis', clim=[0.9, 1.1])
plt.title('Master flat for the r band')
plt.colorbar()
plt.show()


#%%
### PRODUCING IMAGES ###########################################################################################

### R band ###
imagesR = [((r - master_dark300) / master_flatR) for r in R]

#save
out = fits.PrimaryHDU(imagesR[0])
out.writeto(path+'20241029/20241029_R01.fit', overwrite=True)
#plot
plt.figure()
plt.imshow(imagesR[0], cmap='viridis', clim=[8000, 8400])
plt.title('M74 r band - 1')
plt.colorbar()
plt.show()

#save
out = fits.PrimaryHDU(imagesR[1])
out.writeto(path+'20241029/20241029_R02.fit', overwrite=True)
#plot
plt.figure()
plt.imshow(imagesR[1], cmap='viridis', clim=[7800, 8200])
plt.title('M74 r band - 2')
plt.colorbar()
plt.show()

#save
out = fits.PrimaryHDU(imagesR[2])
out.writeto(path+'20241029/20241029_R03.fit', overwrite=True)
#plot
plt.figure()
plt.imshow(imagesR[2], cmap='viridis', clim=[7800, 8200])
plt.title('M74 r band - 3')
plt.colorbar()
plt.show()

#save
out = fits.PrimaryHDU(imagesR[3])
out.writeto(path+'20241029/20241029_R04.fit', overwrite=True)
#plot
plt.figure()
plt.imshow(imagesR[3], cmap='viridis', clim=[7600, 8100])
plt.title('M74 r band - 4')
plt.colorbar()
plt.show()























