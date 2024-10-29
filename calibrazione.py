import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

### DATA #######################################################################################################

path='/Volumes/Noemi USB/Lab data acquisition/Calibrazione/'

# collecting all the bias images
imbias = []
imbias.append(fits.open(path+'20231011/image_052_bias.fit')[0].data)
imbias.append(fits.open(path+'20231011/image_053_bias.fit')[0].data)
imbias.append(fits.open(path+'20231011/image_054_bias.fit')[0].data)
imbias.append(fits.open(path+'20231011/image_055_bias.fit')[0].data)
imbias.append(fits.open(path+'20231011/image_056_bias.fit')[0].data)
imbias = np.array(imbias)
print('Biases imported')


# collecting all the dark images
imdark_1 = []
imdark_1.append(fits.open(path+'20231011/image_052_dark1.fit')[0].data)
imdark_1.append(fits.open(path+'20231011/image_053_dark1.fit')[0].data)
imdark_1.append(fits.open(path+'20231011/image_054_dark1.fit')[0].data)
imdark_1.append(fits.open(path+'20231011/image_055_dark1.fit')[0].data)
imdark_1.append(fits.open(path+'20231011/image_056_dark1.fit')[0].data)
imdark_1 = np.array(imdark_1)

imdark_10 = []
imdark_10.append(fits.open(path+'20231011/image_052_dark10.fit')[0].data)
imdark_10.append(fits.open(path+'20231011/image_053_dark10.fit')[0].data)
imdark_10.append(fits.open(path+'20231011/image_054_dark10.fit')[0].data)
imdark_10.append(fits.open(path+'20231011/image_055_dark10.fit')[0].data)
imdark_10.append(fits.open(path+'20231011/image_056_dark10.fit')[0].data)
imdark_10 = np.array(imdark_10)

imdark_100 = []
imdark_100.append(fits.open(path+'20231011/image_052_dark100.fit')[0].data)
imdark_100.append(fits.open(path+'20231011/image_053_dark100.fit')[0].data)
imdark_100.append(fits.open(path+'20231011/image_054_dark100.fit')[0].data)
imdark_100.append(fits.open(path+'20231011/image_055_dark100.fit')[0].data)
imdark_100.append(fits.open(path+'20231011/image_056_dark100.fit')[0].data)
imdark_100 = np.array(imdark_100)

imdark_1000 = []
imdark_1000.append(fits.open(path+'20231011/image__107_dark1000.fit')[0].data)
imdark_1000.append(fits.open(path+'20231011/image__108_dark1000.fit')[0].data)
imdark_1000.append(fits.open(path+'20231011/image__109_dark1000.fit')[0].data)
imdark_1000.append(fits.open(path+'20231011/image__110_dark1000.fit')[0].data)
imdark_1000.append(fits.open(path+'20231011/image__111_dark1000.fit')[0].data)
imdark_1000 = np.array(imdark_1000)
print('Darks imported')


# collecting all the flat for Ha images
imflatHa = []
imflatHa.append(fits.open(path+'20231011/image_036_flat_Ha.fit')[0].data)
imflatHa.append(fits.open(path+'20231011/image_037_flat_Ha.fit')[0].data)
imflatHa.append(fits.open(path+'20231011/image_038_flat_Ha.fit')[0].data)
imflatHa.append(fits.open(path+'20231011/image_039_flat_Ha.fit')[0].data)
imflatHa.append(fits.open(path+'20231011/image_040_flat_Ha.fit')[0].data)
imflatHa = np.array(imflatHa)
print('Ha flats imported')


# collecting all the flat for r band images
imflatR = []
imflatR.append(fits.open(path+'20231011/image_020_flat_r.fit')[0].data)
imflatR.append(fits.open(path+'20231011/image_021_flat_r.fit')[0].data)
imflatR.append(fits.open(path+'20231011/image_022_flat_r.fit')[0].data)
imflatR.append(fits.open(path+'20231011/image_023_flat_r.fit')[0].data)
imflatR.append(fits.open(path+'20231011/image_024_flat_r.fit')[0].data)
imflatR = np.array(imflatR)
print('r band flats imported')


#collecting all the Ha images
imHa = []
imHa.append(fits.open(path+'20231011/ngc6946_ha_084.fit')[0].data)
imHa.append(fits.open(path+'20231011/ngc6946_ha_085.fit')[0].data)
imHa.append(fits.open(path+'20231011/ngc6946_ha_086.fit')[0].data)
imHa = np.array(imHa)
print('Ha images imported')


#collecting all the r band images
imR = []
imR.append(fits.open(path+'20231011/ngc6946_r_083.fit')[0].data)
imR.append(fits.open(path+'20231011/ngc6946_r_087.fit')[0].data)
imR.append(fits.open(path+'20231011/ngc6946_r_088.fit')[0].data)
imR = np.array(imR)
print('r band images imported')


#%%
### STACKING CALIBRATION FRAMES ################################################################################

bias_mean = np.mean(imbias, axis=0)
plt.figure()
plt.imshow(bias_mean, cmap='viridis', clim=[300, 450])
plt.title('Bias')
plt.colorbar()
plt.show()

dark_1_mean = np.mean(imdark_1, axis=0)
plt.figure()
plt.imshow(dark_1_mean, cmap='viridis', clim=[300, 450])
plt.title('Dark current $t_{exp}=1$')
plt.colorbar()
plt.show()

dark_10_mean = np.mean(imdark_10, axis=0)
plt.figure()
plt.imshow(dark_10_mean, cmap='viridis', clim=[1100, 1300])
plt.title('Dark current $t_{exp}=10$')
plt.colorbar()
plt.show()

dark_100_mean = np.mean(imdark_100, axis=0)
plt.figure()
plt.imshow(dark_100_mean, cmap='viridis', clim=[300, 500])
plt.title('Dark current $t_{exp}=100$')
plt.colorbar()
plt.show()

dark_1000_mean = np.mean(imdark_1000, axis=0)
plt.figure()
plt.imshow(dark_100_mean, cmap='viridis', clim=[300, 500])
plt.title('Dark current $t_{exp}=1000$')
plt.colorbar()
plt.show()

flatHa_mean = np.mean(imflatHa, axis=0)
normHa = np.max(flatHa_mean)
flatHa_norm = (flatHa_mean - bias_mean) / normHa
plt.figure()
plt.imshow(flatHa_norm, cmap='viridis', clim=[0.55, 0.65])
plt.title('Flat field $H\\alpha$')
plt.colorbar()
plt.show()

flatR_mean = np.mean(imflatR, axis=0)
normR = np.max(flatR_mean)
flatR_norm = (flatR_mean - bias_mean) / normR
plt.figure()
plt.imshow(flatR_norm, cmap='viridis', clim=[0.75, 0.85])
plt.title('Flat field r band')
plt.colorbar()
plt.show()


#%%
### H ALPHA ####################################################################################################

Ha_mean = np.mean(imHa, axis=0)

plt.figure()
plt.imshow(Ha_mean, cmap='viridis', clim=[930, 1100])
plt.title('$H\\alpha$ - mean')
plt.colorbar()
plt.show()

Ha = (Ha_mean - dark_100_mean) / flatHa_norm

plt.figure()
plt.imshow(Ha, cmap='viridis', clim=[930, 1100])
plt.title('$H\\alpha$')
plt.colorbar()
plt.show()

plt.figure()
plt.imshow(Ha, cmap='viridis', clim=[930, 1100])
plt.title('$H\\alpha$')
plt.xlim(1200, 3000)
plt.ylim(2500, 1000)
plt.colorbar()
plt.show()

#%%
### R BAND #####################################################################################################

R_mean = np.mean(imR, axis=0)

plt.figure()
plt.imshow(R_mean, cmap='viridis', clim=[19000, 25000])
plt.title('r band - mean')
plt.colorbar()
plt.show()

R = (R_mean - dark_100_mean) / flatR_norm

plt.figure()
plt.imshow(R, cmap='viridis', clim=[25500, 26500])
plt.title('$H\\alpha$')
plt.colorbar()
plt.show()

plt.figure()
plt.imshow(R, cmap='viridis', clim=[25500, 26500])
plt.title('$H\\alpha$')
plt.xlim(1200, 3000)
plt.ylim(2500, 1000)
plt.colorbar()
plt.show()









