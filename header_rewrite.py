from astropy.io import fits

#%%
### DATA #######################################################################################################

path='/Volumes/NOEMI USB/Lab data acquisition/'

#collecting the final images and their header
R = fits.open(path+'M74_Rband.fits')
Ha = fits.open(path+'M74_Ha.fits')
head_R = R[0].header
head_Ha = Ha[0].header

#collecting the wcs
wcs_R = fits.open(path+'wcs/wcs_Rband.fits')[0].header
wcs_Ha = fits.open(path+'wcs/wcs_Ha.fits')[0].header

head_R = wcs_R
head_R['EXPTIME'] = 300.
head_R['FILTER'] = 'r'
head_R['CCDTEMP'] = -20.

head_Ha = wcs_R
head_Ha['EXPTIME'] = 1200.
head_Ha['FILTER'] = 'Ha'
head_Ha['CCDTEMP'] = -20.

R[0].header = head_R
Ha[0].header = head_Ha

R.writeto(path+'M74_wcs_Rband.fits', overwrite=True)
Ha.writeto(path+'M74_wcs_Ha.fits', overwrite=True)


#%%
### TEST #######################################################################################################

testR = fits.open(path+'M74_wcs_Rband.fits')
testHa = fits.open(path+'M74_wcs_Ha.fits')

#print(testR[0].header)
#print(testHa[0].header)