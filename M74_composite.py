import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import photutils.aperture as pha
from tqdm import tqdm

#%%
### DATA #######################################################################################################

path='/Volumes/NOEMI USB/Lab data acquisition/'

imR = fits.open(path+'M74_Rband.fits')[0].data
imHa = fits.open(path+'M74_Ha.fits')[0].data

F0 = 3631   #reference flux in Junsky
lambda_r = 6225.0   #central wavelenght of the filter in A
lambda_Ha = 6568.894958496094
largHa = 36.64978027734375   #widht of the filter in A
largR = 1290.0

C_r = 17.125059673553615   #calibration coefficients with errors per band
C_r_err = 0.040374811601885774
C_Ha = 15.539801701451697
C_Ha_err = 0.01580691822059472


#%%
### FUNCTIONS ##################################################################################################

def plotimage(data, minclim, maxclim, title):
    plt.figure(dpi=150)
    plt.imshow(data, cmap='viridis', clim=[minclim, maxclim])
    plt.title(title)
    plt.colorbar()
    plt.show()

def m_to_flux(mag, lambda_c):
    f0 = F0 * 10**(-23)   #reference flux in erg/s/cm2/Hz
    f = f0 * 10**(-mag/2.5)   #flux in erg/s/cm2/Hz
    return f * 3*10**(18) / lambda_c**2

def err_median(x):
    return 0.7413 * (np.percentile(x, 75) - np.percentile(x, 25))

def f_norm(fR):
    return fR * largHa / largR

def net_electrons(fHa, fR, k):
    fR_norm = f_norm(fR)
    return fHa - fR_norm * k

def net(fHa, fR, k):
    fR_norm = f_norm(fR)
    return (fHa - fR_norm * k) * largHa * 10**(-C_Ha)


#%%
### CORRECTION FACTOR ESTIMATION ###############################################################################

f = open(path+'stars_coord.txt')   #same stars used in calibration
stars = np.genfromtxt(f, delimiter=';', skip_header=1)

star_flux_net = []
star_flux_Ha = []
kk = np.arange(0.9, 1.3, 0.05)

for n in tqdm(range(len(stars))):
    aperture = pha.CircularAperture((stars[n,0], stars[n,1]), 15)   #define the aperture on the stars
    
    ff_Ha = pha.aperture_photometry(imHa, aperture)   #collect flux from Ha
    star_flux_Ha.append(ff_Ha['aperture_sum'][0])
    
    ff_net = []
    for i in range(len(kk)):
        im = net_electrons(imHa, imR, kk[i])   #computing net Ha varying k
        f = pha.aperture_photometry(im, aperture)   #collect flux from netHa
        ff_net.append(f['aperture_sum'][0])
    star_flux_net.append(ff_net)
    
K = 1.04   #chosen correction factor

plt.figure()
[plt.plot(kk, star_flux_net[n] / star_flux_Ha[n], color='royalblue', lw=1) for n in range(len(stars))]
plt.axhline(0, color='black')
plt.axvline(K, color='black', ls='dashed')   #chosen k
plt.title('Correction factor estimation')
plt.xlabel('correction factor')
plt.ylabel('$F_{net}$ / $F_{H_{\\alpha}}$ ($e^{-}/s$)')
plt.show()

imnetHa = net(imHa, imR, K)
plotimage(imnetHa, 0, 1.7e-16, 'Net $H_{\\alpha}$')











