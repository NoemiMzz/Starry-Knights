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

K = 7.9e-42   #SFR factor in [solar mass * s / yr / erg]
d = 9.5   #distance of the galaxy in Mpc
dist = d * 3e24   #convert distance from Mpc to cm


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

def net_electrons(fHa, fR, a):
    fR_norm = f_norm(fR)
    return fHa - fR_norm * a

def net(fHa, fR, a):
    fR_norm = f_norm(fR)
    return (fHa - fR_norm * a) * largHa * 10**(-C_Ha)

def SFR(flux):
    return ((4 * np.pi * dist**2) * flux) * K

def rms(x):
    return np.sqrt(np.mean((x - np.mean(x))**2))


#%%
### CORRECTION FACTOR ESTIMATION ###############################################################################

f = open(path+'stars_coord.txt')   #same stars used in calibration
stars = np.genfromtxt(f, delimiter=';', skip_header=1)

star_flux_net = []
star_flux_Ha = []
aa = np.arange(0.9, 1.3, 0.05)

for n in tqdm(range(len(stars))):
    aperture = pha.CircularAperture((stars[n,0], stars[n,1]), 15)   #define the aperture on the stars
    
    ff_Ha = pha.aperture_photometry(imHa, aperture)   #collect flux from Ha
    star_flux_Ha.append(ff_Ha['aperture_sum'][0])
    
    ff_net = []
    for i in range(len(aa)):
        im = net_electrons(imHa, imR, aa[i])   #computing net Ha varying k
        f = pha.aperture_photometry(im, aperture)   #collect flux from netHa
        ff_net.append(f['aperture_sum'][0])
    star_flux_net.append(ff_net)
    
A = 1.04   #chosen correction factor

plt.figure()
[plt.plot(aa, star_flux_net[n] / star_flux_Ha[n], color='royalblue', lw=1) for n in range(len(stars))]
plt.axhline(0, color='black')
plt.axvline(A, color='black', ls='dashed')   #chosen a
plt.title('Correction factor estimation')
plt.xlabel('correction factor')
plt.ylabel('$F_{net}$ / $F_{H_{\\alpha}}$ ($e^{-}/s$)')
plt.show()

imnetHa = net(imHa, imR, A)
plotimage(imnetHa, 0, 1.7e-16, 'Net $H_{\\alpha}$')


#%%
### GALAXY APERTURE ############################################################################################

gal_cx = 2287
gal_cy = 1941
radius = 600

gal_aperture = pha.CircularAperture((gal_cx, gal_cy), radius)

plt.figure(dpi=150)
plt.imshow(imR, cmap='viridis', clim=[0, 0.5])
plt.scatter(gal_cx, gal_cy, marker='x', s=70, color='r')
gal_aperture.plot(color='red')
plt.title('Zoom r band')
plt.xlim(1500, 3100)
plt.ylim(2800, 1200)
plt.colorbar()
plt.show()

plt.figure(dpi=150)
plt.imshow(imHa, cmap='viridis', clim=[0, 0.03])
plt.scatter(gal_cx, gal_cy, marker='x', s=70, color='r')
gal_aperture.plot(color='red')
plt.title('Zoom $H_{\\alpha}$')
plt.xlim(1500, 3100)
plt.ylim(2800, 1200)
plt.colorbar()
plt.show()

plt.figure(dpi=150)
plt.imshow(imnetHa, cmap='viridis', clim=[0, 1.7e-16])
plt.scatter(gal_cx, gal_cy, marker='x', s=70, color='r')
gal_aperture.plot(color='red')
plt.title('Zoom net $H_{\\alpha}$')
plt.xlim(1500, 3100)
plt.ylim(2800, 1200)
plt.colorbar()
plt.show()


#%%
### RMS COMPUTATION ############################################################################################

rms_cx = 1000
rms_cy = 2700
rms_aperture = pha.CircularAperture((rms_cx, rms_cy), radius)

plt.figure(dpi=150)
plt.imshow(imnetHa, cmap='viridis', clim=[0, 1.7e-16])
plt.scatter(gal_cx, gal_cy, marker='x', s=70, color='r')
gal_aperture.plot(color='red')
rms_aperture.plot(color='cyan')
plt.title('Net $H_{\\alpha}$')
plt.colorbar()
plt.show()

py, px = np.ogrid[:3599, :4499]
circle = (px - rms_cx)**2 + (py - rms_cy)**2 <= radius**2

rms_pix = np.std(imnetHa[circle])

n_pixels = len(imnetHa[circle])
rms_sky = rms_pix * np.sqrt(n_pixels)

print('\nRMS per pixel:', rms_pix)
print('RMS taken from the sky:', rms_sky)


#%%
### SFR ESTIMATION #############################################################################################

fl = pha.aperture_photometry(imnetHa, gal_aperture)
gal_flux = fl['aperture_sum'][0]

gal_SFR_S = SFR(gal_flux)
gal_SFR_C = gal_SFR_S / 1.57
err_SFR_S = rms_sky * (gal_SFR_S / gal_flux)
err_SFR_C = rms_sky * (gal_SFR_C / gal_flux)
print('\nSFR of the whole galaxy (Salpeter):', np.round(gal_SFR_S, 3), '[solar mass / yr]')
print('with error:', np.round(err_SFR_S, 3))
print('SFR of the whole galaxy (Chabrier):', np.round(gal_SFR_C, 3), '[solar mass / yr]')
print('with error:', np.round(err_SFR_C, 3))



