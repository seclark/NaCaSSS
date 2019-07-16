import numpy as np
import healpy as hp
from astropy.io import fits
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord

import sys
sys.path.insert(0, '../../GalfaCuber/code')
import galfa_vel_helpers as gvh

    
DR2_slice_root = "/data/seclark/GALFADR2/Wide_maps/"
outhdr= fits.getheader(DR2_slice_root+"GALFA_HI_W_S1024_V0000.4kms.fits")

ny = 2432
nx = 21600
gradmag_tot = np.zeros((ny, nx), np.float_)   

count = 0 

for _vel in np.arange(958, 1093):
    
    galfafn = galfaroot+gvh.get_galfa_W_name(_vel)
    galfaslice = fits.getdata(galfafn)
    
    gradslice = np.gradient(galfaslice)
    gradmag = np.sqrt(gradslice[0, :, :]**2 + gradslice[1, :, :]**2)
    
    gradmag_tot += gradmag
    
    if count == 0:
        start_vel = copy.copy(_vel)
    if count == 6:
        stop_vel = copy.copy(_vel)
        outfn = DR2_slice_root+"GALFA_HI_W_sumgradmags_vels{}_to_{}.fits".format(start_vel, stop_vel)
        fits.writeto(outfn, gradmag_tot, outhdr)
        count = 0
        gradmag_tot = np.zeros((ny, nx), np.float_)
    else:
        count += 1
    