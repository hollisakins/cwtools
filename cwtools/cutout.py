'''This module handles cutout generation from mosaics hosted on CANDIDE.'''

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
import aiohttp, pwinput
import numpy as np
import tqdm, warnings, time, os
from astropy.wcs import FITSFixedWarning
warnings.simplefilter('ignore', FITSFixedWarning)
from . import utils

def get_url(band, tile, ext):
    base = 'https://exchg.calet.org/cosmosweb/.misc/'
    if band in ['f115w','f150w','f277w','f444w']:
        assert ext in ['sci']
        return f'{base}/CW_{band}_60mas_tot_v8.fits'
    elif band in ['cfht_u','uvista_Y','uvista_J','uvista_H','uvista_Ks']:
        return f'{base}/{band}_{ext}.fits'
    elif band in ['hsc_g','hsc_r','hsc_i','hsc_z','hsc_y']:
        return f'{base}/{band}_{ext}_{tile}.fits'
    else:
        raise ValueError(f"Value of band '{band}' not understood. Accepted values are: 'f115w','f150w','f277w','f444w',"+
                          "'uvista_Y','uvista_J','uvista_H','uvista_Ks','hsc_g','hsc_r','hsc_i','hsc_z','hsc_y','cfht_u'")

def cutout(band, positions, size, 
           ext='sci',
           pixscale='native',
           outdir=None,
           names=None,
           user=None,
           password=None,
           verbose=False):
    """ Generate cutouts from CANDIDE-hosted mosaics.

    Parameters
    ----------

    band : str
        Desired band for cutouts. Options currently implemented are: 
        f115w, f150w, f277w, f444w, uvista_Y, uvista_J, uvista_H, uvista_Ks, hsc_g, hsc_r, hsc_i, hsc_z, hsc_y, cfht_u

    positions : astropy.coordinates.SkyCoord
        Position or positions of object(s)

    size : astropy.Quantity or float
        Desired cutout size, provided as astropy Quantity with units or else assumed to be in arcsec

    ext : str - optional
        Extension (sci, err, wht). Defaults to 'sci'. Currently, only 'sci' images are available for NIRCam. 

    pixscale : str or float - optional
        Desired pixel scale for output cutouts. Defaults to 'native', i.e. adopt the native pixel scale
        for the target images: 60mas for NIRCam, 150mas UVISTA/CFHT, 168mas for HSC. If specified, the 
        cutouts will be reprojected into the requested pixel scale before output, using the reproject package. 

    outdir : str - optional
        Output directory for cutouts. Defaults to None, i.e. no file output. If the directory doesn't exist, 
        the code will generate it.

    names : list(str) - optional
        Names of source(s). Input should be a list of the same length as positions. If provided, output files, 
        will be named according to the input names. Otherwise, output files will be named with IAU-style 
        coordinate designations. 

    user : str - optional
        Username for CANDIDE COSMOS-Web server. If not provided, the code will check the CW_USER environment 
        variable first, then ask you to input the username. 
    
    password : str - optional
        Password for CANDIDE COSMOS-Web server. If not provided, the code will check the CW_USER environment 
        variable first, then ask you to input the password. 

    verbose : bool - optional
        Whether to print extra info (default: False)

    """
    

    if user and password:
        user, pswd = user, password
    else:
        user = os.getenv('CW_USER')
        pswd = os.getenv('CW_PSWD')
    if (user is None) or (pswd is None):
        print('No username/password found in kwargs or environment variables')
        user = input('Username: ')
        pswd = pwinput.pwinput('Password: ') 
    fsspec_kwargs = dict(client_kwargs={'auth': aiohttp.BasicAuth(user, pswd)})

    if hasattr(size, 'unit'):
        if size.unit != u.arcsec:
            if verbose: print('Converting cutout size to arcsec')
            size = size.to(u.arcsec)
    else:
        if verbose: print('Assuming cutout size in arcsec')
        size = size * u.arcsec

    if pixscale != 'native': # if you're going to reproject afterwards, pad the cutouts a bit 
        size_cut = size*1.2
    else:
        size_cut = size

    cutouts = []
    hdrs = []

    tiles = utils.get_tile(positions)

    for i in tqdm.tqdm(range(len(positions))):
        url = get_url(band, tiles[i], ext)

        with fits.open(url, use_fsspec=True, fsspec_kwargs=fsspec_kwargs) as hdul:  
            hdr = hdul[0].header
            wcs = WCS(hdr)
            cut = Cutout2D(hdul[0].section,
                            position=positions[i],
                            size=size_cut,
                            wcs=wcs)
        hdrs.append(hdr)
        cutouts.append(cut)


    ### If requested, reproject all cutouts onto the same pixel scale
    if pixscale != 'native':
        if verbose: print(f'Reprojecting cutouts to requested pixel scale of {pixscale} arcsec/pix')
        from reproject import reproject_interp
        for i in range(len(positions)):
            pos, cut = positions[i], cutouts[i]
            if not round(np.abs(cut.wcs.proj_plane_pixel_scales()[0]).to(u.arcsec).value*3600,2) == pixscale:

                wcs = WCS(naxis=2)
                wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
                wcs.wcs.crval = [pos.ra.value, pos.dec.value]
                wcs.wcs.cdelt = [-pixscale/3600,pixscale/3600]
                wcs.wcs.crpix = [size.value/pixscale/2, size.value/pixscale/2]

                cut_reproj, _  = reproject_interp((cut.data, cut.wcs), wcs, shape_out=(int(round(size.value/pixscale,0)),int(round(size.value/pixscale,0))))
                cut.data = cut_reproj
                cut.wcs = wcs


    for cut in cutouts:
        ps = cut.wcs.proj_plane_pixel_scales()[0].to(u.arcsec).value
        size = np.shape(cut.data)[0]
        extent = [-size*ps/2, size*ps/2, -size*ps/2, size*ps/2]
        setattr(cut, 'extent', extent)


    ### If requested, write cutouts to file(s) 
    if outdir:
        if verbose: print(f'Writing cutouts to {outdir}')
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
        if names:
            assert len(names)==len(positions)
        else:
            if verbose: print('No object names provided, writing files with IAU-style coordinate designations')
            names = [f'J{pos.ra.to_string(unit=u.hourangle, sep="", precision=2, pad=True)}{pos.dec.to_string(sep="", precision=2, alwayssign=True, pad=True)}' for pos in positions]
        for i in range(len(positions)):
            cut = cutouts[i]
            name = names[i]
            hdr = hdrs[i]
            hdr.update(cut.wcs.to_header())
            fits.writeto(f'{outdir}/cutout_{name}_{band}_{ext}.fits', cut.data, hdr, overwrite=True)

    return cutouts


