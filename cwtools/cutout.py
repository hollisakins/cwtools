from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
import time, os
import aiohttp
import pwinput
import numpy as np
import tqdm


def get_url(band, tile, ext, ps):
    base_nircam = 'https://exchg.calet.org/cosmosweb/COSMOS-Web_Jan24/NIRCam/v0.8/sci_mosaics/'
    base_acsA = 'https://exchg.calet.org/cosmosweb/COSMOS-Web_Apr23/ACS/'
    base_acsB = 'https://exchg.calet.org/cosmosweb/COSMOS-Web_Jan24/ACS/'
    base_misc = 'https://exchg.calet.org/cosmosweb/.misc/'
    if band in ['f115w','f150w','f277w','f444w']:
        assert ps in ['30mas','60mas']
        assert ext in ['sci','err','wht']
        return f'{base_nircam}/mosaic_nircam_{band}_COSMOS-Web_{ps}_{tile}_v0_8_{ext}.fits.gz' 
    if band in ['Y','J','H','Ks']:
        return f'{base_misc}/uvista_{band}_{ext}_{tile}.fits'
    if band in ['g','r','i','z','y']:
        return f'{base_misc}/hsc_{band}_{ext}_{tile}.fits'
    if band in ['u']:
        return f'{base_misc}/cfht_{band}_{ext}_{tile}.fits'
    if band in ['f814w']:
        if ext=='sci': ext='drz'
        if tile.startswith('A'): return f'{base_acsA}/mosaic_cosmos_web_2023apr_{ps}_tile_{tile}_hst_acs_wfc_f814w_{ext}.fits'
        if tile.startswith('B'): return f'{base_acsB}/mosaic_cosmos_web_2024jan_{ps}_tile_{tile}_hst_acs_wfc_f814w_{ext}.fits'

def cutout(band, positions, size, ext='sci', orient='CW', verbose=False, **kwargs):
    """ Generate cutouts from CANDIDE-hosted mosaics.

    Parameters
    ----------

    band : str
        What band you want

    positions : astropy.coordinates.SkyCoord or list
        Position or positions of object(s) 

    size : astropy.Quantity or float
        cutout size, either as astropy Quantity with units or assumed to be in arcsec

    ext : str - optional
        extension (sci, err, wht)
        default: sci

    orient : str - optional
        WCS orientation 
        'N' = default north-up orientation
        'CW' = cosmos-web orientation, 20 deg offset from N-up

    verbose : bool - option
        verbosity (default: False)
    """
    
    if ('user' in kwargs) and ('password' in kwargs):
        user = kwargs['user']
        pswd = kwargs['password']
    else:
        user = os.getenv('CW_USER')
        pswd = os.getenv('CW_PSWD')
    if (user is None) or (pswd is None):
        print('No username/password found in kwargs or environment variables')
        user = input('Username: ')
        pswd = pwinput.pwinput('Password: ') #getpass.getpass('Password: ')
    fsspec_kwargs = dict(client_kwargs={'auth': aiohttp.BasicAuth(user, pswd)})

    tiles = get_tile(positions)
    unique_tiles, indices = np.unique(tiles, return_inverse=True)
    
    start = time.time()

    cutouts = []
    
    for i in tqdm.tqdm(range(len(positions))):
    # for tile in unique_tiles:
        url = get_url(band, tiles[i], ext, '60mas')
        print(url)


        with fits.open(url, use_fsspec=True, fsspec_kwargs=fsspec_kwargs) as hdul:  
            wcs = WCS(hdul[0].header)
            cut = Cutout2D(hdul[0].section,
                            position=positions[i],
                            size=size,
                            wcs=wcs)

    end = time.time()
    print(f'Completed in {end-start:.2f}s')

    ##### reproject CW orient NIRCam/MIRI/ACS bands to N-up if requested
    ##### reproject N-up orient ground-based bands to CW if requested

    return cut

