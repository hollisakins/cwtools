from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
import time, os
import aiohttp
import pwinput
import numpy as np


def get_tile(coords_in):
    """ Get the tile for a given coordinate.
    """
    coords_in = SkyCoord(coords_in)
    tiles = np.array(['A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','B1','B2','B3','B4','B5','B6','B7','B8','B9','B10'])
    ra = np.array([149.83060849224526,149.96618456679428,150.10175329130632,150.23731326578584,150.37286289023743,149.7647850185184,149.9003454431693,150.03589926778474,150.17144509237895,150.30698131694916,149.96220693940455,150.09747486414722,150.23328673819674,150.36919886197347,150.50453028662352,149.89657956583684,150.03225868964083,150.16766126474636,150.30327548908838,150.4387048134958])
    dec = np.array([2.2105308808257207,2.1612435187006245,2.111943456573772,2.062631644446241,2.013309032319107,2.0296854374792526,1.9804026753817368,1.9311082632823562,1.8818031511845486,1.8324881890870612,2.5720372673973584,2.5230015052794137,2.473525843055812,2.423885930793452,2.374833668674315,2.3912872741259417,2.341569861870867,2.292723449827574,2.2433903376662596,2.1941086255207374])
    coords_cen = SkyCoord(ra, dec, unit='deg')
    idx, d2d, d3d = coords_in.match_to_catalog_sky(coords_cen)
    return tiles[idx]



def get_url(band, tile, ext, ps):
    base_nircam = 'https://exchg.calet.org/cosmosweb/COSMOS-Web_Jan24/NIRCam/v0.8/sci_mosaics/'
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


    tiles = get_tile(positions)
    unique_tiles, indices = np.unique(tiles, return_inverse=True)
    
    start = time.time()

    cutouts = []
    
    for tile in unique_tiles:
        url = get_url(band, tile, ext, '')

        with fits.open(url, use_fsspec=True, fsspec_kwargs=dict(client_kwargs={'auth': aiohttp.BasicAuth(user, pswd)})) as hdul:  
            wcs = WCS(hdul[0].header)

            for position in positions:
                cut = Cutout2D(hdul[0].section,
                                position=position,
                                size=size,
                                wcs=wcs)

    end = time.time()
    print(f'Completed in {end-start:.2f}s')

    ##### reproject CW orient NIRCam/MIRI/ACS bands to N-up if requested
    ##### reproject N-up orient ground-based bands to CW if requested

    return cut




positions = [SkyCoord('10h00m26.358s','+02d15m26.92300s'),SkyCoord('10h00m26.358s','+02d15m26.92300s')]
size = 5*u.arcsec
cutouts = cutout('f444w',positions,size)

import matplotlib.pyplot as plt
plt.style.use('hba_sans')
import numpy as np

plt.imshow(cut.data)

plt.show()