'''This module includes useful helper functions for COSMOS-Web'''

from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from regions import Regions
import os
import numpy as np

install_dir = os.path.dirname(os.path.realpath(__file__))
'''Stores the install directory for easy reference'''

def get_tile(coords):
    """
    Get the tile for a given coordinate. Finds the closest match 
    to the center coordinates of the tiles. Returns None for 
    coordinates outside the COSMOS-Web footprint. 

    Parameters
    ----------

    coords : astropy.coordinates.SkyCoord or np.ndarray(N,2)
        Coordinates of input sources. Can be specified as either 
        a single SkyCoord object, a list of SkyCoord objects, or a 
        numpy array with 2 columns giving RA and Dec. 
    """
    
    if type(coords)==np.ndarray:
        coords = SkyCoords(ra=coords[:,0], dec=coords[:,1], unit='deg')
    elif type(coords)==list:
        coords = SkyCoord(coords)

    tiles = np.array(['A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','B1','B2','B3','B4','B5','B6','B7','B8','B9','B10'])
    ra = np.array([149.83060849224526,149.96618456679428,150.10175329130632,150.23731326578584,150.37286289023743,149.7647850185184,149.9003454431693,150.03589926778474,150.17144509237895,150.30698131694916,149.96220693940455,150.09747486414722,150.23328673819674,150.36919886197347,150.50453028662352,149.89657956583684,150.03225868964083,150.16766126474636,150.30327548908838,150.4387048134958])
    dec = np.array([2.2105308808257207,2.1612435187006245,2.111943456573772,2.062631644446241,2.013309032319107,2.0296854374792526,1.9804026753817368,1.9311082632823562,1.8818031511845486,1.8324881890870612,2.5720372673973584,2.5230015052794137,2.473525843055812,2.423885930793452,2.374833668674315,2.3912872741259417,2.341569861870867,2.292723449827574,2.2433903376662596,2.1941086255207374])
    coords_cen = SkyCoord(ra, dec, unit='deg')
    idx, d2d, d3d = coords.match_to_catalog_sky(coords_cen)
    t = list(tiles[idx])
    for i in np.where(~check_in_cw(coords))[0]:
        t[i] = None

    return t

def check_in_cw(coords):
    """
    Check whether the given coordinate(s) fall within the 
    COSMOS-Web NIRCam footprint. 
    
    Parameters
    ----------

    coords : astropy.coordinates.SkyCoord
    """
    
    region = Regions.read(os.path.join(install_dir,'data/lw_outline.reg'))[0]
    # region.contains requires a WCS, so we'll create a simple 180mas WCS
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['RA---TAN','DEC--TAN']
    wcs.wcs.cdelt = [-180e-3/3600, 180e-3/3600]
    wcs.wcs.crval = [150.11632752531, 2.2009681511549] # COSMOS tangent point?
    wcs.wcs.crpix = [19785/2,19212/2]
    return region.contains(coords, wcs)


## more to come...