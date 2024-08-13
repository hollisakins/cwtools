from cwtools.cutout import cutout

import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

size = 5*u.arcsec
positions = SkyCoord(*tuple(np.array([['10h00m26.358s','+02d15m26.92300s'],
                                      ['10h02m26.358s','+02d15m26.92300s']]).T))

cuts = cutout('f814w', positions, size)

import matplotlib.pyplot as plt
plt.style.use('hba_sans')
plt.imshow(cuts.data)
plt.show()