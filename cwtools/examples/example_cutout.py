import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord

from cwtools.cutout import cutout

size = 10*u.arcsec
positions = SkyCoord(*tuple(np.array([[150.2613374, 2.3162508],
                                      [149.8540891, 2.1146467],
                                      [150.0282134, 2.5518383]]).T), unit='deg')

cuts = cutout('f277w', positions, size, 
              pixscale=0.06, 
              outdir='./cutouts/',
              names=['exA','exB','exC'])

import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(12,4), sharex=True, sharey=True, constrained_layout=True)

from astropy.visualization import simple_norm
ax[0].imshow(cuts[0].data, extent=cuts[0].extent, norm=simple_norm(cuts[0].data, 'log'), cmap='Greys')
ax[1].imshow(cuts[1].data, extent=cuts[1].extent, norm=simple_norm(cuts[1].data, 'log'), cmap='Greys')
ax[2].imshow(cuts[2].data, extent=cuts[2].extent, norm=simple_norm(cuts[2].data, 'log'), cmap='Greys')

for axi in ax:
    axi.set_xlim(-5, 5)
    axi.set_ylim(-5, 5)

plt.show()