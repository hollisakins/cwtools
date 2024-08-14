import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord

# import the cutout function
from cwtools.cutout import cutout

# three positions I picked at random
positions = SkyCoord(*tuple(np.array([[150.2613374, 2.3162508],
                                      [149.8540891, 2.1146467],
                                      [150.0282134, 2.5518383]]).T), unit='deg')
size = 10*u.arcsec

# generate cutouts in f277w, saving them to the cutouts directory 
cuts = cutout('f277w', positions, size, 
              outdir='./cutouts/',
              names=['ex_A','ex_B','ex_C'], verbose=True)

import matplotlib.pyplot as plt
from astropy.visualization import simple_norm

fig, ax = plt.subplots(1, 3, figsize=(12,4), sharex=True, sharey=True, constrained_layout=True)
ax[0].imshow(cuts[0].data, extent=cuts[0].extent, norm=simple_norm(cuts[0].data, 'log'), cmap='Greys')
ax[1].imshow(cuts[1].data, extent=cuts[1].extent, norm=simple_norm(cuts[1].data, 'log'), cmap='Greys')
ax[2].imshow(cuts[2].data, extent=cuts[2].extent, norm=simple_norm(cuts[2].data, 'log'), cmap='Greys')
[a.set_xlim(-5, 5) for a in ax]
[a.set_ylim(-5, 5) for a in ax]
plt.show()