# Python tools for COSMOS-Web

This package provides handy python tools specific to the 
COSMOS-Web JWST program. 

The `cutout` module provides a method for programmatic generation 
of cutouts of mosaics hosted on an external server, in this case 
CANDIDE. See `examples/example_cutout.py` for a simple use-case. 

## Installation

To install, I recommend cloning this github repository onto your 
local machine and saving it somewhere you'll remember. Then, install 
using `pip` with the `-e` flag, to make it editable---that way, if 
the code updates, you can update your installation with `git pull`. 

```
git clone https://github.com/hollisakins/cwtools.git
cd cwtools
pip install -e . 
```