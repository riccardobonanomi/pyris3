#!/usr/bin/env python

# .:==============:.
# .:    PyRIS     :.
# .: -------------:.
# .: Setup Script :.
# .:==============:.

# Python - RIvers from Satellites

from distutils.core import setup
import platform

pyris_script = 'pyris/bin/pyris'

# Proper Windows installation
if platform.system() == 'Windows':
    import shutil
    pyris_win_script = 'pyris/bin/pyrisw.py'
    shutil.copyfile(pyris_script, pyris_win_script)
    pyris_script = pyris_win_script

description = 'PyRIS :: Python - RIvers by Satellites'
long_description = '\n'.join((
    description,
    '''

    See Publication:
    
    --------------------------------------------------------------------------------------------------------
    Monegaglia et al., 2018
    "Automated extraction of meandering river morphodynamics from multitemporal remotely sensed data",
    Environmental Modeling & Software [https://www.sciencedirect.com/science/article/pii/S1364815217309118]
    --------------------------------------------------------------------------------------------------------

    Requires: NumPy, SciPy, MatPlotLib, Scikits-Image, GDAL, imagecodecs
    
    Updated to python3 and with a first gee introduction by Riccardo Bonanomi

    '''
))

setup(
    name = 'pyris',
    version = '3.1.0',
    author = 'Federico Monegaglia and Riccardo Bonanomi',
    author_email = 'f.monegaglia@gmail.com',
    maintainer = 'Riccardo Bonanomi',
    maintainer_email = 'riccardo.bonanomi@unitn.it',
    description = description,
    long_description = long_description,
    url = 'https://github.com/riccardobonanomi/pyris3',
    install_requires = [ 'numpy', 'scipy', 'matplotlib', 'scikit-image', 'gdal', 'imagecodecs' ],
    packages = [ 'pyris', 'pyris.config', 'pyris.misc', 'pyris.raster', 'pyris.vector' ],
    scripts = [pyris_script],
    
)
    
