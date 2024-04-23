#!/usr/bin/env python

# .:==============:.
# .:    PyRIS     :.
# .: -------------:.
# .: Setup Script :.
# .:==============:.

# Python - RIvers from Satellites

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
import platform
from pyris.info import __name__, __fullname__, __version__, __url__, __author__, __author_email__
from pyris.info import __maintainer__, __maintainer_email__, __license__, __long_description__

pyris_script = 'pyris/bin/pyris'

# Proper Windows installation
if platform.system() == 'Windows':
    import shutil
    pyris_win_script = 'pyris/bin/pyrisw.py'
    shutil.copyfile(pyris_script, pyris_win_script)
    pyris_script = pyris_win_script

description = ('%s v%s :: %s' % (__name__, __version__, __fullname__))
long_description = __long_description__

setup(
    name = __name__,
    version = __version__,
    author = __author__,
    author_email = __author_email__,
    maintainer = __maintainer__,
    maintainer_email = __maintainer_email__,
    license = __license__,
    description = description,
    long_description = long_description,
    url = __url__,
    install_requires = [ 'numpy', 'scipy', 'matplotlib', 'scikit-image', 'gdal', 'imagecodecs' ],
    packages = [ 'pyris', 'pyris.config', 'pyris.misc', 'pyris.raster', 'pyris.vector' ],
    scripts = [pyris_script],
    
)
    
