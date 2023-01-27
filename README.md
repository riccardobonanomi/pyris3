# README #

PyRIS - Python RIvers by Satellite

### Extract River Features from Landsat Multispectral Data ###

* Contains: pyris (module), pyris (cli script)
* Version 3.0 (development)

### Who do I talk to? ###

* Original author: Federico Monegaglia (f.monegaglia@gmail.com)
* Original repo: https://github.com/fmonegaglia/pyris
* Repo owner or admin: Riccardo Bonanomi (riccardo.bonanomi@unitn.it)


### Install Instructions ###

On Unix/Linux
-------------

1) Install pip:
     sudo apt-get install python-pip

3) Install gdal (https://mothergeo-py.readthedocs.io/en/latest/development/how-to/gdal-ubuntu-pkg.html)

2) Install dependencies for PyRIS:
     sudo pip install numpy scipy matplotlib scikit-image gdal

4) Execute setup script to install PyRIS:
     sudo python setup.py install

On Windows
----------

1) Install Anaconda (https://www.continuum.io/downloads)

2) Add Python to your Path (https://docs.python.org/2/using/windows.html)

3) Install GDAL 
     - through Conda Prompt
       conda install gdal
     - (https://sandbox.oarc.ucla.edu/tutorials/installing-gdal-for-windows)

### Run Instruction ###
* Call pyris from the command line:
      pyris [args]
* 'pyris --help' provides information on the required arguments

### Updated to python 3 ###
Beta version by Riccardo Bonanomi (riccardo.bonanomi@unitn.it)