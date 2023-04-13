# README

PyRIS - Python RIvers by Satellite

### Extract River Features from Landsat Multispectral Data

* Contains: pyris (module), pyris (cli script)
* Version 3.1.1 (development)
* updated 2023/04/13

### Who do I talk to?

* Original author: Federico Monegaglia (f.monegaglia@gmail.com)
* Original repo: https://github.com/fmonegaglia/pyris
* Current repo owner: Riccardo Bonanomi (riccardo.bonanomi@unitn.it)
* Current repo: https://github.com/riccardobonanomi/pyris3

### Repository branch description
- the master branch contains the up-to-date working code
- the original branch contains the original code converted to python 3

### Version tracking
- 1.0   (2018) 	 -> PyRIS deployment
- 3.0.0 (2022)       -> Converted to python 3
- 3.1.0 (2022)       -> Added GEE mask analysis capability
- 3.1.1 (2023/03/20) -> Added ability to skip masks when segmenting Landsat or importing GEE mask

-------------
### Install Instructions

#### On Unix/Linux

1) Install pip:
     sudo apt-get install python-pip

2) Install gdal (https://mothergeo-py.readthedocs.io/en/latest/development/how-to/gdal-ubuntu-pkg.html)

3) Install dependencies for PyRIS:
     sudo pip install numpy scipy matplotlib scikit-image gdal

4) Execute setup script to install PyRIS:
     sudo python setup.py install

#### On Windows

1) Install Anaconda (https://www.continuum.io/downloads)

2) Add Python to your Path (https://docs.python.org/2/using/windows.html)

3) Install GDAL 
     - through Conda Prompt
          conda install gdal
     - (https://sandbox.oarc.ucla.edu/tutorials/installing-gdal-for-windows)

### Run Instruction
* Call pyris from the command line:
      pyris [args]
* 'pyris --help' provides information on the required arguments