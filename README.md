# README #

PyRIS - Python RIvers by Satellite

### Extract River Features from Landsat Multispectral Data ###

* Contains: pyris (module), pyris (cli script)
* Version 2.0 (development)

### Who do I talk to? ###

* Repo owner or admin: Federico Monegaglia (f.monegaglia@gmail.com)

### Install Instructions ###

On Unix/Linux
-------------

1) Install pip:
     sudo apt-get install python-pip

2) Install dependencies for PyRIS:
     sudo pip install numpy scipy matplotlib scikit-image gdal

3) Execute setup script to install PyRIS:
     sudo python setup.py install

On Windows
----------

1) Install Anaconda (https://www.continuum.io/downloads)

2) Add Python to your Path (https://docs.python.org/2/using/windows.html)

3) Install GDAL through Conda Prompt
       conda install gdal

### Run Instruction ###
* Call pyris from the command line:
      pyris [args]
* 'pyris --help' provides information on the required arguments
