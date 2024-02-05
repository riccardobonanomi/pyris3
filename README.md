# README

PyRIS - Python RIvers by Satellite

### Extract River Features from Landsat Multispectral Data

* Contains: pyris (module), pyris (cli script)
* Version 3.1.3 (development)
* README updated 2024/01/18

### Who do I talk to?

* Original author: Federico Monegaglia (f.monegaglia@gmail.com)
* Original repo: https://github.com/fmonegaglia/pyris - no longer available
* Current repo owner: Riccardo Bonanomi (riccardo.bonanomi@unitn.it)
* Current repo: https://github.com/riccardobonanomi/pyris3

### Repository branch description
- the master branch contains the up-to-date working code
- the original branch contains the original code converted to python 3

### Version tracking
- 1.0   (2018) 	 by Federico Monegaglia -> PyRIS deployment
- 3.0.0 (2022)       by Riccardo Bonanomi   -> Converted to python 3
- 3.1.0 (2022)       by Riccardo Bonanomi   -> Added GEE mask analysis capability
- 3.1.1 (2023/03/20) by Riccardo Bonanomi   -> Added ability to skip masks when segmenting Landsat or importing GEE mask
- 3.1.2 (2023/04/19) by Riccardo Bonanomi   -> Added ability to draw black masks on gee mask
- 3.1.3 (2024/01/18) by Riccardo Bonanomi   -> Fixed bug in migration evaluation
- 3.1.4 (2024/02/05) by Riccardo Bonanomi   -> Added skipping capability of non .tif files in external masks

-------------
## Fast guide
This is just a brief installation and usage guide, a more comprehensive documentation can be found in the **PyRISmanual.pdf** file at https://github.com/riccardobonanomi/pyris3/blob/main/PyRISmanual.pdf (note that it is still under development).

### Install Instructions

#### On Unix/Linux

1) Install pip:
     sudo apt-get install python-pip

2) Install gdal (https://mothergeo-py.readthedocs.io/en/latest/development/how-to/gdal-ubuntu-pkg.html)

3) Install dependencies for PyRIS:
     sudo pip install numpy scipy matplotlib scikit-image gdal

4) To install PyRIS from the pyris folder run:
     pip install .

#### On Windows - refer to the manual

### Run Instruction
* Call pyris from the command line:
      pyris [args]
* 'pyris --help' provides information on the required arguments
