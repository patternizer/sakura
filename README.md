![image](https://github.com/patternizer/glosat-sakura/blob/master/glosat-sakura.png)

# glosat-sakura

Python code to plot the phenological timeseries for Prunus jamasakura in Kyoto City and overlay PAGES2k (& HadCRUT 2001+) temperature anomalies as part of ongoing work for the [GloSAT](https://www.glosat.org) project: www.glosat.org. 

## Contents

* `glosat-sakura.py` - python code to extract the timeseries and perform the paleoclimate timeseries plotting
* `glosat-sakura.png` - data graphic

The first step is to clone the latest glosat-sakura code and step into the check out directory: 

    $ git clone https://github.com/patternizer/glosat-sakura.git
    $ cd glosat-sakura

### Using Standard Python

The code should run with the [standard CPython](https://www.python.org/downloads/) installation and was tested in a conda virtual environment running a 64-bit version of Python 3.8+.

glosat-sakura scripts can be run from sources directly, once the required data dependencies are satisfied.

Run with:

    $ python glosat-sakura.py

## License

The code is distributed under terms and conditions of the [Open Government License](http://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/).

## Contact information

* [Michael Taylor](michael.a.taylor@uea.ac.uk)

