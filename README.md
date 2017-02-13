This file gives instructions which will hopefully enable use of the Metbot
methodology described in:
Neil C. G. Hart, Chris J. C. Reason, and Nicolas Fauchereau, 2012: 
Building a Tropical–Extratropical Cloud Band Metbot.  Mon. Wea. Rev., 140,
4005–4016. doi: http://dx.doi.org/10.1175/MWR-D-12-00127.1

This metbot identifies cloud bands over southern Africa and possibly elsewhere
by using image segmentation software.

The python modules implementing the method are available in the MetBot
package. Changes have been made to use scikit-image instead of the originally
used cvblob C++ code.

Dependencies:
1- Numpy, Scipy, matplotlib, basemap
2- scikit-image
3- Also netCDF4 is needed for opening nc files
4- pyclimate
5- windspharm
OPT (don't think needed)- Python Image Processing Library (PIL)
NOTE: This software has been developed in Python 2.7 and is therefore unlikely
to work in Python 3.\*.

Recommended to use with Anaconda python distribution for easy installation of
dependencies. https://www.continuum.io/downloads
This software is known to work with the Linux 64-bit installer for anaconda.
You might however need an academic license (free) or fully paid 
license in order to install all dependencies (eg. scikit-image, opencv). The
academic license is simply available by registering for an anaconda account
at anaconda.org, and then placing your new license file in the directory 
/home/USERNAME/.continuum/. But try without this first.


Follow anaconda install instructions then install dependencies:
conda update
conda install scikit-image
conda install basemap
conda install netcdf4

Finally the best likely way to get your own copy of this code is to fork it to your github account. Alternatively you can clone the directory to your local machine. Forking is preferable as it will then be possible for your improvements to be incorporated into the main project.

Finally you should add path to  add these lines to your .bashrc file (or just execute them at the command line)
Now run test.py in test folder in order to check all is working. Do this as follows:
1- You will need to download the NOAA olr file if don't already have this. Running the getOLRdata.sh file will sort you out.
2- Then open python using "ipython --pylab" is probably easiest
3- Ensure you are within the satttmetbot/test directory
4- Run test.py in ipython with 
          --> run test.py

To see what is actually going on by plotting imagery, change:
showblb=False to showblb=True
To have control of how fast the algorithm proceeds, change the keyword argument
on line 57 of test.py to interact=True.
