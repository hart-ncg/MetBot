# This is the initiation script to import the Automated Meterologist routines
'''MetBot package consists of these modules:
   MetBlobs
   SynopticAnatomy
   EventStats

   These are dependant on:
   MY OWN PYTHON MODULES:   mytools
                            mynetcdf
   OTHER MODULES:           matplotlib
                            mpl_toolkits.basemap
                            ScientificPython
                            numpy
                            scipy

        IMPORTANT: Information about data formats:
                   Two main formats of input data are used
                   1: Gridded data sets (eg. interpolated observations, reanalysis, model output)
                      DIMENSIONS ARE ASSUMED AND TREATED AS
                      TIME X LAT X LON
                   2: Rainfall station data sets
                      DIMENSIONS ARE ASSUMED AND TREATED AS
                      STATIONNO. X TIME (yes, make sense that time should be first as above,
                                        but I started it this way and haven't had a chance
                                        yet to got through the code and change this)

                   Many routines search for the time slice of interest and
                   so it is vital that data as TIME dimension in appropriate place'''


# BUILT-INS
#import os
#import datetime
#from time import time as timer
#import cPickle
#import glob
#import shutil
## AVAILABLE IN REPOSITORIES
#import numpy as np
#import Image
#import matplotlib.pyplot as plt
#from matplotlib.nxutils import points_inside_poly
#import cv
## BUILT LOCALLY
#import mpl_toolkits.basemap as bm
#from mpl_toolkits.basemap import num2date
#import cvblob as cvb
## MY OWN
#import mytools as my
#import mynetcdf as mync
#
## MAIN MODULES IN THIS PACKAGE
#import MetBlobs as blb
#import SynopticAnatomy as sy
#import EventStats as stats
