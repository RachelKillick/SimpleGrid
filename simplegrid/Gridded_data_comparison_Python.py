""" Code to test whether my gridded data match between the original Python code
only going down to one depth and my updated Python code designed to calculate
OHC for multiple depths at once.

"""


import glob
import iris
import matplotlib.pyplot as plt
import numpy as np
import os


# The locations of the files to be compared:
iris_fpath1 = '/scratch/rkillick/BAMS/Simple_grid_profs/Python3/Full_take2/'
iris_fpath2 = '/scratch/rkillick/BAMS/Simple_grid_profs/Python3/Testing_depths200/'

# Get lists of both sets of files:
ifiles1 = sorted(os.listdir(iris_fpath1))
ifiles2 = sorted(os.listdir(iris_fpath2))

# Now compare the files - printing a statemtent if masks don't 
# match or if the biggest difference between them is > 1e-6

for f in range(len(ifiles2)):
    pycons1 = iris.load(iris_fpath1 + ifiles1[f], 'CONS_TEMP')[0][0][0]
    pycons2 = iris.load(iris_fpath2 + ifiles2[f], 'CONS_TEMP')[0][0][0]
    pypcnt1 = iris.load(iris_fpath1 + ifiles1[f], 'pcount')[0][0][0]
    pypcnt2 = iris.load(iris_fpath2 + ifiles2[f], 'pcount')[0][0][0]
   
    if (pycons1.data.mask == pycons2.data.mask).all() != True:
        print('Masks do not match for CONS {0}'.format(ifiles1[f][24:30]))
        print(np.shape(pycons1.data.mask[pycons1.data.mask == False]))
        print(np.shape(pycons2.data.mask[pycons2.data.mask == False]))
    if abs(pycons1.data.data[pycons1.data.mask == False] - 
      pycons2.data.data[pycons1.data.mask == False]).max() > 0:
        print('Data discrepancy greater than 2e-6 for CONS {0}'.format(ifiles1[f][24:30]))
        print(abs(pycons1.data.data[pycons1.data.mask == False] - 
          pycons2.data.data[pycons1.data.mask == False]).max())
    if abs(pypcnt1.data - pypcnt2.data).max() > 0:
        print('Data discrepancy greater than 0 for PCOUNT {0}'.format(ifiles1[f][24:30]))
        print(abs(pypcnt1.data.data[pypcnt1.data.mask == False] - 
          pypcnt2.data.data[pypcnt1.data.mask == False]).max())

# Data currently match exactly which is good news!
