"""
Module containing main routines to execute SimpleGrid.

"""


import argparse
import configparser
import glob
from subprocess import call

import profiles_fixedgap
import tools
import climatology


def get_args():
    """ Get arguments from command line. """
    
    parser = argparse.ArgumentParser(
        description='A package for simple aggregation of ocean profile data onto a regular grid')
    parser.add_argument(
        'namelist', type=str, help='Path to namelist.ini')
    args = parser.parse_args()

    return args


def get_namelist(args):
    """ Read config options into <ConfigParser> object """
    config = configparser.ConfigParser()
    config.read(args.namelist)
    
    return config


def main():
    """ Read config and run """
    
    # Read config
    args = get_args()
    config = get_namelist(args)
    clim = climatology.GridClim(config)
    minyr = config.getint('profiles', 'minyr')
    maxyr = config.getint('profiles', 'maxyr')
    maxmonth = config.getint('profiles', 'maxmonth')
    minmonth = config.getint('profiles', 'minmonth')
    mupdate = config.getboolean('profiles', 'monthly_update')
    gridfiles = []
    dts, fnames = tools.get_dt_files(config, minyr, maxyr, minmonth, maxmonth,
      mupdate)
    #cminyr = config.getint('climatology', 'minyr')
    #cmaxyr = config.geting('climatology', 'maxyr')
    #cdts, cfnames = tools.get_dt_files(config, cminyr, cmaxyr)
    
    # Grid data
    for dt, fname in zip(dts, fnames):
        prof = profiles_fixedgap.Profiles(config, fname, dt)
        prof.grid_data()
        prof.write_grid()
        gridfiles.append(prof.fout)
        clim.accumulate_profiles(prof)
     
    # Calculate monthly climatology
    if config.getboolean('climatology', 'calc_climatology'):
        print('Got up to here')
        clim.calc_clim()
        clim.write_clim()
     
    # Calculate anomalies
    if config.getboolean('anomalies', 'calc_anomalies'):
        for dt, gridfile in zip(dts, gridfiles):
            tools.calc_anom(gridfile, dt, clim)
       
        # Move the anomalies to the anomaly folder:
        anomdir = config.get('anomalies','anomdir')
        profdir = config.get('grid','dir')
        files = glob.glob(profdir + '*_anom_*')
        for f in files:
            call('mv ' + f + ' ' + anomdir, shell = True)

    # Finished
    print('\nFinished!\n')

        
