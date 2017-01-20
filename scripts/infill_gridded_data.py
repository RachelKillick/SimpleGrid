#!/usr/bin/env python2.7
"""
Infill gridded data using an inverse-distance weighting approach.

"""

from netCDF4 import Dataset 
import shutil
import numpy as np
import argparse
import ConfigParser
import iris.analysis.cartography as cart


class ShapeError(Exception):
    pass


def print_progress(task_name, nmax, n, nbar=20):
    """ Print progress to standard out. """
    done = nbar * '|'
    todo = nbar * '.'
    flt_pct = 100. * np.float(n)/nmax
    progind = np.int(flt_pct)/(100/nbar)
    progbar = done[:progind] + todo[progind:]
    
    print ('\r%25s: %s %6.2f%%' %
          (task_name, progbar, flt_pct)),
    
    if np.int(flt_pct) == 100:
        print ''


def get_args():
    """ Return arguments from command line """
    parser = argparse.ArgumentParser(
        description='Infill missing data in gridded data using an inverse-distance weighting approach.')
    parser.add_argument('namelist', type=str, help='Path to namelist.ini file used to create gridded data.')
    parser.add_argument('gridf', type=str, help='Path to netcdf file containing gridded data.')
    parser.add_argument('infillf', type=str, help='Output file containing infilled data.')
    parser.add_argument('--bathyf', type=str, help='Path to netcdf file containing bathymetry data for estimation of land-sea mask', default='None')
    parser.add_argument('--bathyvar', type=str, help='Variable name for bathymetry data.', default='None')  
    parser.add_argument('-c', '--cutoff', type=float, help='Maximum distance (km) used for interpolation.', default=1000.) 
    parser.add_argument('-p', '--power', type=float, help='Power used for inverse-distance weighting.', default=2.)
    parser.add_argument('-s', '--smooth', type=str, help='If True, apply inverse-distance weighting everywhere to smooth field.', default='False')
    args = parser.parse_args()

    return args


def get_namelist(args):
    """ Read config options into <ConfigParser> object """
    config = ConfigParser.ConfigParser()
    config.read(args.namelist)
    
    return config


def read_bathymetry(args):
    """ Read bathymetry data """
    nc = Dataset(args.bathyf)
    bathy = nc.variables[args.bathyvar][:].squeeze()    
    
    if bathy.ndim != 2:
        raise ShapeError('Bathymetry data must be 2D')
    
    nc.close()
    
    return bathy
    
    
def create_land_sea_mask(args, config):
    """ Create land sea mask from specified bathymetry data set """
    bathy = read_bathymetry(args)
    zbounds = np.array(config.get('grid', 'zbounds').split(','), dtype=np.float64)[:-1]
    zbounds = np.abs(zbounds) * -1.
    ny, nx = bathy.shape
    nz = zbounds.size
    lsmask = np.ones((nz,ny,nx)) == 1
    
    for n, zb in enumerate(zbounds):
        lsmask[n, bathy < zb] = False
    
    return lsmask


def reshape_latlon(lat, lon):
    """ Reshape latitude and longitude into 2D coordinates """ 
    lat_2d = np.ones((len(lat), len(lon))) * lat[:, np.newaxis]
    lon_2d = np.ones((len(lat), len(lon))) * lon[np.newaxis, :]  
    
    return lat_2d, lon_2d


def calc_distances(lons, lats, lon, lat):
    """ Calculate distances between all lats/lons to specified lat/lon """
    lons_rot, lats_rot = cart.rotate_pole(lons, lats, lon, lat)
    dists = 111. * (90 - lats_rot)
    
    return dists


def interp_idw(dat, lats, lons, lsmask, cutoff,
                power=2, smooth=False):
    """ Interpolate field using inverse-distance weighting """
    
    datmask = dat.mask
    filled = np.ma.MaskedArray(dat)
    nt, nz, ny, nx = filled.shape
    lats_2d, lons_2d = reshape_latlon(lats, lons)
    
    for ii in range(nx):
        print_progress('Interpolating', nx, ii+1, nbar=20)
        for ij in xrange(ny):
            if lsmask[0, ij, ii] == False:            
                lon = lons[ii]
                lat = lats[ij]
                dists = calc_distances(lons_2d, lats_2d, lon, lat)
                dist_ind = dists < cutoff
                for ik in xrange(nz):                    
                    for it in xrange(nt):
                        if (lsmask[ik, ij, ii] == False):
                            if (datmask[it, ik, ij, ii] == True) or (smooth == True):
                                loc_dat = dat[it, ik, dist_ind]
                                loc_dists = dists[dist_ind]
                                dat_ind = (loc_dat.mask == False)
                                loc_dat = loc_dat[dat_ind]
                                loc_dists = loc_dists[dat_ind]
                                inv_dists = 1. / (loc_dists ** power)
                                
                                if len(loc_dat) > 0:
                                    filled[it, ik, ij, ii] = (loc_dat * inv_dists).sum() / inv_dists.sum()
    # Mask filled field
    lsmask = (lsmask[np.newaxis] * np.ones(filled.shape)) == 1
    filled = np.ma.MaskedArray(filled, mask = (filled.mask | lsmask) )      
                            
    return filled


def infill(args, config):
    """ Infill missing data """
    
    # Create file for infill data
    cutoff = args.cutoff
    fname = args.gridf
    fout = args.infillf
    print 'Writing: %s' % fout
    shutil.copy(fname, fout)
    
    # Open file and read vars
    nc = Dataset(fout, 'r+')
    dat = nc.variables[config.get('profiles', 'datavar')]
    lons = nc.variables[config.get('profiles', 'xvar')]
    lats = nc.variables[config.get('profiles', 'yvar')]
    
    # Create land sea mask
    if args.bathyf != 'None':
        lsmask = create_land_sea_mask(args, config)
        if dat[0].shape != lsmask.shape:
            raise ShapeError('Data and land-sea mask must have the same shape!') 
    else:
        lsmask = np.zeros(dat[0].shape) == 1
               
    # Apply interpolation
    if args.smooth.lower() == 'false':
        smooth = False
    else:
        smooth = True

    dat[:] = interp_idw(dat[:], lats[:], lons[:], lsmask, cutoff, 
                        power=args.power, smooth=smooth)

    # Close file
    nc.close()


def main():
    """
    Infill missing data in gridded analysis using specified method.
        
    """ 
    args = get_args()
    config = get_namelist(args)
    infill(args, config)
    print '\nFinished\n'
        

if __name__ == '__main__':
    main()
