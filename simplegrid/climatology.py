"""
Module containing class to create a gridded climatology.

"""
from netCDF4 import Dataset
import numpy as np
import time


class GridClim(object):
    """ Class containing methods to generate a monthly gridded climatology """ 
    
    def __init__(self, config):
        """ Initialize using configuration options""" 
        self.config = config
        self.minyr = self.config.getint('climatology', 'minyr')
        self.maxyr = self.config.getint('climatology', 'maxyr')
        self.initialized = False
        
    def init_grids(self, profiles):
        """ Initialize climatology grid using existing <Profile> class object """ 
        self.xvar = profiles.xvar
        self.yvar = profiles.yvar
        self.zvar = profiles.zvar
        self.datavar = profiles.datavar
        self.xgrid = profiles.xgrid
        self.ygrid = profiles.ygrid
        self.zgrid = profiles.zgrid
        self.xminbounds = profiles.xminbounds
        self.xmaxbounds = profiles.xmaxbounds
        self.yminbounds = profiles.yminbounds
        self.ymaxbounds = profiles.ymaxbounds
        self.zminbounds = profiles.zminbounds
        self.zmaxbounds = profiles.zmaxbounds
        self.nx = profiles.nx
        self.ny = profiles.ny
        self.nz = profiles.nz
        self.grid_sum = np.zeros((12, self.nz, self.ny, self.nx))
        self.grid_count = np.zeros((12, self.nz, self.ny, self.nx))
        self.grid_meant = np.zeros((12, self.nz, self.ny, self.nx))
        self.grid_meansum = np.zeros((12, self.nz, self.ny, self.nx))
        self.grid_pcount = np.zeros((12, self.nz, self.ny, self.nx))
        self.grid_nonzerocount = np.zeros((12, self.nz, self.ny, self.nx))
        
    def accumulate_profiles(self, profiles):
        """ Accumulate data from profiles that match date constraints """
        dt = profiles.dt
        mon, yr = dt.month, dt.year
    
        if (yr >= self.minyr) & (yr <= self.maxyr):
            if not self.initialized:
                self.init_grids(profiles)
                self.initialized=True

            self.grid_meansum[mon-1] += profiles.grid_meansum
            self.grid_pcount[mon-1] += profiles.grid_pcount
            self.grid_nonzerocount[mon-1][np.where(profiles.grid_pcount != 0.)] += 1
            self.grid_count[mon-1] += profiles.grid_count
            self.grid_sum[mon-1] += profiles.grid_sum
            self.grid_meant[mon-1][np.where(profiles.grid_meantmean.mask == False)] += profiles.grid_meantmean[np.where(profiles.grid_meantmean.mask == False)]

    def calc_clim(self):
        """ Calculate monthly climatological means"""
        self.grid_mean = self.grid_sum / self.grid_count
        self.grid_mean = np.ma.MaskedArray(self.grid_mean, mask = (self.grid_count == 0))
        self.grid_pmean = self.grid_meant / self.grid_nonzerocount
        self.grid_pmean = np.ma.MaskedArray(self.grid_pmean, mask = (self.grid_nonzerocount == 0))
        self.grid_meansumo = self.grid_meansum / self.grid_pcount
        self.grid_meansumo = np.ma.MaskedArray(self.grid_meansumo, mask = (self.grid_pcount == 0))

    def create_savename(self):
        """ Generate file name based on file name and grid specification """
        climf = self.config.get('climatology', 'cldir') + self.config.get('profiles', 'fpattern')
        climf = climf.replace('${YYYY}${MM}', 'climatology_%i-%i' % (self.minyr, self.maxyr))
        newsuffix = '_gridded_%ix%ix%i.nc' % (self.nx, self.ny, self.nz)
        climf = climf.replace('.nc', newsuffix)
        
        return climf

    def write_clim(self):
        """ Write gridded climatology to netcdf  """
        self.fout = self.create_savename()
        ncout = Dataset(self.fout, 'w')
        print('Writing: {0}'.format(self.fout))
        
        # Create dimensions
        xdim = ncout.createDimension(self.xvar, self.nx)
        ydim = ncout.createDimension(self.yvar, self.ny)
        zdim = ncout.createDimension(self.zvar, self.nz) 
        mdim = ncout.createDimension('month', 12)
        bndsDim = ncout.createDimension('bnds', 2)

        # Create variables
        varx = ncout.createVariable(self.xvar, 'float64', (self.xvar,))
        vary = ncout.createVariable(self.yvar, 'float64', (self.yvar,))
        varz = ncout.createVariable(self.zvar, 'float64', (self.zvar,))
        varmeansum = ncout.createVariable('meansum', 'float32', ('month',self.zvar,self.yvar,self.xvar))
        varpcount = ncout.createVariable('pcount', 'float32', ('month',self.zvar,self.yvar,self.xvar))
        varcount = ncout.createVariable('count', 'float32', ('month',self.zvar,self.yvar,self.xvar))
        varsum = ncout.createVariable('sum', 'float32', ('month',self.zvar,self.yvar,self.xvar))
        varmeant = ncout.createVariable('meant', 'float32', ('month',self.zvar,self.yvar,self.xvar))
        varmean = ncout.createVariable('mean', 'float32', ('month',self.zvar,self.yvar,self.xvar))
        varpmean = ncout.createVariable(self.datavar, 'float32', ('month',self.zvar,self.yvar,self.xvar))
        varmeansumo = ncout.createVariable('meansumo', 'float32', ('month',self.zvar,self.yvar,self.xvar))
        varmon = ncout.createVariable('month', 'int32', ('month',))
        
        varx.standard_name = 'longitude'
        varx.units = 'degrees'
        ncout.variables['LONGITUDE'].bounds = 'lon_bnds'
        lonBndsVar = ncout.createVariable('lon_bnds', 'float64', (self.xvar, 'bnds'))
        xboundaries = np.concatenate([self.xminbounds, np.reshape(self.xmaxbounds[-1],(1,1))[0]])
        lonBndsVar[:,:] = np.array([xboundaries[:-1], xboundaries[1:]]).T

        vary.standard_name = 'latitude'
        vary.units = 'degrees'
        ncout.variables['LATITUDE'].bounds = 'lat_bnds'
        latBndsVar = ncout.createVariable('lat_bnds', 'float64', (self.yvar, 'bnds'))
        yboundaries = np.concatenate([self.yminbounds, np.reshape(self.ymaxbounds[-1],(1,1))[0]])
        latBndsVar[:,:] = np.array([yboundaries[:-1], yboundaries[1:]]).T
        
        varz.standard_name = 'depth'
        varz.units = 'metres'
        ncout.variables['DEPH_CORRECTED'].bounds = 'depth_bnds'
        depthBndsVar = ncout.createVariable('depth_bnds', 'float64', (self.zvar, 'bnds'))
        zboundaries = np.concatenate([self.zminbounds, np.reshape(self.zmaxbounds[-1],(1,1))[0]])
        depthBndsVar[:,:] = np.array([zboundaries[:-1], zboundaries[1:]]).T

        # Write to variables
        varx[:] = self.xgrid
        vary[:] = self.ygrid
        varz[:] = self.zgrid
        varmeansum[:] = self.grid_meansum
        varpcount[:] = self.grid_pcount
        varcount[:] = self.grid_count
        varsum[:] = self.grid_sum
        varmeant[:] = self.grid_meant
        varmean[:] = self.grid_mean
        varpmean[:] = self.grid_pmean
        varmeansumo[:] = self.grid_meansumo
        varmon[:] = np.arange(12) + 1
        
        # Add  global attributes
        ncout.history = 'Created ' + time.ctime(time.time())
        
        # Save
        ncout.close()  
    
