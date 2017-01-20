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
        
    def accumulate_profiles(self, profiles):
        """ Accumulate data from profiles that match date constraints """
        dt = profiles.dt
        mon, yr = dt.month, dt.year
    
        if (yr >= self.minyr) & (yr <= self.maxyr):
            if not self.initialized:
                self.init_grids(profiles)
                self.initialized=True
                
            self.grid_count[mon-1] += profiles.grid_count
            self.grid_sum[mon-1] += profiles.grid_sum
            
    def calc_clim(self):
        """ Calculate monthly climatological means"""
        self.grid_mean = self.grid_sum / self.grid_count
        self.grid_mean = np.ma.MaskedArray(self.grid_mean, mask = (self.grid_count == 0))

    def create_savename(self):
        """ Generate file name based on file name and grid specification """
        climf = self.config.get('grid', 'dir') + self.config.get('profiles', 'fpattern')
        climf = climf.replace('${YYYY}${MM}', 'climatology_%i-%i' % (self.minyr, self.maxyr))
        newsuffix = '_gridded_%ix%ix%i.nc' % (self.nx, self.ny, self.nz)
        climf = climf.replace('.nc', newsuffix)
        
        return climf

    def write_clim(self):
        """ Write gridded climatology to netcdf  """
        self.fout = self.create_savename()
        ncout = Dataset(self.fout, 'w')
        print 'Writing: %s' % self.fout
        
        # Create dimensions
        xdim = ncout.createDimension(self.xvar, self.nx)
        ydim = ncout.createDimension(self.yvar, self.ny)
        zdim = ncout.createDimension(self.zvar, self.nz) 
        mdim = ncout.createDimension('month', 12) 

        # Create variables
        varx = ncout.createVariable(self.xvar, 'float64', (self.xvar,))
        varxmin = ncout.createVariable('%s_minbounds' % (self.xvar), 'float32', (self.xvar,))
        varxmax = ncout.createVariable('%s_maxbounds' % (self.xvar), 'float32', (self.xvar,))
        vary = ncout.createVariable(self.yvar, 'float64', (self.yvar,))
        varymin = ncout.createVariable('%s_minbounds' % (self.yvar), 'float32', (self.yvar,))
        varymax = ncout.createVariable('%s_maxbounds' % (self.yvar), 'float32', (self.yvar,))
        varz = ncout.createVariable(self.zvar, 'float64', (self.zvar,))
        varzmin = ncout.createVariable('%s_minbounds' % (self.zvar), 'float32', (self.zvar,))
        varzmax = ncout.createVariable('%s_maxbounds' % (self.zvar), 'float32', (self.zvar,))
        varmean = ncout.createVariable(self.datavar, 'float32', ('month',self.zvar,self.yvar,self.xvar))
        varsum = ncout.createVariable('sum', 'float32', ('month',self.zvar,self.yvar,self.xvar))
        varcount = ncout.createVariable('count', 'float32', ('month',self.zvar,self.yvar,self.xvar))
        varmon = ncout.createVariable('month', 'int32', ('month',))

        # Write to variables
        varx[:] = self.xgrid
        varxmin[:] = self.xminbounds
        varxmax[:] = self.xmaxbounds
        vary[:] = self.ygrid
        varymin[:] = self.yminbounds
        varymax[:] = self.ymaxbounds
        varz[:] = self.zgrid
        varzmin[:] = self.zminbounds
        varzmax[:] = self.zmaxbounds
        varmean[:] = self.grid_mean
        varsum[:] = self.grid_sum
        varcount[:] = self.grid_count
        varmon[:] = np.arange(12) + 1
        
        # Add  global attributes
        ncout.history = 'Created ' + time.ctime(time.time())
        
        # Save
        ncout.close()  
    