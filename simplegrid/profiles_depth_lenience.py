""" Module containing class to load and aggregate profile data """ 

from netCDF4 import Dataset, date2num
import numpy as np
import scipy.stats
import inst_type
import time
import warnings

class ShapeError(Exception):
    pass

class Profiles(object):
    """ Class containing methods to read and manipulate profile data """ 
    
    def __init__(self, config, fname, dt, preload=True):
        """ Initialize using configuration options """ 
        
        self.config = config
        self.fname = fname
        self.dt = dt
        self.xvar = config.get('profiles', 'xvar')
        self.yvar = config.get('profiles', 'yvar')
        self.zvar = config.get('profiles', 'zvar')
        self.zbounds = np.array(self.config.get('grid', 'zbounds').split(','), dtype=np.float64)
        self.pvar = config.get('profiles', 'pvar')
        self.pnvar = config.get('profiles', 'pnvar')
        self.irvar = config.get('profiles', 'irvar')
        self.psvar = config.get('profiles', 'psvar')
        self.psalqcvar = config.get('profiles', 'psalqcvar')
        self.qcvar = config.get('profiles', 'qcvar')
        self.posqcvar = config.get('profiles', 'posqcvar')
        self.datavar = config.get('profiles', 'datavar')
        self.fixedgap = config.get('profiles', 'fixedgap')
        
        if preload: 
            self.load_data()
            self.load_x()
            self.load_y()
            self.load_z()
            self.load_p()
            self.load_pn()
            self.load_ir()
            self.load_ps()
            self.load_psalqc()
            self.load_qc()
            self.load_posqc()
        
    def read_var(self, ncvar):
        """ Read data from specified variable """
        ncf = Dataset(self.fname)
        dat = ncf.variables[ncvar][:]
        ncf.close()
        return dat   
                       
    def load_x(self):
        """ Load x-coordinate data as <np.array> """
        self.x = self.read_var(self.xvar)
        self.test_shape(self.xvar, self.x.shape, 1)
    
    def load_y(self):
        """ Load y-coordinate data as <np.array> """
        self.y = self.read_var(self.yvar)
        self.test_shape(self.yvar, self.y.shape, 1)

    def load_z(self):
        """ Load z-coordinate data as <np.array> """
        self.z = self.read_var(self.zvar)
        self.test_shape(self.zvar, self.z.shape, 2)

    def load_p(self):
        """ Load p-coordinate data as <np.array> """
        self.p = self.read_var(self.pvar)
        new_arr = []
        for p in range(np.shape(self.p)[0]):
            new_arr.append(p)
        self.p = new_arr
        self.p = np.array(self.p)
        self.test_shape(self.pvar, self.p.shape, 1)

    def load_pn(self):
        """ Load project-name-coordinate data as <np.array> """
        self.pn = self.read_var(self.pnvar)
        self.pn = self.pn.astype('unicode')
        new_arr = []
        for p in range(np.shape(self.pn)[0]):
            new_arr.append(''.join(self.pn[p]))
        self.pn = new_arr
        self.pn = np.array(self.pn)
        self.test_shape(self.pnvar, self.pn.shape, 1)

    def load_ir(self):
        """ Load instrument-reference-coordinate data as <np.array> """
        self.ir = self.read_var(self.irvar)
        self.ir = self.ir.astype('unicode')
        new_arr = []
        for p in range(np.shape(self.ir)[0]):
            new_arr.append(''.join(self.ir[p]))
        self.ir = new_arr
        self.ir = np.array(self.ir)
        self.test_shape(self.irvar, self.ir.shape, 1)

    def load_data(self):
        """ Load profile data as <np.array> """
        self.data = self.read_var(self.datavar)
        self.test_shape(self.datavar, self.data.shape, 2)
    
    def load_ps(self):
        """ Load salinity data as <np.array> """
        self.ps = self.read_var(self.psvar)
        self.test_shape(self.psvar, self.ps.shape, 2)    
    
    def load_psalqc(self):
        """ Load salinity QC data as boolean <np.array> """
        rejectval = self.config.getint('profiles', 'qcreject')
        self.psalqc = self.read_var(self.psalqcvar)
        self.psalqc = self.psalqc.astype('unicode')
        self.psalqc[np.where(self.psalqc!='4')]= 1
        self.psalqc[np.where(self.psalqc=='4')]= 0
        self.psalqc = self.psalqc.astype(bool)
        self.test_shape(self.psalqcvar, self.psalqc.shape, 2)
    
    def load_qc(self):
        """ Load data quality control flags as <np.array> """
        rejectval = self.config.getint('profiles', 'qcreject')        
        self.qc = self.read_var(self.qcvar)
        self.qc = self.qc.astype('unicode')
        self.qc[self.qc == '1'] = 1
        self.qc[np.logical_or(self.qc=='4', self.qc == '0')] = 0
        self.qc = self.qc.astype(bool)
        self.test_shape(self.qcvar, self.qc.shape, 2)
        
    def load_posqc(self):
        """ Load position quality control flags as <np.array> """
        rejectval = self.config.getint('profiles', 'posqcreject')        
        self.posqc = self.read_var(self.posqcvar)
        self.posqc = self.posqc.astype('unicode')
        self.posqc[np.where(self.posqc != '4')] = 1
        self.posqc[np.where(self.posqc == '4')] = 0
        self.posqc = self.posqc.astype(bool)
        self.test_shape(self.posqcvar, self.posqc.shape, 1)          

    def test_shape(self, varname, varshape, ndim):
        """ Raise error if shape unexpected """ 
        if len(varshape) != ndim:
            raise ShapeError('Shape=%s. Expected %i-D array for %s' %
                              (repr(varshape), ndim, varname))
    
    def reshape_1d(self, dat):
        """ Reshape data into 1d arrays """
        if dat.ndim == 1:
            dat = np.ones_like(self.data) * dat[:, np.newaxis]
            
        if (dat.shape != self.data.shape):
            raise ShapeError('%s != %s: reshaped variables must have same shape as data.' 
                             % (repr(dat.shape), repr(self.data.shape)))

        try: 
            mask = self.data.mask
            dat = dat[mask == False]
        except AttributeError:
            pass
        
        dat = np.reshape(dat, dat.size)
        
        return dat
    
    def grid_data(self, method='mean'):
        """ Grid data using specifications in config attribute """ 
        
        # Add a few extra variables in that I print out each month as I can then
        # compare these with the IDL output to see if there are discrepancies in
        # the number of profiles being used:
        rej_profiles = 0.
        nWOD9 = 0.
        nGTSPP0 = 0.
        nGTSPP999 = 0.
        index = np.where(self.qc == False)
        print('No. of rejected temperature values', np.shape(index))
        
        # Make the salinity values missing where the salinity QC flag is bad:
        fv = self.z.fill_value
        self.ps[np.where(self.psalqc == False)] = fv
        index = np.where(self.psalqc == False)
        print('No. of rejected salinity values', np.shape(index))
        
        # Restrict to profiles with good posqc:
        self.data = self.data[self.posqc]        
        self.x = self.x[self.posqc]
        self.y = self.y[self.posqc]
        self.p = self.p[self.posqc]
        self.pn = self.pn[self.posqc]
        self.ir = self.ir[self.posqc]
        self.ps = self.ps[self.posqc]
        self.z = self.z[self.posqc]
        self.qc = self.qc[self.posqc]
        self.posqc = self.posqc[self.posqc]
        
        # Store the maximum depth (with a good quality flag) of each profile:
        zvar1 = np.ma.masked_where(self.qc == False, self.z)
        maxdepth = np.amax(zvar1, axis = 1)
        # Get only the profiles that aren't entirely bad:
        self.data = self.data[maxdepth.mask == False]
        self.x = self.x[maxdepth.mask == False]
        self.y = self.y[maxdepth.mask == False]
        self.p = self.p[maxdepth.mask == False]
        self.pn = self.pn[maxdepth.mask == False]
        self.ir = self.ir[maxdepth.mask == False]
        self.ps = self.ps[maxdepth.mask == False]
        self.z = self.z[maxdepth.mask == False]
        self.qc = self.qc[maxdepth.mask == False]
        self.posqc = self.posqc[maxdepth.mask == False]      
        maxdepth = maxdepth[maxdepth.mask == False]
        
        # I could thin to the second element of zbounds as if profiles don't go
        # that deep then I won't be able to use them (this would be true even if
        # the first element of zbounds wasn't zero) - it also removes the need 
        # for an OHC dep value in my configuration file:
        # In this _depth_lenience version, allow a 1% tolerance, i.e for my 
        # 0 - 2000m level, I will now accept profiles going down to 1980m.
        self.data = self.data[maxdepth >= self.zbounds[1] * 0.99]
        self.x = self.x[maxdepth >= self.zbounds[1] * 0.99]
        self.y = self.y[maxdepth >= self.zbounds[1] * 0.99]
        self.p = self.p[maxdepth >= self.zbounds[1] * 0.99]
        self.pn = self.pn[maxdepth >= self.zbounds[1] * 0.99]
        self.ir = self.ir[maxdepth >= self.zbounds[1] * 0.99]
        self.ps = self.ps[maxdepth >= self.zbounds[1] * 0.99]
        self.z = self.z[maxdepth >= self.zbounds[1] * 0.99]
        self.qc = self.qc[maxdepth >= self.zbounds[1] * 0.99]
        self.posqc = self.posqc[maxdepth >= self.zbounds[1] * 0.99]
        year = int(self.fname[-9:-5])
        
        # Filter out low quality XBTs:
        rem = []
        for p in range(len(self.p)):
            xbt = inst_type.is_xbt(self.pn[p], self.ir[p], self.ps[p], fv, 
              self.z[p], fv)
            if xbt[0][0] >= 0:
                # Remove any XBTs sourced from WOD where the fall rate equation
                # is unknown:
                if xbt[3] == 9:
                    rem.append(p)
                    nWOD9 += 1
                # Remove any GTSPP XBTs where the type is unknown and year is 
                # >= 1995. 
                # Or if type is unknown and it may not be a T4/T6/T7/DB because 
                # the depth it reaches is too deep. Some of these will have been
                # given the Hanawa correction and so be inaccurate - this 
                # happens in EN processing if probe code is zero:
                projectName = ''.join(self.pn[p])
                if projectName[0:5] == 'GTSPP':
                    if (xbt[4] == 0 or xbt[4] == 99 or xbt[4] == 999) and year >= 1995:
                        rem.append(p)
                        if xbt[4] == 0:
                            nGTSPP0 += 1
                        else:
                            nGTSPP999 += 1
                    if (xbt[4] == 0 and xbt[1] > 900):
                        rem.append(p)
                
        # Get rid of the low quality XBTs:
        nolowxbt = np.array(np.setdiff1d(range(len(self.p)),rem))
        self.data = self.data[nolowxbt]
        self.x = self.x[nolowxbt]
        self.y = self.y[nolowxbt]
        self.p = self.p[nolowxbt]
        self.pn = self.pn[nolowxbt]
        self.ps = self.ps[nolowxbt]
        self.ir = self.ir[nolowxbt]
        self.z = self.z[nolowxbt]
        self.qc = self.qc[nolowxbt]
        self.posqc = self.posqc[nolowxbt]
                
        # Do the vertical averaging:
        self.p = np.array(range(len(self.p)))

        # Loop over these profiles and do vertical averages:
        # Storage vectors need to start off filled with fill values otherwise
        # you get zero values being included in the averages instead of being
        # discarded. Make these storage matrices, not vectors so you can store 
        # avergae values for multiple depth levels in them.
        all_mT = np.zeros((len(self.p),np.shape(self.zbounds)[0]-1))
        all_mT.fill(fv)
        all_lT = np.zeros((len(self.p),np.shape(self.zbounds)[0]-1))
        all_lT.fill(fv)
        all_dep = np.zeros((len(self.p),np.shape(self.zbounds)[0]-1))
        all_dep.fill(fv)
        all_x = np.zeros((len(self.p),np.shape(self.zbounds)[0]-1))
        all_x.fill(fv)
        all_y = np.zeros((len(self.p),np.shape(self.zbounds)[0]-1))
        all_y.fill(fv)
        for p in range(len(self.p)):
            # 1. Select the profile of interest and make sure no tar_t1 values
            # are carried over from a previous profile:
            x_p = self.x[p]
            y_p = self.y[p]
            qc_p = np.where(np.logical_and(self.qc[p] == True, self.z[p].mask == False))
            data_p = self.data[p][qc_p]
            z_p = self.z[p][qc_p].data
            tar_t1 = fv
            # 1a. Sanity check to make sure there are no missing data going into
            # the averaging process:
            tempanddeppres = np.where(np.logical_and(data_p != fv, z_p != fv))[0]
            data_p = data_p[tempanddeppres]
            z_p = z_p[tempanddeppres]
            # 2. Make sure that depths are in correct order, but sorting takes time
            # therefore only sort if I've identified non-ascending depths:
            if np.any((z_p[1:] - z_p[:-1]) < 0):
                sortz = np.argsort(z_p)
                data_p = data_p[sortz]
                z_p = z_p[sortz]
            # To be in line with the IDL code also need to check for very wrong 
            # depth values that might have slipped through:
            udep = np.where(np.logical_and(z_p > -99.9, z_p < 10000))
            z_p = z_p[udep]
            data_p = data_p[udep]
            # 3. Find the temperature at the exact depth level - this is now 
            # more in depth as there are multiple depth levels to look at - AT 
            # THE MOMENT THIS STILL ASSUMES THAT THE FIRST LEVEL IS 0 TO SOME
            # DEPTH - THIS IS STILL A SLIGHT SIMPLIFICATION:
            dval = 0
            for dep in self.zbounds[1:]:
                # Have added in this switch to allow you to specify whether
                # you're using a fixed gap or not:
                if self.fixedgap == True:
                    maxgap = 200
                else:
                    maxgap = max(0.3*(self.zbounds[dval+1]),100)
                # Get only the levels of the profile in the depth range of 
                # interest:
                LTi1 = np.where(np.logical_and(z_p < dep, 
                  z_p >= self.zbounds[dval]))
                GEi1 = np.where(z_p >= dep)
                # Add a tolerance of 1% if profiles don't quite go down to the 
                # target depth, but almost do. If this is the case then repeat
                # the last temperature value at a new depth level assigned to be
                # the target depth:
                if (np.shape(GEi1)[1] == 0 and z_p[LTi1[0][-1]] >= dep * 0.99):
                    z_p = np.concatenate((z_p, dep), axis = None)
                    data_p = np.concatenate((data_p, data_p[-1]), axis = None)
                    LTi1 = np.where(np.logical_and(z_p < dep, 
                      z_p >= self.zbounds[dval]))
                    GEi1 = np.where(z_p >= dep)                    
                if (np.shape(LTi1)[1] != 0 and np.shape(GEi1)[1] != 0):
                    # Get the depth differences between layers and the mean 
                    # temps across layers:
                    nk = np.shape(LTi1)[1] + 1
                    dz = np.zeros(nk)
                    mt = np.zeros(nk)
                    if dval == 0:
                        dz[0] = z_p[0]
                        mt[0] = data_p[0]
                        for kk in range(1, nk):
                            dz[kk] = z_p[kk] - z_p[kk-1]
                            mt[kk] = 0.5 * (data_p[kk] + data_p[kk-1])
                    else:
                        # Effectively missing the first layer as dz will be zero
                        # there if you've calculated a temperature at that depth
                        # for tar_t1 on the previous loop, but it won't exist if
                        # you haven't been able to calculate a tar_t1 value, so
                        # then you'll have to do what you do when you're at the 
                        # first depth level and aren't garunteed a value at 0m.
                        if tar_t1 != fv:
                            dz[0] = z_p[LTi1[0][0]] - self.zbounds[dval]
                            mt[0] = 0.5 * (data_p[LTi1[0][0]] + tar_t1)
                        else:
                            dz[0] = z_p[LTi1[0][0]] - self.zbounds[dval]
                            mt[0] = data_p[LTi1[0][0]]
                        for kk in range(0, nk-1):
                            dz[kk+1] = z_p[LTi1[0][kk]+1] - z_p[LTi1[0][kk]]
                            mt[kk+1] = 0.5 * (data_p[LTi1[0][kk]+1] + data_p[LTi1[0][kk]])
                    
                    # Work out the temp at the target depth:
                    if z_p[GEi1[0][0]] == dep:
                        #print('A sampled depth is equal to the desired level')
                        tar_t1 = data_p[GEi1[0][0]]
                    else:
                        deltaT = data_p[GEi1[0][0]] - data_p[GEi1[0][0] -1]
                        deltaZ = z_p[GEi1[0][0]] - z_p[GEi1[0][0] -1]
                        tar_t1 = (dep - z_p[GEi1[0][0] -1])*(deltaT/deltaZ) + data_p[GEi1[0][0] -1]
                        dz[nk -1] = dep - z_p[GEi1[0][0] -1]
                        mt[nk -1] = 0.5*(data_p[GEi1[0][0] -1] + tar_t1)
                    
                    # Check if there are unacceptable gaps between layers:
                    test_gap = np.where(dz > maxgap)
                    if np.shape(test_gap)[1] != 0:
                        mean_t1 = fv
                        mean_t2 = fv
                        tar_t1 = fv
                    else:
                        mean_t1 = sum(np.multiply(mt,dz))/(self.zbounds[dval+1] - self.zbounds[dval])
                    
                    # Make sure there are no crazy mean values:
                    if (abs(mean_t1) > 100 and mean_t1 != fv):
                        raise ValueError('Extreme values found')
                    # Save the mean temperature at that depth and the temp at the target 
                    # depth, also save the depth and profile number:
                    all_mT[p,dval] = mean_t1
                    all_lT[p,dval] = tar_t1
                    all_dep[p,dval] = dep # Lower bound of depth
                    all_x[p,dval] = x_p
                    all_y[p,dval] = y_p
                else:
                    tar_t1 = fv # Make sure that if you have no data in a 
                    # specific depth range you don't carry an old tar_t1 value
                    # over.
                
                dval +=1

        # Reshape
        self.data_1d = self.reshape_1d(self.data)
        self.x_1d = self.reshape_1d(self.x)
        self.y_1d = self.reshape_1d(self.y)
        self.p_1d = self.reshape_1d(self.p)
        self.z_1d = self.reshape_1d(self.z)
        self.qc_1d = self.reshape_1d(self.qc)
        self.posqc_1d = self.reshape_1d(self.posqc)
        
        # Apply QC - Technically the self.posqc_1d step shouldn't be needed as 
        # these profiles will have got filtered out with the earlier QC, but 
        # it's good to leave it in there as then if I remove the depth 
        # restriction step, this step will still catch badly positioned profiles:
        #qcind = (self.qc_1d == True) & (self.posqc_1d == True)
        qcind = (self.qc_1d == True)
        self.qc_1d = self.qc_1d[qcind]
        self.posqc_1d = self.posqc_1d[qcind]
        self.data_1d = self.data_1d[qcind] # Still seems to have 99999 values in
        # it, and I can't see where they get filtered out, but they must get 
        # filtered out somewhere or the mean values wouldn't be sensible.
        self.x_1d = self.x_1d[qcind]
        self.y_1d = self.y_1d[qcind]
        self.p_1d = self.p_1d[qcind]
        self.z_1d = self.z_1d[qcind]
        # Trying to work out why I use self.posqc here, not qcind, I think it's 
        # because qcind will be one value for each level for each profile,
        # whereas posqc is one value for a whole profile - though I don't think
        # there should be any values left with a bad self.posqc - so this step 
        # is probably redundant:
        #all_mTqc = all_mT[self.posqc]
        # Make sure I'm not including values where a mean_t1 value couldn't be 
        # calculated:
        #all_mTqc1 = all_mTqc[np.where(all_mTqc != 99999.0)]
        # Need to remember that all_mT may well now be an array not a vector...
        all_mTqc1 = all_mT[all_mT != 99999.0]
        all_depqc = all_dep[all_mT != 99999.0]
        all_xqc = all_x[all_mT != 99999.0]
        all_yqc = all_y[all_mT != 99999.0]
               
        # Prepare data for gridding
        self.init_xgrid()
        self.init_ygrid()
        self.init_zgrid()
        
        # Getting the unique profile references:
        punique = np.unique(self.p_1d, return_index = True)[1]
        
        # From the unique profile references selecting only those profiles that 
        # have a mean temperature over the depth of interest:
        #puniqueqc = punique[all_mT != 99999.0]
        
        # Getting references for all the points - so every temp value will have 
        # a depth, lat and long i.e. you will still have multiple points per 
        # profile:
        points = np.vstack([self.z_1d, self.y_1d, self.x_1d]).transpose()
        
        # Getting a single reference for each profile, this will just take the 
        # first value in points, for each profile => points2 will have a 
        # latitude, a longitude and the shallowest accepted depth of the profile
        # => This point will always (pretty much) be put into the top set of 
        # grid boxes, would therefore need to consider this further if I wanted
        # to populate multiple depth level grid boxes at once:
        points2 = np.vstack([self.z_1d[punique], self.y_1d[punique], 
          self.x_1d[punique]]).transpose()
        
        # points3 is like points 2, but gets coordinates for only profiles that
        # have a mean temperature down to the depth of interest:
        #points3 = np.vstack([np.ones(len(puniqueqc)), self.y_1d[puniqueqc], 
        #  self.x_1d[puniqueqc]]).transpose()
        points3 = np.vstack([all_depqc - 1, all_yqc, all_xqc]).transpose()
        
        # Pretty self explanatory - the boundaries of the grid boxes:
        bins = [self.zbounds, self.ybounds, self.xbounds]
                
        # Grid data:
        grid_count, binedges, binno = scipy.stats.binned_statistic_dd(
            points, self.data_1d, statistic='count', bins=bins)
        grid_sum, binedges, binno = scipy.stats.binned_statistic_dd(
            points, self.data_1d, statistic='sum', bins=bins)
        grid_meansum, binedges, binno = scipy.stats.binned_statistic_dd(
            points3, all_mTqc1, statistic = 'sum', bins = bins)
        grid_pcount, binedges, binno = scipy.stats.binned_statistic_dd(
            points3, all_mTqc1, statistic='count', bins=bins)
#        grid_pcount, binedges, binno = scipy.stats.binned_statistic_dd(
#            points3, self.data_1d[puniqueqc], statistic='count', bins=bins)
#        grid_max, binedges, binno = scipy.stats.binned_statistic_dd(
#            points, self.data_1d, statistic = 'max', bins=bins)
#        grid_min, binedges, binno = scipy.stats.binned_statistic_dd(
#            points, self.data_1d, statistic = 'min', bins=bins)
#        grid_med, binedges, binno = scipy.stats.binned_statistic_dd(
#            points, all_mTqc1, statistic = 'median', bins = bins)
        
        # Sum of valid temps/ number of valid obs:
        grid_tmean = grid_sum / grid_count
        grid_tmean = np.ma.MaskedArray(grid_tmean, mask = (grid_count == 0))
        # Sum of valid mean temps/ number of valid profiles:
        grid_meantmean = grid_meansum / grid_pcount
        grid_meantmean = np.ma.MaskedArray(grid_meantmean, mask = (grid_pcount == 0))
        self.grid_tmean = grid_tmean
        self.grid_count = grid_count
        self.grid_sum = grid_sum
        self.grid_meansum = grid_meansum
        self.grid_pcount = grid_pcount
        self.grid_meantmean = grid_meantmean
#        self.grid_max = grid_max
#        self.grid_min = grid_min
#        self.grid_med = grid_med

    def create_savename(self):
        """ Generate file name based on file name and grid specification """
        
        savename = self.config.get('grid', 'dir') + self.fname.split('/')[-1]
        newsuffix = '_gridded_%ix%ix%i.nc' % (self.nx, self.ny, self.nz)
        savename = savename.replace('.nc', newsuffix)
        
        return savename

    def write_grid(self):
        """ Write gridded data to netcdf  """
        
        self.fout = self.create_savename()
        ncout = Dataset(self.fout, 'w')
        print('Writing: %s' % self.fout)
        
        # Create dimensions
        lon = ncout.createDimension(self.xvar, self.nx)
        lat = ncout.createDimension(self.yvar, self.ny)
        depth = ncout.createDimension(self.zvar, self.nz)
        tdim = ncout.createDimension('time', None)
        bndsDim = ncout.createDimension('bnds', 2)

        # Create variables
        varx = ncout.createVariable(self.xvar, 'float64', (self.xvar,))
        vary = ncout.createVariable(self.yvar, 'float64', (self.yvar,))
        varz = ncout.createVariable(self.zvar, 'float64', (self.zvar,))

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

        vartmean = ncout.createVariable('tmean', 'float32', ('time',self.zvar,self.yvar,self.xvar))
        varmtmean = ncout.createVariable(self.datavar, 'float32', ('time',self.zvar,self.yvar,self.xvar))
        varsum = ncout.createVariable('sum', 'float32', ('time',self.zvar,self.yvar,self.xvar))
        varmsum = ncout.createVariable('meansum', 'float32', ('time',self.zvar,self.yvar,self.xvar))
        varcount = ncout.createVariable('count', 'float32', ('time',self.zvar,self.yvar,self.xvar))
#        varmax = ncout.createVariable('gmax', 'float32', ('time', self.zvar, self.yvar, self.xvar))
#        varmin = ncout.createVariable('gmin', 'float32', ('time', self.zvar, self.yvar, self.xvar))
#        varmed = ncout.createVariable('median', 'float32', ('time', self.zvar, self.yvar, self.xvar))
        varpcount = ncout.createVariable('pcount', 'float32', ('time', self.zvar, self.yvar, self.xvar))
        vartime = ncout.createVariable('time', 'float64', ('time',))
        vartime.units = 'hours since 0001-01-01 00:00:00'
        vartime.calendar = 'gregorian'

        # Write to variables
        varx[:] = self.xgrid
        vary[:] = self.ygrid
        varz[:] = self.zgrid
        vartmean[:] = self.grid_tmean[np.newaxis]
        varmtmean[:] = self.grid_meantmean[np.newaxis]
        varsum[:] = self.grid_sum[np.newaxis]
        varmsum[:] = self.grid_meansum[np.newaxis]
        varcount[:] = self.grid_count[np.newaxis]
        varpcount[:] = self.grid_pcount[np.newaxis]
#        varmax[:] = self.grid_max[np.newaxis]
#        varmin[:] = self.grid_min[np.newaxis]
#        varmed[:] = self.grid_med[np.newaxis]
        vartime[:] = date2num(self.dt, units=vartime.units, calendar=vartime.calendar)
        
        # Add  global attributes
        ncout.history = 'Created ' + time.ctime(time.time())
        
        # Save
        ncout.close()  
    
    def init_xgrid(self):
        """ Initialize x dimension """
        dx = self.config.getfloat('grid', 'dx')
        xmin = self.config.getfloat('grid', 'xmin')
        xmax = self.config.getfloat('grid', 'xmax')
        self.xbounds = np.arange(xmin, xmax+dx, dx)
        self.xminbounds = self.xbounds[:-1]
        self.xmaxbounds = self.xbounds[1:]
        self.xgrid = 0.5 * (self.xminbounds + self.xmaxbounds)
        self.nx = len(self.xgrid)

    def init_ygrid(self):
        """ Initialize y dimension """
        dy = self.config.getfloat('grid', 'dy')
        ymin = self.config.getfloat('grid', 'ymin')
        ymax = self.config.getfloat('grid', 'ymax')
        self.ybounds = np.arange(ymin, ymax + dy, dy) 
        self.yminbounds = self.ybounds[:-1]
        self.ymaxbounds = self.ybounds[1:]
        self.ygrid = 0.5 * (self.yminbounds + self.ymaxbounds)
        self.ny = len(self.ygrid)

    def init_zgrid(self):
        """ Initialize y dimension """
        self.zbounds = np.array(self.config.get('grid', 'zbounds').split(','), dtype=np.float64)
        self.zminbounds = self.zbounds[:-1]
        self.zmaxbounds = self.zbounds[1:]
        self.zgrid = 0.5 * (self.zminbounds + self.zmaxbounds)
        self.nz = len(self.zgrid)


