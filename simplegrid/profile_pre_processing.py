""" Code to try and test the pre-processing of EN4 files so that they can be 
passed through Chris's SimpleGrid code in a similar way to how they would be 
passed through the IDL code. Initial focus is on obtaining only profiles that go
down to 700m."""

import inst_type
import iris
import numpy as np
from netCDF4 import Dataset


# Load the profiles:
fname = '/scratch/rkillick/BAMS/g10_profiles/EN.4.2.1.f.profiles.g10.200001.nc'
year = int(fname[-9:-5])
minyr = 1950
maxyr = 1950
qcreject = 4
posqcreject = 4
ncf = Dataset(fname,'r')
xvar = ncf.variables['LONGITUDE'][:]
yvar = ncf.variables['LATITUDE'][:]
zvar = ncf.variables['DEPH_CORRECTED'][:]
pvar = ncf.variables['DC_REFERENCE'][:]
pnvar = ncf.variables['PROJECT_NAME'][:]
irvar = ncf.variables['INST_REFERENCE'][:]
psvar = ncf.variables['PSAL_CORRECTED'][:]
datavar = ncf.variables['CONS_TEMP'][:]
qcvar = ncf.variables['POTM_CORRECTED_QC'][:]
qcvarps = ncf.variables['PSAL_CORRECTED_QC'][:]
posqcvar = ncf.variables['POSITION_QC'][:]
juldvar = ncf.variables['JULD'][:]
OHCdep = 700.
maxgap = 200.

# Make sure the variables are the right shape to start with:

new_arr = []
for p in range(np.shape(pvar)[0]):
    new_arr.append(p)

pvar = new_arr
pvar = np.array(pvar)

rejectval = '4'
posqc = posqcvar != rejectval
qc = qcvar != rejectval

# Mimicking what is done in IDL - set PSAL values where the PSAL_CORRECTED_QC 
# value is '4' to be missing:
fv = zvar.fill_value
psvar[np.where(qcvarps == '4')] = fv

# Try thinning to only deep profiles first and those without a bad pos QC flag:
posqcset = np.where(posqcvar != '1')[0] # Bad Pos QC flag
deepset = np.unique(np.where(np.logical_and(np.logical_and(zvar.data >= 700, 
  zvar.mask == False),qcvar != '4'))[0]) # Profiles that go down to 700m with
# acceptable measurements.
combset = np.array(np.setdiff1d(deepset,posqcset)) # Profiles that are deep and 
# aren't bad at depth - this hasn't removed any levels of profiles, just whole 
# profiles that either don't go reliably to 700m or have a bad pos QC flag.

# Get only the good, deep profiles:
datavar = datavar[combset]
xvar = xvar[combset]
yvar = yvar[combset]
pvar = pvar[combset]
pnvar = pnvar[combset]
psvar = psvar[combset]
irvar = irvar[combset]
zvar = zvar[combset]
qcvar = qc[combset]
posqcvar = posqc[combset] # Only have profiles with good posqc.
fv = zvar.fill_value

# I don't think this thins out bad levels from 'good' profiles. So this step should
# have got rid of whole profiles where posqc is bad and profiles where they either
# don't go down to 700m or all depths >= 700m are bad. Values that are less than 
# 700m (or greater than or equal to, but not the only one) and are bad, won't 
# yet have been filtered out, hence the need still for the next QC thinning step.

# Now filter out low quality XBTs (this is worth doing here as then even getting 
# rid of the depth restriction we'll still be filtering out doubtful XBTs).
# Have already only got good data => Want to loop over all the profiles - not 
# just the numbers in pvar (as some of these will now be out of range):

# Set up an empty list to store the values that need to be removed:
rem = []

for p in range(len(pvar)):
    projname = ''.join(pnvar[p])
    instref = ''.join(irvar[p])
    xbt = inst_type.is_xbt(projname, instref, psvar[p,:], fv, zvar[p,:], fv)
    if xbt[0][0] >= 0:
        # Remove any XBTs sourced from WOD where the fall rate equation is
        # unknown:
        if xbt[3] == 9:
            rem.append(p)
        # Remove any GTSPP XBTs where the type is unknown and year is >= 1995. 
        # Or if type is unknown and it may not be a T4/T6/T7/DB because the 
        # depth it reaches is too deep. Some of these will have been given the 
        # Hanawa correction and so be inaccurate - thus happens in EN processing
        # if probe code is zero:
        projectName = ''.join(pnvar[p])
        if projectName[0:5] == 'GTSPP':
            if (xbt[4] == 0 or xbt[4] == 99 or xbt[4] == 999) and year >= 1995:
                rem.append(p)
            if (xbt[4] == 0 and xbt[1] > 900):
                rem.append(p)

nolowxbt = np.array(np.setdiff1d(range(len(pvar)),rem)) # Profiles that are deep and aren't bad

# Get only the good XBT profiles:
datavar = datavar[nolowxbt]
xvar = xvar[nolowxbt]
yvar = yvar[nolowxbt]
pvar = pvar[nolowxbt]
pnvar = pnvar[nolowxbt]
psvar = psvar[nolowxbt]
irvar = irvar[nolowxbt]
zvar = zvar[nolowxbt]
qcvar = qcvar[nolowxbt]
posqcvar = posqcvar[nolowxbt]

# QC currently happens once the data have been made one dimensional - is this 
# necessary? I don't think so, if I QC the data before making them one 
# dimensional I can also do a vertical average of each profile before making
# them one dimensional. Then the averaging step will be able to do simple 
# averaging.

# Actually, no, the QC step has to happen after everything has been made 1d to 
# make sure the right levels get removed - unless I did a QC step inside the 
# vertical averaging?

pvar = np.array(range(len(pvar)))

# Loop over these profiles and do vertical averages:
all_mT = np.zeros(len(pvar))
all_lT = np.zeros(len(pvar))
for p in pvar:
    # 1. Select the profile of interest:
    x_p = xvar[p]
    y_p = yvar[p]
    qc_p = np.where(np.logical_and(qcvar[p] == True, zvar[p].mask == False))
    data_p = datavar[p][qc_p]
    z_p = zvar[p][qc_p].data
    # 1a. Sanity check to make sure no depth or temp data are missing:
    tempanddeppres = np.where(np.logical_and(data_p != fv, z_p != fv))[0]
    data_p = data_p[tempanddeppres]
    z_p = z_p[tempanddeppres]
    # 2. Make sure that depths are in correct order, but sorting takes time
    # therefore only sort if I've identified non-ascending depths:
    if np.any((z_p[1:] - z_p[:-1]) < -1):
        sortz = np.argsort(z_p)
        data_p = data_p[sortz]
        z_p = z_p[sortz]
    # 2a. To conform with what IDL does I should check at this stage if there 
    # are any outrageous depth values that have thus far slipped through:
    udep = np.where(np.logical_and(z_p > -99.9, z_p < 10000))[0]
    z_p = z_p[udep]
    data_p = data_p[udep]
    # 3. Find the temperature at the exact depth level:
    LTi1 = np.where(z_p < OHCdep)
    GEi1 = np.where(z_p >= OHCdep)
    if (np.shape(LTi1)[1] != 0 and np.shape(GEi1)[1] != 0):
        # Get the depth differences between layers and the mean temps across
        # layers:
        nk = np.shape(LTi1)[1] + 1
        dz = np.zeros(nk)
        mt = np.zeros(nk)
        dz[0] = z_p[0]
        mt[0] = data_p[0]
        for kk in (np.array(range(nk -1)) + 1):
            dz[kk] = z_p[kk] - z_p[kk-1]
            mt[kk] = 0.5 * (data_p[kk] + data_p[kk-1])
        # Work out the temp at the target depth:
        if z_p[GEi1[0][0]] == OHCdep:
            print('A sampled depth is equal to the desired level')
            tar_t1 = data_p[GEi1[0][0]]
        else:
            deltaT = data_p[GEi1[0][0]] - data_p[GEi1[0][0] -1]
            deltaZ = z_p[GEi1[0][0]] - z_p[GEi1[0][0] -1]
            tar_t1 = (OHCdep - z_p[GEi1[0][0] -1])*(deltaT/deltaZ) + data_p[GEi1[0][0] -1]
            dz[nk -1] = OHCdep - z_p[GEi1[0][0] -1]
            mt[nk -1] = 0.5*(data_p[GEi1[0][0] -1] + tar_t1)
        # Check if there are unacceptable gaps between layers:
        test_gap = np.where(dz > maxgap)
        if np.shape(test_gap)[1] != 0:
            mean_t1 = fv
            mean_t2 = fv
        else:
            mean_t1 = sum(np.multiply(mt,dz))/OHCdep
        # Make sure there are no crazy mean values:
        if (abs(mean_t1) > 100 and mean_t1 != fv):
            raise ValueError, 'Extreme values found'
        # Save the mean temperature at that depth and the temp at the target 
        # depth:
        all_mT[p] = mean_t1
        all_lT[p] = tar_t1
        # Now what I want to do (as a weighted average fix for the time being)
        # is to make sure that this weighted mean gets used instead of the 
        # simple mean - I think the simplest (VERY messy) way of doing this is
        # to replace all the actual profile values with this mean:
        if mean_t1 != fv:
            datavar[p][np.where(zvar[p] <= OHCdep)] = mean_t1
        else:
            qcvar[p][np.where(zvar[p] <= OHCdep)] = False
        # THINK OF A BETTER WAY OF DOING THIS - OVERWRITING ISN'T A GREAT IDEA!
    
# Make them all one-dimensional:
datavar_1d = reshape_1d(datavar)
xvar_1d = reshape_1d(xvar)
yvar_1d = reshape_1d(yvar)
pvar_1d = reshape_1d(pvar)
zvar_1d = reshape_1d(zvar)
qcvar_1d = reshape_1d(qcvar)
posqcvar_1d = reshape_1d(posqcvar)

# Apply QC:
qcind = (qcvar_1d == True) & (posqcvar_1d == True)
qc_1d = qcvar_1d[qcind]
posqc_1d = posqcvar_1d[qcind]
data_1d = datavar_1d[qcind]
x_1d = xvar_1d[qcind]
y_1d = yvar_1d[qcind]
p_1d = pvar_1d[qcind]
z_1d = zvar_1d[qcind]

# Get only the profiles that go down to at least 700m (with a good QC flag):
#pnew = np.array([])
#for p in set(p_1d):
#    pind = np.where(np.logical_and(p_1d == p,z_1d.mask == False))
#    try:
#        if z_1d[pind].data[-1] >= 700.:
#            pnew = np.concatenate((pnew, pind[0]))
#    except IndexError:
#        pass


# Defining functions to do the same as in Chris' script:

def reshape_1d(dat):
    """ Reshape data into 1d arrays """
    if dat.ndim == 1:
        dat = np.ones_like(datavar) * dat[:, np.newaxis]
    if dat.shape != datavar.shape :
        raise ShapeError('%s != %s: reshaped variables must have same shape as data.' 
                             % (repr(dat.shape), repr(data.shape)))
        
    try:
        mask = datavar.mask
        dat = dat[mask == False]
    except AttributeError:
        pass
    
    dat = np.reshape(dat, dat.size)
    return dat



