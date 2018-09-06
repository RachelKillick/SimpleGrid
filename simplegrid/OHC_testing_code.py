""" This code is designed to test the OHC calculation code from the gridding 
step right through to the final output. It's impractical to test every month of 
data (although if this code proves effective it could be expanded) and so five 
months expected to have relatively different characteristics were chosen:
- April 1956 - pre XBTs - and largely shallower profiles.
- January 1977 - after XBTs, but pre-Argo.
- July 1997 - WOCE era => likely more deep profiles.
- October 2005 - beginning(ish) of Argo.
- March 2016 - Argo in full swing, could make this test month even later in the 
series when DeepArgo are being deployed more.

The profile files that are used for this data are those used for the BAMS time 
series at the beginning of 2018 - these are stored in /project/EN/OHC_testing - 
the same profile files should always be used (unless it is absolutely necessary 
to change them) as then it will be known that if any tests fail it is because of 
changes to the code and not changes in the profiles.

This code is currently designed to test 2 by 2 latitude-longitude gridded data
and should be run in Data/fcm/SimpleGrid/simplegrid so that the correct version 
of the inst_type function is used.

"""


import inst_type
import iris
import matplotlib.pyplot as plt
import numpy as np
import os


# Load the test files:
fpath = '/project/EN/OHC_testing/Profile_files/'
files = sorted(os.listdir(fpath))

ftpath = '/project/EN/OHC_testing/SimpleGridFiles/'
ftfiles = sorted(os.listdir(ftpath))

# Define the test depths:
tdeps = np.array([100.,300.,700.,2000.,5500.])

# Define the fill value:
fv = 99999.

# Loop through the files to test:
fno = -1
for f in files:
    fno += 1
    # Define a numpy array - initially with all values set to 0 where 
    # locations of profiles can be stored and define another numpy array with
    # all the values initially made missing where the mean temperatures can be
    # summed:
    plocs = np.zeros((len(tdeps),90,180))
    meansum = np.zeros((len(tdeps),90,180))
    # Get the year, which will be useful when filtering out unreliable XBTs:
    yr = int(f[-9:-5])
    # Get the position QC flag:
    posqc = iris.load(fpath + f, 'quality on position (latitude and longitude)')[0]
    posqc = posqc.data.astype('unicode')
    # Get the temperature values:
    potm = iris.load(fpath + f, 'sea_water_conservative_temperature')[0]
    # Get the temperature QC values:
    potmqc = iris.load(fpath + f, 'quality on pot. temperature')[0]
    potmqc = potmqc.data.astype('unicode')
    # Get the depths:
    depth = iris.load(fpath + f, 'corrected depth')[0]
    # Values with bad qc aren't used, therefore mask all these values out:
    potm = np.ma.masked_where(potmqc != '1', potm.data)
    depth = np.ma.masked_where(potmqc != '1', depth.data)
    # Get the latitudes and longitudes:
    lats = iris.load(fpath + f, 'latitude of the station, best estimated value')[0]
    lons = iris.load(fpath + f, 'longitude of the station, best estimated value')[0]
    # Get the other variables that will be needed when removing XBTs so you can 
    # remove the bad profile references from these too:
    projname = iris.load(fpath + f, 'name of the project')[0]
    projname = projname.data.astype('unicode')
    psal = iris.load(fpath + f, 'corrected practical salinity')[0]
    instref = iris.load(fpath + f, 'instrument type')[0] # If you were to 
    # accidentally include an underscore here then you'd change the variable!
    # Perhaps something to consider if I do any work on the profile netcdf files.
    instref = instref.data.astype('unicode')
    # Remove the profiles with bad position QC:
    lats = lats[posqc == '1']
    lons = lons[posqc == '1']
    potm = potm[posqc == '1']
    depth = depth[posqc == '1']
    projname = projname[posqc == '1']
    psal = psal[posqc == '1']
    instref = instref[posqc == '1']
    # Remove the profiles that don't even go down to the first depth of interest
    # as this means they won't be able to be used in any of the averaging (this 
    # is tricky as the best way to get maximum depths is probably how I did it 
    # in my main code, but if I use the same code then it's not really a test!).
    # A test needs to be comparing to known values - which is why I wasn't sure
    # if in addition to these broad checks I should choose a few grid boxes to 
    # investigate more thoroughly:
    lats = lats[np.ma.any(np.logical_and(depth.data >= tdeps[0], 
      depth.mask == False), axis = 1)]
    lons = lons[np.ma.any(np.logical_and(depth.data >= tdeps[0], 
      depth.mask == False), axis = 1)]
    potm = potm[np.ma.any(np.logical_and(depth.data >= tdeps[0], 
      depth.mask == False), axis = 1)]
    projname = projname[np.ma.any(np.logical_and(depth.data >= tdeps[0], 
      depth.mask == False), axis = 1)]
    psal = psal[np.ma.any(np.logical_and(depth.data >= tdeps[0], 
      depth.mask == False), axis = 1)]
    instref = instref[np.ma.any(np.logical_and(depth.data >= tdeps[0], 
      depth.mask == False), axis = 1)]
    depth = depth[np.ma.any(np.logical_and(depth.data >= tdeps[0], 
      depth.mask == False), axis = 1)]
    # Remove questionable XBT profiles - this is any GTSPP XBTs where the type 
    # is unknown and the year is >= 1995 or if the type is unknown and it can't 
    # be a T4/T6/T7 or DB because the depth it reaches is too deep. Also remove
    # any XBTs sourced from WOD where the fall rate equation is unknown.
    badxbt = []
    rejset = ['  0', ' 99', '999']
    for p in range(np.shape(lats)[0]):
        xbt = inst_type.is_xbt(''.join(projname[p]), ''.join(instref[p]), 
          psal.data[p], fv, depth[p].data, fv)
        if xbt[0][0] >= 0:
            if np.logical_and(projname[p][0] == 'W', instref[p][-1] == '9'):
                badxbt.append(p)
            elif projname[p][0] == 'G':
                if np.logical_and(yr >= 1995, ''.join(instref[p][60:63]) in rejset):
                    badxbt.append(p)
                if np.logical_and(''.join(instref[p][60:63]) == '  0', xbt[1] >= 900):
                    badxbt.append(p)
    goodprofs = list(set(np.arange(np.shape(lats)[0])) - set(badxbt))
    # Now only keep the non-questionable profiles:
    lats = lats.data[goodprofs]
    lons = lons.data[goodprofs]
    potm = potm[goodprofs]
    projname = projname[goodprofs]
    psal = psal.data[goodprofs]
    instref = instref[goodprofs]
    depth = depth[goodprofs]
    # Set up the arrays for storing the average temps and temps at target depth
    # for each profile:
    all_mT = np.full((len(lats),len(tdeps)), fv)
    all_lT = np.full((len(lats), len(tdeps)), fv)
    # Now loop over the profiles:
    for p in range(len(lats)):
        #if (lats[p] < 30 or lats[p] >= 32):
        #    continue
        #if (lons[p] < -140 or lons[p] >= -138):
        #    continue
        #print('p',p)
        # Make sure the profile contains only good data, no missing data and is 
        # in depth order (values with bad potm QC have been masked further up):
        dep = depth[p]
        temp = potm[p]
        nonmiss = np.logical_and(dep.mask == False, temp.mask == False)
        dep = dep[nonmiss]
        temp = temp[nonmiss]
        # Get rid of any clearly wrong values:
        temp = temp[np.logical_and(dep > -99.9, dep < 10000)]
        dep = dep[np.logical_and(dep > -99.9, dep < 10000)]
        if np.any((dep[:-1] - dep[1:]) > 0):
            sor = np.argsort(dep)
            dep = dep[sor]
            temp = temp[sor]
        # Find the temperature at the depth level of interest:
        tar_t1 = fv
        for d in tdeps:
            i = np.where(tdeps == d)[0][0]
            # If you've WORKED OUT (i.e. you didn't already have it) a value at 
            # the target depth previously, start your prof segment with this:
            if (i != 0 and tar_t1 != fv and tset == False):
                dep1 = np.append(tdeps[i-1], dep)
                temp1 = np.append(tar_t1, temp)
            else:
                dep1 = dep
                temp1 = temp
            # Work out which levels of the profile are in the range of interest:
            if i == 0:
                inside = np.logical_and(dep1 >= 0, dep1 < d)
            else:
                inside = np.logical_and(dep1 >= tdeps[i-1], dep1 < d)
            over = (dep1 >= d)
            # Select these levels:
            indep = dep1[inside]
            odep = dep1[over]
            intemp = temp1[inside]
            otemp = temp1[over]
            # As long as inside and over both have data in you can work out the temp
            # at the target depth:
            if (inside[inside == True].size != 0 and over[over == True].size != 0):
                if odep[0] == d:
                    tset = True
                    tar_t1 = otemp[0]
                else:
                    tset = False
                    tar_t1 = (d - indep[-1])*((otemp[0] - intemp[-1])/(odep[0] - 
                      indep[-1])) + intemp[-1]
                # Once you have the temperature at the target depth you can work out 
                # the average temperature over the depth window:
                indep = np.append(indep, d)
                intemp = np.append(intemp, tar_t1)
                dz = np.full(len(indep), fv)
                dt = np.full(len(intemp), fv)
                if i == 0:
                    dz[0] = indep[0]
                else:
                    dz[0] = indep[0] - tdeps[i-1]
                dt[0] = intemp[0] # Assuming the first level is at the same temp as 
                # the first measurement is an assumption that might want to be 
                # worked on in the future.
                for l in range(1, len(dz)):
                    dz[l] = indep[l] - indep[l-1]
                    dt[l] = (intemp[l] + intemp[l-1])/2
                # Now check if maxgap is exceeded at any point in time:
                maxgap = max(100, d*0.3)
                if np.any(dz > maxgap):
                    mean_t1 = fv
                    tar_t1 = fv
                else:
                    if i == 0:
                        mean_t1 = sum(dz * dt)/d
                    else:
                        mean_t1 = sum(dz * dt)/(d - tdeps[i-1])
                    lat1 = int(abs(np.ceil((-90 - lats[p])/2)))
                    lon1 = int(abs(np.ceil((-180 - lons[p])/2)))
                    plocs[i, lat1, lon1] += 1
                    meansum[i, lat1, lon1] += mean_t1
            else:
                tar_t1 = fv
                mean_t1 = fv
            all_mT[p,i] = mean_t1
            all_lT[p,i] = tar_t1
    
    # Once all my profiles and the appropriate grid boxes have been incremented 
    # to indicate where averages have been performed I need to average the 
    # averages:
    output = (meansum/plocs)
    output = np.ma.masked_where(plocs == 0, output)
    
    # Now check if these output arrays match those from the SimpleGrid code I'm
    # testing:
    pcount = iris.load(ftpath + ftfiles[fno], 'pcount')[0].data
    msum = iris.load(ftpath + ftfiles[fno], 'meansum')[0].data
    constemp = iris.load(ftpath + ftfiles[fno], 'CONS_TEMP')[0].data
    
    if np.any(plocs != pcount[0]):
        print('Discrepancies in pcount for {0}'.format(f))
        print((plocs - pcount[0])[plocs != pcount[0]])
    if np.any(meansum != msum[0]):
        print('Discrepancies in meansum for {0}'.format(f))
        print((meansum - msum[0])[meansum != msum[0]])
    if np.any(output != constemp[0]):
        print('Discrepancies in constemp for {0}'.format(f))
        print((output - constemp[0])[output != constemp[0]])
        plt.plot((output - constemp[0])[output != constemp[0]])
        plt.show()
        plt.scatter(pcount[0][output != constemp[0]],
          (output - constemp[0])[output != constemp[0]])
        plt.show()

# The areas where I could see the difference was smaller seemed to be clustered 
# around the same latitude boxes (10-15 (for March 2016 at least) - note this is BOX reference not a 
# position on the globe, though this can be worked out from the box reference), 
# but when I did a contourf plot of the differences (this is mind bending and 
# eye-hurting to look at) there didn't seem to be any systematic differences in 
# errors in different locations - so I looked at pcount versus the difference 
# this showed quite a nice normal-like distribution with bigger errors for fewer
# profiles in a grid box (potentially when there's more the differences cancel 
# each other out) so I think that my code is fine, but I'm not ruling out 
# further investigation at a later date - it is worth noting though that the 
# differences are small! Could run this test code on any new months too though, 
# then if I find larger differences I'll know that there's a mistake and that I 
# shouldn't publish anything until I've got to the bottom of it.
# The normal distributions from my test months are saved in the Documents/OHC
# folder and take the form pcount_[year][month].png. I think that October 2005 
# might have the largest discrepancies so if I wanted to do some further 
# investigation this would probably be a good place to start with January 1977 
# being the second best candidate.
