""" Following on from OHC_testing_code.py this code is seeking to work out why
there are discrepancies between my test code and my run code - so I can work out
which one is doing the right thing and which one is not! For this I'm using 
April 1956 as my test month as there are fewer profiles to work with then and 
I'm investigating 4 grid boxes:
[0][0][18][170] - lats: [-54, -52), lons: [160, 162); profs 4653
[0][0][60][20] - lats: [30, 32), lons: [-140, -138); profs 0, 8, 22, 33, 42, 45,
49, 56, 66, 90, 102, 109, 116, 122, 135, 158, 183,  208,  224,  261,  295,  317,
353,  369,  388,  431,  498,  530,  561,  584,  608,  654,  725, 750,  777, 807,
833,  884,  927,  951,  976,  995, 1019, 1061, 1109, 1124, 1144, 1165, 1182,
1213, 1252, 1276, 1294, 1315, 1334, 1367, 1403, 1419, 1436, 1456, 1471, 1501, 
1536, 1558, 1574, 1594, 1609, 1683, 1706, 1732, 1760, 1782, 1810, 1837, 1856, 
1870, 1887, 1904, 1989, 2004, 2083, 2114, 2142, 2145, 2159, 2161, 2189, 2192,
2202, 2208, 2250, 2282, 2293, 2350, 2367, 2463, 2499, 2606, 2631, 2679, 2833, 
2905, 2948, 3064, 3129, 3172, 3245, 3269, 3279, 3442, 3459, 3639, 3687, 3734,
3750, 3770, 3786, 3823, 3914, 3972, 4043, 4063, 4140, 4216, 4317, 4325, 4339,
4368, 4376, 4388, 4402, 4414, 4434, 4456 - this is 134 profs, although only 130
get used according to pcount.
[0][1][16][169] - lats: [-58, -56), lons: [158, 160]; profs 4483, 4516
[0][2][66][66] - lats: [42, 44), lons: [-48, -46); profs 4752, 4823, 4824, 4829,
4884, 4885, 4886 (although only six get used according to pcount)
[0][3][71][80] - lats: [52, 54), lons: [-20, -18); profs 87,  403,  512,  611, 
735,  830,  926, 1018, 1255, 1331, 1542, 1679, 1987, 2283, 2421, 2593, 2688, 
6336, 6341, 6344,6350 (although according to pcount only one of these gets used)
The count of profiles contributing to gridboxes is the same for all grid boxes 
in the five months I tested, so I think it's the code to get the average over a 
profile that is causing the discrepancies.
"""


import iris
import numpy as np
import matplotlib.pyplot as plt


# The profile file:

filen = '/project/EN/OHC_testing/Profile_files/EN.4.2.1.f.profiles.g10.195604.nc'

# Look at the profiles in the first grid box under investigation - this is for 
# the 0 - 100m layer:

potm = iris.load(filen, 'sea_water_conservative_temperature')[0]
# Get the temperature QC values:
potmqc = iris.load(filen, 'quality on pot. temperature')[0]
potmqc = potmqc.data.astype('unicode')
# Get the depths:
depth = iris.load(filen, 'corrected depth')[0]
# Values with bad qc aren't used, therefore mask all these values out:
potm = np.ma.masked_where(potmqc != '1', potm.data)
depth = np.ma.masked_where(potmqc != '1', depth.data)
fv = 99999.0
tdeps = np.array([100.,300.,700.,2000.,5500.])

# Select the profile(s) of interest:

profs = np.array([4653])
profs = np.array([0,8,22,33,42,45,49,56,66,90,102,109,116,122,135,158,183,208,
224,261,295,317,353,369,388,431,498,530,561,584,608,  654,  725, 750,  777, 807,
833,  884,  927,  951,  976,  995, 1019, 1061, 1109, 1124, 1144, 1165, 1182,
1213, 1252, 1276, 1294, 1315, 1334, 1367, 1403, 1419, 1436, 1456, 1471, 1501, 
1536, 1558, 1574, 1594, 1609, 1683, 1706, 1732, 1760, 1782, 1810, 1837, 1856, 
1870, 1887, 1904, 1989, 2004, 2083, 2114, 2142, 2145, 2159, 2161, 2189, 2192,
2202, 2208, 2250, 2282, 2293, 2350, 2367, 2463, 2499, 2606, 2631, 2679, 2833, 
2905, 2948, 3064, 3129, 3172, 3245, 3269, 3279, 3442, 3459, 3639, 3687, 3734,
3750, 3770, 3786, 3823, 3914, 3972, 4043, 4063, 4140, 4216, 4317, 4325, 4339,
4368, 4376, 4388, 4402, 4414, 4434, 4456])

posqc = iris.load(filen, 'quality on position (latitude and longitude)')[0]
posqc = posqc.data.astype('unicode')
pno = np.arange(np.shape(depth)[0])
pno = pno[np.logical_and(posqc == '1', np.ma.any(np.logical_and(depth.data >= tdeps[0], 
      depth.mask == False), axis = 1))]
profs1 = sorted((set(profs).intersection(set(pno))))

tot = 0
for p in profs1:
    print(p)
    dep = depth[p]
    temp = potm[p]
    nonmiss = np.logical_and(dep.mask == False, temp.mask == False)
    dep = dep[nonmiss]
    temp = temp[nonmiss]
    # Get rid of any clearly wrong values and ensure the depths are descending:
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
            #print('d, indep', d, indep)
            #print('intemp', intemp)
            #print('odep',odep)
            #print('otemp',otemp)
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
            #print('dz', dz)
            #print('dt', dt)
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
                #lat1 = int(abs(np.ceil((-90 - lats[p])/2)))
                #lon1 = int(abs(np.ceil((-180 - lons[p])/2)))
                #plocs[i, lat1, lon1] += 1
                #meansum[i, lat1, lon1] += mean_t1
        else:
            tar_t1 = fv
            mean_t1 = fv
        print('mean_t1', mean_t1)
        print('tar_t1', tar_t1)
    tot += mean_t1
    print('meansum', tot)
