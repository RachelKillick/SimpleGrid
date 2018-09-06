""" Focusing on the same month and year (April 1956) and the same profiles (4653
initially) as in OHC_advanced_testing.py, but this is reproducing the 
profiles.py code not the test code so I can try and work out where the 
discrepancies are arising.

"""

import iris
import matplotlib.pyplot as plt
import numpy as np

# The profile file:

filen = '/project/EN/OHC_testing/Profile_files/EN.4.2.1.f.profiles.g10.195604.nc'

# Look at the profiles in the first grid box under investigation - this is for 
# the 0 - 100m layer:

# Get the temperature:
potm = iris.load(filen, 'sea_water_conservative_temperature')[0].data

# Get the temperature QC values:
potmqc = iris.load(filen, 'quality on pot. temperature')[0]
potmqc = potmqc.data.astype('unicode')
potmqc[potmqc == '1'] = 1
potmqc[np.logical_or(potmqc=='4', potmqc == '0')] = 0
potmqc = potmqc.astype(bool)

# Get the depths:
depth = iris.load(filen, 'corrected depth')[0].data

# Select the profiles to investigate and set some necessary constants:
#profs = np.array([4653]) # For first grid box under investigation
profs = np.array([0,8,22,33,42,45,49,56,66,90,102,109,116,122,135,158,183,208,
224,261,295,317,353,369,388,431,498,530,561,584,608,  654,  725, 750,  777, 807,
833,  884,  927,  951,  976,  995, 1019, 1061, 1109, 1124, 1144, 1165, 1182,
1213, 1252, 1276, 1294, 1315, 1334, 1367, 1403, 1419, 1436, 1456, 1471, 1501, 
1536, 1558, 1574, 1594, 1609, 1683, 1706, 1732, 1760, 1782, 1810, 1837, 1856, 
1870, 1887, 1904, 1989, 2004, 2083, 2114, 2142, 2145, 2159, 2161, 2189, 2192,
2202, 2208, 2250, 2282, 2293, 2350, 2367, 2463, 2499, 2606, 2631, 2679, 2833, 
2905, 2948, 3064, 3129, 3172, 3245, 3269, 3279, 3442, 3459, 3639, 3687, 3734,
3750, 3770, 3786, 3823, 3914, 3972, 4043, 4063, 4140, 4216, 4317, 4325, 4339,
4368, 4376, 4388, 4402, 4414, 4434, 4456]) # For second grid box under 
# investigation - some of these profiles don't get used though - work out which 
# ones (the code for this is in OHC_advanced_testing.py):
profs1 = np.array([0, 8, 22, 33, 42, 45, 49, 56, 66, 90, 102, 109, 116, 122,135,
158, 183, 208, 224, 261, 295, 317, 353, 369, 388, 431, 498, 530, 561, 584, 608, 
654, 725, 750, 777, 807, 833, 884, 927, 951, 976, 995, 1019, 1061, 1109, 1124, 
1144, 1165, 1182, 1213, 1252, 1276, 1294, 1315, 1334, 1367, 1403, 1419, 1436, 
1456, 1471, 1501, 1536, 1558, 1574, 1594, 1609, 1683, 1706, 1732, 1760, 1782, 
1810, 1856, 1870, 1887, 1904, 1989, 2004, 2083, 2114, 2142, 2159, 2161, 2189, 
2202, 2208, 2250, 2282, 2293, 2350, 2367, 2463, 2499, 2606, 2631, 2679, 2833, 
2905, 2948, 3064, 3129, 3172, 3245, 3269, 3279, 3442, 3459, 3639, 3687, 3750,
3770, 3786, 3823, 3914, 3972, 4043, 4063, 4140, 4216, 4317, 4325, 4339, 4368,
4376, 4388, 4402, 4414, 4434, 4456])

fv = 99999.0
zbounds = np.array([0.0,100.0,300.0,700.0,2000.0,5500.0])

tot = 0
for p in profs1:
    print(p)
    # Values with bad qc aren't used, therefore mask all these values out:
    qc_p = np.where(np.logical_and(potmqc[p] == True, depth[p].mask == False))
    data_p = potm[p][qc_p]
    z_p = depth[p][qc_p].data
    tar_t1 = fv
    # Check no missing data are going to be processed:
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
    for dep in zbounds[1:]:
        maxgap = max(0.3*(zbounds[dval+1]),100)
        # Get only the levels of the profile in the depth range of 
        # interest:
        LTi1 = np.where(np.logical_and(z_p < dep, z_p >= zbounds[dval]))
        GEi1 = np.where(z_p >= dep)
        if (np.shape(LTi1)[1] != 0 and np.shape(GEi1)[1] != 0):
            # Get the depth differences between layers and the mean temps across
            # layers:
            nk = np.shape(LTi1)[1] + 1
            dz = np.zeros(nk)
            mt = np.zeros(nk)
            if dval == 0:
                dz[0] = z_p[0]
                mt[0] = data_p[0]
                for kk in range(1, nk):
                    dz[kk] = z_p[kk] - z_p[kk-1]
                    mt[kk] = 0.5 * (data_p[kk] + data_p[kk-1])
                #print('dval', dval)
                #print('dz', dz)
                #print('mt', mt)
            else:
                # Effectively missing the first layer as dz will be zero
                # there if you've calculated a temperature at that depth
                # for tar_t1 on the previous loop, but it won't exist if
                # you haven't been able to calculate a tar_t1 value, so
                # then you'll have to do what you do when you're at the 
                # first depth level and aren't garunteed a value at 0m.
                if tar_t1 != fv:
                    dz[0] = z_p[LTi1[0][0]] - zbounds[dval]
                    mt[0] = 0.5 * (data_p[LTi1[0][0]] + tar_t1)
                    #print('dz1', dz)
                    #print('mt1', mt)
                else:
                    dz[0] = z_p[LTi1[0][0]] - zbounds[dval]
                    mt[0] = data_p[LTi1[0][0]]
                for kk in range(0, nk-1):
                    dz[kk+1] = z_p[LTi1[0][kk]+1] - z_p[LTi1[0][kk]]
                    mt[kk+1] = 0.5 * (data_p[LTi1[0][kk]+1] + data_p[LTi1[0][kk]])
                    #print('dz2', dz)
                    #print('mt2', mt)
            # Work out the temp at the target depth:
            if z_p[GEi1[0][0]] == dep:
                print('A sampled depth is equal to the desired level')
                tar_t1 = data_p[GEi1[0][0]]
            else:
                deltaT = data_p[GEi1[0][0]] - data_p[GEi1[0][0] -1]
                deltaZ = z_p[GEi1[0][0]] - z_p[GEi1[0][0] -1]
                tar_t1 = (dep - z_p[GEi1[0][0] -1])*(deltaT/deltaZ) + data_p[GEi1[0][0] -1]
                dz[nk -1] = dep - z_p[GEi1[0][0] -1]
                mt[nk -1] = 0.5*(data_p[GEi1[0][0] -1] + tar_t1)
                #print('dz3', dz)
                #print('mt3', mt)
                
            # Check if there are unacceptable gaps between layers:
            test_gap = np.where(dz > maxgap)
            if np.shape(test_gap)[1] != 0:
                mean_t1 = fv
                mean_t2 = fv
                tar_t1 = fv
            else:
                mean_t1 = sum(np.multiply(mt,dz))/(zbounds[dval+1] - zbounds[dval])
            
            # Make sure there are no crazy mean values:
            if (abs(mean_t1) > 100 and mean_t1 != fv):
                raise ValueError('Extreme values found')
            # Save the mean temperature at that depth and the temp at the target 
            # depth, also save the depth and profile number:
            #all_mT[p,dval] = mean_t1
            #all_lT[p,dval] = tar_t1
            #all_dep[p,dval] = dep # Lower bound of depth
            #all_x[p,dval] = x_p
            #all_y[p,dval] = y_p
        else:
            tar_t1 = fv # Make sure that if you have no data in a 
            # specific depth range you don't carry an old tar_t1 value
            # over.
        
        dval +=1
        print('mean_t1', mean_t1)
        print('tar_t1', tar_t1)
    tot += mean_t1
    print('meansum', tot)
