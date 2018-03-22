"""Progs to say the type of instrument.

IMPORTANT NOTE: This version of inst_type.py is different from the other 
versions (e.g. those used in the classification of XBTs for correction). This 
code has been changed slightly so that the output can be used as the output
from is_xbt.pro is used when filtering out low quality XBTs in the IDL OHC code.
The changes are that probeCode is now output and that 0 or 99 or 999 are treated
as missing probeCode types instead of just 0 or 99. The restrictions for OHC are
you shouldn't use the data if probeCode is 0 or 999 AND the year is greater than
or equal to 1995 or the maxDepth is greater than 900.

"""

from netCDF4 import Dataset
import netCDF4
import numpy as np
# Have removed the call to import sag_utilities as sag as it didn't seem to be 
# used.

def is_xbt(projectName,
           instRef,
           salinities,
           salinityFV,
           depths,
           depthFV):

    """Say if something is an XBT and the type."""
    maxDepth = -99
    WODCountryCode = '-99'
    WODCorrectionCode = -99
    validDepths = depths != depthFV
    if np.count_nonzero(validDepths) > 0:
        maxDepth = np.max(depths[validDepths])

    # Definition of types.
    codes = [[2  ,  -1 ,         0 ,    0], # Unknown
             [-1 ,  101,         1 ,    0], # Unknown
             [-1 ,  102,         2 ,    0], # Unknown
             [201,  -1 ,        70 ,    0], # T7 Brand unknown   
             [202,  -1 ,        40 ,    0], # T4 Brand unknown   
             [203,  -1 ,        60 ,    0], # T6 Brand unknown   
             [204,  -1 ,        50 ,    0], # T5 Brand unknown   
             [205,  -1 ,        100,    0], # T10 Brand unknown  
             [206,  -1 ,        110,    0], # T11 Brand unknown  
             [207,   41,        71 ,    1], # T7 Sippican        
             [207,   42,        72 ,    1], # T7 Sippican        
             [208,    1,        41 ,    1], # T4 Sippican             
             [208,    2,        42 ,    1], # T4 Sippican        
             [209,   31,        60 ,    1], # T6 Sippican        
             [209,   32,        60 ,    1], # T6 Sippican        
             [210,   11,        50 ,    1], # T5 Sippican        
             [211,   61,        100,    1], # T10 Sippican       
             [212,   71,        110,    1], # T11 Sippican       
             [-1 ,  900,        120,    1], # T12 Sippican       
             [213,   21,        130,    1], # Fast Deep Sippican 
             [214,   51,        141,    1], # Deep Blue Sippican 
             [214,   52,        142,    1], # Deep Blue Sippican 
             [-1 ,   81,        150,    1], # AXBT Sippican      
             [215,  201,        41 ,    2], # T4 TSK             
             [215,  202,        42 ,    2], # T4 TSK             
             [216,  211,        60 ,    2], # T6 TSK             
             [216,  212,        60 ,    2], # T6 TSK             
             [217,  221,        71 ,    2], # T7 TSK             
             [217,  222,        72 ,    2], # T7 TSK             
             [218,   -1,         0 ,    0], # MDI; Academy of Sc 
             [219,  231,        50 ,    2], # T5 TSK             
             [220,  241,        100,    2], # T10 TSK            
             [221,  401,        10 ,    3], # XBT-1 Sparton      
             [222,  411,        30 ,    3], # XBT-3 Sparton      
             [223,  421,        40 ,    3], # XBT-4 Sparton      
             [224,  431,        50 ,    3], # XBT-5 Sparton      
             [225,  441,        51,     3], # XBT-5DB Sparton    
             [226,  451,        60 ,    3], # XBT-6 Sparton      
             [227,  461,        70 ,    3], # XBT-7 Sparton      
             [-1 ,  462,        71 ,    3], # XBT-7 Sparton      
             [228,  471,        72 ,    3], # XBT-7DB Sparton    
             [229,  481,        100,    3], # XBT-10 Sparton     
             [230,  491,        200,    3], # XBT-20 Sparton     
             [231,  501,        201,    3], # XBT-20DB Sparton   
             [232,  251,        140,    2], # Deep Blue TSK      
             [232,  252,        140,    2], # Deep Blue TSK       
             [233,  261,        150,    2], # AXBT TSK           
             [234,  -1 ,        150,    0], # AXBT Unknown       
             [235,  -1 ,        140,    0], # Deep Blue Unknown  
             [236,  -1 ,        130,    0], # Fast Deep Unknown  
             [237,  -1 ,        160,    0], # SSXBT Sippican     
             [238,  -1 ,        150,    0]] # AXBT Sparton           
    codes = np.array(codes)

    notXBT = [np.array([-1, -1])]
    unkXBT = [np.array([0, 0])] 

    try:
        probeCode = int(instRef[60:63])
    except ValueError:
        probeCode = ''

    if np.any(salinities != salinityFV):
        return notXBT + [maxDepth, WODCountryCode, WODCorrectionCode, probeCode]
    
    if projectName[0:3] == 'WOD':
        # Set the column in the table to look at.
        codeIndex = 0

        # Check the instrument type code and return if can go no further.
        instCode = instRef[0:60]
        instCode = instCode.lstrip()
        instCode = instCode.rstrip()

        # Get some other information to return to user.
        WODCorrectionCode = int(instRef[63])
        WODCountryCode = projectName[8:10]

        # Extract the number for the probe type.
        probeCode = int(instRef[60:63])

        if instCode == '':
            if projectName[5:8] == 'XBT':
                return unkXBT + [maxDepth, WODCountryCode, WODCorrectionCode, probeCode]
            else:
                return notXBT + [maxDepth, WODCountryCode, WODCorrectionCode, probeCode]
        if instCode != '2':
            return notXBT + [maxDepth, WODCountryCode, WODCorrectionCode, probeCode]

    elif projectName[0:5] == 'GTSPP':
        # Set the column in the table to look at.
        codeIndex = 1

        # Get the number for the probe type.
        probeCode = int(instRef)

        # Check for situations where can go no further.
        if probeCode == 0 or probeCode == 99 or probeCode == 999:
            if projectName[5:7] == 'XB':
                return unkXBT + [maxDepth, WODCountryCode, WODCorrectionCode, probeCode]
            else:
                return notXBT + [maxDepth, WODCountryCode, WODCorrectionCode, probeCode]
    else: 
        return notXBT + [maxDepth, WODCountryCode, WODCorrectionCode, probeCode]

    matches = codes[:, codeIndex] == probeCode
    nMatches = np.count_nonzero(matches)
    if nMatches == 0:
        return notXBT + [maxDepth, WODCountryCode, WODCorrectionCode, probeCode]

    iMatches = np.argwhere(matches)
    iMatches = iMatches[0] # Repeated codes for WOD data.
    result = np.reshape(codes[iMatches, 2:4], 2)

    # Finally, follow prescription in EN3 processing to
    # pick out mislabelled T5s. These will have 5000 added
    # to their code. 840 m max depth is used as depth 
    # criterion as this is as used in the EN3 processing and
    # is also max depth of T7s (Tim Boyer; personal communication).
    # Deep Blues can now go to 920 m (Tim Boyer; personal 
    # communication) and so these are not tested here.
    typeCode = result[0] // 10
    notSparton = result[1] < 3
    if notSparton:
        if (typeCode == 4 or   # T4 max depth = 460m
            typeCode == 6 or   # T6 max depth = 460m
            typeCode == 7 or   # T7 max depth = 760m
            typeCode == 10 or  # T10 max depth = 200m
            typeCode == 11):   # T11 max depth = 460m
            if maxDepth > 840.0:
                result[0] += 5000

    return [result] + [maxDepth, WODCountryCode, WODCorrectionCode, probeCode]

def is_mbt(projectName,
           instRef):
    """Identify MBTs."""
    notMBT = np.array([-1, -1])
    unkMBT = np.array([0, 0])

    codes = [[  1,  800,        00,    00], # MBT type/make unknown.
             [101,   -1,        01,    01]] # GM39 (Russia).     
    codes = np.array(codes)   

    if projectName[0:3] == 'WOD':
        codeIndex = 0
        instCode = instRef[0:60]
        instCode = instCode.lstrip()
        instCode = instCode.rstrip()
        if instCode == '':
            if projectName[5:8] == 'MBT':
                return unkMBT
            else:
                return notMBT
        elif instCode != '1':
            return notMBT
        probeCode = int(instRef[60:63])     
    elif projectName[0:5] == 'GTSPP':
        codeIndex = 1
        probeCode = int(instRef)
        if probeCode == 0 or probeCode == 999:
            if projectName[5:7] == 'MB':
                return unkMBT
            else:
                return notMBT
    else:
        return notMBT  
     
    matches = codes[:, codeIndex] == probeCode
    nMatches = np.count_nonzero(matches)
    if nMatches == 0:
        return notMBT

    iMatches = np.argwhere(matches)
    iMatches = iMatches[0] 
    result = np.reshape(codes[iMatches, 2:4], 2)

    return result

def test():
    file = '/data/local/hadgs/Data/EN3_v2a_NoCWT/Profiles/EN3_v2a_NoCWT_Profiles_199501.nc'
    ncid = Dataset(file)
    pn = netCDF4.chartostring(ncid.variables['PROJECT_NAME'][:])
    ir = netCDF4.chartostring(ncid.variables['INST_REFERENCE'][:])
    ps = ncid.variables['PSAL_CORRECTED'][:]
    psfv = ncid.variables['PSAL_CORRECTED']._fillvalue
    de = ncid.variables['DEPH_CORRECTED'][:]
    defv = ncid.variables['DEPH_CORRECTED']._fillvalue
    ncid.close()

    f = open('test.txt', 'w')
    WODCountryCode = -1
    for i in np.arange(pn.size):
        vals = is_xbt(pn[i], ir[i], ps[i, :], 
                      psfv, de[i, :], defv)
        f.write('%i %i %f %s %i\n' % (vals[0][0], vals[0][1], vals[1], vals[2], vals[3]))
        vals = is_mbt(pn[i], ir[i])
        f.write('%i %i\n' % (vals[0], vals[1]))

    return
