# SimpleGrid
SimpleGrid is a package for efficient aggregation of ocean profile data onto a regular x-y grid. This tool was designed to work with profiles in the netcdf format used by the [EN4 database][EN4-ref] or output from the package SynthPro.


## Using SimpleGrid
#### Required python libraries
SimpleGrid was developed using Python 2.7 and requires the installation of the following python libraries and their associated dependencies: `numpy`, `netCDF4`, `ConfigParser`, `argparse`, `sys`, `shutil`, `os`, `time`, and `scipy`.

#### Cloning the git repository
To retrieve a copy of the SimpleGrid source code and create a working directory, run the following on the command line: 

```> git clone git@github.com:cdr30/SimpleGrid.git```

or 

```> git clone https://github.com/cdr30/SimpleGrid.git```


#### Running SimpleGrid
SimpleGrid is invoked from the command line using the `run_simplegrid.py` script. Executing this script without any arguments will return an error message that demonstrates the correct usage:

```
> python2.7 run_simplegrid.py 
usage: run_simplegrid.py [-h] namelist
run_simplegrid.py: error: too few arguments
```

The namelist argument is a path to an [INI][INI-ref] configuration file containing data paths and other options. An example namelist is available in the `./config` subdirectory and an annotated version is described the section below. 


### Configuring a `namelist.ini` file
This section provides an annotated examples of a SimpleGrid namelist congiguration file. 

###### `[profiles]`
```
dir = /profile/dir/  						### Directory containing profile data
fpattern = ${YYYY}${MM}profiles.nc 			### File pattern for profile netcdf files.
datavar = POTM_CORRECTED 					### Variable name of data to be gridded.
zvar = DEPH_CORRECTED 						### Variable name for z-coord.
yvar = LATITUDE 							### Variable name for y-coord.
xvar = LONGITUDE 							### Variable name for x-coord.
minyr = 2000 								### Start year for analysis.
maxyr = 2015 								### End year of analysis.
qcvar = POTM_CORRECTED_QC 					### Variable name for profile quality control flag.
posqcvar = POSITION_QC 						### Variable name for position quality control flag.
qcreject = 4 								### Reject value for profile quality control flag.
posqcreject = 4 							### Reject value for position quality control flag.
```

##### `[climatology]`
```
minyr=2000 		  							### First year to use for calculation of climatology.
maxyr=2015									### Last year to use for calculation of climatology.

##### `[grid]`
dir = /output/dir/							### Output directory
xmin=-180									### Minimum boundary for x-coord. Used to generate grid.
xmax=180									### Maximum boundary for x-coord. Used to generate grid.
ymin=-90									### Minimum boundary for y-coord. Used to generate grid.
ymax=90										### Maximum boundary for y-coord. Used to generate grid.
dx=2										### Grid spacing in x direction.
dy=2										### Grid spacing in y direction.
zbounds=0.0,10.0,20.1,30.2,40.4,50.6,60.9,71.3,81.9,92.7,104.0,115.8,128.4,142.0,157.3,174.8,195.5,220.6,251.9,291.6,342.4,407.6,490.9,595.9,725.2,880.5,1061.8,1267.5,1494.7,1739.9,1999.4,2270.0,2548.9,2833.9,3123.3,3415.9,3710.7,4007.0,4304.5,4602.7,4901.5,5200.6,5500										### Comma-separated list of boundaries for z-coordinate.
```

##### `[anomalies]`
```
calc_anomalies = True					    ### If True, calculate anomalies relative to the specified climatology.
```

## Infilling missing data 

An example script for the interpolation of missing data using an inverse-distance weighting approach is provided in the ```scripts``` directory. Help is available from the command line: 

 ```
>python2.7 infill_gridded_data.py -h 
usage: infill_gridded_data.py [-h] [--bathyf BATHYF] [--bathyvar BATHYVAR]
                              [-c CUTOFF] [-p POWER] [-s SMOOTH]
                              namelist gridf infillf
Infill missing data in gridded data using an inverse-distance weighting
approach.
positional arguments:
  namelist              Path to namelist.ini file used to create gridded data.
  gridf                 Path to netcdf file containing gridded data.
  infillf               Output file containing infilled data.
optional arguments:
  -h, --help            show this help message and exit
  --bathyf BATHYF       Path to netcdf file containing bathymetry data for
                        estimation of land-sea mask
  --bathyvar BATHYVAR   Variable name for bathymetry data.
  -c CUTOFF, --cutoff CUTOFF
                        Maximum distance (km) used for interpolation.
  -p POWER, --power POWER
                        Power used for inverse-distance weighting.
  -s SMOOTH, --smooth SMOOTH
                        If True, apply inverse-distance weighting everywhere
                        to smooth field.
```
