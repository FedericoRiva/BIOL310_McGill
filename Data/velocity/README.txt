Release date:
6/17/2016

Last update:
6/17/2016

Authors: 
Sean Parks
Rocky Mountain Research Station
US Forest Service
sean_parks@fs.fed.us

Solomon Dobrowski
College of Forestry and Conservation
University of Montana
solomon.dobrowski@umontana.edu

File name: 
Dobrowski.and.Parks.Data.7z

Brief description: 
This 7z file contains the raster products that were used in the paper cited below.

Citation:
Dobrowski and Parks. 2016. Climate change velocity underestimates climate change exposure in mountainous regions. Nature Communications.

This 7z file contains gridded datasets depicting:
1. Minimum exposure distance (MED) (filename = ’MED.tif’)
2. MED-based climate velocity (filename = ‘velocity.MED.tif’)
3. Minimum cumulative exposure (MCE) (filename = ‘MCE.tif’)
4. ED-based velocity (ED = Euclidean distance) (filename = ‘velocity.ED.tif’)
5. Ratio of MED:ED (filename = ‘ratio.MED.to.ED.tof’)

More in-depth description of rater datasets (please see Dobrowski and Parks, 2016 for further details):
1. Minimum exposure distance (MED) – this dataset represents the distance (km) from each pixel to its nearest climate analog. The values do not represent Euclidean distance (ED), but instead represent the distance that minimizes exposure to dissimilar climates based on a resistance surface that penalized dissimilar climates. Pixels on islands and whose future climate analog is on an island are assigned a value = -1. The reference climate is represented my mean annual temperature for the years 1981-2010 and the future climate is represented by mean annual temperature for the years 2071-2100. Future climate is an ensemble of 15 CMIP5 GCMs and use the RCP 8.5 emissions scenario. All climate datasets were obtained from AdaptWest (https://adaptwest.databasin.org/). 

2. MED-based climate velocity – this dataset is simply MED/90; the units are km/year. MED is the previously described dataset and 90 represents the number of years between the reference period (1995) and the future climatic conditions (2085). Simply put, this dataset represents climate change velocity. Pixels on islands and whose future climate analog is on an island are assigned a value = -1. See figure 2a in Dobrowski and Parks. 

3. Minimum cumulative exposure (MCE) – this dataset represents the exposure to dissimilar climates (in degrees C) along trajectories between each pixel and is nearest climate analog. These trajectories were built on the assumption that organisms would minimize exposure to dissimilar climates. Pixels on islands and whose future climate analog is on an island are assigned a value = -1. See figure 3 in Dobrowski and Parks.

4. ED-based velocity – this dataset represents the velocity (km/year) based on Euclidean distance (ED). ED-based velocity minimizes distance between each pixel and its future climate analog. This dataset is provided to contrast MED- and ED-based estimates of climate change velocity. Pixels on islands and whose future climate analog is on an island are assigned a value = -1. For more information on ED-based velocity, please see Hamman et al. 2015. Velocity of climate change algorithms for guiding conservation and management. Global Change Biology. 21: 997-1004.

5. Ratio of MED to ED – This dataset represents the ratio between MED and ED; it also represent the ratio between MED-based velocity and ED-based velocity. The two comparisons are by definition identical. This dataset is provided to compare MED- to ED-based approaches for characterizing climate change velocity. Pixels on islands and whose future climate analog is on an island are assigned a value = -1. See figures 2b and 2c in Dobrowski and Parks.

Gridded dataset properties:
File format – geotiffs; each *.tif also contains a world file (*.tfw)
Resolution – 5km
Projection – LAMBERT
Datum – WGS84
Spheroid – WGS84
Units – METERS
Zunits – NO
Xshift – 0.0
Yshift – 0.0
Parameters
49  0  0.0 /* 1st standard parallel
77  0  0.0 /* 2nd standard parallel
-95  0  0.0 /* central meridian
0  0  0.0 /* latitude of projection's origin
0.0 /* false easting (meters)
0.0 /* false northing (meters)


Other pertinent information:
1. Sample code and data are also provided (http://adaptwest.databasin.org/pages/adaptwest-velocitymed). 

2. Pixels on islands and whose future climate analog is on an island are assigned a value = -1 in all datasets
