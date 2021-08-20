#!/usr/bin/env python
# coding: utf-8

# Author: Noemi Vergopolan - Princeton University
# Contact: noemi@princeton.edu, noemi.v.rocha@gmail.com
# Date: 18th August 2021
# 
# Summary:
#
#      This script post-process the SMAP-HydroBlocks dataset from the Hydrological Response Unit (HRU) space 
#      into the geographically gridded space. The data output is at 30-m 6-h resolution, using a Plate Carrée 
#      projection and stored in netCDF4 format. The SMAP-HydroBlocks dataset reports the top 5-cm surface
#      soil moisture in volumetric units (m3/m3).
#
# Usage: 
#   
#      1. Download the SMAP-HydroBlocks Database at https://doi.org/10.5281/zenodo.5206725
#      2. Update database folder path, and desired data extent, period, compression level
#      3. Run: python ./SMAPHB_hru2grid.py
#
# Note: 
#      SMAP-HydroBlocks is a very big dataset. If remaped entirely it comprises of 22TB of data in 
#      maximum compression (option 9) or 600 TB with no compression (option 0). This script allows 
#      for subsetting and postprocessing this dataset according to user's needs and resources. 
#      For example, SMAP-HydroBlocks at 30-meter 6-hour resolution at a 10-km by 10-km box extent 
#      over a 1 year period is expected to output 55 MB of data using maximum compression (option 9) 
#      or 1.5 GB of data with no compression (option 0). As such, please keep in mind that the domain
#      extent, time period, and compression option selected will determine the output file size and 
#      running time.
#
#      If you run out of memmory, please consider:
#          1. Reduce the domain extent
#          2. Reduce the time period
#          3. Perform data post-processing in batches
#          4. Deploy this script in a HPC system
#
# Reference:
# 
#       Please cite the following paper when using the dataset in any publication:
#       Vergopolan, N., Chaney, N. W., Beck, H. E., Pan, M., Sheffield, J., Chan, S., & Wood, E. F. (2020).
#             Combining hyper-resolution land surface modeling with SMAP brightness temperatures to 
#             obtain 30-m soil moisture estimates. Remote Sensing of Environment, 242, 111740.
#             https://doi.org/10.1016/j.rse.2020.111740
#



# Load libraries needed
import os
import numpy as np
import xarray as xr
import datetime
import SMAPHB_functions 



# 1. Download the SMAP-HydroBlocks database from:
# Database path - Edit the workspace path according to where in your computer you save the database 
database_path = 'database'

# 2. Output file name
output_file = 'SMAPHB_sample.nc'

# 3. Define the domain of interest in geographic coordinates
# If you have limited computing memmory, you may want to keep the domain small

# In this example, we are going to look at the Little River Basin
minlat = 31.465 
maxlat = 31.76 
minlon = -83.78 
maxlon = -83.51     

# 4. Define dates of interest. It must be between 31 March 2015 and 31 December 2019 (SMAP retrieval dates)
start_date = "2016-06-01"
end_date = "2016-06-31"

# Define compression level [0-9]. If ouput file is too big, increase compression
compression_level = 9

# 5. Output dataset
SMAPHB_functions.save_data(database_path, 
                           minlon, maxlon, minlat, maxlat, 
                           start_date, end_date,
                           compression_level,
                           output_file
                          )



