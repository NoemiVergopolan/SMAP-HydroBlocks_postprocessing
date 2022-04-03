#!/usr/bin/env python
# coding: utf-8

# Author: Noemi Vergopolan - Princeton University
# Contact: noemi@princeton.edu, noemi.v.rocha@gmail.com
# Last Update: Arpil 1st, 2022
#
# Summary:
#
#      This script post-process the SMAP-HydroBlocks dataset from the Hydrological Response Unit (HRU) space
#      into the geographically gridded space. The final domain extent, temporal coverage, and temporal resolution
#      are user-defined parameters. The data uses a Plate Carrée projection and can be output as Zarr storage
#      (recommended) or NetCDF files. Both can be open with python libraries such as xarrays.
#
# Usage:
#
#      1. Download the SMAP-HydroBlocks Database: wget https://zenodo.org/record/5206725/files/SMAP-HB_hru_6h.zip
#      2. Unzip the file: unzip SMAP-HB_hru_6h.zip
#      3. Create environemnt: conda env create --name remapping -f yml/default/environment.yml
#      4. Load environment: source activate remapping
#      3. Run: python ./SMAPHB_hru2grid.py
#
# Note:
#      SMAP-HydroBlocks is a very big dataset. If remapped entirely it comprises>122 TB of NetCDF data. To reduce
#      data storage and improve data access, outputting the data in Zarr format is recommended.
#      This script allows for subsetting and postprocessing this dataset according to the user's needs and resources.
#      For example, SMAP-HydroBlocks mapped to a 30-meter daily resolution, at ~ 50km x 50km domain, for a 1 year
#      period, with maximum compression [option 9], is expected to output 786 MB of NetCDF files or 656 MB of Zarr
#      storage. As such, please keep in mind that the domain extent, time period, temporal resolution, and compression
#      option selected will determine runtime and data storage
#
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
import datetime
import SMAPHB_functions as SMAPHB_functions

# 1. Define workspace folder
workspace = os.getcwd() # current folder

# Database path
database_path = '%s/database' % workspace
variable = 'SMAPHB_SM'

# 2. Output folder name
output_folder = '%s/SMAPHB_sample' % workspace

# 3. Output file format. Options: 'zarr', 'netcdf', or 'both'
# Netcdf is only recommended for small subsets of the data. Alternatively, please use zarr. Both netcdf and zarr files can be open with xarrays.
output_data_format = 'netcdf'

# 4. Set the domain of interest in geographic coordinates
minlat = 37.0
maxlat = 37.5
minlon = -91.7
maxlon = -91.2

# 5. Set the output time period of interest. It must be between 1 April 2015 and 31 December 2019 (SMAP-HB retrieval dates)
start_date = datetime.date(2016,1,1)
end_date = datetime.date(2016,12,31)

# 6. Set the output time resolution of interest. Original retrievals are at '6h' temporal resolution, but data can also be output at 'daily', 'monthly', 'annual' temporal resolution.
final_temporal_resolution = 'daily' 

# 7. Set the output spatial resolution. This is an experimental feature that regrids data to spatial resolutions > 30 meters. This uses xESMF, but only works without MPI support (mpi_run = False)
final_spatial_resolution = 30   # meters

# 8. Set the output data compression level [0-9]. More compression saves storage but will slow down data procressing.
compression_level = 9

# The remapping can be performed in parallel with mpi_run = True and:
# 'mpirun -np <number of processes> python SMAPHB_hru2grid.py'
# In the future, it would be nice to have dask parallelization instead.
mpi_run = False

# Retrive data
SMAPHB_functions.retrieve_data(database_path, variable,
                  output_folder, output_data_format,
                  final_temporal_resolution,
                  start_date, end_date,
                  minlat, maxlat, minlon, maxlon,
                  final_spatial_resolution,
                  compression_level,
                  mpi_run)

