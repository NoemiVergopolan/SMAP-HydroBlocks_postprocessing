# Summary
SMAP-HydroBlocks (SMAP-HB) is a hyper-resolution satellite-based surface soil moisture product that combines NASA's Soil Moisture Active-Passive (SMAP) L3 Enhance product, hyper-resolution land surface modeling, radiative transfer modeling, machine learning, and in-situ observations. This dataset was developed over the continental United States at 30-m 6-hourly resolution (2015–2019). This script post-process the SMAP-HydroBlocks dataset from the Hydrological Response Unit (HRU) space into the geographically gridded space. The data output is at 30-m 6-h resolution, using a Plate Carrée projection and stored in netCDF4 format. The SMAP-HydroBlocks dataset reports the top 5-cm surface soil moisture in volumetric units (m3/m3).

# Usage
1. Download the SMAP-HydroBlocks Database at X
2. Update database folder path, desired data extent, period, compression level in the SMAPHB_hru2grid.py file
3. Make sure to have installed python and the following libraries: numpy, xarray, rasterio, datetime
4. Run: ```python ./SMAPHB_hru2grid.py```

### Notes
SMAP-HydroBlocks is a very big dataset. If remaped entirely it comprises of 22TB of data in maximum compression (option 9) or 600 TB with no compression (option 0). This script allows for subsetting and postprocessing this dataset according to user's needs and resources. For example, SMAP-HydroBlocks at 30-meter 6-hour resolution at a 10-km by 10-km box extent over a 1 year period is expected to output 55 MB of data using maximum compression (option 9) or 1.5 GB of data with no compression (option 0). As such, please keep in mind that the domain extent, time period, and compression option selected will determine the output file size and running time.

### If you run out of memmory, please consider:
 - Reduce the domain extent
 - Reduce the time period
 - Perform data post-processing in batches
 - Deploy this script in a HPC system

# Reference

Please cite the following paper when using the dataset in any publication:

Vergopolan, N., Chaney, N. W., Beck, H. E., Pan, M., Sheffield, J., Chan, S., & Wood, E. F. (2020). Combining hyper-resolution land surface modeling with SMAP brightness temperatures to obtain 30-m soil moisture estimates. Remote Sensing of Environment, 242, 111740. https://doi.org/10.1016/j.rse.2020.111740

# Contact
 - Noemi Vergopolan, Princeton University
 - Website: www.waterai.earth/smaphb
 - Email: noemi@princeton.edu
