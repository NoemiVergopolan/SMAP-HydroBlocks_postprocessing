# Summary
[SMAP-HydroBlocks (SMAP-HB)](https://waterai.earth/smaphb/) is a hyper-resolution satellite-based surface soil moisture product that combines NASA's Soil Moisture Active-Passive (SMAP) L3 Enhance product, hyper-resolution land surface modeling, radiative transfer modeling, machine learning, and in-situ observations. This dataset reports the surface soil moisture (top 5-cm of the soil layer) in volumetric unit (m3/m3). Data is available across the continental United States at an effective 30-m spatial resolution (2015–2019). Retrievals are organized 6h, but are only reported at the time of SMAP retrieval. This script post-process the SMAP-HydroBlocks dataset from the Hydrological Response Unit (HRU) space (time, hru) into the geographically gridded space (time, latitude, longitude) using the Plate Carrée projection. Data can be output in Zarr and/or netCDF4 format. 

# Usage

### 1. Clone this repository
```
git clone https://github.com/NoemiVergopolan/SMAP-HydroBlocks_postprocessing.git
cd SMAP-HydroBlocks_postprocessing
```

### 2. Download the SMAP-HydroBlocks Database from https://doi.org/10.5281/zenodo.5206725 and unzip it
```
wget https://zenodo.org/record/5206725/files/SMAP-HB_hru_6h.zip
unzip SMAP-HB_hru_6h.zip
```

### 3. Create an anaconda environment

Option 1 – Create a conda environment from the yml file
```
conda env create --name SMAPHB -f yml/github/environment.yml
source activate SMAPHB
```

Option 2 – Create a conda environment yourself:
```
conda create -n SMAPHB -y
source activate SMAPHB
conda install -c conda-forge numpy xarray rioxarray pandas dask netcdf4 zarr mpi4py xesmf esmpy
```

### 4. Run the script
```
python ./SMAPHB_hru2grid.py
```

Gridded SMAP-HydroBlocks data will be saved on ```SMAPHB_sample``` folder. Please edit the ```SMAPHB_hru2grid.py``` script to change the desired spatial extent, time period, time resolution, data format, compression level, output file and folder name, etc.



### Support for MPI-Parallel via mpirun
The remapping can be performed in a multi-core setup with ```mpi_run = True``` in ```SMAPHB_hru2grid.py``` and run as:
```
mpirun -np <number of processes> python ./SMAPHB_hru2grid.py
```


### Support for remapping to coarser spatial resolutions
There is an experimental feature that allows remapping the SMAP-HydroBlocks data at 30-m resolution to a coarser spatial resolution (e.g., 100m, 500m, 1000m, etc.). This feature relies on the xESMF library, which currently does not provide MPI parallel support. Thus, to remap SMAP-HydroBlocks data to a coarse spatial resolution, please edit ```SMAPHB_hru2grid.py``` and set ```mpi_run = False``` and  ```final_spatial_resolution = <desire resolution in meters>```. Run the script as:
```
python ./SMAPHB_hru2grid.py
```

### Notes
SMAP-HydroBlocks is a very big dataset. If remapped entirely at a daily scale comprises ~120TB of data in maximum compression (option 9). This script allows for post-processing subsets of the dataset according to the user's needs and resources. For example, SMAP-HydroBlocks mapped to a 30-meter daily resolution, at ~ 50km x 50km domain, for a 1 year period, with maximum compression [option 9], is expected to output 786 MB of NetCDF files or 656 MB of Zarr storage. As such, please keep in mind that the domain extent, time period, temporal resolution, and compression option selected will determine runtime and data storage requirements.


# Reference

Please cite the following paper when using the dataset in any publication:

Vergopolan, N. et al. SMAP-HydroBlocks, a 30-m satellite-based soil moisture dataset for the conterminous US. Sci Data 8, 264 (2021). https://doi.org/10.1038/s41597-021-01050-2



# Contact
 - Dr. Noemi Vergopolan, Rice University
 - Website: www.waterai.earth/smaphb
 - Email: vergopolan@rice.edu
