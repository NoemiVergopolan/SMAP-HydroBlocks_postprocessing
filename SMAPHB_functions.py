import os
import sys
import numpy as np
import xarray as xr
from pandas import Series as pd_Series
import datetime
from dateutil.relativedelta import relativedelta
import warnings
import gc
import dask.array as da
import zarr
from itertools import product

# Hopes that it helps minimize memory leak when open and reading many netcdf files
# It limits the number of files that can be simultaneously opened
xr.set_options(file_cache_maxsize=1)


def retrive_time_stepping(database_path, start_date, end_date, final_temporal_resolution):

    # check if variable names are correct
    if final_temporal_resolution not in ['6h', 'daily', 'monthly', 'annual']:
        sys.exit("Error: date_time_step = %s not available. Please use one: '6h', 'daily', 'monthly', 'annual'" % final_temporal_resolution)

    ic = 0
    infile = '%s/SMAPHB_hru_6hr/%i.nc' % (database_path, ic)
    while os.path.exists(infile) == False:
        infile = '%s/SMAPHB_hru_6hr/%i.nc' % (database_path, ic)
        ic = ic +1
    ds = xr.open_dataset(infile, cache=False)
    ds = ds.sel(time=slice(start_date,end_date))
    times = ds.time.values

    ds = xarray_time_agregation(ds, final_temporal_resolution)
    ag_times = ds.time.values
    ds.close()
    del ds

    return times, ag_times

def update_latlon_to_new_spatial_resolution(final_spatial_resolution, lats, lons):

    data_res = 27.7777777 # meters
    lats = lats[::-1]
     
    if final_spatial_resolution == 30:
        ag_lats = lats
        ag_lons = lons
        ag_delta_lat = (lats[-1]-lats[0])/(lats.size-1)
        ag_delta_lon = (lons[-1]-lons[0])/(lons.size-1)

    else: 
        lat_delta = (lats[-1]-lats[0])/(lats.size-1)
        lat_delta_half = lat_delta/2.0
        lon_delta = (lons[-1]-lons[0])/(lons.size-1)
        lon_delta_half = lon_delta/2.0

        ag_delta_lat = lat_delta*final_spatial_resolution/data_res
        ag_delta_lat_half = ag_delta_lat/2.0
        ag_delta_lon = lon_delta*final_spatial_resolution/data_res
        ag_delta_lon_half = ag_delta_lon/2.0

        ilat = lats[0]  -lat_delta_half
        flat = lats[-1] +lat_delta_half
        ilon = lons[0]  -lon_delta_half
        flon = lons[-1] +lon_delta_half

        n_ag_lats = int(np.ceil((flat-ilat)/ag_delta_lat))
        n_ag_lons = int(np.ceil((flon-ilon)/ag_delta_lon))

        # Not pythonic elegant, but minimize truncation error for small deltas
        ag_lats = [ilat+ag_delta_lat_half]
        for i in range(n_ag_lats-1): ag_lats.append(ag_lats[-1]+ag_delta_lat)
        ag_lats = np.array(ag_lats)

        ag_lons = [ilon+ag_delta_lon_half]
        for i in range(n_ag_lons-1): ag_lons.append(ag_lons[-1]+ag_delta_lon)
        ag_lons = np.array(ag_lons) 

    return ag_lats, ag_lons, ag_delta_lat, ag_delta_lon

def update_chunks_to_new_spatial_resolution(final_spatial_resolution,chunk):
    data_res = 27.7777777 # meters
    if final_spatial_resolution == 30:
        pass
    else:
        chunk_factor = final_spatial_resolution/data_res
        chunk['lat'] = int(np.ceil(chunk['lat']/chunk_factor))
        chunk['lon'] = int(np.ceil(chunk['lon']/chunk_factor))
    return chunk

def create_zarr_template(final_path, variable, lats, lons, final_spatial_resolution, chunk, 
                        ag_times, final_temporal_resolution, compression_level):

    nlon = len(lons)
    nlat = len(lats)
    nt = len(ag_times)

    # Create a uninitialize xarray template
    template = xr.DataArray(da.empty((nt,nlat,nlon), dtype=np.float32, compute=False),
                            coords=[ag_times, lats, lons],
                            dims=["time", "lat", "lon"])
    template = template.to_dataset(name=variable)

    # Defin e and rechunk data
    template = template.chunk(chunk)

    # Grab chunk ranges
    clats = template.chunks['lat']
    clons = template.chunks['lon']
    ctimes = template.chunks['time']

    # Double chunk in time
    chunk['time'] = int(np.ceil(chunk['time']/10.))
    template = template.chunk(chunk)
    
    # Define atributes
    attrs = dict(unit = 'm3/m3',
                 title = "SMAP-HydroBlocks Surface Soil Moisture Data (m3/m3)",
                 description = 'SMAP-HydroBlocks (SMAP-HB) is a 30-m hyper-resolution satellite-based surface soil moisture product (2015-2019). The dataset combines NASA Soil Moisture Active-Passive (SMAP) L3 Enhance product, hyper-resolution land surface modeling, radiative transfer modeling, machine learning, and in-situ observations. This subset was mapped to %i-m %s resolution using geographic coordinates and Plate Carrée projection.' % (final_spatial_resolution, final_temporal_resolution),
                 creator_name = "Noemi Vergopolan (noemi@princeton.edu)",
                 institution = 'Princeton University',
                 citation = "Vergopolan et al. (2021). SMAP-HydroBlocks, a 30-m satellite-based soil moisture dataset for the conterminous US. Scientific Data, 8, 264. https://doi.org/10.1038/s41597-021-01050-2",
                 )
    template = template.assign_attrs(attrs)

    # Create zarr template on the disk
    os.system('rm -rf %s' % (final_path))
    zarr_compressor = zarr.Blosc(cname="zstd", clevel=compression_level, shuffle=1)
    zarr_encoding = { variable    : {'_FillValue': -9999,
                                    'compressor': zarr_compressor,
                                    'chunks': (chunk['time'],chunk['lat'],chunk['lon'])},
                                    'time': {'dtype': 'i4'},
                    }
    template.to_zarr(final_path, encoding=zarr_encoding, compute=False, consolidated=True, mode='w')
    print(datetime.datetime.now(), 'Data template is ready', flush=True)

    template.close()
    del template

    regions = {}

    # define slicing regions
    regions['lon_slice'] = [(i*clons[0],i*clons[0]+clons[i]-1) for i in range(len(clons))]
    regions['lat_slice'] = [(i*clats[0],i*clats[0]+clats[i]-1) for i in range(len(clats))]
    regions['time_slice'] = [(i*ctimes[0],i*ctimes[0]+ctimes[i]-1) for i in range(len(ctimes))]

    # defining regions ranges
    regions['lat_range'] = [(lats[region[0]],lats[region[1]])
                             for region in regions['lat_slice']]
    regions['lon_range'] = [(lons[region[0]],lons[region[1]])
                             for region in regions['lon_slice']]
    regions['time_range'] = [( ag_times[region[0]],
                               get_final_time( ag_times[region[1]], final_temporal_resolution)
                             ) for region in regions['time_slice']]

    return regions

def get_final_time(date, final_temporal_resolution):
    if final_temporal_resolution == '6h':
        pass
    else:
        date = pd_Series(date)[0]
        if final_temporal_resolution == 'daily':
            date = date + relativedelta(days=1) - relativedelta(seconds=1)
        if final_temporal_resolution == 'monthly':
            date = date + relativedelta(months=1) - relativedelta(seconds=1)
        if final_temporal_resolution == 'annual':
            date = date + relativedelta(years=1) - relativedelta(seconds=1)
        date = np.datetime64( date, 'ns')
    return date

def define_output_file_name_zarr(final_spatial_resolution, final_temporal_resolution):
    core_name = 'SMAP-HB_surface-soil-moisture'
    output_file = '%s_%im_%s' % (core_name, final_spatial_resolution, final_temporal_resolution)
    return output_file

def open_mosaic_object(database_path, minlat, maxlat, minlon, maxlon):

    # check the lat lon inputs
    if maxlat < minlat: sys.exit('minlat and maxlat inputs may be inverted!')
    if maxlon < minlon: sys.exit('minlon and maxlon inputs may be inverted!')

    # Read the catchments and HRU maps
    catch_map = xr.open_rasterio(filename = '%s/mapping_catchments/all_catchments.vrt' % database_path,
                                 chunks   = {'x': 1000, 'y': 1000})
    catch_map = catch_map.sel(x=slice(minlon,maxlon),y=slice(maxlat,minlat))
    catch_map = catch_map.sel(band=1)

    hru_map = xr.open_rasterio(filename ='%s/mapping_hrus/all_hrus.vrt' % database_path,
                               chunks   = {'x': 1000, 'y': 1000})
    hru_map = hru_map.sel(x=slice(minlon,maxlon),y=slice(maxlat,minlat))
    hru_map = hru_map.sel(band=1)

    return catch_map, hru_map

def retrieve_mosaic(catch_map, hru_map, lat_range, lon_range, lat_delta, lon_delta):

    minlat = lat_range[0]-lat_delta/2
    maxlat = lat_range[1]+lat_delta/2
    minlon = lon_range[0]-lon_delta/2
    maxlon = lon_range[1]+lon_delta/2

    # Subset the catchments and HRU maps
    catch_map = catch_map.sel(x=slice(minlon,maxlon),y=slice(maxlat, minlat))
    catch_map = catch_map.values

    hru_map = hru_map.sel(x=slice(minlon,maxlon),y=slice(maxlat, minlat))
    xlons = hru_map.x.values
    ylats = hru_map.y.values
    hru_map = hru_map.values

    # retrieve unique catchments ID
    unique_catchments = np.unique(catch_map.astype(np.int16).ravel())
    unique_catchments = unique_catchments[unique_catchments != -9999]

    return catch_map, hru_map, unique_catchments, ylats, xlons

def xarray_time_agregation(ds, final_temporal_resolution):
    if final_temporal_resolution == '6h':
        pass
    else:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if final_temporal_resolution == 'daily':
                ds = ds.resample(time='1D').reduce(np.nanmean)
            if final_temporal_resolution == 'monthly':
                ds = ds.resample(time='1MS').reduce(np.nanmean)
            if final_temporal_resolution == 'annual':
                ds = ds.resample(time='1YS').reduce(np.nanmean)
    return ds

def map_data(database_path, variable,
              region_lats, region_lons,
              #lat_range, lat_delta,
              #lon_range, lon_delta,
              spatial_regrid,
              time_range, final_temporal_resolution,
              catch_map, hru_map, unique_catchments,
              final_spatial_resolution):
    
    start_date = time_range[0]
    end_date = time_range[1]

    # Define time array
    infile = '%s/SMAPHB_hru_6hr/%i.nc' % (database_path, unique_catchments[0])
    ds = xr.open_dataset(infile, cache=False)
    ds = ds.sel(time=slice(start_date,end_date))
    times = ds.time.values

    # Define interpolated time array
    ds = xarray_time_agregation(ds, final_temporal_resolution)
    ag_times = ds.time.values

    ds.close()
    del ds

    # Create and allocate data
    nt = len(ag_times)
    nlat = len(region_lats)
    nlon = len(region_lons)
    total_size = nt*nlat*nlon
    new_data = np.full((nt,nlat,nlon), fill_value=np.nan, dtype=np.float32)

    # Loop over catchments to remap data
    for i, icatch in enumerate(unique_catchments):

        # read data at the HRU-space
        infile = '%s/SMAPHB_hru_6hr/%i.nc' % (database_path,icatch)
        ds = xr.open_dataset(infile, cache=False)
        ds = ds.sel(time=slice(start_date,end_date)).load()
        ds = xarray_time_agregation(ds, final_temporal_resolution)
        sm_data = ds[variable].values
        hrus = ds.hru.values
        ds.close()
        del ds
   
        # if no data in this time-step, continue the loop
        if np.all(np.isnan(sm_data)):
            continue

        # Let's subset the catchment and hru masks, so we can work more efficiently on smaller areas
        mask_position = np.where(catch_map == icatch)
        cimin = np.min(mask_position[0])
        cimax = np.max(mask_position[0])
        cjmin = np.min(mask_position[1])
        cjmax = np.max(mask_position[1])
        mini_mask_catchment = (catch_map[cimin:cimax+1,cjmin:cjmax+1] == icatch)
        mini_hru_map = hru_map[cimin:cimax+1,cjmin:cjmax+1]
        mini_new_data = new_data[:,cimin:cimax+1,cjmin:cjmax+1]

        #print('Loop over HRU')
        # Loop through the HRUs and remap the data
        for ihru in hrus:
            mini_mask_hrus = mini_mask_catchment & (mini_hru_map == ihru)
            if np.any(mini_mask_hrus):
                position = np.where(mini_mask_hrus)
                mini_new_data[:, position[0], position[1]] = sm_data[:,ihru][:, np.newaxis]            

    # Create xarray data 
    final_map_xr = xr.DataArray(new_data, 
                                coords=[ag_times, region_lats, region_lons], 
                                dims=["time", "lat", "lon"])
    final_map_xr.attrs["units"]="m3/m3"
    final_map_xr = final_map_xr.to_dataset(name=variable)

    # Flip map in the y direction
    final_map_xr = final_map_xr.reindex(lat=final_map_xr.lat[::-1])

    # Regrid in space
    if final_spatial_resolution > 30:
        with warnings.catch_warnings():
                warnings.simplefilter("ignore") 
                regridder = spatial_regrid['regridder']
                final_map_xr = regridder(final_map_xr)
                final_map_xr = final_map_xr.reindex(lat=final_map_xr.lat[::-1])

    elif final_spatial_resolution < 30:
        final_map_xr = final_map_xr.interp(lat=spatial_regrid['ag_region_lats_center'], 
                                           lon=spatial_regrid['ag_region_lons_center'], method='linear')

    return final_map_xr

def define_regrid_information(final_spatial_resolution, r, region_lats, region_lons,
                              lat_range, lon_range, lat_delta, lon_delta):

    if final_spatial_resolution == 30:
        spatial_regrid = None

    if final_spatial_resolution != 30:

        # Create local range
        # Not pythonic elegant, but consistent with update_latlon_to_new_spatial_resolution
        # This implementation minimizes truncation error for small deltas

        ag_region_lats_center = [lat_range[0]]
        while ag_region_lats_center[-1] < lat_range[1]:
            ag_region_lats_center.append(ag_region_lats_center[-1]+lat_delta)
        ag_region_lats_center = np.array(ag_region_lats_center)

        ag_region_lons_center = [lon_range[0]]
        while ag_region_lons_center[-1] < lon_range[1]:
            ag_region_lons_center.append(ag_region_lons_center[-1]+lon_delta)
        ag_region_lons_center = np.array(ag_region_lons_center)

        spatial_regrid = {}
        if final_spatial_resolution < 30: 
 
            spatial_regrid['ag_region_lats_center'] = ag_region_lats_center
            spatial_regrid['ag_region_lons_center'] = ag_region_lons_center

        if final_spatial_resolution > 30:

            # Regrid with ESMF - grid settings
            grid_in  = {'lat':   region_lats,
                        'lon':   region_lons,    # center
                        'lat_b': get_bounds_from_centers(region_lats),
                        'lon_b': get_bounds_from_centers(region_lons)}  # bounds
            grid_out = {'lat':   ag_region_lats_center,
                        'lon':   ag_region_lons_center,
                        'lat_b': get_bounds_from_centers(ag_region_lats_center,lat_delta),
                        'lon_b': get_bounds_from_centers(ag_region_lons_center,lon_delta)}
            ds_out = xr.Dataset({
                    "lat": (["lat"], ag_region_lats_center),
                    "lon": (["lon"], ag_region_lons_center),})


            # Spatial regridding uses the xESMF package. However, it will only work if installed TOGETHER with esmpy:
            # source activate SMAPHB
            # conda install -c conda-forge esmpy xesmf
            # This will include a mpi installation into the working environment
            # https://github.com/JiaweiZhuang/xESMF/issues/102
            import xesmf as xe
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                regridder = xe.Regridder(grid_in, grid_out, 'conservative', reuse_weights=False)
                spatial_regrid['regridder'] = regridder

    return spatial_regrid

def get_bounds_from_centers(array_centers, delta=None):
    if delta == None: delta = array_centers[1]-array_centers[0]
    array_bounds = np.concatenate((array_centers, [array_centers[-1]+delta]), axis=0)
    array_bounds = array_bounds-delta/2.
    return array_bounds


def subset_zarr_into_netcdf_groups(final_path,final_temporal_resolution):
    ds = xr.open_zarr(final_path)
    #if final_temporal_resolution in ['6h']:
    #    index = ds.time.dt.strftime('%Y-%m-%d')
    if final_temporal_resolution in ['6h','daily','monthly']:
        index = ds.time.dt.strftime('%Y-%m')
    elif final_temporal_resolution in ['annual']:
        index = ds.time.dt.strftime('%Y')
    dates, sub_datasets = zip(*ds.groupby(index))
    paths = ['%s_%s.nc' % (final_path, d) for d in dates]
    return sub_datasets, paths

def write_groups_into_netcdf_files(sub_datasets, paths, compression_level, chunk):
    
    for dataset, path in zip(sub_datasets, paths):
        chunk['time'] = len(dataset.time)
        netcdf_encoding = {'SMAPHB_SM' : {'_FillValue': -9999,
                                          'complevel': compression_level, # Output data compression level:
                                                #  [0] No compression (fast)
                                                #  [9] max compression (slow)
                                         'zlib': True,
                                         #'chunksizes': (chunk['time'],chunk['lat'],chunk['lon']),
                                        },
                        'time': {'dtype': 'i4'},
                        }
        for var in dataset:
            del dataset[var].encoding['chunks']
            del dataset[var].encoding['preferred_chunks']
        dataset = dataset.chunk(chunk)
        dataset.to_netcdf(path, encoding=netcdf_encoding, format='NETCDF4', engine="netcdf4")

    return

def retrieve_data(database_path, 
                  variable,
                  output_folder, output_format,
                  final_temporal_resolution,
                  start_date, end_date,
                  minlat, maxlat, minlon, maxlon,
                  final_spatial_resolution,
                  compression_level,
                  mpi_run):

    if final_spatial_resolution != 30 and mpi_run == True:
        sys.exit('Data regridding to spatial resolution > 30 meters only works without MPI, please set mpi_run = False')
    if final_spatial_resolution > 30 and mpi_run == False:
        import xesmf as xe

    if mpi_run == True:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.Get_rank()
    else:
        rank=0
        size=1

    if rank == 0: print(datetime.datetime.now(), 'Start', flush=True)

    # Retrive time arrays
    times, ag_times = retrive_time_stepping(database_path, start_date, end_date, final_temporal_resolution)

    # Retrievel HRU mapping and lat lon arrays
    catch_map, hru_map = open_mosaic_object(database_path, minlat, maxlat, minlon, maxlon)
    lats = catch_map.y.values
    lons = catch_map.x.values

    # get lat lons and deltas
    lats, lons, lat_delta, lon_delta = update_latlon_to_new_spatial_resolution(final_spatial_resolution, lats, lons)

    # If data is output at corser spatial resolution, update the lat lon arrays
    if  final_spatial_resolution != 30 :
        catch_map, hru_map = open_mosaic_object(database_path, lats[0] -lat_delta/2.0, 
                                                               lats[-1]+lat_delta/2.0,
                                                               lons[0] -lon_delta/2.0,
                                                               lons[-1]+lon_delta/2.0)

    # Define final data chunk for zarr
    if final_temporal_resolution == '6h':      chunk = {"time":140, "lat": 500,  "lon": 500}
    if final_temporal_resolution == 'daily':   chunk = {"time":70,  "lat": 900,  "lon": 900}
    if final_temporal_resolution == 'monthly': chunk = {"time":30, "lat": 1400, "lon": 1400}
    if final_temporal_resolution == 'annual':  chunk = {"time":10, "lat": 2500, "lon": 2500}
    original_chunks = chunk.copy()

    # update chunk in lat lon to match final_spatial_resolution
    chunk = update_chunks_to_new_spatial_resolution(final_spatial_resolution,chunk)

    # create zarr template
    file_name = define_output_file_name_zarr(final_spatial_resolution, final_temporal_resolution)
    final_path = '%s/%s' % (output_folder, file_name)
    if rank == 0:
        regions = create_zarr_template(final_path, variable,
                                       lats, lons,
                                       final_spatial_resolution, chunk,
                                       ag_times, final_temporal_resolution, 
                                       compression_level)
    else:
        regions = None
    
    if mpi_run: 
        comm.Barrier()
        # retrive master information
        regions = comm.bcast(regions, root = 0)

    # Create a combination of sub-domains (regions) based on chunk sizes
    regions_range_elements = list(product(regions['lat_range'], regions['lon_range']))
    regions_slice_elements = list(product(regions['lat_slice'], regions['lon_slice']))
    nregions = len(regions_slice_elements)
    if rank == 0: print(datetime.datetime.now(), 'total sub-domains to work on:', nregions, flush=True) 

    # Loop over the regions
    for r in np.arange(nregions)[rank::size]:

        lat_range, lon_range = regions_range_elements[r]
        lat_slice, lon_slice = regions_slice_elements[r]
            
        # Retrieve mosaic data
        region_catch_map, region_hru_map, unique_catchments, region_lats, region_lons = retrieve_mosaic(
                                                                  catch_map, hru_map,
                                                                  lat_range, lon_range,
                                                                  lat_delta, lon_delta)

        if len(unique_catchments) == 0:
            #print(datetime.datetime.now(), 'domain:', r, 'empty', flush=True)
            continue

        # If data change changes spatial resolution use regridder
        spatial_regrid = define_regrid_information(final_spatial_resolution, r,
                                                   region_lats, region_lons,
                                                   lat_range, lon_range, lat_delta, lon_delta)

        # Loop through dates and save final file
        print(datetime.datetime.now(), 'domain:', r, 'computing...',flush=True)
        for time_slice, time_range in zip(regions['time_slice'],regions['time_range']):


                # grab data, interpolate in time, and space
                region_map = map_data(database_path, variable,
                                      region_lats, region_lons,
                                      spatial_regrid,
                                      time_range, final_temporal_resolution,
                                      region_catch_map, region_hru_map, unique_catchments,
                                      final_spatial_resolution)

                # send data to master zarr
                region_map.to_zarr(final_path, region={
                                               "lat": slice(lat_slice[0], lat_slice[1]+1),
                                               "lon": slice(lon_slice[0], lon_slice[1]+1),
                                               "time": slice(time_slice[0], time_slice[1]+1)},
                                  )
                region_map.close()
                del region_map
    
        if final_spatial_resolution > 30:
            os.system('rm -rf conservative_*.nc') 

    hru_map.close()
    del catch_map, hru_map
    if mpi_run: comm.Barrier()

    # Save netcdf output
    if output_format == 'netcdf' or output_format == 'both':
        sub_datasets, paths = subset_zarr_into_netcdf_groups(final_path,final_temporal_resolution)
        if rank == 0: print(datetime.datetime.now(),'Writting %i NetCDF files' % len(paths),flush=True)
        write_groups_into_netcdf_files(sub_datasets[rank::size], paths[rank::size], 
                                       compression_level, original_chunks)
        del sub_datasets, paths
        if mpi_run: comm.Barrier()
        if rank == 0:
            os.system('mkdir %s_netcdf' % final_path)
            os.system('mv %s/*.nc %s_netcdf/' % (output_folder,final_path))
            if output_format == 'netcdf': 
                os.system('rm -rf %s' % final_path)

    # Save zarr output
    if output_format == 'zarr' or output_format == 'both':
        if rank == 0:
            final_file = '%s.zarr' % final_path
            os.system('rm -rf %s' % (final_file))
            os.system('mv %s %s' % (final_path, final_file))

    if mpi_run: comm.Barrier()
    if rank == 0: print(datetime.datetime.now(),'Completed',flush=True)

    return
