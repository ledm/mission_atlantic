import matplotlib as mpl
mpl.use('Agg')

from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot
from scipy import interpolate

import os


# interpolate to this grid:





def regrid(arr, lats, lons, new_lats, new_lons, method='linear',fill_value=1.E20):
    """
     Regrid into a higher res image.

    arr is input data
    lons is input longitude
    lats in input longitube
    It's expecting a
    """
    if lons.ndim == 2 and lats.ndim == 2:
        in_points = np.array([[la, lo] for la, lo in zip(lats.flatten(), lons.flatten())])
    elif lons.ndim == 1 and lats.ndim == 1:
        #mg_lats, mg_longs = np.meshgrid(lats, lons)
        mg_longs, mg_lats = np.meshgrid(lons, lats)

        in_points = np.array([[la, lo] for la, lo in zip(mg_lats.flatten(), mg_longs.flatten())])
    else: assert 0

    if new_lons.ndim == 2 and new_lats.ndim == 2:
        out_points = np.array([[la, lo] for la, lo in zip(new_lats.flatten(), new_lons.flatten())])
        mg_nlats = new_lats
    elif new_lons.ndim == 1 and new_lats.ndim == 1:
#        mg_nlats, mg_nlongs = np.meshgrid(new_lats, new_lons)
        mg_nlongs, mg_nlats = np.meshgrid(new_lons, new_lats)

        out_points = np.array([[la, lo] for la, lo in zip(mg_nlats.flatten(), mg_nlongs.flatten())])

    else: assert 0
    values = np.array(arr.data.flatten())
    print('regrid:', len(in_points), len(out_points), len(values), (arr.min(), arr.max()))

    new_arr = interpolate.griddata(in_points, values, out_points, method=method, )
    new_arr = new_arr.reshape(mg_nlats.shape)
    new_arr= np.ma.masked_where(new_arr==fill_value, new_arr)
    print('regridded:',(new_arr.min(), new_arr.max()))
    print(arr.shape, new_arr.shape, mg_lats.shape, mg_nlats.shape)


    return new_arr




def run(fni, fng, fno, datasetFormat='NETCDF4'):
    """
    fni: input data
    fng: input grid
    fno: output path
    """
    debug = True
    Twelve_months_only=False

    if debug:
        print('Opening', fni)

    if os.path.exists(fno):
        print('THIS file already exists:', fno)
        return
    nci = Dataset(fni,'r')
    ncg = Dataset(fng,'r')

    # create dataset and netcdf attributes.
    nco = Dataset(fno,'w',format=datasetFormat)
    attributes = nci.ncattrs()
    attributes = {a:nci.getncattr(a) for a in attributes}
    attributes['Regrid'] = 'Regridded from original grid.'
    print(attributes)
    for att,attribute in attributes.items():
        nco.setncattr(att, attribute)
        if debug: print ('Adding attribute: ',att,'\t(',attribute,')')

    # copy dimensions from input grid, fng/ncg
    for dim, dimension in ncg.dimensions.items():
        dimSize=len(ncg.dimensions[dim])
        if dim in ['time', 'TIME']: 
            dimSize = None
        if debug: print ('Adding dimension: ',dim,'\t(',dimSize,')')
        nco.createDimension(dim, dimSize)

    # list of variables to save
    keys = ncg.variables.keys()

    removeVar = 'pCO2a' # from ncg
    if fni.find('ndep-nhx')>-1: 
        newVars = ['nhx', ] # from nci
    if fni.find('ndep-noy') > -2:
        newVars = ['noy', ] # from nci



    # create Variables:
    for var in keys:
        if var in newVars: continue
        if var == removeVar: continue
        newname = var
        dimensions = list(ncg.variables[var].dimensions)
        vartype = ncg.variables[var].dtype
        nco.createVariable(newname, vartype, tuple(dimensions), zlib=True,complevel=5)
        if debug: print('adding variable: ',var,'-->', newname, '\t(',dimensions,')')

    for var in newVars:
        dimensions = list(ncg.variables[removeVar].dimensions)
        vartype = ncg.variables[removeVar].dtype
        nco.createVariable(var, vartype, tuple(dimensions),zlib=True,complevel=5)
        if debug: print('adding new variable: ',var, '\t(',dimensions,')')

    # long names, standard_name, units, fill values and data for dimensions
    for var in keys:
        if var in newVars: continue
        if var == removeVar: continue

        try: nco.variables[var].units = ncg.variables[var].units
        except: pass
        try: nco.variables[var].standard_name = ncg.variables[var].standard_name
        except: pass
        try: nco.variables[var].long_name = ncg.variables[var].long_name
        except: pass
        try: nco.variables[var].axis = ncg.variables[var].axis
        except: pass
        try: nco.variables[var].calendar = ncg.variables[var].calendar
        except: pass
        try: nco.variables[var].bounds = ncg.variables[var].bounds
        except: pass
        try: nco.variables[var].missing_value = ncg.variables[var].missing_value
        except: pass
        try: nco.variables[var]._FillValue = ncg.variables[var]._FillValue
        except: pass
        arr = ncg.variables[var][:]
        if Twelve_months_only and var == 'time': arr=arr[:12]

        nco.variables[var][:] = arr

    for var in newVars:
        try: nco.variables[var].units = nci.variables[var].units
        except: pass
        try: nco.variables[var].standard_name = nci.variables[var].standard_name
        except: pass
        try: nco.variables[var].long_name = nci.variables[var].long_name
        except: pass
        try: nco.variables[var].axis = nci.variables[var].axis
        except: pass
        try: nco.variables[var].calendar = nci.variables[var].calendar
        except: pass
        try: nco.variables[var].bounds = nci.variables[var].bounds
        except: pass
        try: nco.variables[var].missing_value = nci.variables[var].missing_value
        except: pass
        try: nco.variables[var]._FillValue = nci.variables[var]._FillValue
        except: pass

        old_time = nci.variables['time'][:]
        old_lats = nci.variables['lat'][:]
        old_lons = nci.variables['lon'][:]

        new_lats = nco.variables['LAT'][:]
        new_lons = nco.variables['LON'][:]

        if Twelve_months_only:
            out_arr = np.ma.empty((12, len(new_lats), len(new_lons)))
        else:
            out_arr = np.ma.empty((len(old_time), len(new_lats), len(new_lons)))

        print('making out array like:', out_arr.shape)
        for t in np.arange(old_time.size)[:]:
            if Twelve_months_only and t >=12: continue
            if t%6==0: 
                print('iterating', t, 'of', len(old_time))
            out_arr[t] = regrid(nci.variables[var][t,:,:], old_lats, old_lons, new_lats,new_lons, method='linear',)

        nco.variables[var][:] = np.ma.array(out_arr)
        print('successfully populated:\t', var)


    # Close netcdfs:
    nco.close()
    nci.close()
    ncg.close()
    print('successfully created:\t', fno)




outgrid = '/users/modellers/gig/Documents/MissionAtantic/eORCA1_Stuff/pCO2a_y1958.nc'
# outgrid2 = '/users/modellers/gig/Documents/MissionAtantic/eORCA1_Stuff/gelbstoff_absorption_satellite_y1958.nc'

in_data_fn  = '/data/sthenno1/scratch/yuti/MA_MissionAtlantic/N_deposition/ndep-nhx_histsoc_monthly_1901_2018.nc'
in_data_fn2 = '/data/sthenno1/scratch/yuti/MA_MissionAtlantic/N_deposition/ndep-noy_histsoc_monthly_1901_2018.nc'

outpath   = '/data/sthenno1/scratch/ledm/missionatlantic/n_dep/ndep-nhx_histsoc_monthly_1901_2018_atmosgrid_v3.nc'
outpath_2 = '/data/sthenno1/scratch/ledm/missionatlantic/n_dep/ndep-noy_histsoc_monthly_1901_2018_atmosgrid_v1.nc'


#run(in_data_fn, outgrid,outpath) # nhx
run(in_data_fn2, outgrid,outpath_2)

