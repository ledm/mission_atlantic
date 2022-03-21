
# In this script we look at the old restart files
# Compare them against the new ones
# Copy over a bunch of data while adding new pixels to the grid.
# save the file - while keeping track of data versions.
# There may be some regridding as well.

import matplotlib as mpl
mpl.use('Agg')

from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot
import os


imarnet_forcing_path = "/data/sthenno1/to_archive/ledm/Projects/iMarNet/OUTPUT_v2/iMarNet_forcing/"
iMarNet_inputs = "/data/sthenno1/to_archive/ledm/Projects/iMarNet/inputNetCDFs/"

imarnet_restart = 'input/iMARNET_fields_1890_UEA_ERSEM_newFeR3c_20140919.nc'
eorca_restart = 'input/eORCA1_00000992_restart.nc'
ersem_restart = 'input/sthenno1/restart_trc.HC.1986-1990.12.nc'
grid_nc = 'input/eORCA1_mesh_mask.nc'

def regrid(arr, lons, lats, new_lons, new_lats, method='linear',):
    """
     Regrid into a higher res image.

    arr is input data
    lons is input longitude
    lats in input longitube

    """
    in_points = np.array([[la, lo] for la, lo in zip(lats.flatten(), lons.flatten())])
    out_points = np.array([[la, lo] for la, lo in zip(new_lats.flatten(), new_lons.flatten())])
    values = np.array(arr.data.flatten())

    new_arr = interpolate.griddata(in_points, values, out_points, method=method, fill_value=1.E20)
    new_arr = new_arr.reshape(new_lats.shape)
    return new_arr

def print_nc(nc, key):
    print('\n------------------\nNetcdf:', key)

    #print('keys', sorted(nc.variables.keys()))
    for dim in nc.dimensions.values():
            print(dim, dim.name)

    for key in  sorted(nc.variables.keys()):
        continue
        v = nc.variables[key]
        #print(v.long_name)
        #ln = v.long_name
        #units = v.units
        try: ln = v.long_name
        except: ln='long_name'

        try:units = v.units
        except: units = 'units'

        print(key,':\t', ln, '\t',units, '\t',v.shape)

master_dict = {
        'TRBP1_Chl' : 'TRBChl1',
        'TRBP2_Chl' : 'TRBChl2',
        'TRBP3_Chl' : 'TRBChl3',
        'TRBP4_Chl' : 'TRBChl4',
        'TRNP1_Chl' : 'TRNChl1',
        'TRNP2_Chl' : 'TRNChl2',
        'TRNP3_Chl' : 'TRNChl3',
        'TRNP4_Chl' : 'TRNChl4',

        'TRBO3_bioalk' : 'TRBbioalk',
        'TRNO3_bioalk' : 'TRNbioalk',
        'TRBO3_TA': 'TRBbioalk',
        'TRNO3_TA': 'TRNbioalk',

        'TRNlight_ADY' : 'ZERO', # Need to write to James
        'TRBlight_ADY' : 'ZERO',
        'fabm_st2DbG2_o_deep' : 'ZERO', # ask Gennadi, put zero for now
        'fabm_st2DnG2_o_deep' : 'ZERO',
        'fabm_st2DbbL2_c' : 'TR2Dbl2c',
        'fabm_st2DnbL2_c' : 'TR2Dbl2c',
        'fabm_st2Dbben_nit_G4n' : 'TR2DG4n',# ask Gennadi, put zero for now
        'fabm_st2Dnben_nit_G4n' : 'TR2DG4n',
        'rdttrc1' : 'ZERO',

        # Iron:
        'fabm_st2DnQ6_pen_depth_f': 'TR2DQ6f',
        'fabm_st2DbQ6_pen_depth_f': 'TR2DQ6f',
        'TRBN7_f': 'TRBN7f',
        'TRBP1_f': 'TRBP1f',
        'TRBP2_f': 'TRBP2f',
        'TRBP3_f': 'TRBP3f',
        'TRBP4_f': 'TRBP4f',
        'TRBR4_f': 'TRBR4f',
        'TRBR6_f': 'TRBR6f',
        'TRNN7_f': 'TRNN7f',
        'TRNP1_f': 'TRNP1f',
        'TRNP2_f': 'TRNP2f',
        'TRNP3_f': 'TRNP3f',
        'TRNP4_f': 'TRNP4f',
        'TRNR4_f': 'TRNR4f',
        'TRNR6_f': 'TRNR6f',

        'fabm_st2DnQ6_pen_depth_c' : 'TR2DD6m',
        'fabm_st2DnQ6_pen_depth_n' : 'TR2DD7m',
        'fabm_st2DnQ6_pen_depth_p' : 'TR2DD8m',
        'fabm_st2DnQ6_pen_depth_s' : 'TR2DD9m',
        'fabm_st2DnQ7_pen_depth_c' : 'TR2DD3m',
        'fabm_st2DnQ7_pen_depth_n' : 'TR2DD4m',
        'fabm_st2DnQ7_pen_depth_p' : 'TR2DD5m',
        'fabm_st2Dnben_col_D1m' : 'TR2DD1m',
        'fabm_st2Dnben_col_D2m' : 'TR2DD2m',

        'fabm_st2DbQ6_pen_depth_c': 'TR2DD6m',
        'fabm_st2DbQ6_pen_depth_n': 'TR2DD7m',
        'fabm_st2DbQ6_pen_depth_p': 'TR2DD8m',
        'fabm_st2DbQ6_pen_depth_s': 'TR2DD9m',
        'fabm_st2DbQ7_pen_depth_c': 'TR2DD3m',
        'fabm_st2DbQ7_pen_depth_n': 'TR2DD4m',
        'fabm_st2DbQ7_pen_depth_p': 'TR2DD5m',
        'fabm_st2Dbben_col_D1m': 'TR2DD1m',
        'fabm_st2Dbben_col_D2m': 'TR2DD2m',
        }


def link_old_to_new_name(new_name, nc1, nc2):
    """
    ie convert TRNB1n to TRNB1_n
    new_name is key from recent ERSEM netcdf restart
    nc1 is recent ERSEM restart
    nc2 is old iMarNEt ERSEM restart.
    """

    # Bespoke list:
    name = master_dict.get(new_name, "")
    if name != "":
        return name

    # Identical bewteen lists:
    nc2vars = list(nc2.variables.keys())
    if new_name in nc2vars:
        return new_name

    # Yuri said to set these to zero:
    if new_name[:3] in ['sbc', 'atf', 'rdb']:
        return 'ZERO'

    # Create list of iMarNet keys without preficeSx.
#    nc2_real_names = {}
#    for n in nc2vars:
#        nn = n.replace('TRN','').replace('TRB', '')
#        try: nc2_real_names[nn].append(n)
#        except: nc2_real_names[nn] = [n,]

    # Replace some bits in modern ERSEM:
    replacements = {'_': '',
        'fabm': '',
        'st2Db':'',
        'st2Dn': '',
            }
    appends = ['', 'TRN', 'TRB', 'TR2D',]

    new_name_copy = new_name[:]
    for r,a in replacements.items():
        new_name_copy = new_name_copy.replace(r,a)

        for app in appends:
            print(new_name, 'trying:', app+new_name_copy)
            if app+new_name_copy in nc2vars:
                return app+new_name_copy

    print('not found:', new_name, '->', new_name_copy)
    return False

def run():
    orcnc = Dataset(eorca_restart, 'r')
    print_nc(orcnc, 'eORCA')
    """Netcdf: eORCA
    keys ['a_fwb', 'a_fwb_b', 'adatrj', 'avm_k', 'avt_k', 'dissl', 'e3t_b', 'e3t_n', 'emp_b', 'en', 'fraqsr_1lev', 'fwf_isf_b', 'isf_hc_b', 'isf_sc_b', 'kt', 'nav_lat', 'nav_lev', 'nav_lon', 'ndastp', 'ntime', 'q
    ns_b', 'qsr_hc_b', 'rdt', 'rhop', 'rnf_b', 'rnf_hc_b', 'rnf_sc_b', 'sb', 'sbc_hc_b', 'sbc_sc_b', 'sfx_b', 'sn', 'sshb', 'sshn', 'tb', 'time_counter', 'tn', 'ub', 'ub2_b', 'un', 'un_bf', 'utau_b', 'vb', 'vb2_b
    ', 'vn', 'vn_bf', 'vtau_b']
    <class 'netCDF4._netCDF4.Dimension'>: name = 'x', size = 362
    <class 'netCDF4._netCDF4.Dimension'>: name = 'y', size = 332
    <class 'netCDF4._netCDF4.Dimension'>: name = 'nav_lev', size = 75
    <class 'netCDF4._netCDF4.Dimension'> (unlimited): name = 'time_counter', size = 1"""
    #temp = orcnc.variables[

    ersnc = Dataset(ersem_restart, 'r')
    print_nc(ersnc, 'ERSEM')
    """Netcdf: ERSEM
    keys ['TRBB1_c', 'TRBB1_n', 'TRBB1_p', 'TRBL2_c', 'TRBN1_p', 'TRBN3_n', 'TRBN4_n', 'TRBN5_s', 'TRBO2_o', 'TRBO3_bioalk', 'TRBO3_c', 'TRBP1_Chl', 'TRBP1_c', 'TRBP1_n', 'TRBP1_p', 'TRBP1_s', 'TRBP2_Chl', 'TRBP2
    _c', 'TRBP2_n', 'TRBP2_p', 'TRBP3_Chl', 'TRBP3_c', 'TRBP3_n', 'TRBP3_p', 'TRBP4_Chl', 'TRBP4_c', 'TRBP4_n', 'TRBP4_p', 'TRBR1_c', 'TRBR1_n', 'TRBR1_p', 'TRBR2_c', 'TRBR3_c', 'TRBR4_c', 'TRBR4_n', 'TRBR4_p', '
    TRBR6_c', 'TRBR6_n', 'TRBR6_p', 'TRBR6_s', 'TRBR8_c', 'TRBR8_n', 'TRBR8_p', 'TRBR8_s', 'TRBZ4_c', 'TRBZ5_c', 'TRBZ5_n', 'TRBZ5_p', 'TRBZ6_c', 'TRBZ6_n', 'TRBZ6_p', 'TRBlight_ADY', 'TRNB1_c', 'TRNB1_n', 'TRNB1
    _p', 'TRNL2_c', 'TRNN1_p', 'TRNN3_n', 'TRNN4_n', 'TRNN5_s', 'TRNO2_o', 'TRNO3_bioalk', 'TRNO3_c', 'TRNP1_Chl', 'TRNP1_c', 'TRNP1_n', 'TRNP1_p', 'TRNP1_s', 'TRNP2_Chl', 'TRNP2_c', 'TRNP2_n', 'TRNP2_p', 'TRNP3_
    Chl', 'TRNP3_c', 'TRNP3_n', 'TRNP3_p', 'TRNP4_Chl', 'TRNP4_c', 'TRNP4_n', 'TRNP4_p', 'TRNR1_c', 'TRNR1_n', 'TRNR1_p', 'TRNR2_c', 'TRNR3_c', 'TRNR4_c', 'TRNR4_n', 'TRNR4_p', 'TRNR6_c', 'TRNR6_n', 'TRNR6_p', 'T
    RNR6_s', 'TRNR8_c', 'TRNR8_n', 'TRNR8_p', 'TRNR8_s', 'TRNZ4_c', 'TRNZ5_c', 'TRNZ5_n', 'TRNZ5_p', 'TRNZ6_c', 'TRNZ6_n', 'TRNZ6_p', 'TRNlight_ADY', 'adatrj', 'atf_trend_B1_c', 'atf_trend_B1_n', 'atf_trend_B1_p'
    , 'atf_trend_L2_c', 'atf_trend_N1_p', 'atf_trend_N3_n', 'atf_trend_N4_n', 'atf_trend_N5_s', 'atf_trend_O2_o', 'atf_trend_O3_bioalk', 'atf_trend_O3_c', 'atf_trend_P1_Chl', 'atf_trend_P1_c', 'atf_trend_P1_n', '
    atf_trend_P1_p', 'atf_trend_P1_s', 'atf_trend_P2_Chl', 'atf_trend_P2_c', 'atf_trend_P2_n', 'atf_trend_P2_p', 'atf_trend_P3_Chl', 'atf_trend_P3_c', 'atf_trend_P3_n', 'atf_trend_P3_p', 'atf_trend_P4_Chl', 'atf_
    trend_P4_c', 'atf_trend_P4_n', 'atf_trend_P4_p', 'atf_trend_R1_c', 'atf_trend_R1_n', 'atf_trend_R1_p', 'atf_trend_R2_c', 'atf_trend_R3_c', 'atf_trend_R4_c', 'atf_trend_R4_n', 'atf_trend_R4_p', 'atf_trend_R6_c
    ', 'atf_trend_R6_n', 'atf_trend_R6_p', 'atf_trend_R6_s', 'atf_trend_R8_c', 'atf_trend_R8_n', 'atf_trend_R8_p', 'atf_trend_R8_s', 'atf_trend_Z4_c', 'atf_trend_Z5_c', 'atf_trend_Z5_n', 'atf_trend_Z5_p', 'atf_tr
    end_Z6_c', 'atf_trend_Z6_n', 'atf_trend_Z6_p', 'atf_trend_light_ADY', 'fabm_st2DbG2_o', 'fabm_st2DbG2_o_deep', 'fabm_st2DbG3_c', 'fabm_st2DbH1_c', 'fabm_st2DbH2_c', 'fabm_st2DbK1_p', 'fabm_st2DbK3_n', 'fabm_s
    t2DbK4_n', 'fabm_st2DbK5_s', 'fabm_st2DbQ17_c', 'fabm_st2DbQ17_n', 'fabm_st2DbQ17_p', 'fabm_st2DbQ1_c', 'fabm_st2DbQ1_n', 'fabm_st2DbQ1_p', 'fabm_st2DbQ6_c', 'fabm_st2DbQ6_n', 'fabm_st2DbQ6_p', 'fabm_st2DbQ6_
    pen_depth_c', 'fabm_st2DbQ6_pen_depth_n', 'fabm_st2DbQ6_pen_depth_p', 'fabm_st2DbQ6_pen_depth_s', 'fabm_st2DbQ6_s', 'fabm_st2DbQ7_c', 'fabm_st2DbQ7_n', 'fabm_st2DbQ7_p', 'fabm_st2DbQ7_pen_depth_c', 'fabm_st2D
    bQ7_pen_depth_n', 'fabm_st2DbQ7_pen_depth_p', 'fabm_st2DbY2_c', 'fabm_st2DbY3_c', 'fabm_st2DbY4_c', 'fabm_st2DbbL2_c', 'fabm_st2Dbben_col_D1m', 'fabm_st2Dbben_col_D2m', 'fabm_st2Dbben_nit_G4n', 'fabm_st2DnG2_
    o', 'fabm_st2DnG2_o_deep', 'fabm_st2DnG3_c', 'fabm_st2DnH1_c', 'fabm_st2DnH2_c', 'fabm_st2DnK1_p', 'fabm_st2DnK3_n', 'fabm_st2DnK4_n', 'fabm_st2DnK5_s', 'fabm_st2DnQ17_c', 'fabm_st2DnQ17_n', 'fabm_st2DnQ17_p'
    , 'fabm_st2DnQ1_c', 'fabm_st2DnQ1_n', 'fabm_st2DnQ1_p', 'fabm_st2DnQ6_c', 'fabm_st2DnQ6_n', 'fabm_st2DnQ6_p', 'fabm_st2DnQ6_pen_depth_c', 'fabm_st2DnQ6_pen_depth_n', 'fabm_st2DnQ6_pen_depth_p', 'fabm_st2DnQ6_
    pen_depth_s', 'fabm_st2DnQ6_s', 'fabm_st2DnQ7_c', 'fabm_st2DnQ7_n', 'fabm_st2DnQ7_p', 'fabm_st2DnQ7_pen_depth_c', 'fabm_st2DnQ7_pen_depth_n', 'fabm_st2DnQ7_pen_depth_p', 'fabm_st2DnY2_c', 'fabm_st2DnY3_c', 'f
    abm_st2DnY4_c', 'fabm_st2DnbL2_c', 'fabm_st2Dnben_col_D1m', 'fabm_st2Dnben_col_D2m', 'fabm_st2Dnben_nit_G4n', 'kt', 'nav_lat', 'nav_lev', 'nav_lon', 'ndastp', 'rdb_trend_B1_c', 'rdb_trend_B1_n', 'rdb_trend_B1
    _p', 'rdb_trend_L2_c', 'rdb_trend_N1_p', 'rdb_trend_N3_n', 'rdb_trend_N4_n', 'rdb_trend_N5_s', 'rdb_trend_O3_c', 'rdb_trend_P1_Chl', 'rdb_trend_P1_c', 'rdb_trend_P1_n', 'rdb_trend_P1_p', 'rdb_trend_P1_s', 'rd
    b_trend_P2_Chl', 'rdb_trend_P2_c', 'rdb_trend_P2_n', 'rdb_trend_P2_p', 'rdb_trend_P3_Chl', 'rdb_trend_P3_c', 'rdb_trend_P3_n', 'rdb_trend_P3_p', 'rdb_trend_P4_Chl', 'rdb_trend_P4_c', 'rdb_trend_P4_n', 'rdb_tr
    end_P4_p', 'rdb_trend_R1_c', 'rdb_trend_R1_n', 'rdb_trend_R1_p', 'rdb_trend_R2_c', 'rdb_trend_R3_c', 'rdb_trend_R4_c', 'rdb_trend_R4_n', 'rdb_trend_R4_p', 'rdb_trend_R6_c', 'rdb_trend_R6_n', 'rdb_trend_R6_p',
     'rdb_trend_R6_s', 'rdb_trend_R8_c', 'rdb_trend_R8_n', 'rdb_trend_R8_p', 'rdb_trend_R8_s', 'rdb_trend_Z4_c', 'rdb_trend_Z5_c', 'rdb_trend_Z5_n', 'rdb_trend_Z5_p', 'rdb_trend_Z6_c', 'rdb_trend_Z6_n', 'rdb_tren
    d_Z6_p', 'rdb_trend_light_ADY', 'rdttrc1', 'sbc_B1_c_b', 'sbc_B1_n_b', 'sbc_B1_p_b', 'sbc_L2_c_b', 'sbc_N1_p_b', 'sbc_N3_n_b', 'sbc_N4_n_b', 'sbc_N5_s_b', 'sbc_O2_o_b', 'sbc_O3_bioalk_b', 'sbc_O3_c_b', 'sbc_P
    1_Chl_b', 'sbc_P1_c_b', 'sbc_P1_n_b', 'sbc_P1_p_b', 'sbc_P1_s_b', 'sbc_P2_Chl_b', 'sbc_P2_c_b', 'sbc_P2_n_b', 'sbc_P2_p_b', 'sbc_P3_Chl_b', 'sbc_P3_c_b', 'sbc_P3_n_b', 'sbc_P3_p_b', 'sbc_P4_Chl_b', 'sbc_P4_c_
    b', 'sbc_P4_n_b', 'sbc_P4_p_b', 'sbc_R1_c_b', 'sbc_R1_n_b', 'sbc_R1_p_b', 'sbc_R2_c_b', 'sbc_R3_c_b', 'sbc_R4_c_b', 'sbc_R4_n_b', 'sbc_R4_p_b', 'sbc_R6_c_b', 'sbc_R6_n_b', 'sbc_R6_p_b', 'sbc_R6_s_b', 'sbc_R8_
    c_b', 'sbc_R8_n_b', 'sbc_R8_p_b', 'sbc_R8_s_b', 'sbc_Z4_c_b', 'sbc_Z5_c_b', 'sbc_Z5_n_b', 'sbc_Z5_p_b', 'sbc_Z6_c_b', 'sbc_Z6_n_b', 'sbc_Z6_p_b', 'sbc_light_ADY_b', 'time_counter']
    <class 'netCDF4._netCDF4.Dimension'>: name = 'y', size = 375
    <class 'netCDF4._netCDF4.Dimension'>: name = 'x', size = 297
    <class 'netCDF4._netCDF4.Dimension'>: name = 'z', size = 51
    <class 'netCDF4._netCDF4.Dimension'> (unlimited): name = 't', size = 1"""


    imnnc = Dataset(imarnet_restart, 'r')
    print_nc(imnnc, 'iMarNet')

    """Netcdf: iMarNet
    keys ['TR2DD1m', 'TR2DD2m', 'TR2DD3m', 'TR2DD4m', 'TR2DD5m', 'TR2DD6m', 'TR2DD7m', 'TR2DD8m', 'TR2DD9m', 'TR2DG2o', 'TR2DG3c', 'TR2DG4n', 'TR2DH1c', 'TR2DH2c', 'TR2DK1p', 'TR2DK3n', 'TR2DK4n', 'TR2DK5s', 'TR2
    DQ17c', 'TR2DQ17n', 'TR2DQ17p', 'TR2DQ1c', 'TR2DQ1n', 'TR2DQ1p', 'TR2DQ6c', 'TR2DQ6f', 'TR2DQ6n', 'TR2DQ6p', 'TR2DQ6s', 'TR2DQ7c', 'TR2DQ7n', 'TR2DQ7p', 'TR2DY2c', 'TR2DY3c', 'TR2DY4c', 'TR2Dbl2c', 'TRBB1c',
    'TRBB1n', 'TRBB1p', 'TRBChl1', 'TRBChl2', 'TRBChl3', 'TRBChl4', 'TRBL2c', 'TRBN1p', 'TRBN3n', 'TRBN4n', 'TRBN5s', 'TRBN7f', 'TRBO2o', 'TRBO3c', 'TRBP1c', 'TRBP1f', 'TRBP1n', 'TRBP1p', 'TRBP1s', 'TRBP2c', 'TRB
    P2f', 'TRBP2n', 'TRBP2p', 'TRBP3c', 'TRBP3f', 'TRBP3n', 'TRBP3p', 'TRBP4c', 'TRBP4f', 'TRBP4n', 'TRBP4p', 'TRBR1c', 'TRBR1n', 'TRBR1p', 'TRBR2c', 'TRBR3c', 'TRBR4c', 'TRBR4f', 'TRBR4n', 'TRBR4p', 'TRBR6c', 'T
    RBR6f', 'TRBR6n', 'TRBR6p', 'TRBR6s', 'TRBR8c', 'TRBR8n', 'TRBR8p', 'TRBR8s', 'TRBZ4c', 'TRBZ5c', 'TRBZ5n', 'TRBZ5p', 'TRBZ6c', 'TRBZ6n', 'TRBZ6p', 'TRBbioalk', 'TRNB1c', 'TRNB1n', 'TRNB1p', 'TRNChl1', 'TRNCh
    l2', 'TRNChl3', 'TRNChl4', 'TRNL2c', 'TRNN1p', 'TRNN3n', 'TRNN4n', 'TRNN5s', 'TRNN7f', 'TRNO2o', 'TRNO3c', 'TRNP1c', 'TRNP1f', 'TRNP1n', 'TRNP1p', 'TRNP1s', 'TRNP2c', 'TRNP2f', 'TRNP2n', 'TRNP2p', 'TRNP3c', '
    TRNP3f', 'TRNP3n', 'TRNP3p', 'TRNP4c', 'TRNP4f', 'TRNP4n', 'TRNP4p', 'TRNR1c', 'TRNR1n', 'TRNR1p', 'TRNR2c', 'TRNR3c', 'TRNR4c', 'TRNR4f', 'TRNR4n', 'TRNR4p', 'TRNR6c', 'TRNR6f', 'TRNR6n', 'TRNR6p', 'TRNR6s',
     'TRNR8c', 'TRNR8n', 'TRNR8p', 'TRNR8s', 'TRNZ4c', 'TRNZ5c', 'TRNZ5n', 'TRNZ5p', 'TRNZ6c', 'TRNZ6n', 'TRNZ6p', 'TRNbioalk', 'adatrj', 'arak0', 'kt', 'nav_lat', 'nav_lev', 'nav_lon', 'ndastp', 'time_counter']

    <class 'netCDF4._netCDF4.Dimension'>: name = 'x', size = 362
    <class 'netCDF4._netCDF4.Dimension'>: name = 'y', size = 292
    <class 'netCDF4._netCDF4.Dimension'>: name = 'z', size = 75
    <class 'netCDF4._netCDF4.Dimension'> (unlimited): name = 't', size = 1 """
    outnc = Dataset(new_restart, 'w')

    # tasks:
    # load data
    # for each variable in the ersem restart file
    #    find the relevant imarnet field
    #    load the iMarNEt data
    #    expand the iMarNet field to the same size as the new eORCA1 grid.
    #
    # save output

def plot_map(arr, key,fold='images/'):
    fn = fold+key+'.png'
    fn = fn.replace(' ', '_')
    if os.path.exists(fn):
        return
    pyplot.pcolormesh(arr)
    pyplot.colorbar()
    pyplot.title(key)
    print('Saving figure:', fn)
    try:os.mkdir(os.path.dirname(fn))
    except: pass
    pyplot.savefig(fn)
    pyplot.close()


def extend_ORCA(arr):
    """Converts old ORCA array and returns the same data in the new eORCA array"""
    if arr.ndim == 2:
        b = np.ma.ones((332, 362)) * 1E20
        b[-292:,:] = arr
    elif arr.ndim == 3:
        b = np.ma.ones((arr.shape[0], 332, 362)) * 1E20
        b[:,-292:,:] = arr

    elif arr.ndim == 4:
        b = np.ma.ones((arr.shape[0],arr.shape[1], 332, 362)) * 1E20
        b[:,:,-292:,:] = arr
    else: assert 0
    b = np.ma.masked_where(b.mask + (b==1E20), b)
    return b


#run()
def test_grid_convert():
    #nc = Dataset(imarnet_restart, 'r') # eORCA1 #
    nc = Dataset(eorca_restart, 'r')
    lat = nc.variables['nav_lat'][:]
    lon = nc.variables['nav_lon'][:]
    plot_map(lat, 'eORCA lat')
    plot_map(lon, 'eORCA lon')
    print_nc(nc, 'eORCA1')
    sb = nc.variables['sb'][0,0]
    mask = np.ma.masked_where(sb==0. + sb.mask,sb).mask
    plot_map( nc.variables['sb'][0,0], 'sb')
    plot_map( nc.variables['sb'][0,0].mask, 'sb mask')

    nc2 = Dataset(imarnet_restart, 'r') # ORCA1 grid
    lati = nc2.variables['nav_lat'][:]
    loni = nc2.variables['nav_lon'][:]

    plot_map(lati, 'ORCA lat')
    plot_map(loni, 'ORCA lon')

    late = extend_ORCA(lati)
    lone = extend_ORCA(loni)

    plot_map(late, 'extended ORCA lat')
    plot_map(lone, 'extended ORCA lon')

    lat_diff = late - lat
    lon_diff = lone - lon

    lat_diff = np.ma.masked_where(mask, lat_diff)
    lon_diff =  np.ma.masked_where(mask, lon_diff)

    plot_map(lat_diff, 'Difference lat')
    plot_map(lon_diff, 'Difference lon')
    plot_map(np.ma.clip(lon_diff,-2  ,2  ), 'Difference lon clipped')

# test_grid_convert()
# Test the link between new data nad old data.

#    nc = Dataset(eorca_restart, 'r')

def test_linking():
    """
    Test the link bewteen the new ERSEM restart and the iMArNEt one.,
    """
    nc = Dataset(ersem_restart, 'r')
    nc2 = Dataset(imarnet_restart, 'r') # ORCA1 grid
    successes = {}
    failures = []
    missings = []
    iMArNEt_successes = {key:[] for key in nc2.variables.keys()}

    for key in sorted(nc.variables.keys()):
        oldkey = link_old_to_new_name(key, nc, nc2)
        if oldkey:
            print('success', key, '-->', oldkey)
            successes[key] = oldkey
            try:iMArNEt_successes[oldkey].append(key)
            except:iMArNEt_successes[oldkey] = [key, ]
        else:
            failures.append(key)

    # Look at iMarNEt file for missing variables.
    for ikey in sorted(nc2.variables.keys()):
        if len(iMArNEt_successes[ikey]): continue
        missings.append(ikey)

    for i, missing in enumerate(missings):
        if nc2.variables[missing][:].min() == nc2.variables[missing][:].max():
            print('missing:', i,missing, 'flat:', nc2.variables[missing][:].min())
        else:
            print('missing:', i,missing,nc2.variables[missing][:].min(), nc2.variables[missing][:].max())

    for i, missing in enumerate(missings):
        print('        "'+missing+'",',nc2.variables[missing].dtype, tuple([dim.name for dim in nc2.variables[missing].get_dims()]))

    for i,fail in enumerate(failures):
        print('not found:', i, fail)
#test_linking()


def calc_bathymetry():
    """
    Use thge grid file to calculate a bathymetry field.
    .
    """
    ncgrid = Dataset(grid_nc, 'r')
    bathy_ints = ncgrid.variables['mbathy'][:]
    nav_lev = ncgrid.variables['nav_lev'][:]

    bathy_m = np.zeros_like(bathy_ints[0]) # 2d field
    for (z,y,x), bathy in np.ndenumerate(bathy_ints):
        if np.ma.is_masked(bathy): continue
        if bathy == -2147483647: continue

        bathy_m[y,x] = nav_lev[bathy]
    plot_map(bathy_m, 'bathy', fold = 'images/')

    return bathy_m





###
# Generate a new netCDF.
def generate_ncdf(new_restart, fnkey, debug=False):

    if os.path.exists(new_restart):
        print('file exists', generate_ncdf)
        return
    nc = Dataset(eorca_restart, 'r')
    ncersem = Dataset(ersem_restart, 'r')
    nciMarNet = Dataset(imarnet_restart, 'r')
    print_nc(nc, 'eORCA physics')
    print_nc(ncersem, 'ERSEM restart')

    # Generate new netcdf.
    newnc = Dataset(new_restart, 'w', format='NETCDF4')
    newnc.description = 'New ERSEM RESTART file for Mission Atlantic'

    # Copy over Attributes?

    #  copy dimensions from physics eORCA
    dimensions = {}
    dim_sizes = {}
    for dim in nc.dimensions.keys():
        dimname = nc.dimensions[dim].name
        dimsize = nc.dimensions[dim].size
        if debug:
            dimsize = 1
        if dimname == 'nav_lev': # physics uses nav_lev for this, but ERSEM uses z.
            dimname = 'z'
        if dimname == 'time_counter': # physics uses time_counter for this, but ERSEM uses t.
            dimname = 't'
            dimsize = None
        print('Adding dimension:', dim, dimname, dimsize)
        dimensions[dimname] = newnc.createDimension(dimname, dimsize)
        if dimsize == None: dimsize = 1
        dim_sizes[dimname] = dimsize

    # copy variables from ERSEM
    variables = {}
    All_Created_variables = {}
    for key, var in ncersem.variables.items():
        dims = tuple([dim.name for dim in var.get_dims()])
        print('Adding variable:', key, var.dtype, dims)
        if dims is None: assert 0
        variables[key] = newnc.createVariable(key, var.dtype, dims)
        All_Created_variables[key] = 'unfilled (copy from ersem)'

    # Add "missing" data from iMarNet which isn't in the AMM domain ERSEM (ie Iron)
    imatnettype = nciMarNet.variables['TRBP1f'].dtype
    ersemtype = ncersem.variables['TRBP1_c'].dtype

    bathy_m = calc_bathymetry()
    print(var.dtype, imatnettype, ersemtype)

    missing_vars = [
        ["TR2DQ6_f", imatnettype, ('t', 'y', 'x')],
        ["TRBN7_f", imatnettype, ('t', 'z', 'y', 'x')],
        ["TRBP1_f", imatnettype, ('t', 'z', 'y', 'x')],
        ["TRBP2_f", imatnettype, ('t', 'z', 'y', 'x')],
        ["TRBP3_f", imatnettype, ('t', 'z', 'y', 'x')],
        ["TRBP4_f", imatnettype, ('t', 'z', 'y', 'x')],
        ["TRBR4_f", imatnettype, ('t', 'z', 'y', 'x')],
        ["TRBR6_f", imatnettype, ('t', 'z', 'y', 'x')],
        ["TRNN7_f", imatnettype, ('t', 'z', 'y', 'x')],
        ["TRNP1_f", imatnettype, ('t', 'z', 'y', 'x')],
        ["TRNP2_f", imatnettype, ('t', 'z', 'y', 'x')],
        ["TRNP3_f", imatnettype, ('t', 'z', 'y', 'x')],
        ["TRNP4_f", imatnettype, ('t', 'z', 'y', 'x')],
        ["TRNR4_f", imatnettype, ('t', 'z', 'y', 'x')],
        ["TRNR6_f", imatnettype, ('t', 'z', 'y', 'x')],]

    new_vars = [
        ['TRBO3_TA', imatnettype, ('t', 'z', 'y', 'x')],
        ['TRNO3_TA', imatnettype, ('t', 'z', 'y', 'x')],
        ['fabm_st2DbQ6_f', imatnettype, ('t', 'y', 'x')],
        ['fabm_st2DbQ6_pen_depth_f', imatnettype, ('t', 'y', 'x')],
        ['fabm_st2DbK7_f', imatnettype, ('t', 'y', 'x')],
        ['fabm_st2DnQ6_f', imatnettype, ('t', 'y', 'x')],
        ['fabm_st2DnQ6_pen_depth_f', imatnettype, ('t', 'y', 'x')],
        ['fabm_st2DnK7_f', imatnettype, ('t', 'y', 'x')],
            ]
    new_vars_flats = {
        'fabm_st2DbQ6_f':1.695 ,
        'fabm_st2DbQ6_pen_depth_f': 0.02,
        'fabm_st2DbK7_f': 20.1,
        'fabm_st2DnQ6_f': 1.695,
        'fabm_st2DnQ6_pen_depth_f': 0.02,
        'fabm_st2DnK7_f': 20.1,
        }

    new_benthic_lambdas = {
        'Y2_c': lambda h: 8834.4 * h ** (-0.1651),
        'Y3_c': lambda h: -1300*np.log(h) + 8830.8,
        'Y4_c': lambda h: 0.1689*h + 55.342,

        'H2_c': lambda h: -156.98*np.log(h) + 1139.8,
        'H1_c': lambda h: 0.61*h + 23.98,

        'Q1_c': lambda h: 828.74*h**(-0.6246),
        'Q1_n': lambda h: 21.103*h**(-0.6882),
        'Q1_p' : lambda h: 1.60278*h**(-0.6725),

        'Q6_c': lambda h: -1069*np.log(h) + 10900,
        'Q6_n': lambda h: -7.6368*np.log(h) + 78.564,
        'Q6_p': lambda h: -0.545*np.log(h) + 6.0114,
        'Q6_s': lambda h: -64.598*np.log(h) + 391.61,
        'Q6_pen_depth_c': lambda h: 0.0486*h**0.103,
        'Q6_pen_depth_n': lambda h: 0.0486*h**0.104,
        'Q6_pen_depth_p': lambda h: 0.0484*h**0.1042,
        'Q6_pen_depth_s': lambda h: 1e-5*h + 0.0232,

        'Q7_c': lambda h: 50*(-1069*np.log(h) + 10900),
        'Q7_n': lambda h: 50*(-7.6368*np.log(h) + 78.564),
        'Q7_p': lambda h: 50*(-0.545*np.log(h) + 6.0114),
        'Q7_pen_depth_c': lambda h: 1e-5*h + 0.0155,
        'Q7_pen_depth_n': lambda h: 8e-6*h + 0.0155,
        'Q7_pen_depth_p': lambda h: 0.0495*h**0.0965,

        'Q17_c': lambda h: h*0.,
        'Q17_n': lambda h: h*0.,
        'Q17_p': lambda h: h*0.,
        'G2_o_deep': lambda h: h*0.,
    }
    new_benthic_lambdas_keys = new_benthic_lambdas.keys()
    for key in new_benthic_lambdas_keys:
        new_benthic_lambdas['fabm_st2Db'+key] = new_benthic_lambdas[key]
        new_benthic_lambdas['fabm_st2Dn'+key] = new_benthic_lambdas[key]


    miss_vars = []
    for (key,dtype, dims) in missing_vars:
         print('Adding missing variable:', key, dtype, dims)
         variables[key] = newnc.createVariable(key, dtype, dims)
         All_Created_variables[key] = 'unfilled (missing)'
         miss_vars.append(key)

    for (key,dtype, dims) in new_vars:
         print('Adding new variable:', key, dtype, dims)
         variables[key] = newnc.createVariable(key, dtype, dims)
         All_Created_variables[key] = 'unfilled (new)'

    #fill variables with data:

    for (key,dtype, dims) in new_vars:
        print('Adding new variable data:', key)
        if key in ['TRBO3_TA', 'TRNO3_TA']:
            iMNkey = link_old_to_new_name(key, ncersem, nciMarNet)
            arr = nciMarNet.variables[iMNkey][:]
            arr = extend_ORCA(arr)
            All_Created_variables[key] = 'filled from imarnet'

        elif key in new_vars_flats:
            shape = tuple([dim_sizes[dimname] for dimname in dims])
            arr = np.ones(shape) * new_vars_flats[key]
            All_Created_variables[key] = 'filled_with '+str(new_vars_flats[key])

        else: assert 0
        variables[key][:] = arr

    """
    for (key,dtype, dims) in missing_vars:
        print('Adding missing variable data:', key, 'gets data from iMarNet:', dims)
        iMNkey = link_old_to_new_name(key, ncersem, nciMarNet)
        arr = nciMarNet.variables[iMNkey][:]
        #arr = extend_ORCA(arr)
        if np.ma.is_masked(arr.max()) or arr.max()>1.E10:
            print(key, '->', iMNkey, arr.min(),arr.max())
            assert 0
        arr = extend_ORCA(arr)
        if np.ma.is_masked(arr.max()) or arr.max()>1.E10:
            print(key, '->', iMNkey, arr.min(),arr.max())
            assert 0
        variables[key] = arr
        print(key, '->', iMNkey, arr.min(),arr.max())

#        assert 0
        All_Created_variables[key] = 'filled'
        print(arr)
        # assert 0
    """

    #fill variables with data:
    for key, var in ncersem.variables.items():
        dims = tuple([dim.name for dim in var.get_dims()])
        shape = tuple([dim_sizes[dimname] for dimname in dims])
        if key in miss_vars: assert 0
        if debug:
            print('Adding variable data:', key, shape)
            variables[key][:] = np.ones(shape)
            All_Created_variables[key] = 'filled_with_ones'

        else:
            # values where we take data from the physics model:
            if key in ['nav_lon', 'nav_lat' , 'nav_lev', 'time_counter', 'kt', 'ndastp', 'adatrj']:
                print('Adding variable data:', key, 'gets data from physics!', shape, nc.variables[key].shape, dims)
                variables[key][:] = nc.variables[key][:]
                All_Created_variables[key] = 'filled_from_physics'
                continue

            iMNkey = link_old_to_new_name(key, ncersem, nciMarNet)
            if iMNkey.upper() == 'ZERO':
                print('Adding variable data:', key, 'gets set to zero:',iMNkey, shape, dims)
                variables[key][:] = np.zeros(shape)
                All_Created_variables[key] = 'filled_with_zeroes'

            else:
                print('Adding variable data:', key, 'gets data from iMarNet:',iMNkey, nciMarNet.variables[iMNkey].shape, '->', shape, dims)
                arr = extend_ORCA(nciMarNet.variables[iMNkey][:])
                if np.ma.is_masked(arr.max()):
                    print('Fully masked:', key, var, iMNkey)
                variables[key][:] = arr
                All_Created_variables[key] = 'filled_with_iMNkey'

    for (key,dtype, dims) in missing_vars:
        print('Adding missing variable data:', key, 'gets data from iMarNet:', dims)
        iMNkey = link_old_to_new_name(key, ncersem, nciMarNet)
        arr = nciMarNet.variables[iMNkey][:]
        #arr = extend_ORCA(arr)
        if np.ma.is_masked(arr.max()) or arr.max()>1.E10:
            print(key, '->', iMNkey, arr.min(),arr.max())
            assert 0
        arr = extend_ORCA(arr)
        if np.ma.is_masked(arr.max()) or arr.max()>1.E10:
            print(key, '->', iMNkey, arr.min(),arr.max())
            assert 0
        variables[key][:] = arr
        print(key, '->', iMNkey, arr.min(),arr.max())

#        assert 0
        All_Created_variables[key] = 'filled at last'
        print('iMarNet iron', arr.min(), arr.max(), variables[key][:].min(), variables[key][:].max() )

#    assert 0


    for key in sorted(All_Created_variables.keys()):
        print(key, All_Created_variables[key])
    #assert 0
    # data
    print('Saving:', new_restart)
    newnc.close()


def plot_all(fn, fnkey):
    nc = Dataset(fn, 'r')
    for key, var in nc.variables.items():
        ndim = var.ndim
        if ndim == 2:
            data = var[:]
        elif ndim == 3:
            data = var[0]
        elif ndim == 4:
            data = var[0,0]
        else:
            #rint('skip plot:', key)
            continue
        if data.min()==data.max() and data.min()!=0.:
            print(key,':', data.min())
        if np.ma.is_masked(data.max()):
            print('Fully masked:', key, var, fnkey)
            #assert 0
            fold = 'images/'+fnkey+'/fully_masked/'
        fold = 'images/'+fnkey+'/'
        try: os.mkdir(fold)
        except: pass
        plot_map(data, fnkey + ' '+key, fold = fold)
    nc.close()


#plot_all(new_restart, 'newrestart')
#plot_all(imarnet_restart, 'iMarNet')
#generate_ncdf()
# plot_all(new_restart, 'newrestart_v3')
# ncersem = Dataset(ersem_restart, 'r')

def main():
    calc_bathymetry()
    assert 0
    new_restart = 'input/sthenno1/restart_trc_v9.nc'

    generate_ncdf(new_restart, 'restart_trc_v9')
    plot_all(new_restart, 'restart_trc_v9')

main()
