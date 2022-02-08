
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


imarnet_forcing_path = "/data/sthenno1/to_archive/ledm/Projects/iMarNet/OUTPUT_v2/iMarNet_forcing/"
iMarNet_inputs = "/data/sthenno1/to_archive/ledm/Projects/iMarNet/inputNetCDFs/"

imarnet_restart = 'input/iMARNET_fields_1890_UEA_ERSEM_newFeR3c_20140919.nc'
eorca_restart = 'input/eORCA1_00000992_restart.nc'
ersem_restart = 'input/sthenno1/restart_trc.HC.1986-1990.12.nc'

new_restart = 'input/sthenno1/new_restart_v2.nc'

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

# def print_nc(nc, key):
#     print('\n------------------\nNetcdf:', key)
#
#     #print('keys', sorted(nc.variables.keys()))
#     for dim in nc.dimensions.values():
#             print(dim, dim.name)
#
#     for key in  sorted(nc.variables.keys()):
#         continue
#         v = nc.variables[key]
#         #print(v.long_name)
#         #ln = v.long_name
#         #units = v.units
#         try: ln = v.long_name
#         except: ln='long_name'
#
#         try:units = v.units
#         except: units = 'units'
#
#         print(key,':\t', ln, '\t',units, '\t',v.shape)
#
# master_dict = {
#         'TRBP1_Chl' : 'TRBChl1',
#         'TRBP2_Chl' : 'TRBChl2',
#         'TRBP3_Chl' : 'TRBChl3',
#         'TRBP4_Chl' : 'TRBChl4',
#         'TRNP1_Chl' : 'TRNChl1',
#         'TRNP2_Chl' : 'TRNChl2',
#         'TRNP3_Chl' : 'TRNChl3',
#         'TRNP4_Chl' : 'TRNChl4',
#
#         'TRBO3_bioalk' : 'TRBbioalk',
#         'TRNO3_bioalk' : 'TRNbioalk',
#
#         'TRNlight_ADY' : 'ZERO', # Need to write to James
#         'TRBlight_ADY' : 'ZERO',
#         'fabm_st2DbG2_o_deep' : 'ZERO', # ask Gennadi, put zero for now
#         'fabm_st2DnG2_o_deep' : 'ZERO',
#         'fabm_st2DbbL2_c' : 'TR2Dbl2c',
#         'fabm_st2DnbL2_c' : 'TR2Dbl2c',
#         'fabm_st2Dbben_nit_G4n' : 'TR2DG4n',# ask Gennadi, put zero for now
#         'fabm_st2Dnben_nit_G4n' : 'TR2DG4n',
#         'rdttrc1' : 'ZERO',
#
#         # Iron:
# 		'fabm_st2DnQ6_pen_depth_f': 'TR2DQ6f',
# 		'fabm_st2DbQ6_pen_depth_f': 'TR2DQ6f',
# 		'TRBN7_f': 'TRBN7f',
# 		'TRBP1_f': 'TRBP1f',
# 		'TRBP2_f': 'TRBP2f',
# 		'TRBP3_f': 'TRBP3f',
# 		'TRBP4_f': 'TRBP4f',
# 		'TRBR4_f': 'TRBR4f',
# 		'TRBR6_f': 'TRBR6f',
# 		'TRNN7_f': 'TRNN7f',
# 		'TRNP1_f': 'TRNP1f',
# 		'TRNP2_f': 'TRNP2f',
# 		'TRNP3_f': 'TRNP3f',
# 		'TRNP4_f': 'TRNP4f',
# 		'TRNR4_f': 'TRNR4f',
# 		'TRNR6_f': 'TRNR6f',
#
#         'fabm_st2DnQ6_pen_depth_c' : 'TR2DD6m',
#         'fabm_st2DnQ6_pen_depth_n' : 'TR2DD7m',
#         'fabm_st2DnQ6_pen_depth_p' : 'TR2DD8m',
#         'fabm_st2DnQ6_pen_depth_s' : 'TR2DD9m',
#         'fabm_st2DnQ7_pen_depth_c' : 'TR2DD3m',
#         'fabm_st2DnQ7_pen_depth_n' : 'TR2DD4m',
#         'fabm_st2DnQ7_pen_depth_p' : 'TR2DD5m',
#         'fabm_st2Dnben_col_D1m' : 'TR2DD1m',
#         'fabm_st2Dnben_col_D2m' : 'TR2DD2m',
#
#         'fabm_st2DbQ6_pen_depth_c': 'TR2DD6m',
#         'fabm_st2DbQ6_pen_depth_n': 'TR2DD7m',
#         'fabm_st2DbQ6_pen_depth_p': 'TR2DD8m',
#         'fabm_st2DbQ6_pen_depth_s': 'TR2DD9m',
#         'fabm_st2DbQ7_pen_depth_c': 'TR2DD3m',
#         'fabm_st2DbQ7_pen_depth_n': 'TR2DD4m',
#         'fabm_st2DbQ7_pen_depth_p': 'TR2DD5m',
#         'fabm_st2Dbben_col_D1m': 'TR2DD1m',
#         'fabm_st2Dbben_col_D2m': 'TR2DD2m',
#         }

# def link_old_to_new_name(new_name, nc1, nc2):
#     """
#     ie convert TRNB1n to TRNB1_n
#     new_name is key from recent ERSEM netcdf restart
#     nc1 is recent ERSEM restart
#     nc2 is old iMarNEt ERSEM restart.
#     """
#
#     # Bespoke list:
#     name = master_dict.get(new_name, "")
#     if name != "":
#         return name
#
#     # Identical bewteen lists:
#     nc2vars = list(nc2.variables.keys())
#     if new_name in nc2vars:
#         return new_name
#
#     # Yuri said to set these to zero:
#     if new_name[:3] in ['sbc', 'atf', 'rdb']:
#         return 'ZERO'
#
#     # Create list of iMarNet keys without preficeSx.
# #    nc2_real_names = {}
# #    for n in nc2vars:
# #        nn = n.replace('TRN','').replace('TRB', '')
# #        try: nc2_real_names[nn].append(n)
# #        except: nc2_real_names[nn] = [n,]
#
#     # Replace some bits in modern ERSEM:
#     replacements = {'_': '',
#         'fabm': '',
#         'st2Db':'',
#         'st2Dn': '',
#             }
#     appends = ['', 'TRN', 'TRB', 'TR2D',]
#
#     new_name_copy = new_name[:]
#     for r,a in replacements.items():
#         new_name_copy = new_name_copy.replace(r,a)
#
#         for app in appends:
#             print(new_name, 'trying:', app+new_name_copy)
#             if app+new_name_copy in nc2vars:
#                 return app+new_name_copy
#
#     print('not found:', new_name, '->', new_name_copy)
#     return False

# def run():
#     orcnc = Dataset(eorca_restart, 'r')
#     print_nc(orcnc, 'eORCA')
#     """Netcdf: eORCA
#     keys ['a_fwb', 'a_fwb_b', 'adatrj', 'avm_k', 'avt_k', 'dissl', 'e3t_b', 'e3t_n', 'emp_b', 'en', 'fraqsr_1lev', 'fwf_isf_b', 'isf_hc_b', 'isf_sc_b', 'kt', 'nav_lat', 'nav_lev', 'nav_lon', 'ndastp', 'ntime', 'q
#     ns_b', 'qsr_hc_b', 'rdt', 'rhop', 'rnf_b', 'rnf_hc_b', 'rnf_sc_b', 'sb', 'sbc_hc_b', 'sbc_sc_b', 'sfx_b', 'sn', 'sshb', 'sshn', 'tb', 'time_counter', 'tn', 'ub', 'ub2_b', 'un', 'un_bf', 'utau_b', 'vb', 'vb2_b
#     ', 'vn', 'vn_bf', 'vtau_b']
#     <class 'netCDF4._netCDF4.Dimension'>: name = 'x', size = 362
#     <class 'netCDF4._netCDF4.Dimension'>: name = 'y', size = 332
#     <class 'netCDF4._netCDF4.Dimension'>: name = 'nav_lev', size = 75
#     <class 'netCDF4._netCDF4.Dimension'> (unlimited): name = 'time_counter', size = 1"""
#     #temp = orcnc.variables[
#
#     ersnc = Dataset(ersem_restart, 'r')
#     print_nc(ersnc, 'ERSEM')
#     """Netcdf: ERSEM
#     keys ['TRBB1_c', 'TRBB1_n', 'TRBB1_p', 'TRBL2_c', 'TRBN1_p', 'TRBN3_n', 'TRBN4_n', 'TRBN5_s', 'TRBO2_o', 'TRBO3_bioalk', 'TRBO3_c', 'TRBP1_Chl', 'TRBP1_c', 'TRBP1_n', 'TRBP1_p', 'TRBP1_s', 'TRBP2_Chl', 'TRBP2
#     _c', 'TRBP2_n', 'TRBP2_p', 'TRBP3_Chl', 'TRBP3_c', 'TRBP3_n', 'TRBP3_p', 'TRBP4_Chl', 'TRBP4_c', 'TRBP4_n', 'TRBP4_p', 'TRBR1_c', 'TRBR1_n', 'TRBR1_p', 'TRBR2_c', 'TRBR3_c', 'TRBR4_c', 'TRBR4_n', 'TRBR4_p', '
#     TRBR6_c', 'TRBR6_n', 'TRBR6_p', 'TRBR6_s', 'TRBR8_c', 'TRBR8_n', 'TRBR8_p', 'TRBR8_s', 'TRBZ4_c', 'TRBZ5_c', 'TRBZ5_n', 'TRBZ5_p', 'TRBZ6_c', 'TRBZ6_n', 'TRBZ6_p', 'TRBlight_ADY', 'TRNB1_c', 'TRNB1_n', 'TRNB1
#     _p', 'TRNL2_c', 'TRNN1_p', 'TRNN3_n', 'TRNN4_n', 'TRNN5_s', 'TRNO2_o', 'TRNO3_bioalk', 'TRNO3_c', 'TRNP1_Chl', 'TRNP1_c', 'TRNP1_n', 'TRNP1_p', 'TRNP1_s', 'TRNP2_Chl', 'TRNP2_c', 'TRNP2_n', 'TRNP2_p', 'TRNP3_
#     Chl', 'TRNP3_c', 'TRNP3_n', 'TRNP3_p', 'TRNP4_Chl', 'TRNP4_c', 'TRNP4_n', 'TRNP4_p', 'TRNR1_c', 'TRNR1_n', 'TRNR1_p', 'TRNR2_c', 'TRNR3_c', 'TRNR4_c', 'TRNR4_n', 'TRNR4_p', 'TRNR6_c', 'TRNR6_n', 'TRNR6_p', 'T
#     RNR6_s', 'TRNR8_c', 'TRNR8_n', 'TRNR8_p', 'TRNR8_s', 'TRNZ4_c', 'TRNZ5_c', 'TRNZ5_n', 'TRNZ5_p', 'TRNZ6_c', 'TRNZ6_n', 'TRNZ6_p', 'TRNlight_ADY', 'adatrj', 'atf_trend_B1_c', 'atf_trend_B1_n', 'atf_trend_B1_p'
#     , 'atf_trend_L2_c', 'atf_trend_N1_p', 'atf_trend_N3_n', 'atf_trend_N4_n', 'atf_trend_N5_s', 'atf_trend_O2_o', 'atf_trend_O3_bioalk', 'atf_trend_O3_c', 'atf_trend_P1_Chl', 'atf_trend_P1_c', 'atf_trend_P1_n', '
#     atf_trend_P1_p', 'atf_trend_P1_s', 'atf_trend_P2_Chl', 'atf_trend_P2_c', 'atf_trend_P2_n', 'atf_trend_P2_p', 'atf_trend_P3_Chl', 'atf_trend_P3_c', 'atf_trend_P3_n', 'atf_trend_P3_p', 'atf_trend_P4_Chl', 'atf_
#     trend_P4_c', 'atf_trend_P4_n', 'atf_trend_P4_p', 'atf_trend_R1_c', 'atf_trend_R1_n', 'atf_trend_R1_p', 'atf_trend_R2_c', 'atf_trend_R3_c', 'atf_trend_R4_c', 'atf_trend_R4_n', 'atf_trend_R4_p', 'atf_trend_R6_c
#     ', 'atf_trend_R6_n', 'atf_trend_R6_p', 'atf_trend_R6_s', 'atf_trend_R8_c', 'atf_trend_R8_n', 'atf_trend_R8_p', 'atf_trend_R8_s', 'atf_trend_Z4_c', 'atf_trend_Z5_c', 'atf_trend_Z5_n', 'atf_trend_Z5_p', 'atf_tr
#     end_Z6_c', 'atf_trend_Z6_n', 'atf_trend_Z6_p', 'atf_trend_light_ADY', 'fabm_st2DbG2_o', 'fabm_st2DbG2_o_deep', 'fabm_st2DbG3_c', 'fabm_st2DbH1_c', 'fabm_st2DbH2_c', 'fabm_st2DbK1_p', 'fabm_st2DbK3_n', 'fabm_s
#     t2DbK4_n', 'fabm_st2DbK5_s', 'fabm_st2DbQ17_c', 'fabm_st2DbQ17_n', 'fabm_st2DbQ17_p', 'fabm_st2DbQ1_c', 'fabm_st2DbQ1_n', 'fabm_st2DbQ1_p', 'fabm_st2DbQ6_c', 'fabm_st2DbQ6_n', 'fabm_st2DbQ6_p', 'fabm_st2DbQ6_
#     pen_depth_c', 'fabm_st2DbQ6_pen_depth_n', 'fabm_st2DbQ6_pen_depth_p', 'fabm_st2DbQ6_pen_depth_s', 'fabm_st2DbQ6_s', 'fabm_st2DbQ7_c', 'fabm_st2DbQ7_n', 'fabm_st2DbQ7_p', 'fabm_st2DbQ7_pen_depth_c', 'fabm_st2D
#     bQ7_pen_depth_n', 'fabm_st2DbQ7_pen_depth_p', 'fabm_st2DbY2_c', 'fabm_st2DbY3_c', 'fabm_st2DbY4_c', 'fabm_st2DbbL2_c', 'fabm_st2Dbben_col_D1m', 'fabm_st2Dbben_col_D2m', 'fabm_st2Dbben_nit_G4n', 'fabm_st2DnG2_
#     o', 'fabm_st2DnG2_o_deep', 'fabm_st2DnG3_c', 'fabm_st2DnH1_c', 'fabm_st2DnH2_c', 'fabm_st2DnK1_p', 'fabm_st2DnK3_n', 'fabm_st2DnK4_n', 'fabm_st2DnK5_s', 'fabm_st2DnQ17_c', 'fabm_st2DnQ17_n', 'fabm_st2DnQ17_p'
#     , 'fabm_st2DnQ1_c', 'fabm_st2DnQ1_n', 'fabm_st2DnQ1_p', 'fabm_st2DnQ6_c', 'fabm_st2DnQ6_n', 'fabm_st2DnQ6_p', 'fabm_st2DnQ6_pen_depth_c', 'fabm_st2DnQ6_pen_depth_n', 'fabm_st2DnQ6_pen_depth_p', 'fabm_st2DnQ6_
#     pen_depth_s', 'fabm_st2DnQ6_s', 'fabm_st2DnQ7_c', 'fabm_st2DnQ7_n', 'fabm_st2DnQ7_p', 'fabm_st2DnQ7_pen_depth_c', 'fabm_st2DnQ7_pen_depth_n', 'fabm_st2DnQ7_pen_depth_p', 'fabm_st2DnY2_c', 'fabm_st2DnY3_c', 'f
#     abm_st2DnY4_c', 'fabm_st2DnbL2_c', 'fabm_st2Dnben_col_D1m', 'fabm_st2Dnben_col_D2m', 'fabm_st2Dnben_nit_G4n', 'kt', 'nav_lat', 'nav_lev', 'nav_lon', 'ndastp', 'rdb_trend_B1_c', 'rdb_trend_B1_n', 'rdb_trend_B1
#     _p', 'rdb_trend_L2_c', 'rdb_trend_N1_p', 'rdb_trend_N3_n', 'rdb_trend_N4_n', 'rdb_trend_N5_s', 'rdb_trend_O3_c', 'rdb_trend_P1_Chl', 'rdb_trend_P1_c', 'rdb_trend_P1_n', 'rdb_trend_P1_p', 'rdb_trend_P1_s', 'rd
#     b_trend_P2_Chl', 'rdb_trend_P2_c', 'rdb_trend_P2_n', 'rdb_trend_P2_p', 'rdb_trend_P3_Chl', 'rdb_trend_P3_c', 'rdb_trend_P3_n', 'rdb_trend_P3_p', 'rdb_trend_P4_Chl', 'rdb_trend_P4_c', 'rdb_trend_P4_n', 'rdb_tr
#     end_P4_p', 'rdb_trend_R1_c', 'rdb_trend_R1_n', 'rdb_trend_R1_p', 'rdb_trend_R2_c', 'rdb_trend_R3_c', 'rdb_trend_R4_c', 'rdb_trend_R4_n', 'rdb_trend_R4_p', 'rdb_trend_R6_c', 'rdb_trend_R6_n', 'rdb_trend_R6_p',
#      'rdb_trend_R6_s', 'rdb_trend_R8_c', 'rdb_trend_R8_n', 'rdb_trend_R8_p', 'rdb_trend_R8_s', 'rdb_trend_Z4_c', 'rdb_trend_Z5_c', 'rdb_trend_Z5_n', 'rdb_trend_Z5_p', 'rdb_trend_Z6_c', 'rdb_trend_Z6_n', 'rdb_tren
#     d_Z6_p', 'rdb_trend_light_ADY', 'rdttrc1', 'sbc_B1_c_b', 'sbc_B1_n_b', 'sbc_B1_p_b', 'sbc_L2_c_b', 'sbc_N1_p_b', 'sbc_N3_n_b', 'sbc_N4_n_b', 'sbc_N5_s_b', 'sbc_O2_o_b', 'sbc_O3_bioalk_b', 'sbc_O3_c_b', 'sbc_P
#     1_Chl_b', 'sbc_P1_c_b', 'sbc_P1_n_b', 'sbc_P1_p_b', 'sbc_P1_s_b', 'sbc_P2_Chl_b', 'sbc_P2_c_b', 'sbc_P2_n_b', 'sbc_P2_p_b', 'sbc_P3_Chl_b', 'sbc_P3_c_b', 'sbc_P3_n_b', 'sbc_P3_p_b', 'sbc_P4_Chl_b', 'sbc_P4_c_
#     b', 'sbc_P4_n_b', 'sbc_P4_p_b', 'sbc_R1_c_b', 'sbc_R1_n_b', 'sbc_R1_p_b', 'sbc_R2_c_b', 'sbc_R3_c_b', 'sbc_R4_c_b', 'sbc_R4_n_b', 'sbc_R4_p_b', 'sbc_R6_c_b', 'sbc_R6_n_b', 'sbc_R6_p_b', 'sbc_R6_s_b', 'sbc_R8_
#     c_b', 'sbc_R8_n_b', 'sbc_R8_p_b', 'sbc_R8_s_b', 'sbc_Z4_c_b', 'sbc_Z5_c_b', 'sbc_Z5_n_b', 'sbc_Z5_p_b', 'sbc_Z6_c_b', 'sbc_Z6_n_b', 'sbc_Z6_p_b', 'sbc_light_ADY_b', 'time_counter']
#     <class 'netCDF4._netCDF4.Dimension'>: name = 'y', size = 375
#     <class 'netCDF4._netCDF4.Dimension'>: name = 'x', size = 297
#     <class 'netCDF4._netCDF4.Dimension'>: name = 'z', size = 51
#     <class 'netCDF4._netCDF4.Dimension'> (unlimited): name = 't', size = 1"""
#
#
#     imnnc = Dataset(imarnet_restart, 'r')
#     print_nc(imnnc, 'iMarNet')
#
#     """Netcdf: iMarNet
#     keys ['TR2DD1m', 'TR2DD2m', 'TR2DD3m', 'TR2DD4m', 'TR2DD5m', 'TR2DD6m', 'TR2DD7m', 'TR2DD8m', 'TR2DD9m', 'TR2DG2o', 'TR2DG3c', 'TR2DG4n', 'TR2DH1c', 'TR2DH2c', 'TR2DK1p', 'TR2DK3n', 'TR2DK4n', 'TR2DK5s', 'TR2
#     DQ17c', 'TR2DQ17n', 'TR2DQ17p', 'TR2DQ1c', 'TR2DQ1n', 'TR2DQ1p', 'TR2DQ6c', 'TR2DQ6f', 'TR2DQ6n', 'TR2DQ6p', 'TR2DQ6s', 'TR2DQ7c', 'TR2DQ7n', 'TR2DQ7p', 'TR2DY2c', 'TR2DY3c', 'TR2DY4c', 'TR2Dbl2c', 'TRBB1c',
#     'TRBB1n', 'TRBB1p', 'TRBChl1', 'TRBChl2', 'TRBChl3', 'TRBChl4', 'TRBL2c', 'TRBN1p', 'TRBN3n', 'TRBN4n', 'TRBN5s', 'TRBN7f', 'TRBO2o', 'TRBO3c', 'TRBP1c', 'TRBP1f', 'TRBP1n', 'TRBP1p', 'TRBP1s', 'TRBP2c', 'TRB
#     P2f', 'TRBP2n', 'TRBP2p', 'TRBP3c', 'TRBP3f', 'TRBP3n', 'TRBP3p', 'TRBP4c', 'TRBP4f', 'TRBP4n', 'TRBP4p', 'TRBR1c', 'TRBR1n', 'TRBR1p', 'TRBR2c', 'TRBR3c', 'TRBR4c', 'TRBR4f', 'TRBR4n', 'TRBR4p', 'TRBR6c', 'T
#     RBR6f', 'TRBR6n', 'TRBR6p', 'TRBR6s', 'TRBR8c', 'TRBR8n', 'TRBR8p', 'TRBR8s', 'TRBZ4c', 'TRBZ5c', 'TRBZ5n', 'TRBZ5p', 'TRBZ6c', 'TRBZ6n', 'TRBZ6p', 'TRBbioalk', 'TRNB1c', 'TRNB1n', 'TRNB1p', 'TRNChl1', 'TRNCh
#     l2', 'TRNChl3', 'TRNChl4', 'TRNL2c', 'TRNN1p', 'TRNN3n', 'TRNN4n', 'TRNN5s', 'TRNN7f', 'TRNO2o', 'TRNO3c', 'TRNP1c', 'TRNP1f', 'TRNP1n', 'TRNP1p', 'TRNP1s', 'TRNP2c', 'TRNP2f', 'TRNP2n', 'TRNP2p', 'TRNP3c', '
#     TRNP3f', 'TRNP3n', 'TRNP3p', 'TRNP4c', 'TRNP4f', 'TRNP4n', 'TRNP4p', 'TRNR1c', 'TRNR1n', 'TRNR1p', 'TRNR2c', 'TRNR3c', 'TRNR4c', 'TRNR4f', 'TRNR4n', 'TRNR4p', 'TRNR6c', 'TRNR6f', 'TRNR6n', 'TRNR6p', 'TRNR6s',
#      'TRNR8c', 'TRNR8n', 'TRNR8p', 'TRNR8s', 'TRNZ4c', 'TRNZ5c', 'TRNZ5n', 'TRNZ5p', 'TRNZ6c', 'TRNZ6n', 'TRNZ6p', 'TRNbioalk', 'adatrj', 'arak0', 'kt', 'nav_lat', 'nav_lev', 'nav_lon', 'ndastp', 'time_counter']
#
#     <class 'netCDF4._netCDF4.Dimension'>: name = 'x', size = 362
#     <class 'netCDF4._netCDF4.Dimension'>: name = 'y', size = 292
#     <class 'netCDF4._netCDF4.Dimension'>: name = 'z', size = 75
#     <class 'netCDF4._netCDF4.Dimension'> (unlimited): name = 't', size = 1 """
#     outnc = Dataset(new_restart, 'w')

    # tasks:
    # load data
    # for each variable in the ersem restart file
    #    find the relevant imarnet field
    #    load the iMarNEt data
    #    expand the iMarNet field to the same size as the new eORCA1 grid.
    #
    # save output

def plot_map(arr, key,fold='images/'):
    pyplot.pcolormesh(arr)
    pyplot.colorbar()
    pyplot.title(key)
    fn = fold+key+'.png'
    print('Saving figure:', fn)
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
    attributes['Regrid'] = 'Converted to eORCA1 grid (from ORCA1 grid)'
    print(attributes)
    for att,attribute in attributes.items():
        nco.setncattr(att, attribute)
        if debug: print ('Adding attribute: ',att,'\t(',attribute,')')


    # copy dimensions from input grid, fng/ncg
    for dim, dimension in ncg.dimensions.items():
        dimSize=len(ncg.dimensions[dim])
        if dim in ['time', 'TIME', 'time_counter']:
            dimSize = None
        if debug: print ('Adding dimension: ',dim,'\t(',dimSize,')')
        nco.createDimension(dim, dimSize)

    # list of variables to save
    keys = [ 'nav_lon', 'nav_lat', 'nav_lev', 'time_counter', ] # ncg.variables.keys()

    # removeVar = 'pCO2a' # from ncg
    # if fni.find('ndep-nhx')>-1:
    #     newVars = ['nhx', ] # from nci
    # if fni.find('ndep-noy') > -2:
    #     newVars = ['noy', ] # from nci
    newVars = ['dust', ]
    removeVar = 'sfx_b' # a 2D+t var

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

        # old_time = nci.variables['time'][:]
        # old_lats = nci.variables['lat'][:]
        # old_lons = nci.variables['lon'][:]

        # new_lats = nco.variables['LAT'][:]
        # new_lons = nco.variables['LON'][:]

        # out_arr = np.ma.empty((len(old_time), len(new_lats), len(new_lons)))
        out_arr =  extend_ORCA(nci.variables[var][:,:,:])
        print('making out array like:', out_arr.shape)
        # for t in np.arange(old_time.size)[:]:
        #     if Twelve_months_only and t >=12: continue
        #     if t%6==0:
        #         print('iterating', t, 'of', len(old_time))
        #     out_arr[t] = regrid(nci.variables[var][t,:,:], old_lats, old_lons, new_lats,new_lons, method='linear',)

        nco.variables[var][:] = np.ma.array(out_arr)
        print('successfully populated:\t', var)


    # Close netcdfs:
    nco.close()
    nci.close()
    ncg.close()
    print('successfully created:\t', fno)



def main():
    input_duct_fn = 'input/dust.orca.nM.nc'
    new_grid = 'input/eORCA1_00000992_restart.nc'
    output_dust_fn = 'output/dust_orca_nM_eORCA1_v1.nc'
    run(input_duct_fn, new_grid, output_dust_fn)

main()
