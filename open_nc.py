#!/usr/bin/ipython -i

from netCDF4 import Dataset
from sys import argv

if __name__ == "__main__":
    print('opening:', argv[1])
    nc = Dataset(argv[1], 'r')
    for key in nc.variables.keys():
        #if key.find('_f')==-1: continue
        if key.find('_f')>-1 or key[-1]=='f':
                print('\n', key, nc.variables[key].shape)
                #if key.find('_f')>-1:
                print('range:', nc.variables[key][:].min(), nc.variables[key][:].max())
                print(nc.variables[key])
