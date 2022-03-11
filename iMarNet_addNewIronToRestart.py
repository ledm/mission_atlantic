from changeNC import changeNC, AutoVivification
from ncdfView import ncdfView
from deMoraTools import folder
from numpy import zeros
import numpy as np
oldFile = '/data/euryale7/scratch/ledm/iMarNet/RESTARTS/newSiRestart/iMARNET_fields_1890_UEA_ERSEM_newSi_flatField.nc'

ironNC =ncdfView('/home_nfs/momm/ORCA100-inputs/ORCA100-75-ICs-N-O2-DIC-Fe.nc',Quiet=True)
newFile = folder('/data/euryale7/scratch/ledm/iMarNet/RESTARTS/newFe_NotErics/')+'iMARNET_fields_1890_UEA_ERSEM_newFeR3c_20140919.nc'


# Rebuilding a restart file with Momme's old field., which has identical depth levels and ORCA grid.

iron = ironNC('N7f')[:]
inds = np.where(iron.mask)
iron[inds] = 1.

iron = iron.data

av = AutoVivification()
av['newVar']['TRNN7f']['name'] = 'TRNN7f'
av['newVar']['TRBN7f']['name'] = 'TRBN7f'
av['newVar']['TRNN7f']['long_name'] = 'iron'
av['newVar']['TRBN7f']['long_name'] = 'iron'
av['newVar']['TRNN7f']['units'] = ''
av['newVar']['TRBN7f']['units'] = ''
av['newVar']['TRNN7f']['dtype'] = 'f8'
av['newVar']['TRBN7f']['dtype'] = 'f8'
av['newVar']['TRNN7f']['newDims'] =  [u't', u'z', u'y', u'x'] 
av['newVar']['TRBN7f']['newDims'] =  [u't', u'z', u'y', u'x'] 
av['newVar']['TRNN7f']['newData'] = iron
av['newVar']['TRBN7f']['newData'] = iron

r3c = np.zeros_like(iron) + 0.1
av['newVar']['TRNR3c']['name'] = 'TRNR3c'
av['newVar']['TRBR3c']['name'] = 'TRBR3c'
av['newVar']['TRNR3c']['long_name'] = ''
av['newVar']['TRBR3c']['long_name'] = ''
av['newVar']['TRNR3c']['units'] = ''
av['newVar']['TRBR3c']['units'] = ''
av['newVar']['TRNR3c']['dtype'] = 'f8'
av['newVar']['TRBR3c']['dtype'] = 'f8'
av['newVar']['TRNR3c']['newDims'] =  [u't', u'z', u'y', u'x'] 
av['newVar']['TRBR3c']['newDims'] =  [u't', u'z', u'y', u'x'] 
av['newVar']['TRNR3c']['newData'] = r3c
av['newVar']['TRBR3c']['newData'] = r3c


av['att']['Description'] = 'RestartFile based on '+oldFile+' and /home_nfs/momm/ORCA100-inputs/ORCA100-75-ICs-N-O2-DIC-Fe.nc, includes new R3c tracer. Created 2014-09-18 by ledm'

a = changeNC( oldFile, newFile, av,)

