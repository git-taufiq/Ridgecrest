import numpy as np
from scipy.interpolate import interp1d

node64_y = np.asarray([3850000.0,3936500.0,3946500.0,3950000.0,3951250.0,3956000.0,3964000.0,4052000.0])
node64_z = np.asarray([   2000.0,   2000.0, -10500.0, -11200.0, -11200.0, -10200.0,   2000.0,   2000.0])

f64 = interp1d(node64_y, node64_z, kind='cubic')
grid = 100
node64_y_new = np.arange(3850000.0,4052000.0+grid,grid)
node64_z_new = f64(node64_y_new)
node64_z_new[np.argwhere(node64_y_new<node64_y[1])] = 2000.0
node64_z_new[np.argwhere(node64_y_new>node64_y[-2])] = 2000.0

node71_y = np.asarray([3850000.0,3932000.0,3943000.0,3950000.0,3951250.0,3962000.0,3970000.0,3975000.0,4052000.0])
node71_z = np.asarray([   2000.0,   2000.0,  -9500.0, -11200.0, -11200.0, -10000.0,  -6000.0,   2000.0,   2000.0])


f71 = interp1d(node71_y, node71_z, kind='cubic')
grid = 100
node71_y_new = np.arange(3850000.0,4052000.0+grid,grid)
node71_z_new = f71(node71_y_new)
node71_z_new[np.argwhere(node71_y_new<node71_y[1])] = 2000.0
node71_z_new[np.argwhere(node71_y_new>node71_y[-2])] = 2000.0

ny = len(node64_y_new)
print('writing SeismoDepth64.nc netcdf')
####Creating the netcdf file
from netCDF4 import Dataset
rootgrp = Dataset("output_files/seismoDepth64.nc", "w", format="NETCDF4")
rootgrp.createDimension("y", ny)
vy = rootgrp.createVariable("y","f4",("y",))
vy[:]= node64_y_new
vSeismoDepth =  rootgrp.createVariable("seismoDepth64","f4",("y"))
vSeismoDepth[:] = node64_z_new
rootgrp.close()

ny = len(node71_y_new)
print('writing SeismoDepth71.nc netcdf')
####Creating the netcdf file
from netCDF4 import Dataset
rootgrp = Dataset("output_files/seismoDepth71.nc", "w", format="NETCDF4")
rootgrp.createDimension("y", ny)
vy = rootgrp.createVariable("y","f4",("y",))
vy[:]= node71_y_new
vSeismoDepth =  rootgrp.createVariable("seismoDepth71","f4",("y"))
vSeismoDepth[:] = node71_z_new
rootgrp.close()
