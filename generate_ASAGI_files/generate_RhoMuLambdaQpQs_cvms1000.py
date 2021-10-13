import numpy as np

X_min = 352000
X_max = 556000
Y_min = 3850000
Y_max = 4052000
Z_min =-5000
Z_max = 105000

#Vp = 14; Vs = 15; Rho = 16
VpVsRho = np.loadtxt("input_files/Ridgecrest_cvms1000", usecols=[14,15,16])
print('done reading ascii Ridgecrest_cvms1000')
nx=205
ny=203
nz=111

VpVsRho = VpVsRho.reshape((nz,ny,nx,3))

#compute rho average of the column of earth above.
rho = VpVsRho[:,:,:,2]
rho_col_average = np.cumsum(rho,axis=0)
num_data = np.linspace(1,nz,nz)-5
num_data[num_data<=0] = 1
rho_col_average = rho_col_average/num_data[:,None,None]

#copy the 5th slices to above (where no data)
for i in range(0,5):
   VpVsRho[i,:,:,:] = VpVsRho[5,:,:,:]
   rho_col_average[i,:,:] = rho_col_average[5,:,:]

rho = VpVsRho[:,:,:,2]
mu  = VpVsRho[:,:,:,1] ** 2 * rho
lam = VpVsRho[:,:,:,0] ** 2 * rho - 2 * mu

# CVM-S4.26 attenuation relationship (based on Olsen et al., 2003)
qs = VpVsRho[:,:,:,1] * 1
qs[qs<1500.0] *= 0.02
qs[qs>=1500.0] *= 0.10
qp = qs * 1.5

dgrid = 1000
zgrid = 1000

x = np.arange(X_min,X_max+dgrid,dgrid)
y = np.arange(Y_min,Y_max+dgrid,dgrid)
z = -np.arange(Z_min,Z_max+zgrid,zgrid)

########## paraview file
print('writing netcdf file for paraview')
####Creating the netcdf file
from netCDF4 import Dataset
rootgrp = Dataset("output_files/cvms1000_paraview.nc", "w", format="NETCDF4")

rootgrp.createDimension("x", nx)
rootgrp.createDimension("y", ny)
rootgrp.createDimension("z", nz)


vx = rootgrp.createVariable("x","f4",("x",))
vx[:]=x
vy = rootgrp.createVariable("y","f4",("y",))
vy[:]=y
vz = rootgrp.createVariable("z","f4",("z",))
vz[:]=z
vvp = rootgrp.createVariable("vp","f4",("z","y","x"))
vvp[:,:,:]=VpVsRho[:,:,:,0]
vvs = rootgrp.createVariable("vs","f4",("z","y","x"))
vvs[:,:,:]=VpVsRho[:,:,:,1]
vrho = rootgrp.createVariable("rho","f4",("z","y","x"))
vrho[:,:,:]=VpVsRho[:,:,:,2]
vqp = rootgrp.createVariable("Qp","f4",("z","y","x"))
vqp[:,:,:]=qp
vqs = rootgrp.createVariable("Qs","f4",("z","y","x"))
vqs[:,:,:]=qs

rootgrp.close()


########## rho col average
print('writing rho_col_average.nc netcdf')
####Creating the netcdf file
from netCDF4 import Dataset
rootgrp = Dataset("output_files/rho_col_average_cvms1000.nc", "w", format="NETCDF4")

rootgrp.createDimension("x", nx)
rootgrp.createDimension("y", ny)
rootgrp.createDimension("z", nz)


vx = rootgrp.createVariable("x","f4",("x",))
vx[:]=x
vy = rootgrp.createVariable("y","f4",("y",))
vy[:]=y
vz = rootgrp.createVariable("z","f4",("z",))
vz[:]=z
vrho_col_average= rootgrp.createVariable("rho_col_average","f4",("z","y","x"))
vrho_col_average[:,:,:]=rho_col_average[:,:,:]
rootgrp.close()

########## rho col average
print('writing PlastCo.nc netcdf')
####Creating the netcdf file
from netCDF4 import Dataset
rootgrp = Dataset("output_files/plastCo_cvms1000.nc", "w", format="NETCDF4")

rootgrp.createDimension("x", nx)
rootgrp.createDimension("y", ny)
rootgrp.createDimension("z", nz)


vx = rootgrp.createVariable("x","f4",("x",))
vx[:]=x
vy = rootgrp.createVariable("y","f4",("y",))
vy[:]=y
vz = rootgrp.createVariable("z","f4",("z",))
vz[:]=z
vplastCo =  rootgrp.createVariable("plastCo","f4",("z","y","x"))
vplastCo[:,:,:] = 0.0001*mu[:,:,:]
rootgrp.close()


########## ASAGI file
print('writing netcdf file for SeisSol')
rootgrp = Dataset("output_files/RhoMuLambdaQpQs_cvms1000.nc", "w", format="NETCDF4")

rootgrp.createDimension("x", nx)
rootgrp.createDimension("y", ny)
rootgrp.createDimension("z", nz)


vx = rootgrp.createVariable("x","f4",("x",))
vx[:]=x
vy = rootgrp.createVariable("y","f4",("y",))
vy[:]=y
vz = rootgrp.createVariable("z","f4",("z",))
vz[:]=z

mattype4 = np.dtype([('rho','f4'),('mu','f4'),('lambda','f4'),('Qp','f4'),('Qs','f4')])
mattype8 = np.dtype([('rho','f8'),('mu','f8'),('lambda','f8'),('Qp','f8'),('Qs','f8')])
mat_t = rootgrp.createCompoundType(mattype4,'material')

#this transform the 4 D array into an array of tuples
arr = np.stack((rho,mu,lam,qp,qs), axis=3)
newarr = arr.view(dtype=mattype8)
newarr = newarr.reshape(newarr.shape[:-1])
mat = rootgrp.createVariable("data",mat_t,("z","y","x"))
mat[:] = newarr
rootgrp.close()
