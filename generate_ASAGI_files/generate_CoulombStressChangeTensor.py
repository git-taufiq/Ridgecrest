import utm
import numpy as np
a = np.loadtxt('input_files/Tensor_PreRidgecrest.csv')

nx = 173
ny = 173
nz = 21

c_ij = np.zeros((nz,ny,nx,6))
coords = np.zeros((nx,ny,2))

idx=0
for k in range(0,nz):
   for i in range(0,nx):
      for j in range(0,ny):
         c_ij[k,j,i,:] = a[idx,3:9]
         idx=idx+1

print('print range (MPa)')
for i in range(0,6):
   print(np.amin(c_ij[:,:,:,i])/1e6, np.amax(c_ij[:,:,:,i])/1e6)

tapper=3e6
print('tapper CFS to %fPa' %(tapper))
c_ij[c_ij>tapper]=tapper
c_ij[c_ij<-tapper]=-tapper

idx=0
for i in range(0,nx):
   for j in range(0,ny):
      coords[i,j,:] = a[idx,0:2]
      idx=idx+1

#Define the structured grid
#remove 2 values on both side to avoid bounds_error
x = coords[2:-2,ny//2,0]
y = coords[nx//2,2:-2,1]
Xg,Yg = np.meshgrid(x,y)
#project grid to UTM
glat_grid, glon_grid = utm.to_latlon(Xg.flatten(), Yg.flatten(), 11, 'S')
coords_grid = np.column_stack((glon_grid,glat_grid))

#get the lons and lats array (data are structured in lat/long)
glat_grid, glon_grid = utm.to_latlon(coords[:,:,0], coords[:,:,1], 11, 'S')
lats = glat_grid[0,:]
lons = glon_grid[:,0]

c_ij_grid = np.zeros((c_ij.shape[0], c_ij.shape[1]-4, c_ij.shape[2]-4, c_ij.shape[3]))

from scipy.interpolate import RegularGridInterpolator
for k in range(c_ij.shape[0]):
   print(k)
   for ixx in range(6):
      fc_ij = RegularGridInterpolator((lons, lats), c_ij[k,:,:,ixx].T, fill_value=np.nan, bounds_error=True)
      out = fc_ij(coords_grid)
      #reverse order in the k dimension
      c_ij_grid[c_ij.shape[0] - (k+1),:,:,ixx] = out.reshape(Xg.shape)

#extend for grids higher than depth 0 km
c_ij_grid_ = np.zeros((c_ij.shape[0]+5, c_ij.shape[1]-4, c_ij.shape[2]-4, c_ij.shape[3]))
c_ij_grid_[:-5,:,:,:] = c_ij_grid * 1
c_ij_grid_zero = c_ij_grid[-2,:,:,:] #take interpolated data at depth 2 km
for g in range(6):
   c_ij_grid_[-g-1,:,:,:] = c_ij_grid_zero

z = np.linspace(-20e3,0+5e3,21+5)

#the initial file is in NEU, and this needs to be ported to Seissol ENU
# 0   1  2  3  4  5
# NN EE UU NE NU EU
# yy xx zz yx yz xz
# N=y, E=x, U=z
# 1   0  2  3  5  4
#we need to switch columns 0 and 1 and 4 and 5
c_ij_grid_[:,:,:,[0, 1]] = c_ij_grid_[:,:,:,[1, 0]]
c_ij_grid_[:,:,:,[4, 5]] = c_ij_grid_[:,:,:,[5, 4]]

########## paraview file
print('writing netcdf file for paraview')
####Creating the netcdf file
from netCDF4 import Dataset
rootgrp = Dataset("output_files/c_ij_paraview.nc", "w", format="NETCDF4")

rootgrp.createDimension("x", nx-4)
rootgrp.createDimension("y", ny-4)
rootgrp.createDimension("z", nz+5)

vx = rootgrp.createVariable("x","f4",("x",))
vx[:]=x
vy = rootgrp.createVariable("y","f4",("y",))
vy[:]=y
vz = rootgrp.createVariable("z","f4",("z",))
vz[:]=z
vxx = rootgrp.createVariable("c_xx","f4",("z","y","x"))
vxx[:,:,:]=c_ij_grid_[:,:,:,0]
vyy = rootgrp.createVariable("c_yy","f4",("z","y","x"))
vyy[:,:,:]=c_ij_grid_[:,:,:,1]
vzz = rootgrp.createVariable("c_zz","f4",("z","y","x"))
vzz[:,:,:]=c_ij_grid_[:,:,:,2]
vxy = rootgrp.createVariable("c_xy","f4",("z","y","x"))
vxy[:,:,:]=c_ij_grid_[:,:,:,3]
vxz = rootgrp.createVariable("c_xz","f4",("z","y","x"))
vxz[:,:,:]=c_ij_grid_[:,:,:,4]
vyz = rootgrp.createVariable("c_yz","f4",("z","y","x"))
vyz[:,:,:]=c_ij_grid_[:,:,:,5]

rootgrp.close()


########## ASAGI file
print('writing netcdf file for SeisSol')
rootgrp = Dataset("output_files/c_ij.nc", "w", format="NETCDF4")

rootgrp.createDimension("x", nx-4)
rootgrp.createDimension("y", ny-4)
rootgrp.createDimension("z", nz+5)

vx = rootgrp.createVariable("x","f4",("x",))
vx[:]=x
vy = rootgrp.createVariable("y","f4",("y",))
vy[:]=y
vz = rootgrp.createVariable("z","f4",("z",))
vz[:]=z

mattype4 = np.dtype([('c_xx','f4'),('c_yy','f4'),('c_zz','f4'), ('c_xy','f4'), ('c_xz','f4'), ('c_yz','f4')])
mattype8 = np.dtype([('c_xx','f8'),('c_yy','f8'),('c_zz','f8'), ('c_xy','f8'), ('c_xz','f8'), ('c_yz','f8')])
mat_t = rootgrp.createCompoundType(mattype4,'c_ij_mat')
#this transform the 4 D array into an array of tuples
newarr = c_ij_grid_.view(dtype=mattype8)
newarr = newarr.reshape(newarr.shape[:-1])
mat = rootgrp.createVariable("c_ij",mat_t,("z","y","x"))
mat[:] = newarr
rootgrp.close()
