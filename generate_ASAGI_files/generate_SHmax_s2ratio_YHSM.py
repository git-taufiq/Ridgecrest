import numpy as np
from scipy import ndimage as nd

def fill(data, invalid=None):
    #replace nan by nearest neighbor
    invalid = np.isnan(data)
    ind = nd.distance_transform_edt(invalid, return_distances=False, return_indices=True)
    return data[tuple(ind)]

"""
# here is the format of the file
# LON LAT DEP 1 SHmax_trend SHmax_mag SHmin_mag SV_mag  See Sen Seu Snn Snu Suu                  RATIO DEV  ISO
# LON LAT DEP 2  1-sigma     1-sigma   1-sigma  1-sigma DOT dot-1sig ANG s1-1sig s2-1sig s3-1sig 1sig  1sig 1sig
#
# metedata with Yang and Huaksson (2012) 2D stress inversion:
# LON LAT NaN 1 SHmax_trend NaN NaN NaN  See Sen Seu Snn Snu Suu   1 0  0
# LON LAT NaN 2  1-sigma    NaN NaN NaN  0 NaN 0 NaN NaN NaN NaN  NaN NaN
"""

stress = np.loadtxt('input_files/CSM_YHSM-2013_Yang_Hauksson_2012.txt')
stress = stress[0::2] # take only odd rows
print(stress)
drawplots=False

# crop the data
stress_crop = stress[np.argwhere(stress[:,0]>-118.7),:][:,0,:]
stress_crop = stress_crop[np.argwhere(stress_crop[:,0]<-116.3),:][:,0,:]
stress_crop = stress_crop[np.argwhere(stress_crop[:,1]>34.7),:][:,0,:]
stress_crop = stress_crop[np.argwhere(stress_crop[:,1]<36.7),:][:,0,:]

# shift to have continuous values
stress_crop[:,4][stress_crop[:,4]>90.0] -= 180.0

ids = np.argwhere(stress_crop[:,2]==1).flatten()
stress_crop = stress_crop[:,:][ids]
lonsdata = stress_crop[:,0]
latsdata = stress_crop[:,1]

# Compute s2ratio (we use col 5)
from numpy import linalg
ndata=len(lonsdata)
for i in range(ndata):
    A = np.array([[stress_crop[i, 8],stress_crop[i, 9],stress_crop[i,10]],
                  [stress_crop[i, 9],stress_crop[i,11],stress_crop[i,12]],
                  [stress_crop[i,10],stress_crop[i,12],stress_crop[i,13]]])
    w0, v = linalg.eig(A)
    #s2ratio = (s_2-s_3)/(s_1-s_3) where s_1>s_2>s_3
    w = np.sort(w0)[::-1]
    print(w)
    stress_crop[i,5] = (w[1]-w[2])/(w[0]-w[2])
    #Sv defines which of the principal stresses s_i is vertical
    idv = int(np.argmax(np.absolute(np.dot(v,np.array([0,0,1])))))
    stress_crop[i,6] = np.where(abs(w-w0[idv])<1e-8)[0]+1
print('done computing S2ratio')

#determine the lat lon grid
lons = np.unique(lonsdata)
lats = np.unique(latsdata)
X,Y = np.meshgrid(lons, lats)
nx=len(lons)
ny=len(lats)

#We use the same procedure for SHmax and s2ratio 
#for ll in [4,5]:
asname={4:'SHmax',5:'s2ratio', 6:'Sv', 8:'bxx', 9:'bxy',10:'bxz',11:'byy',12:'byz',13:'bzz'}

for ll in range(4,7):
   #Grid the SHmax data
   SHmaxdata = stress_crop[:,ll]
   SHmax = np.zeros_like(X)
   SHmax[:,:] = np.nan
   for kk in range(ndata):
      i = np.searchsorted(lons, lonsdata[kk])
      j = np.searchsorted(lats, latsdata[kk])
      SHmax[j,i] = SHmaxdata[kk]
   SHmax=fill(SHmax)

   if drawplots:
      #plot the 2d distribution
      import matplotlib.pyplot as plt
      plt.pcolor(X, Y, SHmax)
      plt.colorbar()
      plt.figure()
      smin = np.amin(SHmax)
      smax = np.amax(SHmax)
      SHmax2 = nd.gaussian_filter(SHmax,sigma=2)
      plt.pcolor(X, Y, SHmax2)
      plt.colorbar()
      plt.clim(smin,smax)
      plt.show()

   from scipy.interpolate import RegularGridInterpolator
   fSHmax = RegularGridInterpolator((lons, lats), SHmax.T, fill_value=np.nan, bounds_error=False)
   X_min = 352000
   X_max = 556000
   Y_min = 3850000
   Y_max = 4052000
   grid = 2000

   x=np.arange(X_min, X_max+grid, grid)
   y=np.arange(Y_min, Y_max+grid, grid)
   nxg=len(x)
   nyg=len(y)

   Xg,Yg = np.meshgrid(x,y)
   import utm
   glat_grid, glon_grid = utm.to_latlon(Xg.flatten(), Yg.flatten(), 11, 'S')

   coords_grid = np.column_stack((glon_grid,glat_grid))
   SHmaxg = fSHmax(coords_grid)
   SHmaxg = SHmaxg.reshape(np.shape(Xg))
   SHmaxg=fill(SHmaxg)
   sig=2
   print('applying gaussian filter, sigma=', sig)
   SHmaxg = nd.gaussian_filter(SHmaxg,sigma=sig)

   if drawplots:
      #plot the 2d distribution
      import matplotlib.pyplot as plt
      plt.pcolor(Xg, Yg, SHmaxg)
      plt.colorbar()
      plt.show()

   sname = asname[ll]
   ########## rho col average
   print('writing '+sname+'.nc netcdf')
   ####Creating the netcdf file
   from netCDF4 import Dataset
   rootgrp = Dataset('output_files/'+sname+'.nc', "w", format="NETCDF4")

   rootgrp.createDimension("x", nxg)
   rootgrp.createDimension("y", nyg)

   vx = rootgrp.createVariable("x","f4",("x",))
   vx[:]=x
   vy = rootgrp.createVariable("y","f4",("y",))
   vy[:]=y
   vSHmax =  rootgrp.createVariable(sname,"f4",("y","x"))
   vSHmax[:,:] = SHmaxg[:,:]
   if ll==6: 
      vSHmax[:,:] = np.floor(SHmaxg[:,:])
   rootgrp.close()
