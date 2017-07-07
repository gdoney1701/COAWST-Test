
# coding: utf-8

# # Mplleaflet

# In[46]:

get_ipython().magic('matplotlib inline')
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
import math


# In[47]:

tidx = -1       # just get the final frame, for now.
scale = .05 ## This variable can be played with to determine arrow size
isub = 20
url = 'http://geoport-dev.whoi.edu/thredds/dodsC/coawst_4/use/fmrc/coawst_4_use_best.ncd'


# In[48]:

def shrink(a,b):
    """Return array shrunk to fit a specified shape by triming or averaging.
    
    a = shrink(array, shape)
    
    array is an numpy ndarray, and shape is a tuple (e.g., from
    array.shape). a is the input array shrunk such that its maximum
    dimensions are given by shape. If shape has more dimensions than
    array, the last dimensions of shape are fit.
    
    as, bs = shrink(a, b)
    
    If the second argument is also an array, both a and b are shrunk to
    the dimensions of each other. The input arrays must have the same
    number of dimensions, and the resulting arrays will have the same
    shape.
    Example
    -------
    
    >>> shrink(rand(10, 10), (5, 9, 18)).shape
    (9, 10)
    >>> map(shape, shrink(rand(10, 10, 10), rand(5, 9, 18)))        
    [(5, 9, 10), (5, 9, 10)]   
       
    """

    if isinstance(b, np.ndarray):
        if not len(a.shape) == len(b.shape):
            raise Exception(
                  'input arrays must have the same number of dimensions')
        a = shrink(a,b.shape)
        b = shrink(b,a.shape)
        return (a, b)

    if isinstance(b, int):
        b = (b,)

    if len(a.shape) == 1:                # 1D array is a special case
        dim = b[-1]
        while a.shape[0] > dim:          # only shrink a
            if (dim - a.shape[0]) >= 2:  # trim off edges evenly
                a = a[1:-1]
            else:                        # or average adjacent cells
                a = 0.5*(a[1:] + a[:-1])
    else:
        for dim_idx in range(-(len(a.shape)),0):
            dim = b[dim_idx]
            a = a.swapaxes(0,dim_idx)        # put working dim first
            while a.shape[0] > dim:          # only shrink a
                if (a.shape[0] - dim) >= 2:  # trim off edges evenly
                    a = a[1:-1,:]
                if (a.shape[0] - dim) == 1:  # or average adjacent cells
                    a = 0.5*(a[1:,:] + a[:-1,:])
            a = a.swapaxes(0,dim_idx)        # swap working dim back

    return a


# In[49]:

def rot2d(x, y, ang):
    '''rotate vectors by geometric angle'''
    xr = x*np.cos(ang) - y*np.sin(ang)
    yr = x*np.sin(ang) + y*np.cos(ang)
    return xr, yr


# In[ ]:




# In[50]:

nc = netCDF4.Dataset(url)
mask = nc.variables['mask_rho'][:]
lon_rho = nc.variables['lon_rho'][:]
lat_rho = nc.variables['lat_rho'][:]
anglev = nc.variables['angle'][:]

u = nc.variables['u'][tidx, -1, :, :]
v = nc.variables['v'][tidx, -1, :, :]

u = shrink(u, mask[1:-1, 1:-1].shape)
v = shrink(v, mask[1:-1, 1:-1].shape)

u, v = rot2d(u, v, anglev[1:-1, 1:-1])

#u = u[np.logical_not(np.isnan(u))]
#v = v[np.logical_not(np.isnan(v))]
#print("u =",u)
#print("v =",v)


# In[51]:

lon_c = lon_rho[1:-1, 1:-1]
lat_c = lat_rho[1:-1, 1:-1]


# In[7]:

import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature, COLORS
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


# In[62]:

LAND = NaturalEarthFeature('physical', 'land', '110m', edgecolor='face',
                           facecolor=COLORS['land'])


# In[63]:

fig, ax = plt.subplots(figsize=(12,12),
                       subplot_kw=dict(projection=ccrs.PlateCarree()))

ax.set_extent([lon_c.min(), lon_c.max(), lat_c.min(), lat_c.max()])
ax.add_feature(LAND)
ax.coastlines()
gl = ax.gridlines(draw_labels=True)
gl.xlabels_top = gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

#kw = dict(scale=20, headwidth=2)
#n = 5
#q = ax.quiver(lon.points[::n, ::n], lat.points[::n, ::n],
#              u[::n, ::n], v[::n, ::n], color='black', **kw)
legend_vel=1.0
Q = ax.quiver( lon_c[::isub,::isub], lat_c[::isub,::isub], u[::isub,::isub], v[::isub,::isub], 
        scale=1.0/scale, pivot='middle', zorder=1e35, width=0.003)
legend_str='{} m/s'.format(legend_vel)
qk = ax.quiverkey(Q, 0.92, 0.88, legend_vel, legend_str)


# In[52]:

import mplleaflet


# In[53]:

scale


# Flatten the 2D arrays to 1D

# In[54]:

lon1d = lon_c[::isub,::isub].flatten()
lat1d = lat_c[::isub,::isub].flatten()


# In[55]:

u1d = u[::isub,::isub].flatten()
v1d = v[::isub,::isub].flatten()


# Find only the finite (non NaN) vectors 

# In[56]:

water = np.isfinite(u1d)


# In[ ]:

get_ipython().set_next_input('It seems that mplleaflet chokes on too many vectors.  How many do we have');get_ipython().magic('pinfo have')


# In[64]:

print(len(u1d[water]))


# In[61]:

fig, ax = plt.subplots(figsize=(10,10))
Q = ax.quiver(lon1d[water], lat1d[water], u1d[water], v1d[water],
        scale=1.0/scale, pivot='middle', zorder=1e35, width=0.003, color='white')
mplleaflet.display(fig=fig, tiles='esri_aerial')


# In[ ]:




# In[ ]:



