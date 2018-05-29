##!/uufs/chpc.utah.edu/sys/pkg/python/2.5.4/bin/python2.5

import numpy as np
import matplotlib.pyplot as plt
import array
import os

def myfwrite(fh, formatstring, ndarray):
    arr = array.array(formatstring, ndarray.flatten('F'))
    arr.tofile(fh)
    return 0
def myfread(fh, n):
    ndarray = np.zeros((n))
    arr = np.array(ndarray, dtype='d', order = 'F')
    arr = np.fromfile(file = fh, dtype = 'd', count = n)
    return arr

### Assuming homogeneous field ie uniform roughness  ###
z_o = 0.1

Nx = (512)
Ny = (512)
Nz = (64)


u_star = 0.4154 ## used to create the log profile
u_star_LESin = 0.45  ## specified in LESinputs.txt as a normalization parameter
l_z = 1500.00
# calculating parameters
zH_1 = np.arange(0.0, Nz)/(Nz-1); zH_1[0] = 0.5/(Nz-1);  #For dudz,dvdz,dtdz
zH_2 = np.arange(0.0, Nz)/(Nz-1);             #For w,uw,vw,wt
zH_3 = np.arange(0.5 , Nz)/(Nz-1.0);          #For u,v,T,Pr,Cs
print "zH_1.shape", zH_1.shape
print "zH_2.shape", zH_2.shape
print "zH_3.shape", zH_3.shape
z_1  = zH_1*l_z;
z_w  = zH_2*l_z;
z_u  = zH_3*l_z;

dz = l_z / (Nz - 1)
sscale = 1  # scalar scale from LESinputs.txt
### Allocate velocity fields  ###
u = np.zeros((Nx, Ny))
v = np.zeros((Nx, Ny))
w = np.zeros((Nx, Ny))
scalar = np.zeros((Nx, Ny))
# write to files

fwu = open('u.ini', 'wb')
fwv = open('v.ini', 'wb')
fww = open('w.ini', 'wb')
fws = open('scalar.ini', 'wb')
#f_uvp = open('zarray_uv_p.ini', 'wb')
#f_wpos = open('zarray_w.ini', 'wb')

for k in np.arange(Nz):
    uw = -(l_z - z_w[k]) / l_z * u_star ** 2
    print 'uw', uw
    if(k < Nz - 2) :
        ### Initial perturbation: random noise with an amplitude of 0.25*ustar
        ### on u,v and w
        varu = 0.25 * u_star 
        varv = 0.25 * u_star 
        varw = 0.25 * u_star 
        var_scalar = 0.125 * sscale
    else:
        varu = 0.0
        varv = 0.0
        varw = 0.0
        var_scalar = 0.0
    if(k == 0):
        varw = 0.0
    ur = np.random.randn(Nx, Ny)
    vr = np.random.randn(Nx, Ny)
    wr = np.random.randn(Nx, Ny)
    sr = np.random.randn(Nx, Ny)

    u[:, :] = (u_star/ 0.4 * np.log(z_u[k] / z_o) + ur * varu)/u_star_LESin
    v[:, :] = (vr * varv) / u_star_LESin
    w[:, :] = (wr* varw) / u_star_LESin
    scalar[:, :] = 1.00 + sr * var_scalar
    ### w = 0 on the lowest level due to staggered grid setup
    if(k == 0):
        w[:, :] = 0.0
    print 'mean(u)', np.mean(u)
    ### write each levels into files ###    
    myfwrite(fwu, 'd', u)
    myfwrite(fwv, 'd', v)
    myfwrite(fww, 'd', w)
    myfwrite(fws, 'd', scalar)
#  close files
fwu.close()
fwv.close()
fww.close()
fws.close()
zo=np.ones((Nx,Ny))*z_o;
fz = open('zo.ini','wb')
myfwrite(fz, 'd', zo)
fz.close()

fid = open('u.ini', 'rb')
utemp = np.zeros((Nz, 1))
for i in np.arange(Nz):
    read_data = np.fromfile(fid, dtype=np.float64, count = Nx * Ny)
    utemp[i] = np.mean(read_data) * u_star_LESin
print 'temp.shape', utemp.shape
print utemp
fid.close()
plt.figure()
plt.plot(utemp, z_u, '-.*r')
plt.xlabel('u (m/s)')
plt.ylabel('z(m)')
plt.title('Log profile')
plt.savefig('logProfile.pdf')
plt.show()

fid = open('v.ini', 'rb')
utemp = np.zeros((Nz, 1))
for i in np.arange(Nz):
    read_data = np.fromfile(fid, dtype='d', count = Nx * Ny)
    utemp[i] = np.mean(read_data) * u_star_LESin
fid.close()
plt.figure()
plt.plot(utemp, z_u, '-.*b')
plt.xlabel('v (m/s)')
plt.ylabel('z(m)')
plt.title('v profile')
plt.savefig('vProfile.pdf')
plt.show()

fid = open('w.ini', 'rb')
utemp = np.zeros((Nz, 1))
for i in np.arange(Nz):
    read_data = np.fromfile(fid, dtype=np.float64, count = Nx * Ny)
    utemp[i] = np.mean(read_data) * u_star_LESin
fid.close()
plt.figure()
plt.plot(utemp, z_w, '-.*g')
plt.xlabel('w (m/s)')
plt.ylabel('z(m)')
plt.title('w profile')
plt.savefig('wProfile.pdf')
plt.show()
