## This script generates a grid{xyz}.out from a (3D) grid.out file
import numpy as np

fh = open('grid.out','r')
pos0 = fh.tell()
nl = 0; 
while fh.readline(): nl+=1

fh.seek(pos0); 
nx = int(fh.readline())
gx = np.genfromtxt(fh, skip_footer=nl-nx-1,     usecols=2)

fh.seek(pos0); 
for i in range(1+nx): fh.readline();
ny = int(fh.readline())
gy = np.genfromtxt(fh, skip_footer=nl-nx-ny-2,  usecols=2)

fh.seek(pos0); 
for i in range(2+nx+ny): fh.readline();
nz = int(fh.readline())
gz = np.genfromtxt(fh, usecols=2)
fh.close()

gx.tofile('gridx.out')
gy.tofile('gridy.out')
gz.tofile('gridz.out')



