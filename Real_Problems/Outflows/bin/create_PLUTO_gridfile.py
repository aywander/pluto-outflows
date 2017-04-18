# 
# This script creates a grid.out file for starting
# from external files with PLUTO 4.0.0.
#

import numpy as np
import datetime

#__________________________ 
# Set these 

# The output file
fname = 'grid_in.out'

# All the data needed to write. Ghost cells are ignored, i think.
dim, geom = 2, 'CARTESIAN'
x1_beg, x1_end, nx1, ngh1 = 0, 1.5, 64, 4
x2_beg, x2_end, nx2, ngh2 = -0.5, 0.5, 128, 4
x3_beg, x3_end, nx3, ngh3 = -0.1875, 0.1875, 16, 4



#__________________________ 
# Set automatically

# Date
now = datetime.datetime.now()
date = now.strftime('%a %b %d %H:%M:%S %Y')

# Total number of cells. Not used. Ghost cells are ignored.
nx1_tot = nx1 + 2 * ngh1
nx2_tot = nx2 + 2 * ngh2
nx3_tot = nx3 + 2 * ngh3

# Cell widths
dx1 = (x1_end - x1_beg) / nx1
dx2 = (x2_end - x2_beg) / nx1
dx3 = (x3_end - x3_beg) / nx1

# Cell edge arrays
x1_edges = np.linspace(x1_beg, x1_end, nx1 + 1)
x2_edges = np.linspace(x2_beg, x2_end, nx2 + 1)
x3_edges = np.linspace(x3_beg, x3_end, nx3 + 1)

# Cell numbers
x1_num = np.arange(1, nx1 + 1)
x2_num = np.arange(1, nx2 + 1)
x3_num = np.arange(1, nx3 + 1)

# Combined array for output
x1_arr = np.column_stack((x1_num, x1_edges[:-1], x1_edges[1:]))
x2_arr = np.column_stack((x2_num, x2_edges[:-1], x2_edges[1:]))
x3_arr = np.column_stack((x3_num, x3_edges[:-1], x3_edges[1:]))
fmt = ' %4d %16.8e %16.8e'

#__________________________ 
# Outuput bit

fh = open(fname, 'w')

# The header
outstr = [
'# ******************************************************\n',
'# PLUTO 4.0 Grid File\n',
'# Generated on '+date+'\n',
'#\n',
'# DIMENSIONS: '+str(dim)+'\n',
'# GEOMETRY:   '+geom+'\n',
'# X1: ['+format(x1_beg, '>7')+', '+format(x1_end, '>7')+'], '+format(nx1, '>5')+' point(s), '+format(ngh1, '>3')+' ghosts'+'\n'
'# X2: ['+format(x2_beg, '>7')+', '+format(x2_end, '>7')+'], '+format(nx2, '>5')+' point(s), '+format(ngh2, '>3')+' ghosts'+'\n'
'# X3: ['+format(x3_beg, '>7')+', '+format(x3_end, '>7')+'], '+format(nx3, '>5')+' point(s), '+format(ngh3, '>3')+' ghosts'+'\n'
'# ******************************************************\n',
]
fh.writelines(outstr)
fh.close()

fh = open(fname, 'ab')

# Write and save it all
np.savetxt(fh, np.atleast_1d(nx1), '%d')
np.savetxt(fh, x1_arr, fmt)
np.savetxt(fh, np.atleast_1d(nx2), '%d')
np.savetxt(fh, x2_arr, fmt)
np.savetxt(fh, np.atleast_1d(nx3), '%d')
np.savetxt(fh, x3_arr, fmt)

fh.close()


