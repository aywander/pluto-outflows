# Just a function that computes the number of cells needed
# in the second dimension in a 2D spherical-polar setup
# to make the grid cell shape roughly equilateral.
# Input are xl (left boundary), xr (right boundary), and
# nr (number of radial cells).

import numpy as np
import argparse

def grid_log_equilateral(xl, xr, nx, theta_range=np.pi):
    '''
    Computes the number of cells needed in the second dimension in a 2D spherical-polar setup
    to make the grid cell shape roughly equilateral.

    :param xl:           Left boundary
    :param xr:           Right boundary
    :param nx:           Number of radial cells
    :param theta_range:  Range in theta (in radians, default pi)
    :return:             Number of uniform cells in theta direction across theta_range
    '''

    dx = xl * ((xr / xl) ** (1. / nx) - 1.)
    nt = theta_range / np.arctan(dx / xl)

    return nt


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='computes the number of cells needed in the second '
                    'dimension in a 2D spherical-polar setup to make '
                    'the grid cell shape roughly equilateral.')

    parser.add_argument('xl', type=np.float, help='Left Boundary')
    parser.add_argument('xr', type=np.float, help='Right Boundary')
    parser.add_argument('nx', type=np.int, help='Number of radial cells')
    parser.add_argument('th', type=np.float, nargs='?', help='Range in theta (in radians, default pi)', default=np.pi)

    args = parser.parse_args()

    nt = grid_log_equilateral(args.xl, args.xr, args.nx, args.th)

    print(nt)


