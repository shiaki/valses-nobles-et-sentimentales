#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

from readorbits import ReadOrbits

datapath = '/d0/Shanghai/galex2/dgalex2/qinyj/data/3042_orbsave/'

if __name__ == '__main__':

    # read the 27th file.
    orbits = ReadOrbits(27, rawdata_path=datapath)

    # and get the 2000th orbit from this file. (index from 0 to 99999)
    pos, vel, time = orbits.get_orbit(2000, velocity_method='full', packed=False)
    # here pos and vel are (5001, 3) arrays, and time is a (5001,) array

    # or alternatively,
    orb_t = orbits.get_orbit(2000, velocity_method='full', packed=True)
    # here orb_t is a (5001, 7) array, last dim.: [x, y, z, vx, vy, vz, t]

    # NOTE: the orbit is in inertia frame!!

    plt.plot(pos[:, 0], pos[:, 1])
    plt.show()