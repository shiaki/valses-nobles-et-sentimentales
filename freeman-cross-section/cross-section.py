#!/usr/bin/env python

'''
    Make cross section.
'''

import numpy as np
import matplotlib.pyplot as plt

# our own functions.
from makeinit import makeinit
from freeman import rk4int
from makesos import sos

if __name__ == '__main__':

    # set model parameters.
    A, B, W = 1., 2., 0.25 # slow bar,
    ic_radius = 10. # maximal radius for init condition.
    N_pts = 8
    EJ_t = 3.

    # generate initial condition.
    ic_pts = makeinit(A, B, W, ic_radius, N_pts, EJ_t)

    # integrate to t=100, 8001 steps.
    orb = rk4int(ic_pts, A * A, B * B, W, 0.05, 64001)

    for orb_i in orb:
        k = sos(orb_i, axis=0, direction='+')
        plt.scatter(k[:, 1], k[:, 3], s=2)
    plt.show()