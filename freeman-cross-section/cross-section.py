#!/usr/bin/env python

'''
    Make cross section.
'''

import numpy as np
import matplotlib.pyplot as plt

# our own functions.
from makeinit import make_init
from freeman import rk4int

if __name__ == '__main__':

    # set model parameters.
    A, B, W = 1., 2., 0.5 # slow bar,
    ic_radius = 10. # maximal radius for init condition.
    N_pts = 10000
    EJ_t = 2.

    # generate initial condition.
    ic_pts = make_init(A, B, W, ic_radius, N_pts, EJ_t)

    # integrate to t=100, 8001 steps.
    orb = rk4int(ic_pts, A * A, B * B, W, 0.0125, 8001)