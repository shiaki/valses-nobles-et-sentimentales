#!/usr/bin/env python

'''
    Construct SoS using integrated orbits.
'''

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import brentq

def sos(orb, axis=0):

    ''' Construct SoS for integrated orbits. '''

    # interpolator,
    k_intp = interp1d(np.arange(orb.shape[0]), orb[:, axis], \
                      copy=True, assume_sorted=True)

    # find zeros,
    i_zeros = np.arange(np.arange(orb.shape[0]) - 1)[ \
            orb[:-1, axis] * orb[:-1, axis] < 0.]
    t_zeros = np.array([brentq(k_intp, zi, zi + 1) for zi in i_zeros])

    # interpolate other axes,
    w_intp = interp1d(np.arange(orb.shape[0]), orb, \
                      axis=0, copy=False, assume_sorted=True)

    # return surface if sectuibm
    return w_intp(t_zeros)