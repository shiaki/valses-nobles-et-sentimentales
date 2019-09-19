#!/usr/bin/env python

'''
    Generates initial condition from EJ
'''

import numpy as np
from scipy.optimize import brentq
import matplotlib.pyplot as plt

def makeinit(Vsq, Rsq, Qsq, Omega_p, N_pts, EJ_t):

    # check if params are valid.
    Wsq = Omega_p ** 2
    if Vsq / Wsq - Rsq < 0: # complain if no L_12 saddle points.
        raise ValueError('Model has no L_12 point.')

    # calculate radii and energies of L12, L45 points.
    xsq_L12, ysq_L45 = Vsq / Wsq - Rsq, Vsq / Wsq - Rsq * Qsq
    EJ_L12 = Vsq * np.log(Rsq + xsq_L12      ) / 2. - Wsq * xsq_L12 / 2.
    EJ_L45 = Vsq * np.log(Rsq + ysq_L45 / Qsq) / 2. - Wsq * ysq_L45 / 2.

    # central EJ
    EJ_c = Vsq * np.log(Rsq) / 2.

    # check energies.
    if (EJ_t > EJ_L12) or (EJ_t < EJ_c):
        raise ValueError('Requested EJ outside range.')

    # calculate range of domain,
    xsq_lim = brentq(lambda xsq_t: Vsq * np.log(Rsq + xsq_t      ) / 2. \
            - Wsq * xsq_t / 2. - EJ_t, 0., xsq_L12)
    ysq_lim = brentq(lambda ysq_t: Vsq * np.log(Rsq + ysq_t / Qsq) / 2. \
            - Wsq * ysq_t / 2. - EJ_t, 0., ysq_L45)
    xlim, ylim = np.sqrt(xsq_lim), np.sqrt(ysq_lim)

    # generate random points,
    ic_pts = list()
    while len(ic_pts) < N_pts:

        # draw random point,
        x_i = (2. * np.random.rand() - 1.) * xlim
        y_i = (2. * np.random.rand() - 1.) * ylim

        # positive velocity --> valid point,
        EJ_i = Vsq * np.log(Rsq + x_i ** 2 + y_i ** 2 / Qsq) / 2. \
             - (x_i ** 2 + y_i ** 2) * Wsq / 2.
        if EJ_i > EJ_t: continue # otherwise, skip

        # find velocity, assign isotropic direction
        vel_i = np.sqrt(2. * (EJ_t - EJ_i))
        w_i = np.random.randn(2)
        u_i, v_i = vel_i * w_i / np.linalg.norm(w_i)

        # put into list,
        ic_pts.append((x_i, y_i, u_i, v_i,))

    # convert to numpy
    return np.array(ic_pts)

if __name__ == '__main__' and None:

    # Model parameters.
    A, B, W = 1., 2., 0.5 # slow bar,
    ic_radius = 10. # maximal radius for init condition.
    N_pts = 1000
    EJ_t = 2.

    # try to make initial condition,
    ic_pts = makeinit(A, B, W, ic_radius, N_pts, EJ_t)

if __name__ == '__main__' and None:

    Vsq, Rsq, Qsq = 1., 0.1 ** 2, 0.82 ** 2
    W = 1.

    # draw EJ
    ax = np.linspace(0, 5, 501)
    E_x  = Vsq * np.log(Rsq + ax ** 2      ) / 2 - (ax ** 2) * W ** 2 / 2
    E_y  = Vsq * np.log(Rsq + ax ** 2 / Qsq) / 2 - (ax ** 2) * W ** 2 / 2
    plt.plot(ax, E_x); plt.plot(ax, E_y); plt.show()

    c = makeinit(Vsq, Rsq, Qsq, W, 30000, -1.)
    plt.scatter(c[:, 0], c[:, 1], c=c[:, 2] ** 2 + c[:, 3] ** 2, s=1), plt.show()