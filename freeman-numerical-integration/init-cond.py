#!/usr/bin/env python 

'''
    Generates initial condition from EJ
'''

import numpy as np
import matplotlib.pyplot as plt

# Model parameters.
A, B, W = 1., 2., 1. # slow bar,
ic_radius = 10. # maximal radius for init condition.
N_pts = 1000
EJ_t = 2.

if __name__ == '__main__':

    # calculate EJ at edge of domain,
    EJ_a = (A * ic_radius) ** 2 / 2 - (W * ic_radius) ** 2 / 2
    EJ_b = (B * ic_radius) ** 2 / 2 - (W * ic_radius) ** 2 / 2
    assert (EJ_t > EJ_a) or (EJ_t > EJ_b) or (EJ_t > 0.)
    # EJ at (0, 0) is always 0.

    # generate random points,
    ic_pts = list()
    while len(ic_pts) < N_pts:

        # draw random point,
        r_i = np.sqrt(np.random.rand()) * ic_radius
        w_i = np.random.randn(2)
        x_i, y_i = r_i * w_i / np.linalg.norm(w_i)

        # positive velocity --> valid point,
        EJ_i = ((A * x_i) ** 2 + (B * y_i) ** 2) / 2 \
             - (x_i ** 2 + y_i ** 2) * W ** 2 / 2
        if EJ_i > EJ_t: continue # otherwise, skip

        # find velocity, assign direction
        vel_i = np.sqrt(2. * (EJ_t - EJ_i))
        w_i = np.random.randn(2)
        u_i, v_i = vel_i * w_i / np.linalg.norm(w_i)

        # put into list,
        ic_pts.append((x_i, y_i, u_i, v_i))



if (__name__ == '__main__') and None:

    # draw EJ
    ax = np.linspace(-5, 5, 501)
    xx, yy = np.meshgrid(ax, ax)

    EJ_x  = (A * ax) ** 2 / 2 - (ax ** 2) * W ** 2 / 2
    EJ_y  = (B * ax) ** 2 / 2 - (ax ** 2) * W ** 2 / 2
    EJ_xy = ((A * xx) ** 2 + (B * yy) ** 2) / 2 \
          - (xx ** 2 + yy ** 2) * W ** 2 / 2

    plt.plot(ax, EJ_x)
    plt.plot(ax, EJ_y)
    plt.show()
    # plt.imshow(EJ_xy), plt.colorbar(), plt.show()