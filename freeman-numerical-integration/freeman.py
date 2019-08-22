#!/usr/bin/env python

'''
'''

from tqdm import tqdm, trange
import numpy as np
import matplotlib.pyplot as plt

def rk4int(w_init, Asq, Bsq, Omega, dt, N_step):

    ''' Integrate orbits in Freeman's potential using RK4 integrator. '''

    w0 = np.reshape(w_init, (-1, 4)).T

    N_orb = w0.shape[1]
    orb = np.zeros((N_step, 4, N_orb))
    orb[0, :, :] = w0

    f = lambda x, y, u, v: \
            ( dt * u, \
              dt * v, \
             -dt * (Asq * x + Omega * v), \
             -dt * (Bsq * y - Omega * u))

    xi, yi, ui, vi = w0
    for i_step in trange(1, N_step):
        x1, y1, u1, v1 = f(xi,          yi,          ui,          vi         )
        x2, y2, u2, v2 = f(xi + x1 / 2, yi + y1 / 2, ui + u1 / 2, vi + v1 / 2)
        x3, y3, u3, v3 = f(xi + x2 / 2, yi + y2 / 2, ui + u2 / 2, vi + v2 / 2)
        x4, y4, u4, v4 = f(xi + x3,     yi + y3,     ui + u3,     vi + v3    )
        xi += (x1 + 2 * x2 + 2 * x3 + x4) / 6
        yi += (y1 + 2 * y2 + 2 * y3 + y4) / 6
        ui += (u1 + 2 * u2 + 2 * u3 + u4) / 6
        vi += (v1 + 2 * v2 + 2 * v3 + v4) / 6
        orb[i_step, 0, :] = xi
        orb[i_step, 1, :] = yi
        orb[i_step, 2, :] = ui
        orb[i_step, 3, :] = vi

    return np.swapaxes(np.swapaxes(orb, 1, 2), 0, 1)

if __name__ == '__main__':

    # initial contidion
    w_init = [0., 3., 3., 0]

    orb = rk4int(w_init, 4., 1.2, 1., 1.e-3, 100000)
    plt.plot(orb[0, :, 0], orb[0, :, 1]), plt.show()