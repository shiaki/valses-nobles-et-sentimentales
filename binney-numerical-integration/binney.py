#!/usr/bin/env python

'''
'''

from tqdm import tqdm, trange
import numpy as np
import matplotlib.pyplot as plt

def rk4int(w_init, Vsq, Rsq, Qsq, Omega_p, dt, N_step, desc=''):

    ''' Integrate orbits in Freeman's potential using RK4 integrator. '''

    w0 = np.reshape(w_init, (-1, 4)).T

    N_orb = w0.shape[1]
    orb = np.zeros((N_step, 4, N_orb))
    orb[0, :, :] = w0

    def f(x, y, u, v, dt, Vsq, Rsq, Qsq, Op2, Opsq):
        S = -Vsq / (Rsq + x ** 2 + y ** 2 / Qsq)
        return dt * u, \
               dt * v, \
               dt * (S * x       + Op2 * v + Opsq * x), \
               dt * (S * y / Qsq - Op2 * u + Opsq * y)

    xi, yi, ui, vi = w0
    par = (dt, Vsq, Rsq, Qsq, 2. * Omega_p, Omega_p ** 2)
    for i_step in trange(1, N_step, desc=desc):
        x1, y1, u1, v1 = f(xi,          yi,          ui,          vi         , *par)
        x2, y2, u2, v2 = f(xi + x1 / 2, yi + y1 / 2, ui + u1 / 2, vi + v1 / 2, *par)
        x3, y3, u3, v3 = f(xi + x2 / 2, yi + y2 / 2, ui + u2 / 2, vi + v2 / 2, *par)
        x4, y4, u4, v4 = f(xi + x3,     yi + y3,     ui + u3,     vi + v3    , *par)
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

    V, R, Q, Omega_p = 1.0, 0.1, 0.82, 1. # slow bar,
    x0, y0, u0, v0 = 0., 0.25, -1.25, 0. # from y-axis,

    # initial contidion
    w_init = [x0, y0, u0, v0]

    orb = rk4int(w_init, V ** 2, R ** 2, Q ** 2, Omega_p, 1.e-2, int(200 / 1.e-2))
    plt.plot(orb[0, :, 0], orb[0, :, 1]), plt.show()