#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

from freeman_solution import FreemanLambdified

if __name__ == '__main__':

    A, B, W = 0.99, 2.01, .5 # slow bar,
    x0, y0, u0, v0 = 0., -1., 0., -4. # from y-axis,

    t = np.linspace(0., 800., 8001)
    # t = np.sort(np.random.rand(8001)) * 2000.

    freeman = FreemanLambdified('simplified.pk')
    xt, yt, ut, vt = freeman(A, B, W, x0, y0, u0, v0, t,)
    xp, yp, up, vp = freeman.Wp(A, B, W, x0, y0, u0, v0, t,)
    xq, yq, uq, vq = freeman.Wq(A, B, W, x0, y0, u0, v0, t,)
    vrg = max(max(np.abs(xt)), max(np.abs(yt))) * 1.2

    fig = plt.figure(figsize=(6., 6.))
    ax = fig.add_subplot(1, 1, 1)
    ax.scatter(xt, yt, s=1, marker='.', alpha=0.5)
    ax.scatter(xp, yp, s=1, marker='.', label=r'$\lambda_p$')
    ax.scatter(xq, yq, s=1, marker='.', label=r'$\lambda_q$')

    ax.legend(markerscale=5)
    ax.set_xlim(-vrg, vrg), ax.set_ylim(-vrg, vrg)

    plt.show()