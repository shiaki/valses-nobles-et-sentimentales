#!/usr/bin/env python

'''
    Calculate orbits analytically
'''

import pickle
import numpy as np
from sympy import lambdify

class FreemanLambdified(object):

    def __init__(self, filename='simplified.pk'):
        
        # read simplified equations
        with open(filename, 'rb') as f: s = pickle.load(f)
        up = lambda k: pickle.loads(s[k])

        # define aux variables
        self._V = lambda A_, B_, W_: np.sqrt( \
                  np.sqrt(A_ * A_ + B_ * B_ + W_ * W_ + 2. * A_ * B_) \
                * np.sqrt(A_ * A_ + B_ * B_ + W_ * W_ - 2. * A_ * B_))
        self._P = lambda A_, B_: (A_ * A_ + B_ * B_) / 2.
        self._Q = lambda A_, B_: (A_ * A_ - B_ * B_) / 2.
        self._R = lambda V_, W_: (V_ * V_ + W_ * W_) / 2.
        self._S = lambda V_, W_: (V_ * V_ - W_ * W_) / 2.

        # read symbols
        a, b, c, d = (up(k) for k in 'abcd')
        A, B, P, Q, R, S, W, V, t = (up(k) for k in 'ABPQRSWVt')
        x0, y0, u0, v0 = (up(k + '0') for k in 'xyuv')

        # init condition functions
        vars_t = [A, B, P, Q, R, S, V, W, x0, y0, u0, v0]
        self._a = lambdify(vars_t, up('a_l'))
        self._b = lambdify(vars_t, up('b_l'))
        self._c = lambdify(vars_t, up('c_l'))
        self._d = lambdify(vars_t, up('d_l'))

        # orbit solution (simplified)
        vars_t = [A, B, P, Q, R, S, V, W, a, b, c, d, t]
        self._x = lambdify(vars_t, up('x_l'), modules='numpy')
        self._y = lambdify(vars_t, up('y_l'), modules='numpy')
        self._u = lambdify(vars_t, up('u_l'), modules='numpy')
        self._v = lambdify(vars_t, up('v_l'), modules='numpy')

        # orbit solution (rawexpr)
        vars_t = [A, B, W, x0, y0, u0, v0, t]
        self._xd = lambdify(vars_t, up('x_t'), modules='numpy')
        self._yd = lambdify(vars_t, up('y_t'), modules='numpy')
        self._ud = lambdify(vars_t, up('u_t'), modules='numpy')
        self._vd = lambdify(vars_t, up('v_t'), modules='numpy')

        # tw fundamental frequencies decomposed,
        self._xp = lambdify(vars_t, up('x_p'), modules='numpy')
        self._yp = lambdify(vars_t, up('y_p'), modules='numpy')
        self._up = lambdify(vars_t, up('u_p'), modules='numpy')
        self._vp = lambdify(vars_t, up('v_p'), modules='numpy')

        self._xq = lambdify(vars_t, up('x_q'), modules='numpy')
        self._yq = lambdify(vars_t, up('y_q'), modules='numpy')
        self._uq = lambdify(vars_t, up('u_q'), modules='numpy')
        self._vq = lambdify(vars_t, up('v_q'), modules='numpy')

    def __call__(self, A, B, W, x0, y0, u0, v0, t, rawexpr=False):

        # convert frequency-like variables to complex.
        r2c = lambda x: x * (1. + 0.j)
        A_, B_, W_ = r2c(A), r2c(B), r2c(W)

        if rawexpr: # use directly lambdified

            vars_t = (A_, B_, W_, x0, y0, u0, v0, t)
            x_, y_ = self._xd(*vars_t), self._yd(*vars_t)
            u_, v_ = self._ud(*vars_t), self._vd(*vars_t)

            return np.real(x_), np.real(y_), np.real(u_), np.real(v_)

        # (else) use simplified

        # aux vars
        V_ = self._V(A_, B_, W_)
        P_, Q_ = self._P(A_, B_), self._Q(A_, B_)
        R_, S_ = self._R(V_, W_), self._S(V_, W_)

        # find init consts.
        vars_t = (A_, B_, P_, Q_, R_, S_, V_, W_, x0, y0, u0, v0)
        a_, b_ = self._a(*vars_t), self._b(*vars_t)
        c_, d_ = self._c(*vars_t), self._d(*vars_t)
        
        # find orbits.
        vars_t = (A_, B_, P_, Q_, R_, S_, V_, W_, a_, b_, c_, d_, t)
        x_, y_ = self._x(*vars_t), self._y(*vars_t)
        u_, v_ = self._u(*vars_t), self._v(*vars_t)

        return np.real(x_), np.real(y_), np.real(u_), np.real(v_)

    def W(self, A, B, W, x0, y0, u0, v0, t, rawexpr=False):
        return self.__call__(A, B, W, x0, y0, u0, v0, t, rawexpr)

    def Wp(self, A, B, W, x0, y0, u0, v0, t,):

        # convert frequency-like variables to complex.
        r2c = lambda x: x * (1. + 0.j)
        A_, B_, W_ = r2c(A), r2c(B), r2c(W)

        vars_t = (A_, B_, W_, x0, y0, u0, v0, t)
        x_, y_ = self._xp(*vars_t), self._yp(*vars_t)
        u_, v_ = self._up(*vars_t), self._vp(*vars_t)

        return np.real(x_), np.real(y_), np.real(u_), np.real(v_)

    def Wq(self, A, B, W, x0, y0, u0, v0, t,):

        # convert frequency-like variables to complex.
        r2c = lambda x: x * (1. + 0.j)
        A_, B_, W_ = r2c(A), r2c(B), r2c(W)

        vars_t = (A_, B_, W_, x0, y0, u0, v0, t)
        x_, y_ = self._xq(*vars_t), self._yq(*vars_t)
        u_, v_ = self._uq(*vars_t), self._vq(*vars_t)

        return np.real(x_), np.real(y_), np.real(u_), np.real(v_)
