#!/usr/bin/env python

'''
    Find analytic orbits in Freeman's potential using SYMPY's build-in solver

    *** This one does not work!!
'''

from sympy import Symbol, Function, Derivative, Eq, dsolve

if __name__ == '__main__':

    t = Symbol('t')

    A = Symbol('A', positive=True)
    B = Symbol('B', positive=True)
    W = Symbol('W', positive=True)

    x = Function('x')(t)
    y = Function('y')(t)
    u = Function('u')(t)
    v = Function('v')(t)

    dx_dt = Derivative(x, t)
    dy_dt = Derivative(y, t)
    du_dt = Derivative(u, t)
    dv_dt = Derivative(v, t)

    eqs = (
        Eq(dx_dt, u),
        Eq(dy_dt, v),
        Eq(du_dt, -(A ** 2) * x - W * v),
        Eq(dv_dt, -(B ** 2) * y + W * u)
    )

    sol = dsolve(eqs, hint='system_of_odes_linear_neq_order1_type1')
    print(sol)