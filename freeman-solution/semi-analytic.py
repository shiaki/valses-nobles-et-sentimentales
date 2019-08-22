#!/usr/bin/env python

'''
    Find general solution of orbits in Freeman's bar potential.
'''

from sympy import Symbol, symbols, Matrix, Eq, \
        exp, sqrt, solve, init_printing
init_printing(use_unicode=True)

if __name__ == '__main__':

    # variable
    t = Symbol('t', real=True)

    # constants
    A = Symbol('A', positive=True, real=True)
    B = Symbol('B', positive=True, real=True)
    W = Symbol('W', positive=True, real=True)

    # ODE system
    F = Matrix([
        [     0,      0,      1,      0],
        [     0,      0,      0,      1],
        [-A * A,      0,      0,     -W],
        [     0, -B * B,      W,      0]
    ])

    # eigen vec and vals.
    lam, N_lam, vec = zip(*F.eigenvects())

    # check: should have four eigenvalues
    assert len(N_lam) == 4

    # constants in general solution.
    a = Symbol('a',)
    b = Symbol('b',)
    c = Symbol('c',)
    d = Symbol('d',)

    # general solution (with undetermined constants)
    x =   a * exp(lam[0] * t) * vec[0][0][0] \
        + b * exp(lam[1] * t) * vec[1][0][0] \
        + c * exp(lam[2] * t) * vec[2][0][0] \
        + d * exp(lam[3] * t) * vec[3][0][0] 

    y =   a * exp(lam[0] * t) * vec[0][0][1] \
        + b * exp(lam[1] * t) * vec[1][0][1] \
        + c * exp(lam[2] * t) * vec[2][0][1] \
        + d * exp(lam[3] * t) * vec[3][0][1]

    u =   a * exp(lam[0] * t) * vec[0][0][2] \
        + b * exp(lam[1] * t) * vec[1][0][2] \
        + c * exp(lam[2] * t) * vec[2][0][2] \
        + d * exp(lam[3] * t) * vec[3][0][2] 

    v =   a * exp(lam[0] * t) * vec[0][0][3] \
        + b * exp(lam[1] * t) * vec[1][0][3] \
        + c * exp(lam[2] * t) * vec[2][0][3] \
        + d * exp(lam[3] * t) * vec[3][0][3] 

    # initial conditions --> const in solution.
    x0, y0, u0, v0 = symbols('x0, y0, u0, v0', real=True)
    ic_consts = solve(
        (
            Eq(x.subs(t, 0), x0),
            Eq(y.subs(t, 0), y0),
            Eq(u.subs(t, 0), u0),
            Eq(v.subs(t, 0), v0),
        ), 
        (a, b, c, d),
    )

    # plug back to the general solution.
    x_t = x.subs(ic_consts)
    y_t = y.subs(ic_consts)
    u_t = u.subs(ic_consts)
    v_t = v.subs(ic_consts)

    # define EJ (3.113)
    Ej =  (u_t ** 2 + v_t ** 2) / 2 \
        + ((A * x_t) ** 2 + (B * y_t) ** 2) / 2 \
        - (x_t ** 2 + y_t ** 2) * W ** 2 / 2
    # 't' is not automatically eliminated here

    # define orbital integral Ej0,
    Ej0 = Symbol('Ej0', real=True)

    # eliminate u0^2 from EJ0 and X(t=0, x0=0) to find v0^2
    u0sq_a = solve(
        Eq(   (u0 ** 2 + v0 ** 2) / 2 \
            + (B * y0) ** 2 / 2 \
            - (W * y0) ** 2 / 2 ,
           Ej0
        ),
        u0 ** 2
    )[0] # should have unique solution!!

    # force u0 to be positive
    u0_a = sqrt(u0sq_a)

    # and solve v0^2
    u0sq_b = solve(
        Eq(x_t.subs({x0: 0, t: 0, u0: u0_a}), 0),
        v0 ** 2
    )[0]