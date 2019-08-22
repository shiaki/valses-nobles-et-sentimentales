#!/usr/bin/env python

'''
    Calculate SOS
'''

import pickle
from sympy import Symbol, symbols, Matrix, Eq, \
        exp, solve, sqrt, simplify, init_printing, latex
init_printing(use_unicode=True)

if __name__ == '__main__':

    # read existing variables.
    with open('solution.pk', 'rb') as f:
        shelf = pickle.load(f)
        for k in shelf:
            if k == 'shelf': continue
            globals()[k] = pickle.loads(shelf[k])

    # now let's assume that x0=0 (orbit starts at y-axis)

    # find u0 using EJ_0
    u0sq_a = solve(
        Eq(   (u0 ** 2 + v0 ** 2) / 2 \
            + (B * y0) ** 2 / 2 \
            - (W * y0) ** 2 / 2 ,
           Ej0
        ),
        u0 ** 2
    )[0] # should have unique solution!!

    # assume positive u0
    u0_a = sqrt(u0sq_a)

    # put this u0 into y_t to find relation between y0, v0
    y0_sol = solve(
        Eq( simplify(
                y_t.subs({u0: u0_a, x0: 0, t: 0})),
            y0),
        v0)
    v0_sol = solve(
        Eq( simplify(
                v_t.subs({u0: u0_a, x0: 0, t: 0})),
            v0),
        y0)

    # at least one of them is valid!
    with open('sos-solution.tex', 'w') as f:
        if isinstance(y0_sol, list) and y0_sol:
            f.write('Solutions for $v_0$:\n')
            for s in y0_sol: f.write(latex(s) + '\n')
        if isinstance(v0_sol, list) and v0_sol:
            f.write('Solutions for $y_0$:\n')
            for s in v0_sol: f.write(latex(s) + '\n')
