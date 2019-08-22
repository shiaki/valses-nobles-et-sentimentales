#!/usr/bin/env python

'''
    Find general solution of orbits in Freeman's bar potential.
    (using deconstructed matrix ODE)
'''

import pickle
from sympy import *
init_printing(use_unicode=True)

if __name__ == '__main__':

    # read previous results.
    with open('solution.pk', 'rb') as f:
        sol_f = pickle.load(f)
        sol = {k: pickle.loads(v) for k, v in sol_f.items()}

    # variable
    t = Symbol('t', real=True)

    # constants in general solution.
    a, b, c, d = symbols('a, b, c, d',)
    lam_1, lam_2 = symbols('l1, l2')
    K0x, K1x, K2x, K3x = symbols('K0x, K1x, K2x, K3x')
    K0y, K1y, K2y, K3y = symbols('K0y, K1y, K2y, K3y')
    K0u, K1u, K2u, K3u = symbols('K0u, K1u, K2u, K3u')
    K0v, K1v, K2v, K3v = symbols('K0v, K1v, K2v, K3v')

    # general solution (with undetermined constants)
    x =   a * exp( I * lam_1 * t) * K0x \
        + b * exp(-I * lam_1 * t) * K1x \
        + c * exp( I * lam_2 * t) * K2x \
        + d * exp(-I * lam_2 * t) * K3x

    y =   a * exp( I * lam_1 * t) * K0y \
        + b * exp(-I * lam_1 * t) * K1y \
        + c * exp( I * lam_2 * t) * K2y \
        + d * exp(-I * lam_2 * t) * K3y

    u =   a * exp( I * lam_1 * t) * K0u \
        + b * exp(-I * lam_1 * t) * K1u \
        + c * exp( I * lam_2 * t) * K2u \
        + d * exp(-I * lam_2 * t) * K3u

    v =   a * exp( I * lam_1 * t) * K0v \
        + b * exp(-I * lam_1 * t) * K1v \
        + c * exp( I * lam_2 * t) * K2v \
        + d * exp(-I * lam_2 * t) * K3v

    # angular momentum
    L = x * v - y * u
    L_terms = expand(L).args

    # find time-dependent/independent terms in L
    const_terms = [w for w in L_terms if t not in w.free_symbols]
    tvari_terms = [w for w in L_terms if t     in w.free_symbols]

    # two components: const and variable.
    L_avr, L_res = sum(const_terms), sum(tvari_terms)

    # find central variance
    L_var = sum([w for w in expand(L_res ** 2).args \
                 if t not in w.free_symbols])

    # read a, b, c, d, K
    subst_dict = {
        K0x: simplify(sol['vec'][0][0][0]),
        K1x: simplify(sol['vec'][1][0][0]),
        K2x: simplify(sol['vec'][2][0][0]),
        K3x: simplify(sol['vec'][3][0][0]),
        K0y: simplify(sol['vec'][0][0][1]),
        K1y: simplify(sol['vec'][1][0][1]),
        K2y: simplify(sol['vec'][2][0][1]),
        K3y: simplify(sol['vec'][3][0][1]),
        K0u: simplify(sol['vec'][0][0][2]),
        K1u: simplify(sol['vec'][1][0][2]),
        K2u: simplify(sol['vec'][2][0][2]),
        K3u: simplify(sol['vec'][3][0][2]),
        K0v: simplify(sol['vec'][0][0][3]),
        K1v: simplify(sol['vec'][1][0][3]),
        K2v: simplify(sol['vec'][2][0][3]),
        K3v: simplify(sol['vec'][3][0][3]),
    }

    # full expression:
    L_avr_f = simplify(L_avr.subs(subst_dict))
    L_var_f = simplify(L_var.subs(subst_dict))

    # define new symbols to make it easier,
    A, B, W = sol['A'], sol['B'], sol['W']
    P, Q, R, S = symbols('P, Q, R, S', positive=True, real=True)
    subst_dict = {
        A ** 2 + 2 * A * B + B ** 2 + W ** 2: P ** 2,
        A ** 2 - 2 * A * B + B ** 2 + W ** 2: Q ** 2,
        A ** 4 - 2 * A ** 2 * B ** 2 + 2 * A ** 2 * W ** 2 \
            + B ** 4 + 2 * B ** 2 * W ** 2 + W ** 4: P ** 2 * Q ** 2,
        A ** 2 + B ** 2 + W ** 2: R ** 2,
        B ** 2 - A ** 2: S ** 2
    }
    subst_dict.update({expand(2 * kk): 2 * vv \
                       for kk, vv in subst_dict.items()}) # sympy is stupid

    x0, y0, u0, v0 = sol['x0'], sol['y0'], sol['u0'], sol['v0']

    # consts in general solution with aux variables P, Q, R
    a_p = simplify(sol['ic_consts'][a].subs(subst_dict))
    b_p = simplify(sol['ic_consts'][b].subs(subst_dict))
    c_p = simplify(sol['ic_consts'][c].subs(subst_dict))
    d_p = simplify(sol['ic_consts'][d].subs(subst_dict))

    L_avr_f = L_avr_f.subs(subst_dict)
    L_var_f = L_var_f.subs(subst_dict)

    subst_dict = {
        a * b: simplify(a_p * b_p), 
        c * d: simplify(c_p * d_p), 
    }

    # full expression with x0, y0, u0, v0 as free parameters.
    L_avr_f = cancel(L_avr_f.subs(subst_dict))
    L_var_f = cancel(L_var_f.subs(subst_dict))

    # save global variables.
    shelf = dict()
    for k in ['L_avr', 'L_var', 'L_avr_f', 'L_var_f']:
        try: shelf[k] = pickle.dumps(globals()[k])
        except TypeError: continue # ignore un-pickables
    with open('L-solution.pk', 'wb') as f:
        pickle.dump(shelf, f)
