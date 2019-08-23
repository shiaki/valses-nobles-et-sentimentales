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
    t = sol['t']

    # constants in general solution.
    A, B, W = sol['A'], sol['B'], sol['W']
    a, b, c, d = symbols('a, b, c, d',)

    # eigenvalues and eigenvectors
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
    L_avr = simplify(L_avr.subs(subst_dict))
    L_var = simplify(L_var.subs(subst_dict))

    # define new symbols to make it easier,
    V = symbols('V', positive=True, real=True)
    P, Q, R, S = symbols('P, Q, R, S', real=True)
    subst_dict_aux = {
        sqrt(A ** 2 + 2 * A * B + B ** 2 + W ** 2) \
                * sqrt(A ** 2 - 2 * A * B + B ** 2 + W ** 2): V ** 2,
        (A ** 2 + B ** 2) / 2: P,
        (A ** 2 - B ** 2) / 2: Q,
        (V ** 2 + W ** 2) / 2: R,
        (V ** 2 - W ** 2) / 2: S,
    }
    _subs_0 = {expand(2 * kk): 2 * vv for kk, vv in subst_dict_aux.items()}
    _subs_1 = {simplify(kk): vv for kk, vv in subst_dict_aux.items()}
    _subs_2 = {expand(kk): vv for kk, vv in subst_dict_aux.items()}
    subst_dict_aux.update({**_subs_0, **_subs_1, **_subs_2})

    L_avr = simplify(L_avr.subs(subst_dict_aux))
    L_var = simplify(L_var.subs(subst_dict_aux))

    x0, y0, u0, v0 = sol['x0'], sol['y0'], sol['u0'], sol['v0']

    # consts in general solution with aux variables P, Q, R
    a_p = simplify(sol['ic_consts'][a].subs(subst_dict_aux))
    b_p = simplify(sol['ic_consts'][b].subs(subst_dict_aux))
    c_p = simplify(sol['ic_consts'][c].subs(subst_dict_aux))
    d_p = simplify(sol['ic_consts'][d].subs(subst_dict_aux))

    subst_dict_consts = {
        a * b: simplify(a_p * b_p).subs(subst_dict_aux),
        c * d: simplify(c_p * d_p).subs(subst_dict_aux),
    }

    # full expression with x0, y0, u0, v0 as free parameters.
    L_avr_f = L_avr.subs(subst_dict_consts).subs(subst_dict_aux)
    L_var_f = L_var.subs(subst_dict_consts).subs(subst_dict_aux)

    # save global variables.
    shelf = dict()
    for k in ['L_avr', 'L_var', 'L_avr_f', 'L_var_f']:
        try: shelf[k] = pickle.dumps(globals()[k])
        except TypeError: continue # ignore un-pickables
    with open('L-solution.pk', 'wb') as f:
        pickle.dump(shelf, f)

    # save solution.
    with open('L-solution.tex', 'w') as f:
        for sym_t, expr_t in [('L', latex(L_avr)), \
                ('\mathrm{Var}(L)', latex(L_var))]:
            f.write('\[\n{:}={:}\n\]\n'.format(sym_t, expr_t))
