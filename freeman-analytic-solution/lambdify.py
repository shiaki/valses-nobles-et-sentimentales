#!/usr/bin/env python

import pickle
from sympy import *
init_printing(use_unicode=True)

if __name__ == '__main__':

    # read existing variables.
    with open('solution.pk', 'rb') as f:
        s = pickle.load(f)

    up = lambda k: pickle.loads(s[k])
    A, B, W, t = up('A'), up('B'), up('W'), up('t')
    a, b, c, d = up('a'), up('b'), up('c'), up('d')

    x, x0, x_t = up('x'), up('x0'), up('x_t')
    y, y0, y_t = up('y'), up('y0'), up('y_t')
    u, u0, u_t = up('u'), up('u0'), up('u_t')
    v, v0, v_t = up('v'), up('v0'), up('v_t')

    Ej, Ej0 = up('Ej'), up('Ej0')
    ic_consts = up('ic_consts')

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

    # xyuv(P, Q, R, S, W, a, b, c, d, t)
    x_l = x.subs(subst_dict_aux)
    y_l = y.subs(subst_dict_aux)
    u_l = u.subs(subst_dict_aux)
    v_l = v.subs(subst_dict_aux)

    # {A, B, P, Q, R, S, V, W, u₀, v₀, x₀, y₀}
    a_l = ic_consts[a].subs(subst_dict_aux)
    b_l = ic_consts[b].subs(subst_dict_aux)
    c_l = ic_consts[c].subs(subst_dict_aux)
    d_l = ic_consts[d].subs(subst_dict_aux)

    # save global variables.
    shelf = dict()
    for k in ['x_l', 'y_l', 'u_l', 'v_l', 'a_l', 'b_l', 'c_l', 'd_l']:
        try: shelf[k] = pickle.dumps(globals()[k])
        except TypeError: continue # ignore un-pickables
    with open('lambdified.pk', 'wb') as f:
        pickle.dump(shelf, f)
