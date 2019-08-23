#!/usr/bin/env python

'''
    Simplify things.
'''

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

    x_p, y_p, u_p, v_p = up('x_p'), up('y_p'), up('u_p'), up('v_p')
    x_q, y_q, u_q, v_q = up('x_q'), up('y_q'), up('u_q'), up('v_q')

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

    # {A, B, P, Q, R, S, V, W, u0, v0, x0, y0}
    a_l = ic_consts[a].subs(subst_dict_aux)
    b_l = ic_consts[b].subs(subst_dict_aux)
    c_l = ic_consts[c].subs(subst_dict_aux)
    d_l = ic_consts[d].subs(subst_dict_aux)

    # save global variables.
    shelf = dict()
    for k in   [z + '_l' for z in 'abcdxyuv'] \
             + [z +  '0' for z in 'xyuv'] + [z + '_t' for z in 'xyuv'] \
             + [z + '_p' for z in 'xyuv'] + [z + '_q' for z in 'xyuv'] \
             + list('ABPQRSVWabcdt'):
        try: shelf[k] = pickle.dumps(globals()[k])
        except TypeError: continue # ignore un-pickables
    with open('simplified.pk', 'wb') as f:
        pickle.dump(shelf, f)

    # save solution.
    with open('orbit-solution.tex', 'w') as f:
        for var_t, expr_t in [(sk, latex(k)) for sk, k in \
                [('x(t)',   x_l), ('y(t)',   y_l), \
                 ('V_x(t)', u_l), ('V_y(t)', v_l)]]:
            f.write('\[\n{:}={:}\n\]\n'.format(var_t, expr_t))

    with open('initconsts-solution.tex', 'w') as f:
        for var_t, expr_t in [(sk, latex(k)) for sk, k in \
                [('a', a_l), ('b', b_l), ('c', c_l), ('d', d_l)]]:
            f.write('\[\n{:}={:}\n\]\n'.format(var_t, expr_t))
