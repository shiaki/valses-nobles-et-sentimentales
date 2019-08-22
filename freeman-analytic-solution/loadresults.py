#!/usr/bin/env python

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