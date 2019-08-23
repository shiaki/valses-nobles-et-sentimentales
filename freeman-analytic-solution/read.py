#!/usr/bin/env python

import pickle
from sympy import *
init_printing(use_unicode=True)

if __name__ == '__main__':

    # read existing variables.
    with open('simplified.pk', 'rb') as f:
        for k, v in pickle.load(f).items():
            globals()[k] = pickle.loads(v)