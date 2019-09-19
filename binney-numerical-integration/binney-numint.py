#!/usr/bin/env python

'''
'''

import argparse
import numpy as np
from makeinit import makeinit
from binney import rk4int

if __name__ == '__main__':

    # read parameters,
    ps = argparse.ArgumentParser()
    ps.add_argument('V',       type=float, help='Circular velocity Vc')
    ps.add_argument('R',       type=float, help='Core radius Rc')
    ps.add_argument('Q',       type=float, help='Axis ratio Q')
    ps.add_argument('Omega',   type=float, help='Energy level of orbits')
    ps.add_argument('EJ',      type=float, help='Energy level of orbits')
    ps.add_argument('N_steps', type=int,   help='N. of timesteps per orbit')
    ps.add_argument('N_peri',  type=int,   help='N. of periods per orbit')
    ps.add_argument('N_orb',   type=int,   help='N. of orbits')
    ps.add_argument('Fname',   type=str,   help='Output file name')
    args = ps.parse_args()

    # potential parameters
    Vsq, Rsq, Qsq = args.V ** 2, args.R ** 2, args.Q ** 2

    # generate initial condition,
    w_init = makeinit(Vsq, Rsq, Qsq, args.Omega, args.N_steps, args.EJ)

    # integration trial orbits,
    N_testorb = min(32, args.N_orb)
    orb_test = rk4int(w_init[:N_testorb], Vsq, Rsq, Qsq, args.Omega, 0.05, 1000)

    # determine timestep,

    # integrate orbits,

    # save to file.