#!/usr/bin/env python

'''
    Read orbits data.
'''

import os
import numpy as np

class ReadOrbits(object):

    def __init__(self, I_file,
                 rawdata_path='/dgalex2/qinyj/data/3042_orbsave/'):

        rawdata_fname = rawdata_path + ('orbits%03u' % I_file) + '.res'

        print('\nLoading ', rawdata_fname, '...')
        print('  Huge file, may take a few minutes...')

        raw_file = open(rawdata_fname, 'rb')

        raw_orbits = np.fromfile(raw_file, dtype = np.float32)[4:]
        raw_orbits = np.transpose(np.reshape(raw_orbits, (5001, -1)))

        time_pts   = raw_orbits[1]
        delta_t    = (time_pts[-1] - time_pts[0]) / (time_pts.shape[0] - 1)
        raw_orbits = raw_orbits[2:-1]

        raw_file.close()

        self.raw_orbits = raw_orbits
        self.time_pts = time_pts
        self.delta_t = delta_t

        print('  Imported successfully.')
        print('delta_t is', delta_t)

    def get_orbit(self, I_orb, velocity_method='full', packed=True):

        # kernels for 2nd order numerical differentiation
        ctr_diff = -np.array([1. / 12., -2. / 3., 0., 2. / 3., -1. / 12.])
        fwd_diff = -np.array([-3. / 2., 2., -1. / 2.])
        bwd_diff = -np.array([1. / 2., -2., 3. / 2.])
        # Note the negative sign here!

        # read orbits from the big array
        pos = self.raw_orbits[3 * I_orb: 3 * I_orb + 3]

        if velocity_method == 'full':
            # return the entire orbit.

            if self.time_pts.size < 5:
                raise RuntimeError('Orbit is too short' \
                        + ' for numerical differentiation.')

            vel = np.zeros(shape = (3, pos.shape[1]))

            for i in (0, 1, 2): # velocity components,

                # first two points: interpolate with forward method.
                vel[i][  : 2] = np.convolve(pos[i][  :4], fwd_diff, 
                        mode='valid') / self.delta_t                

                vel[i][ 2:-2] = np.convolve(pos[i]      , ctr_diff,
                        mode='valid') / self.delta_t

                # last two points: using backward method.
                vel[i][-2:  ] = np.convolve(pos[i][-4: ], bwd_diff,
                        mode='valid') / self.delta_t

            if packed:
                return np.vstack((pos, vel, self.time_pts)).T
            else:
                return pos.T, vel.T, self.time_pts

        elif velocity_method == 'central':
            # using central diff. only, ignores first and last two points.

            vel = np.zeros(shape=(3, pos.shape[1] - 4))

            for i in (0, 1, 2):
                vel[i] = np.convolve(pos[i], ctr_diff,
                        mode = 'valid') / self.delta_t

            if packed:
                return np.vstack((pos[:, 2:-2], vel, self.time_pts[2:-2])).T
            else:
                return pos[:,2:-2].T, vel.T, self.time_pts[2:-2]

        else:
            raise RuntimeError('Invalid parameter: velocity_method')