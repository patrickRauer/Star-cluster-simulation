from astropy import constants as cons
from astropy import units as u
from astropy.table import Table
import numpy as np
import pylab as pl


class Star:
    mass = 0

    pos = None
    v = None

    def __init__(self, mass, pos, v):
        # ﻿4.4985021515752855e-15 converting factor for G from SI to pc**3/(yr**2*M_sun)
        self.mass = mass*4.4985021515752855e-15
        # ﻿1.0227121650537078e-06 converting factor from SI to pc/yr
        self.pos = pos*1.0227121650537078e-06
        self.v = v


class Cluster:
    bh = None
    stars = None

    def __init__(self, bh, stars):
        self.bh = bh
        self.stars = stars

    def plot(self):
        fig = pl.figure()
        fig.clf()
        sp = fig.add_subplot(111)
        sp.scatter(self.bh.pos[0], self.bh.pos[1], marker='x', c='k')
        x = self.stars[:, 2]
        y = self.stars[:, 3]
        sp.scatter(x, y, marker='.', c='r')
        pl.show()

    def __step__(self, dt):
        pos = self.stars[:, 2: 5]
        r_sq = np.square(pos)
        r = 1/np.power(np.sum(r_sq, axis=-1), 3/2)
        r = pos*np.tile(r, (3, 1)).T

        a = -self.bh.mass*r
        dv = a*dt+self.stars[:, -3:]
        print(pos,
              dv)
        dx = 0.5*a*dt*dt+dv*dt+pos
        self.stars[:, -3:] = dv
        self.stars[:, 2: 5] = dx


if __name__ == '__main__':
    black_hole = Star(1e4, np.zeros(3,), np.zeros(3,))
    star1 = np.array([(1, 1., 0.5, 0., 0, 0, 0.000001, 0),
                      (1, 1., -0.5, 0., 0, 0, 0.000001, 0)])
    print(np.sum(star1[:, 2: 5], axis=-1))
    cl = Cluster(black_hole, star1)
    cl.__step__(10)
    for i in range(2000):
        cl.__step__(100)
    cl.plot()
