#!/usr/bin/python

import math

def ptcount(particledatalist, ptypelist, ptsums, ypoint=0.0, deltay=2.0):
    for particle in particledatalist:
        if (particle.ptype in ptypelist):
            if (abs(particle.rap - ypoint) < deltay / 2.):
                try:
                    ptsum = ptsums[particle.ptype][0] + particle.pt
                    npart = ptsums[particle.ptype][1] + 1
                    ptsums[particle.ptype] = (ptsum, npart)
                except IndexError:
                    continue
