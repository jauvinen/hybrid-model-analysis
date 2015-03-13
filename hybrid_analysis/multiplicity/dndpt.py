#!/usr/bin/python

import math
from .. import dataobjects

def ptdistr(particledatalist, ptypelist, ypoint, deltay, dndptsums):
    for particle in particledatalist:
        if (particle.ptype in ptypelist):
            if (abs(particle.yrap - ypoint) < deltay / 2.):
                ptbin = int(particle.pt / deltapt)
                try:
                    dndptsums[particle.ptype][ptbin] = (dndptsums[particle.ptype][ptbin]
                                                        + 1.0 / particle.pt)
                except IndexError:
                    continue
