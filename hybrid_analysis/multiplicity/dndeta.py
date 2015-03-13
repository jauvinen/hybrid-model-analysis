#!/usr/bin/python

import math
from .. import dataobjects

def etadistr(particles, deltaeta, etamin, dndetasum):
    np = 0
    for particle in particles:
        if (abs(particle.charge) > 0):
            np += 1
            etabin = int((particle.pseudorap - etamin)/deltaeta + 0.5)
            try:
                dndetasum[etabin] += 1.0
            except IndexError:
                continue

    return dndetasum
