#!/usr/bin/python

import math

def ptdistr(particledatalist, ptypelist, deltapt, ypoint, deltay, dndptsums):
    for particle in particledatalist:
        if (particle.ptype in ptypelist):
            if (abs(particle.rap - ypoint) < deltay / 2.):
                ptbin = int(particle.pt / deltapt)
                try:
                    dndptsums[particle.ptype][ptbin] = (dndptsums[particle.ptype][ptbin]
                                                        + 1.0 / particle.pt)
                except IndexError:
                    continue

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
