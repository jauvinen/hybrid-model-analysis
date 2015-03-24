#!/usr/bin/python

import math

def ptdistr(particledatalist, ptypelist, ptpoints, deltapt,
            ypoint, deltay, dndptsums):
    for particle in particledatalist:
        if (particle.ptype in ptypelist):
            if (abs(particle.rap - ypoint) < deltay / 2.):
                # Note that the default index value is beyond the list
                ptbin = len(ptpoints)
                for i in range(0, len(ptpoints)):
                    if (abs(particle.pt - ptpoints[i]) < deltapt / 2.0):
                        ptbin = i
                        break
                try:
                    dndptsums[particle.ptype][ptbin] = (dndptsums[particle.ptype][ptbin]
                                                        + 1.0 / particle.pt)
                except IndexError:
                    continue

def etadistr(particles, etapoints, deltaeta, dndetasum):
    for particle in particles:
        if (abs(particle.charge) > 0):
            etabin = len(etapoints)
            for i in range(0, len(etapoints)):
                if (abs(particle.pseudorap - etapoints[i]) < deltaeta / 2):
                    etabin = i
                    break
            try:
                dndetasum[etabin] += 1.0
            except IndexError:
                continue
