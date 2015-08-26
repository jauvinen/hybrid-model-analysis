#!/usr/bin/python
""" Functions for calculating distributions (pT, eta)"""

def ptdistr(particles, ptypelist, ptpoints, deltapt,
            ypoint, deltay, dndptsums, pseudorap=False):
    """ Calculate identified particle pT distribution.

    Input:
    particles        -- List of ParticleData objects
    ptypelist        -- List of particle types to be used in analysis
    ptpoints         -- List of pT bin center points
    deltapt          -- Width of pT bin
    ypoint           -- Rapidity bin center point
    deltay           -- Width of rapidity bin
    dndptsums        -- defaultdict of pT distributions for specified particle types
    pseudorap        -- Acceptance is pseudorapidity instead of rapidity
    """
    for particle in particles:
        if particle.ptype in ptypelist:
            if pseudorap:
                yvalue = particle.pseudorap
            else:
                yvalue = particle.rap
            if abs(yvalue - ypoint) < deltay / 2.:
                for ptpoint in ptpoints:
                    if abs(particle.pt - ptpoint) < deltapt / 2.0:
                        dndptsums[particle.ptype][ptpoint] += 1.0 / particle.pt
                        break


def etadistr(particles, etapoints, deltaeta, dndetasum):
    """ Calculate charged particle pseudorapidity distribution.

    Input:
    particles -- List of ParticleData objects
    etapoints -- List of pseudorapidity bin center points
    deltaeta  -- Width of pseudorapidity bin
    dndetasum -- defaultdict of charged particle count in eta bin
    """
    for particle in particles:
        if abs(particle.charge) > 0:
            for etapoint in etapoints:
                if abs(particle.pseudorap - etapoint) < deltaeta / 2:
                    dndetasum[etapoint] += 1.0
                    break
