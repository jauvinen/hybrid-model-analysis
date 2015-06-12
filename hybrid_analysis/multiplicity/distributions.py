#!/usr/bin/python
""" Functions for calculating distributions (pT, eta)"""

def ptdistr(particledatalist, ptypelist, ptpoints, deltapt,
            ypoint, deltay, dndptsums):
    """ Calculate identified particle pT distribution.

    Input:
    particledatalist -- List of ParticleData objects
    ptypelist        -- List of particle types to be used in analysis
    ptpoints         -- List of pT bin center points
    deltapt          -- Width of pT bin
    ypoint           -- Rapidity bin center point
    deltay           -- Width of rapidity bin
    dndptsums        -- Dictionary of pT distributions for specified particle types
    """
    for particle in particledatalist:
        if particle.ptype in ptypelist:
            if abs(particle.rap - ypoint) < deltay / 2.:
                # Note that the default index value is beyond the list
                ptbin = len(ptpoints)
                for i in range(0, len(ptpoints)):
                    if abs(particle.pt - ptpoints[i]) < deltapt / 2.0:
                        ptbin = i
                        break
                try:
                    dndptsums[particle.ptype][ptbin] = (dndptsums[particle.ptype][ptbin]
                                                        + 1.0 / particle.pt)
                except IndexError:
                    continue

def etadistr(particles, etapoints, deltaeta, dndetasum):
    """ Calculate charged particle pseudorapidity distribution.

    Input:
    particles -- List of ParticleData objects
    etapoints -- List of pseudorapidity bin center points
    deltaeta  -- Width of pseudorapidity bin
    dndetasum -- Dictionary of charged particle count in eta bin
    """
    for particle in particles:
        if abs(particle.charge) > 0:
            etabin = len(etapoints)
            for i in range(0, len(etapoints)):
                if abs(particle.pseudorap - etapoints[i]) < deltaeta / 2:
                    etabin = i
                    break
            try:
                dndetasum[etabin] += 1.0
            except IndexError:
                continue
