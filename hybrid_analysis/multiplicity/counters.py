#!/usr/bin/python
""" Custom counting functions """

def ptcount(particledatalist, ptypelist, ptsums, ypoint=0.0, deltay=2.0):
    """ Particle pT counter for mean pT calculation.

    Input:
    particledatalist -- List of ParticleData objects
    ptypelist        -- List of particle types for which to do the count
    ptsums           -- Dictionary of pT sums for given particle type
    ypoint           -- Center point of rapidity bin
    deltay           -- Width of rapidity bin
    """
    for particle in particledatalist:
        if particle.ptype in ptypelist:
            if abs(particle.rap - ypoint) < deltay / 2.:
                try:
                    ptsum = ptsums[particle.ptype][0] + particle.pt
                    npart = ptsums[particle.ptype][1] + 1
                    ptsums[particle.ptype] = (ptsum, npart)
                except IndexError:
                    continue
