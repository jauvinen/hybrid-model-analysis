#!/usr/bin/python

# charged_phis: Count the number of charged particles in an event
# after appropriate transverse momentum and pseudorapidity cuts
# and return list of azimuthal angles of these particles
# Needed for cumulant analysis
#
# Input:
# particlelist: List of particle data for one event
# ptmin: Lower limit of transverse momentum
# ptmax: Upper limit of transverse momentum
# etacut: Pseudorapidity cut (assumed symmetric around 0)
#
# Output:
# ncharges: count of charged particles in the event
# phis_in_event: List of charged particle azimuthal angles
#                (angle in the plane perpendicular to beam)

def charged_phis(particlelist, ptmin=0.2, ptmax=2.0, etacut=1.0):
    ncharges = 0
    phis_in_event = []
    for particle in particlelist:
        if ((particle.charge != 0)
            and particle.pt > ptmin and particle.pt < ptmax
            and abs(particle.pseudorap) < etacut):
            ncharges += 1
            phis_in_event.append(particle.phi)

    return (ncharges, phis_in_event)
