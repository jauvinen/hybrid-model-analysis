#!/usr/bin/python

def charged_phis(particlelist, ptmin=0.2, ptmax=2.0):
    ncharges = 0
    phis_in_event = []
    for particle in particlelist:
        if ((particle.charge != 0) and particle.pt > ptmin
            and particle.pt < ptmax):
            ncharges += 1
            phis_in_event.append(particle.phi)

    return (ncharges, phis_in_event)
