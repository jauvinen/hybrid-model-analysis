#!/usr/bin/python
""" Functions related to ParticleData class """

import math

class ParticleData:
    """ Class for handling data needed by analysis """
    def __init__(self, momentum, ptype, charge):
        self.pt = math.sqrt(momentum[1]**2 + momentum[2]**2)
        self.phi = math.atan2(momentum[2], momentum[1])
        self.pz = momentum[3]
        try:
            self.rap = 0.5 * math.log((momentum[0] + momentum[3])
                                      / (momentum[0] - momentum[3]))
        except ZeroDivisionError:
            self.rap = 1000.0
        except ValueError:
            self.rap = -1000.0

        pabs = math.sqrt(momentum[1]**2 + momentum[2]**2 + momentum[3]**2)
        try:
            self.pseudorap = 0.5 * math.log((pabs + momentum[3])
                                            /(pabs - momentum[3]))
        except ZeroDivisionError:
            self.pseudorap = 1000.0
        except ValueError:
            self.pseudorap = -1000.0

        self.ptype = ptype
        self.charge = charge
