#!/usr/bin/python

import math

class ParticleData:
    def __init__(self, momentum, ptype, charge):
        self.pt = math.sqrt(momentum[1]**2 + momentum[2]**2)
        self.phi = math.atan2(momentum[2], momentum[1])
        self.pz = momentum[3]
        try:
            self.rap = 0.5 * math.log((momentum[0] + momentum[3])
                                      / (momentum[0] - momentum[3]))
        except ZeroDivisionError:
            print "!! zero division error with E", momentum[0],
            print "and pz", momentum[3]
            self.rap =  1000.0
            continue
        except ValueError:
            print "!! log value error with E", momentum[0],
            print "and pz", momentum[3]
            self.rap =  -1000.0
            continue

        pabs = math.sqrt(momentum[1]**2 + momentum[2]**2 + momentum[3]**2)
        try:
            self.pseudorap = 0.5 * math.log((pabs + momentum[3])
                                            /(pabs - momentum[3]))
        except ZeroDivisionError:
            print "!! zero division error with pabs", pabs,
            print "and pz", momentum[3]
            self.pseudorap = 1000.0
            continue
        except ValueError:
            print "!! log value error with pabs", pabs, "and pz",
            print momentum[3]
            self.pseudorap = -1000.0
            continue

        self.ptype = ptype
        self.charge = charge
