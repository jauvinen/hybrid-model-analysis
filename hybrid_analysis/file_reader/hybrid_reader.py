#!/usr/bin/python

# Functions for reading hybrid model output(s)

import os
import copy
from ..dataobjects import particledata as pd

def read_afterburner_output(filename, read_initial=False):
    id_event = 0
    eventlist = []
    if os.path.isfile(filename):
        read_data = False
        initial_list = False
        particlelist = []
        for line in open(filename, "r"):
            inputline = line.split()
            if (inputline and inputline[0] == "#event"):
                id_event += 1
                read_data = False
                if (len(particlelist) > 0):
                    eventlist.append(copy.copy(particlelist))
                    del particlelist[:]
            if (inputline and inputline[0] == "#cols:"):
                if not initial_list:
                    initial_list = True
                    if read_initial:
                        read_data = True
                else:
                    initial_list = False
                    if read_initial:
                        read_data = False
                    else:
                        read_data = True

            if read_data and len(inputline) == 13:
                try:
                    px = float(inputline[4])
                    py = float(inputline[5])
                    pz = float(inputline[6])
                    E = float(inputline[7])
                    ptype = int(inputline[8])
                    charge = int(inputline[10])
                except ValueError:
                    continue
                particledata = pd.ParticleData([E, px, py, pz], ptype, charge)
                particlelist.append(particledata)

        if (len(particlelist) > 0):
            eventlist.append(copy.copy(particlelist))
            del particlelist[:]

    return eventlist
