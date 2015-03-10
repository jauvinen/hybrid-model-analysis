#!/usr/bin/python

import sys
import os
import math
import copy
from hic import flow

ptmin = 0.2
ptmax = 2.0

nevents = 0
events_so_far = 1

charges_list = []
phis_list = []

phis_in_event = []

for arg in sys.argv[1:]:
    if (os.path.isfile(arg)):
        initial_list = False
        read_data = False
        ncharges = 0
        with open(arg,'r') as f:
            for line in f:
                inputline = line.split()
                if (inputline and inputline[0] == "#event"):
                    nevents += 1
                    read_data = False
                    initial_list = False
                if (inputline and inputline[0] == "#cols:"):
                    if not initial_list:
                        initial_list = True
                    else:
                        initial_list = False
                        read_data = True

                if read_data and len(inputline) == 13:
                    try:
                        time = float(inputline[3])
                    except ValueError:
                        continue
                    try:
                        charge = int(inputline[10])
                    except ValueError:
                        continue

                    px = float(inputline[4])
                    py = float(inputline[5])
                    pz = float(inputline[6])
                    ptr = math.sqrt(px**2 + py**2)
                    pabs = math.sqrt(px**2 + py**2 + pz**2)
                    eta = 99
                    if ((pabs - pz > 0) and (pabs + pz > 0)):
                        eta = 0.5 * math.log((pabs + pz) / (pabs - pz))

                    if (charge != 0 and abs(eta) < 1
                        and ptr > ptmin and ptr < ptmax):
                        ncharges += 1
                        phis_in_event.append(math.atan2(py,px))

                elif (nevents > events_so_far):
                    if (ncharges > 0):
                        charges_list.append(ncharges)
                        phis_list.append(copy.copy(phis_in_event))
                    events_so_far = nevents
                    ncharges = 0
                    del phis_in_event[:]

q2_list = []
q4_list = []
for event in phis_list:
    q2_list.append(flow.qn(event, 2))
    q4_list.append(flow.qn(event, 4))

vnk = flow.Cumulant(charges_list, q2=q2_list, q4=q4_list)
v22 = vnk.flow(2, 2)
v24 = vnk.flow(2, 4)

print "v2{2}:", v22, "v2{4}:", v24
