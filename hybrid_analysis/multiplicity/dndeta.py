#!/usr/bin/python

import os
import sys
import math
import matplotlib.pyplot as pplot

deltaeta = 0.2
etamin = -5
etabins = int(-2 * etamin / deltaeta)
dndetasum = [0.0] * etabins

etapoints = []
for i in range(0,etabins):
    etapoints.append(i * deltaeta + etamin + deltaeta / 2)

events = 0
for arg in sys.argv[1:]:
    # if the argument is a file, read the contents
    if (os.path.isfile(arg)):
	initial_list = False
        read_data = False
        np = 0
        # read a line from the file
        for line in open(arg, "r"):
            # split the line into words
            inputline = line.split()
	    if (inputline and inputline[0] == "#event"):
                events += 1
                np = 0
	        read_data = False
            if (inputline and inputline[0] == "#cols:"):
                if not initial_list:
                    initial_list = True
                else:
                    initial_list = False
                    read_data = True

            if read_data:
                try:
                    time = float(inputline[0])
                except ValueError:
                    continue
                try:
                    charge = int(inputline[10])
                except ValueError:
                    continue

                if (abs(charge) > 0):
                    np += 1
                    px = float(inputline[4])
                    py = float(inputline[5])
                    pz = float(inputline[6])
                    pabs = math.sqrt(px**2 + py**2 + pz**2)
                    try:
                        eta = 0.5 * math.log((pabs + pz)/(pabs - pz))
                    except ZeroDivisionError:
                        print "!! zero division error with pabs", pabs,
                        print "and pz", pz
                        continue
                    except ValueError:
                        print "!! log value error with pabs", pabs, "and pz", pz
                        continue

                    etabin = int((eta + etamin)/deltaeta + 0.5)
                    try:
                        dndetasum[etabin] += 1.0
                    except IndexError:
                        continue

print "Average Nch:", sum(dndetasum) * deltaeta / events

dndeta = [ sumbin / events / deltaeta for sumbin in dndetasum ]

for i in range(0, len(etapoints)):
    print etapoints[i], dndeta[i]

pplot.title('Events: {0:3d}'.format(events))
pplot.plot(etapoints,dndeta)
pplot.xlabel("$\eta$")
pplot.ylabel("dN/d$\eta$")
pplot.show()
