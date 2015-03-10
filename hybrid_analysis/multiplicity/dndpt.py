#!/usr/bin/python

import os
import sys
import math
import matplotlib.pyplot as pplot

deltay = 1.0
ypoint = 0.0
deltapt = 0.2
maxpt = 3.0
ptbins = int(maxpt / deltapt)
dndptsumpipl = [0.0] * ptbins
dndptsumpimi = [0.0] * ptbins
dndptsump = [0.0] * ptbins
dndptsumpbar = [0.0] * ptbins
dndptsumKpl = [0.0] * ptbins
dndptsumKmi = [0.0] * ptbins

ptpoints = []
for i in range(0,ptbins):
    ptpoints.append(i * deltapt + deltapt / 2)

events = 0
np = 0
for arg in sys.argv[1:]:
    # if the argument is a file, read the contents
    if (os.path.isfile(arg)):
        initial_list = False
        read_data = False
        # read a line from the file
        for line in open(arg, "r"):
            # split the line into words
            inputline = line.split()
	    if (inputline and inputline[0] == "#event"):
                events += 1
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
                    ptype = abs(int(inputline[8]))
                except ValueError:
                    continue

                # analysis for pi, K, p
                if (ptype == 211 or ptype == 321 or ptype == 2212):
                    E = float(inputline[7])
                    px = float(inputline[4])
                    py = float(inputline[5])
                    pz = float(inputline[6])
                    pt = math.sqrt(px**2 + py**2)
                    try:
                        y = 0.5 * math.log((E + pz)/(E - pz))
                    except ZeroDivisionError:
                        print "!! zero division error with E", E,
                        print "and pz", pz
                        continue
                    except ValueError:
                        print "!! log value error with E", E, "and pz", pz
                        continue

                    if (abs(y - ypoint) < deltay / 2.):
                        ptbin = int(pt / deltapt)
                        try:
                            increment = 1.0 / pt #ptpoints[ptbin]
                            if (ptype == 211):
                                if (charge == 1):
                                    dndptsumpipl[ptbin] += increment
                                elif (charge == -1):
                                    dndptsumpimi[ptbin] += increment
                            elif (ptype == 2212):
                                if (charge == 1):
                                    dndptsump[ptbin] += increment
                                elif (charge == -1):
                                    dndptsumpbar[ptbin] += increment
                            elif (ptype == 321):
                                if (charge == 1):
                                    dndptsumKpl[ptbin] += increment
                                elif (charge == -1):
                                    dndptsumKmi[ptbin] += increment
                        except IndexError:
                            continue

distrpos = [dndptsumpipl, dndptsump, dndptsumKpl]
poslabels = ["$\pi^+$", "$p$", "$K^+$"]

distrneg = [dndptsumpimi, dndptsumpbar, dndptsumKmi]
neglabels = ["$\pi^-$", "$\overline{p}$", "$K^-$"]

print "dn/dpT at y =", ypoint
print "pi-"
for i in range(0, len(ptpoints)):
    print ptpoints[i], dndptsumpimi[i] / events / deltay / deltapt / 2 / math.pi

print "K-"
for i in range(0, len(ptpoints)):
    print ptpoints[i], dndptsumKmi[i] / events / deltay / deltapt / 2 / math.pi

print "p"
for i in range(0, len(ptpoints)):
    print ptpoints[i], dndptsump[i] / events / deltay / deltapt / 2 / math.pi


for i in range(0,3):
    print "Average", poslabels[i], "multiplicity at y=", ypoint, ":",
    print sum(distrpos[i]) * deltay * deltapt / events

    dndpt = [ sumbin / events / deltay / deltapt / 2 / math.pi for sumbin in distrpos[i] ]

    pplot.title('Events: {0:3d}'.format(events))
    pplot.yscale('log')
    pplot.plot(ptpoints,dndpt,linestyle='--',marker='o',label=poslabels[i])
    pplot.xlabel("$p_T$")
    pplot.ylabel("$(1/2\pi p_T)dN/dp_Tdy$")
pplot.legend()
pplot.show()

for i in range(0,3):
    print "Average multiplicity at y=", ypoint, ":",
    print sum(distrneg[i]) * deltay * deltapt / events

    dndpt = [ sumbin / events / deltay / deltapt / 2 / math.pi for sumbin in distrneg[i] ]

    pplot.title('Events: {0:3d}'.format(events))
    pplot.yscale('log')
    pplot.plot(ptpoints,dndpt,linestyle='--',marker='o',label=neglabels[i])
    pplot.xlabel("$p_T$")
    pplot.ylabel("$(1/2\pi p_T)dN/dp_Tdy$")
pplot.legend()
pplot.show()
