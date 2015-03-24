#!/usr/bin/python

import argparse
import sys
import os.path
import math
import array

from hic import flow
from hybrid_analysis.event_selection import centrality_filters as cf
from hybrid_analysis.file_reader import hybrid_reader as reader
from hybrid_analysis.multiplicity import distributions as mlt
from hybrid_analysis.v_n import cumulants as cumu
from hybrid_analysis.v_n import eventplane as ep

# Initialization
events = 0

# pT spectra
# To be compared with PHOBOS data
# PRC75, 024910 (2007)
# Rapidity point
ypoint = 0.8
deltay = 0.1
spectraptpoints = [0.25, 0.30, 0.35, 0.40, 0.50, 0.55, 0.60, 0.70,
                   1.0, 1.2, 1.55, 1.85, 2.2]
spectraptbinw = 0.05
dndptsums = {}
particleidlist = [211, -211, 321, -321, 2212, -2212]
for particletype in particleidlist:
    dndptsum = [0.0] * len(spectraptpoints)
    dndptsums[particletype] = dndptsum

# Pseudorapidity distribution
# To be compared with PHOBOS data
# PRC83, 024913 (2011)
deltaeta = 0.2
etamin = -5.3
etabins = int(-2 * etamin / deltaeta + 0.5)
nchetapoints = [ (etamin + deltaeta * i) for i in range(0, etabins+1) ]
dndetasum = [0.0] * len(nchetapoints)

# Flow analysis
# To be compared with STAR data
# PRC86, 054908 (2012)
flowptpoints = [0.26, 0.44, 0.64, 0.84, 1.04, 1.24, 1.44, 1.64, 1.86, 2.19]
flowptbinw = 0.8
qcharges = {}
qphis = {}
for ptpoint in flowptpoints:
    qcharges[ptpoint] = []
    qphis[ptpoint] = []

vn_event_sums = array.array('d', [0.0]*8)

# parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-b", "--impact", type=float, nargs=2,
                    metavar=('bmin', 'bmax'),
                    help="impact parameter range")
parser.add_argument("-n", "--npart", type=int, nargs=2,
                    metavar=('nmin', 'nmax'),
                    help="participant range")
parser.add_argument("datapath",
                    help="path to datafiles")
args = parser.parse_args()

# Centrality filtering
datafiles = []
if args.impact:
    datafiles = cf.filter_events(args.datapath,
                                 b_min=args.impact[0],
                                 b_max=args.impact[1])
elif args.npart:
    datafiles = cf.filter_events(args.datapath,
                                 npart_min=npart[0],
                                 npart_max=npart[1])

# Data analysis
files = 0
for datafile in datafiles:
    if files%100 == 0:
        print "Files read:" ,files
    files += 1
    eventlist = reader.read_afterburner_output(datafile)
    if eventlist:
        events += len(eventlist)
        for particlelist in eventlist:
            mlt.ptdistr(particlelist, particleidlist, spectraptpoints,
                        spectraptbinw, ypoint, deltay, dndptsums)
            mlt.etadistr(particlelist, nchetapoints, deltaeta, dndetasum)
            for ptpoint in flowptpoints:
                minpt = ptpoint - flowptbinw / 2.0
                maxpt = ptpoint + flowptbinw / 2.0
                (ncharges, phis) = cumu.charged_phis(particlelist, ptmin=minpt,
                                                     ptmax=maxpt)
                qcharges[ptpoint].append(ncharges)
                qphis[ptpoint].append(phis)

            ep.v2v3event(particlelist, vn_event_sums,
                         ptmin=0.2, ptmax=2.0, etacut=1.0)

# Analysis output
print "dn/dpT at y =", ypoint
for ptype in dndptsums:
    print ptype
    for i in range(0, len(spectraptpoints)):
        print spectraptpoints[i],
        print dndptsums[ptype][i] / events / deltay / spectraptbinw / 2/math.pi

print "Average Nch:", sum(dndetasum) * deltaeta / events
dndeta = [ sumbin / events / deltaeta for sumbin in dndetasum ]
for i in range(0, len(nchetapoints)):
    print nchetapoints[i], dndeta[i]

print "Flow cumulant analysis"
print "pT v2{2} v2{4}"
for ptpoint in flowptpoints:
    q2_list = []
    q4_list = []
    for event in qphis[ptpoint]:
        q2_list.append(flow.qn(event, 2))
        q4_list.append(flow.qn(event, 4))

    vnk = flow.Cumulant(qcharges[ptpoint], q2=q2_list, q4=q4_list)
    v22 = vnk.flow(2, 2)
    v24 = vnk.flow(2, 4)
    print ptpoint, v22, v24

print "Flow event plane analysis"
ep.v2v3mean(vn_event_sums, events)
