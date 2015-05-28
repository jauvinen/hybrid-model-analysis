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
from hybrid_analysis.multiplicity import counters
from hybrid_analysis.v_n import cumulants as cumu
from hybrid_analysis.v_n import eventplane as ep

# Initialization
events = 0

# Integrated yields
# To be compared with PHOBOS data
# PRC75, 024910 (2007)
midy_min = -0.1
midy_max = 0.4
integrated_p = 0.0
integrated_pbar = 0.0

# mean pT
# To be compared with STAR data
# PRC79, 034909 (2009)
meanpt_deltay = 0.2
ptsums = {}
particleidlist = [211, -211, 321, -321, 2212, -2212]
for particletype in particleidlist:
    ptsums[particletype] = (0.0, 0)

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
cumulant_etacut = 1.0
flowptpoints = [0.26, 0.44, 0.64, 0.84, 1.04, 1.24, 1.44, 1.64, 1.86]
flowptbinw = 0.2
qcharges = {}
qphis = {}
for ptpoint in flowptpoints:
    qcharges[ptpoint] = []
    qphis[ptpoint] = []

vn_event_etacut = 0.3
vn_event_sums = array.array('d', [0.0]*8)

observables = ["np_integ", "meanpt", "dndpt", "dndeta", "v24", "v2ep"]

# parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-b", "--impact", type=float, nargs=2,
                    metavar=('bmin', 'bmax'),
                    help="impact parameter range")
parser.add_argument("-n", "--npart", type=int, nargs=2,
                    metavar=('nmin', 'nmax'),
                    help="participant range")
parser.add_argument("-o", "--only", nargs=1, action='append',
                    choices=observables,
                    help="run only listed parts of analysis")
parser.add_argument("-x", "--exclude", nargs=1, action='append',
                    choices=observables,
                    help="exclude listed parts of analysis")
parser.add_argument("datapath",
                    help="path to datafiles")
args = parser.parse_args()

# Centrality filtering
datafiles = []
if args.impact:
    bmin = 0.0
    bmax = 0.0
    if args.impact[0] < args.impact[1]:
        bmin = args.impact[0]
        bmax = args.impact[1]
    else:
        bmin = args.impact[1]
        bmax = args.impact[0]
    print "Impact parameter range:", bmin, bmax
    cfilter = cf.CentralityFilter(args.datapath, b_min=bmin, b_max=bmax)
    datafiles = cfilter.filter_events()
elif args.npart:
    npmin = 0.0
    npmax = 0.0
    if args.npart[0] < args.npart[1]:
        npmin = args.npart[0]
        npmax = args.npart[1]
    else:
        npmin = args.npart[1]
        npmax = args.npart[0]
    print "Npart range:", npmin, npmax
    cfilter = cf.CentralityFilter(args.datapath, npart_min=npmin,
                                  npart_max=npmax)
    datafiles = cfilter.filter_events()

analysis = set()
for obs in observables:
    # Because 'append', args.only and args.exclude
    # both return a list of 1-element lists
    if args.only and [obs] in args.only:
        analysis.add(obs)
    elif args.exclude and [obs] not in args.exclude:
        analysis.add(obs)
    elif not args.only and not args.exclude:
        analysis.add(obs)

# Data analysis
files = 0
skipped_files = 0
for datafile in datafiles:
    if files%100 == 0:
        print "Files read:" ,files
    files += 1
    if not os.path.isfile(datafile):
        skipped_files += 1
        continue

    with open(datafile, 'r') as f:
        reading = True
        while reading:
            particlelist = reader.next_text_event(f)
            if not particlelist:
                reading = False
                continue

            events += 1
            if "np_integ" in analysis:
                integrated_p += sum([ 1 for x in particlelist
                                      if (x.ptype == 2212
                                          and x.rap > midy_min
                                          and x.rap < midy_max) ])
                integrated_pbar += sum([ 1 for x in particlelist
                                         if (x.ptype == -2212
                                             and x.rap > midy_min
                                             and x.rap < midy_max) ])
            if "meanpt" in analysis:
                counters.ptcount(particlelist, particleidlist,
                                 ptsums, deltay=meanpt_deltay)
            if "dndpt" in analysis:
                mlt.ptdistr(particlelist, particleidlist, spectraptpoints,
                            spectraptbinw, ypoint, deltay, dndptsums)
            if "dndeta" in analysis:
                mlt.etadistr(particlelist, nchetapoints, deltaeta, dndetasum)
            if "v24" in analysis:
                for ptpoint in flowptpoints:
                    minpt = ptpoint - flowptbinw / 2.0
                    maxpt = ptpoint + flowptbinw / 2.0
                    (ncharges, phis) = cumu.charged_phis(particlelist,
                                                         ptmin=minpt,
                                                         ptmax=maxpt,
                                                         etacut=cumulant_etacut)
                    qcharges[ptpoint].append(ncharges)
                    qphis[ptpoint].append(phis)

            if "v2ep" in analysis:
                ep.v2v3event(particlelist, vn_event_sums,
                             ptmin=0.2, ptmax=2.0, etacut=vn_event_etacut)

print "Attempted to read", files, "files in total, failures:", skipped_files

# Analysis output
if "np_integ" in analysis:
    yrange = midy_max - midy_min
    print "Integrated yields at midrapidity:",
    print midy_min, "< y <", midy_max
    print "Proton Antiproton"
    print integrated_p / yrange / events, integrated_pbar / yrange / events

if "meanpt" in analysis:
    print "<pT> at |y| <", meanpt_deltay / 2
    for ptype in ptsums:
        print ptype
        try:
            print ptsums[ptype][0] / ptsums[ptype][1]
        except ZeroDivisionError:
            print 0.0

if "dndpt" in analysis:
    print "dn/dpT at y =", ypoint
    for ptype in dndptsums:
        print ptype
        for i in range(0, len(spectraptpoints)):
            print spectraptpoints[i],
            print dndptsums[ptype][i] / events / deltay / spectraptbinw / 2/math.pi

if "dndeta" in analysis:
    print "Average Nch:", sum(dndetasum) * deltaeta / events
    dndeta = [ sumbin / events / deltaeta for sumbin in dndetasum ]
    for i in range(0, len(nchetapoints)):
        print nchetapoints[i], dndeta[i]

if "v24" in analysis:
    print "Flow cumulant analysis for pseudorapidity <", cumulant_etacut
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

if "v2ep" in analysis:
    print "Flow event plane analysis for pseudorapidity <", vn_event_etacut
    ep.v2v3mean(vn_event_sums, events)
