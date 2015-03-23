#!/usr/bin/python

import argparse
import sys
import os.path
import math
import array

from hic import flow
from hybrid_analysis.file_reader import hybrid_reader as reader
from hybrid_analysis.multiplicity import distributions as mlt
from hybrid_analysis.v_n import cumulants as cumu
from hybrid_analysis.v_n import eventplane as ep
from hybrid_analysis.event_selection import centrality_filters as cf

ypoint = 0.0
deltay = 1.0

deltapt = 0.2
maxpt = 3.0
ptbins = int(maxpt / deltapt)
dndptsums = {}
particleidlist = [211, -211, 321, -321, 2212, -2212]
for particletype in particleidlist:
    dndptsum = [0.0] * ptbins
    dndptsums[particletype] = dndptsum

deltaeta = 0.2
etamin = -5
etabins = int(-2 * etamin / deltaeta)
dndetasum = [0.0] * etabins

events = 0
files = 0
charges_list = []
phis_list = []

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

datafiles = []
if args.impact:
    datafiles = cf.filter_events(args.datapath,
                                 b_min=args.impact[0],
                                 b_max=args.impact[1])
elif args.npart:
    datafiles = cf.filter_events(args.datapath,
                                 npart_min=npart[0],
                                 npart_max=npart[1])

for datafile in datafiles:
    if files%100 == 0:
        print "Files read:" ,files
    files += 1
    eventlist = reader.read_afterburner_output(datafile)
    if eventlist:
        events += len(eventlist)
        for particlelist in eventlist:
            mlt.ptdistr(particlelist, particleidlist, deltapt, ypoint, deltay,
                        dndptsums)
            mlt.etadistr(particlelist, deltaeta, etamin, dndetasum)
            (ncharges, phis) = cumu.charged_phis(particlelist, ptmin=0.2, ptmax=2.0)
            charges_list.append(ncharges)
            phis_list.append(phis)
            ep.v2v3event(particlelist, vn_event_sums,
                         ptmin=0.2, ptmax=2.0, etacut=1.0)

ptpoints = []
for i in range(0, ptbins):
    ptpoints.append(i * deltapt + deltapt / 2)

print "dn/dpT at y =", ypoint
for ptype in dndptsums:
    print ptype
    for i in range(0, len(ptpoints)):
        print ptpoints[i], dndptsums[ptype][i] / events / deltay / deltapt / 2 / math.pi

print "Average Nch:", sum(dndetasum) * deltaeta / events
dndeta = [ sumbin / events / deltaeta for sumbin in dndetasum ]
etapoints = []
for i in range(0,etabins):
    etapoints.append(i * deltaeta + etamin + deltaeta / 2)
for i in range(0, len(etapoints)):
    print etapoints[i], dndeta[i]

q2_list = []
q4_list = []
for event in phis_list:
    q2_list.append(flow.qn(event, 2))
    q4_list.append(flow.qn(event, 4))

vnk = flow.Cumulant(charges_list, q2=q2_list, q4=q4_list)
v22 = vnk.flow(2, 2)
v24 = vnk.flow(2, 4)

print "v2{2}:", v22, "v2{4}:", v24

ep.v2v3mean(vn_event_sums, events)
