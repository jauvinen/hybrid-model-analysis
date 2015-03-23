#!/usr/bin/python

import os.path
import glob

def filter_events(datapath, b_min=1.0, b_max=-1.0,
                  npart_min=1, npart_max=-1):
    events_by_b = []
    events_by_npart = []
    # find jobs that satisfy the selection criteria
    jobid = -1
    impb = -1.0
    npart = -1
    logpath = datapath+"log/"
    if not os.path.isdir(logpath):
        print "filter_events: Bad logfile path:", logpath
    logfiles = glob.glob(logpath+"hybrid.out.*")
    for logfile in logfiles:
        datafile = ""
        with open(logfile, 'r') as lf:
            for line in lf:
                words = line.split()
                if "rsd" in words:
                    jobid = words[1]
                    datafile = datapath+"/afterburner"+jobid+".output"
                if "bimp:" in words:
                    impb = float(words[1])
                    if impb > b_min and impb < b_max:
                        events_by_b.append(datafile)
                if "participants:" in words:
                    npart = int(words[1])
                    if npart > npart_min and npart < npart_max:
                        events_by_npart.append(datafile)

    if events_by_b:
        print len(events_by_b), "data files remain after filtering."
        return events_by_b
    elif events_by_npart:
        print len(events_by_npart), "data files remain after filtering."
        return events_by_npart
    else:
        print "filter_events: None of the events fulfill the required criteria:"
        print "b range:", b_min, b_max, "Npart range:", npart_min, npart_max
