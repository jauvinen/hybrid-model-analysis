#!/usr/bin/python

import os.path
import glob
import sys

class CentralityFilter:
    def __init__(self, datapath, b_min=1.0, b_max=-1.0,
                 npart_min=1, npart_max=-1,
                 suffix="output", logstring=".out."):
        self._path = datapath
        self._bmin = b_min
        self._bmax = b_max
        self._npmin = npart_min
        self._npmax = npart_max
        self._suffix = suffix.replace('.','')
        self._outprefix = os.path.commonprefix(glob.glob(datapath
                                                         +"/*."+self._suffix))
        self._logstring = logstring

    def outputname(self, identifier):
        return self._outprefix+str(identifier)+"."+self._suffix

    def sort_by_logfolder(self, events_by_b, events_by_npart):
        # find jobs that satisfy the selection criteria
        jobid = -1
        impb = -1.0
        npart = -1

        # Information about b and npart is assumed to be
        # found in output log files, stored in "log" folder
        logpath = self._path+"/log/"
        if not os.path.isdir(logpath):
            print "filter_events: Bad logfile path:", logpath
        logfiles = glob.glob(logpath+"*"+self._logstring+"*")
        processed_files = []
        files_to_remove = []
        for logfile in logfiles:
            datafile = ""
            with open(logfile, 'r') as lf:
                for line in lf:
                    words = line.split()
                    # Each output should be generated by different random seed
                    if "rsd" in words:
                        jobid = words[1]
                        datafile = self.outputname(jobid)
                        # If seed is not unique, output has been overwritten
                        # and information about b and npart has become
                        # ambiguous. Mark such files for deletion.
                        if datafile in processed_files:
                            sys.stderr.write("Warning: Multiple log files for one output.\n")
                            sys.stderr.write(str(datafile)+" will be ignored.\n")
                            files_to_remove.append(datafile)
                            break
                        else:
                            processed_files.append(datafile)
                    if "bimp:" in words:
                        impb = float(words[1])
                        if impb > self._bmin and impb < self._bmax:
                            events_by_b.append(datafile)
                    if "participants:" in words:
                        npart = int(words[1])
                        if npart > self._npmin and npart < self._npmax:
                            events_by_npart.append(datafile)

        # Remove data files which were ambiguously determined
        # due to being referred by multiple log files
        for datafile in files_to_remove:
            if datafile in events_by_b:
                events_by_b.remove(datafile)
            if datafile in events_by_npart:
                events_by_npart.remove(datafile)


    def sort_by_bfile(self, bfiles, events_by_b):
        if len(bfiles) > 1:
            print "Warning: Several .b files detected. Using the first in list:"
            print bfiles[0]
        datapath = os.path.dirname(bfiles[0])
        with open(bfiles[0], 'r') as bf:
            for line in bf:
                data = line.split()
                try:
                    jobid = data[0]
                    impb = float(data[1])
                    if impb >= self._bmin and impb <= self._bmax:
                        events_by_b.append(self.outputname(jobid))
                except ValueError:
                    continue


    def sort_by_npartfile(self, npfiles, events_by_npart):
        if len(npfiles) > 1:
            print "Warning: Several .npart files detected. Using the first in list:"
            print npfiles[0]
        datapath = os.path.dirname(npfiles[0])
        with open(npfiles[0], 'r') as npf:
            for line in npf:
                data = line.split()
                try:
                    jobid = data[0]
                    npart = int(data[1])
                    if npart >= self._npmin and npart <= self._npmax:
                        events_by_npart.append(self.outputname(jobid))
                except ValueError:
                    continue

    # Function for selecting a subset of events
    # based on impact parameter b
    # or number of participants npart
    def filter_events(self):
        events_by_b = []
        events_by_npart = []

        bfiles = [ f for f in glob.glob(self._path+"/*.b") if os.path.isfile(f) ]
        npfiles = [ f for f in glob.glob(self._path+"/*.npart") if os.path.isfile(f) ]

        if bfiles:
            print "Found a .b file, doing impact parameter filtering."
            self.sort_by_bfile(bfiles, events_by_b)
        if npfiles:
            print "Found a .npart file, doing participant number filtering."
            self.sort_by_npartfile(npfiles, events_by_npart)

        if not bfiles and not npfiles:
            self.sort_by_logfolder(events_by_b, events_by_npart)

        # Return the appropriate list of events
        if events_by_b:
            print len(events_by_b), "data files remain after filtering."
            return events_by_b
        elif events_by_npart:
            print len(events_by_npart), "data files remain after filtering."
            return events_by_npart
        else:
            print "filter_events: None of the events fulfill the required criteria:"
            print "b range:", self._bmin, self._bmax, "Npart range:", self._npmin, self._npmax
