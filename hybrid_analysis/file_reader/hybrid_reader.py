#!/usr/bin/python

# Functions for reading hybrid model output(s)

import os

def read_afterburner_output(filename, read_initial=False):
    id_event = 0
    particlelist = []
    if os.path.isfile(filename):
        read_data = False
        initial_list = False
        for line in open(filename, "r"):
            inputline = line.split()
            if (inputline and inputline[0] == "#event"):
                id_event += 1
                read_data = False
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

            if read_data:
                try:
                    x = float(inputline[0])
                    ptype = int(inputline[8])
                except ValueError:
                    continue
                particlelist.append((tuple(inputline), id_event))

    return particlelist
