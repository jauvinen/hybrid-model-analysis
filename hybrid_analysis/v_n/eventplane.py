#!/usr/bin/python
""" Functions for computing event plane v2 and v3 """

import math
import numpy
import scipy.special as bessel

def v2v3event(particlelist, vn_event_sums,
              ptmin=0.2, ptmax=2.0, etacut=1.0,
              particletype=None):
    """ Compute single-event vn quantities.

    Input:
    particlelist  -- List of ParticleData objects
    vn_event_sums -- List of <v2>, <v2^2>, <v3>, <v3^2>, and sub-event versions
    ptmin         -- lower pT cut
    ptmax         -- upper pT cut
    etacut        -- pseudorapidity cut (-etacut < x < etacut)
    particletype  -- type specification for identified particle calculations
    """
    if len(particlelist) > 1:
        # Apply pT, rapidity and charge cuts
        filtered_particles = [x for x in particlelist
                              if (x.pt > ptmin and x.pt < ptmax
                                  and abs(x.pseudorap) < etacut
                                  and x.charge != 0)]

        ncharges = len(filtered_particles)
        # Several particles are required for proper analysis
        if ncharges < 2:
            return

        ptarray = numpy.array([x.pt for x in filtered_particles])
        phiarray = numpy.array([x.phi for x in filtered_particles])

        ptphisin2 = numpy.multiply(ptarray, numpy.sin(2 * phiarray))
        ptphicos2 = numpy.multiply(ptarray, numpy.cos(2 * phiarray))
        ptphisin3 = numpy.multiply(ptarray, numpy.sin(3 * phiarray))
        ptphicos3 = numpy.multiply(ptarray, numpy.cos(3 * phiarray))

        # exclude the analysed particle from event plane definition
        # to avoid autocorrelations
        sin2mean = numpy.subtract(ptphisin2.sum(), ptphisin2) / (ncharges - 1)
        cos2mean = numpy.subtract(ptphicos2.sum(), ptphicos2) / (ncharges - 1)
        sin3mean = numpy.subtract(ptphisin3.sum(), ptphisin3) / (ncharges - 1)
        cos3mean = numpy.subtract(ptphicos3.sum(), ptphicos3) / (ncharges - 1)

        # Subevent sums for resolution correction calculation
        sin2subsum = numpy.zeros(2)
        cos2subsum = numpy.zeros(2)
        sin3subsum = numpy.zeros(2)
        cos3subsum = numpy.zeros(2)
        for i in range(0, 2):
            sin2subsum[i] = numpy.sum(ptphisin2[i::2])
            cos2subsum[i] = numpy.sum(ptphicos2[i::2])
            sin3subsum[i] = numpy.sum(ptphisin3[i::2])
            cos3subsum[i] = numpy.sum(ptphicos3[i::2])

        # full event value
        psitwoarray = numpy.arctan2(sin2mean, cos2mean) / 2
        psithreearray = numpy.arctan2(sin3mean, cos3mean) / 3
        if particletype:
            typearray = numpy.array([1 if x.ptype == particletype else 0
                                     for x in filtered_particles])
            v2eventsum = numpy.sum(numpy.multiply(numpy.cos(2 * (phiarray - psitwoarray)), typearray))
            v3eventsum = numpy.sum(numpy.multiply(numpy.cos(3 * (phiarray - psithreearray)), typearray))
            ncharges = numpy.sum(typearray)

            if ncharges < 2:
                return
        else:
            v2eventsum = numpy.sum(numpy.cos(2 * (phiarray - psitwoarray)))
            v3eventsum = numpy.sum(numpy.cos(3 * (phiarray - psithreearray)))

        vn_event_sums[0] += v2eventsum / ncharges
        vn_event_sums[1] += (v2eventsum / ncharges)**2
        vn_event_sums[2] += v3eventsum / ncharges
        vn_event_sums[3] += (v3eventsum / ncharges)**2

        psi2sub = numpy.arctan2(sin2subsum, cos2subsum) / 2
        psi3sub = numpy.arctan2(sin3subsum, cos3subsum) / 3

        psi2term = math.cos(2 * (psi2sub[0] - psi2sub[1]))
        psi3term = math.cos(3 * (psi3sub[0] - psi3sub[1]))
        vn_event_sums[4] += psi2term
        vn_event_sums[5] += psi2term**2
        vn_event_sums[6] += psi3term
        vn_event_sums[7] += psi3term**2



def v2v3mean(vn_event_sums, nevents):
    """ Calculate resolution-corrected, event-averaged v2 and v3.

    Input:
    vn_event_sums -- List of <v2>, <v2^2>, <v3>, <v3^2>, and sub-event versions
    nevents       -- Number of events
    """
    if nevents > 0:
        print "Events:", nevents
        meanv2 = vn_event_sums[0] / nevents
        meanv2sqr = vn_event_sums[1] / nevents
        meanv3 = vn_event_sums[2] / nevents
        meanv3sqr = vn_event_sums[3] / nevents

        errv2 = math.sqrt((meanv2sqr - meanv2**2) / (nevents - 1))
        errv3 = math.sqrt((meanv3sqr - meanv3**2) / (nevents - 1))

        print "Uncorrected v2:", meanv2, "err:", errv2,
        print "v3:", meanv3, "err:", errv3

        submean = [0.0, 0.0, 0.0, 0.0]

        for i in range(0, len(submean)):
            if vn_event_sums[i+4] > 0.0:
                submean[i] = math.sqrt(vn_event_sums[i+4] / nevents)

        print "Rsub2:", submean[0], "err:", submean[1],
        print "Rsub3:", submean[2], "err:", submean[3]

        v2corr = 0.0
        v2err = 0.0
        v3corr = 0.0
        v3err = 0.0

        print "R2_full chi2 R3_full chi3"
        v2 = True
        for vnsub in [submean[0], submean[2]]:
            chisub = 2.0
            delta = 8.0
            for iteration in range(0, 20):
                arg = chisub**2 / 4
                result = (math.sqrt(math.pi / 2) / 2 * chisub
                          * math.exp(-arg) * (bessel.i0(arg) + bessel.i1(arg)))
                if result < vnsub:
                    chisub += delta
                elif result > vnsub:
                    chisub -= delta
                else:
                    break
                delta /= 2

            chifull = math.sqrt(2) * chisub
            arg = chifull**2 / 4
            reso_corr = (math.sqrt(math.pi / 2) / 2 * chifull
                         * math.exp(-arg) * (bessel.i0(arg) + bessel.i1(arg)))
            print reso_corr, chifull,

            if v2:
                if reso_corr > 0.0:
                    v2corr = meanv2 / reso_corr
                    v2err = errv2 / reso_corr
                v2 = False
            elif reso_corr > 0.0:
                v3corr = meanv3 / reso_corr
                v3err = errv3 / reso_corr

        print ""
        print "Corrected v2 err v3 err"
        print v2corr, v2err, v3corr, v3err

    else:
        print "No events found for event plane v_n determination."
