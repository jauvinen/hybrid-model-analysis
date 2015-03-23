#!/usr/bin/python

# Functions for computing event plane v2 and v3

import math
import numpy
import scipy.special as bessel

# Compute single-event vn quantities
# Input: List of particles, containing their momentum and azimuthal angle
# Return: Tuple of <v2>, <v2^2>, <v3>, <v3^2>, and sub-event versions
def v2v3event(particlelist, vn_event_sums,
              ptmin=0.2, ptmax=2.0, etacut=1.0):

    ncharges = 0
    if (len(particlelist) > 1):
        filtered_particles = [ x for x in particlelist
                               if (x.pt > ptmin and x.pt < ptmax
                                   and abs(x.pseudorap) < etacut) ]

        ptarray = numpy.array([ x.pt for x in filtered_particles ])
        phiarray = numpy.array([ x.phi for x in filtered_particles ])

        sin2array = numpy.sin(2 * phiarray)
        cos2array = numpy.cos(2 * phiarray)
        sin3array = numpy.sin(3 * phiarray)
        cos3array = numpy.cos(3 * phiarray)
        ptphisin2 = (numpy.multiply(ptarray, sin2array))
        ptphicos2 = (numpy.multiply(ptarray, cos2array))
        ptphisin3 = (numpy.multiply(ptarray, sin3array))
        ptphicos3 = (numpy.multiply(ptarray, cos3array))

        sin2subsum = [0.0] * 2
        cos2subsum = [0.0] * 2
        sin3subsum = [0.0] * 2
        cos3subsum = [0.0] * 2
        nsub = [0] * 2

        chgphis = []
        sin2mean = []
        cos2mean = []
        sin3mean = []
        cos3mean = []

        npart = (len(filtered_particles) - 1)
        for i in range (0, len(filtered_particles)):
            if (i%2 == 0):
                subgroup = 1
            else:
                subgroup = 0
            if (filtered_particles[i].charge != 0):
                ncharges += 1
                chgphis.append(filtered_particles[i].phi)

                nsub[subgroup] += 1
                sin2subsum[subgroup] += ptphisin2[i]
                cos2subsum[subgroup] += ptphicos2[i]
                sin3subsum[subgroup] += ptphisin3[i]
                cos3subsum[subgroup] += ptphicos3[i]

                # exclude the analysed particle from event plane definition
                # to avoid autocorrelations

                sin2mean.append((ptphisin2.sum() - ptphisin2[i]) / npart)
                cos2mean.append((ptphicos2.sum() - ptphicos2[i]) / npart)
                sin3mean.append((ptphisin3.sum() - ptphisin3[i]) / npart)
                cos3mean.append((ptphicos3.sum() - ptphicos3[i]) / npart)

        # full event value
        chgphiarray = numpy.array(chgphis)
        sin2meanarray = numpy.array(sin2mean)
        cos2meanarray = numpy.array(cos2mean)
        psitwoarray = numpy.arctan2(sin2meanarray, cos2meanarray) / 2
        sin3meanarray = numpy.array(sin3mean)
        cos3meanarray = numpy.array(cos3mean)
        psithreearray = numpy.arctan2(sin3meanarray, cos3meanarray) / 3
        v2eventsum = numpy.sum(numpy.cos(2 * (chgphiarray - psitwoarray)))
        v3eventsum = numpy.sum(numpy.cos(3 * (chgphiarray - psithreearray)))

        if (ncharges > 0):
            vn_event_sums[0] += v2eventsum / ncharges
            vn_event_sums[1] += (v2eventsum / ncharges)**2
            vn_event_sums[2] += v3eventsum / ncharges
            vn_event_sums[3] += (v3eventsum / ncharges)**2

        psi2sub = [0.0] * 2
        psi3sub = [0.0] * 2
        for i in range(0, 2):
            if (nsub[i] > 0):
                sin2submean = sin2subsum[i] / nsub[i]
                cos2submean = cos2subsum[i] / nsub[i]
                psi2sub[i] = math.atan2(sin2submean, cos2submean) / 2
                sin3submean = sin3subsum[i] / nsub[i]
                cos3submean = cos3subsum[i] / nsub[i]
                psi3sub[i] = math.atan2(sin3submean, cos3submean) / 3

        psi2term = math.cos(2 * (psi2sub[0] - psi2sub[1]))
        psi3term = math.cos(3 * (psi3sub[0] - psi3sub[1]))
        vn_event_sums[4] += psi2term
        vn_event_sums[5] += psi2term**2
        vn_event_sums[6] += psi3term
        vn_event_sums[7] += psi3term**2



def v2v3mean(vn_event_sums, nevents):

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

        meanv2sub = math.sqrt(vn_event_sums[4] / nevents)
        meanv3sub = math.sqrt(vn_event_sums[6] / nevents)
        print "Rsub2:", meanv2sub, "Rsub3:", meanv3sub
        meanv2subsqr = math.sqrt(vn_event_sums[5] / nevents)
        meanv3subsqr = math.sqrt(vn_event_sums[7] / nevents)

        v2corr = 0.0
        v2err = 0.0
        v3corr = 0.0
        v3err = 0.0

        print "R2_full chi2 R3_full chi3"
        v2 = True
        for vnsub in [meanv2sub, meanv3sub]:
            chisub = 2.0
            delta = 8.0
            for iteration in range(0, 20):
                arg = chisub**2 / 4
                result = (math.sqrt(math.pi / 2) / 2 * chisub
                          * math.exp(-arg) * (bessel.i0(arg) + bessel.i1(arg)))
                if (result < vnsub):
                    chisub += delta
                elif (result > vnsub):
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
                v2corr = meanv2 / reso_corr
                v2err = errv2 / reso_corr
                v2 = False
            else:
                v3corr = meanv3 / reso_corr
                v3err = errv3 / reso_corr

        print ""
        print "Corrected v2 err v3 err"
        print v2corr, v2err, v3corr, v3err

    else:
        print "No events found for event plane v_n determination."
