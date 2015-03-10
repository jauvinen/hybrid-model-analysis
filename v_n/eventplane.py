#!/usr/bin/python

# Functions for computing event plane v2 and v3

import math
import scipy.special as bessel

# Compute single-event vn quantities
# Input: List of particles, containing their momentum and azimuthal angle
# Return: Tuple of <v2>, <v2^2>, <v3>, <v3^2>, and sub-event versions
def v2v3(particlelist):
    ncharges = 0
    v2eventsum = 0
    v3eventsum = 0
    v2eventmean = 0
    v2eventmeansqr = 0
    v3eventmean = 0
    v3eventmeansqr = 0

    if (len(particlelist) > 1):
        sin2subsum = [0.0] * 2
        cos2subsum = [0.0] * 2
        sin3subsum = [0.0] * 2
        cos3subsum = [0.0] * 2

        for i in range (0, len(particlelist)):
            if (i%2 == 0):
                subgroup = 1
            else:
                subgroup = 0
            if (particlelist[i].charge != 0
                and abs(particlelist[i].pseudorap) < 1):

                ncharges += 1
                sin2sum = 0
                cos2sum = 0
                sin3sum = 0
                cos3sum = 0

                nparticles[subgroup] += 1
                sin2subsum[subgroup] += pt_j * math.sin(2 * phi_j)
                cos2subsum[subgroup] += pt_j * math.cos(2 * phi_j)
                sin3subsum[subgroup] += pt_j * math.sin(3 * phi_j)
                cos3subsum[subgroup] += pt_j * math.cos(3 * phi_j)

                # vn for particle i
                for j in range (0, len(particlelist)):
                    # exclude the analysed particle from event plane definition
                    # to avoid autocorrelations
                    if (j != i):

                        pt_j = particlelist[j].pt
                        phi_j = particlelist[j].phi
                        sin2sum += pt_j * math.sin(2 * phi_j)
                        cos2sum += pt_j * math.cos(2 * phi_j)
                        sin3sum += pt_j * math.sin(3 * phi_j)
                        cos3sum += pt_j * math.cos(3 * phi_j)


                # full event value
                normalization = (len(particlelist) - 1)
                sin2mean = sin2sum / normalization
                cos2mean = cos2sum / normalization
                psitwo = math.atan2(sin2mean, cos2mean) / 2
                sin3mean = sin3sum / normalization
                cos3mean = cos3sum / normalization
                psithree = math.atan2(sin3mean, cos3mean) / 3
                v2eventsum += math.cos(2 * (particlelist[i].phi - psitwo))
                v3eventsum += math.cos(3 * (particlelist[i].phi - psithree))

                # subevent value
                norm_sub = nparticles[subgroup]
                sin2mean[subgroup] = sin2sum[subgroup] / norm_sub
                cos2mean[subgroup] = cos2sum[subgroup] / norm_sub
                psi2[subgroup] = math.atan2(sin2mean[subgroup],
                                            cos2mean[subgroup]) / 2
                sin3mean[subgroup] = sin3sum[subgroup] / norm_sub
                cos3mean[subgroup] = cos3sum[subgroup] / norm_sub
                psi3[subgroup] = math.atan2(sin3mean[subgroup],
                                            cos3mean[subgroup]) / 3

                psi2term = math.cos(2 * (psi2[0] - psi2[1]))
                psi3term = math.cos(3 * (psi3[0] - psi3[1]))
                v2subsum += psi2term
                v3subsum += psi3term
                v2subsumsqr += psi2term**2
                v3subsumsqr += psi3term**2

        if (ncharges > 0):
            v2eventmean = v2eventsum / ncharges
            v2eventmeansqr = (v2eventsum / ncharges)**2
            v3eventmean = v3eventsum / ncharges
            v3eventmeansqr = (v3eventsum / ncharges)**2

    return (v2eventmean, v2eventmeansqr, v3eventmean, v3eventmeansqr,
            v2subsum, v2subsumsqr, v3)


def v2v3mean(vneventslist):
    meanv2 = 0.0
    meanv2sqr = 0.0
    meanv3 = 0.0
    meanv3sqr = 0.0
    for vntuple in vneventslist:
        meanv2 += vntuple[0]
        meanv2sqr += vntuple[1]
        meanv3 += vntuple[2]
        meanv3sqr += vntuple[3]

    meanv2 /= len(vneventslist)
    meanv2sqr /= len(vneventslist)
    meanv3 /= len(vneventslist)
    meanv3sqr /= len(vneventslist)

    errv2 = math.sqrt((meanv2sqr - meanv2**2) / (len(vneventslist) - 1))
    errv3 = math.sqrt((meanv3sqr - meanv3**2) / (len(vneventslist) - 1))

    return (meanv2, errv2, meanv3, errv3)
