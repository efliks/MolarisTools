#!/usr/bin/python

import numpy, sys


class Distances (object):
    """A class to calculate optimal PA...O3' distances for a PMF simulation by taking average distances from an earlier FEP/US (EVB) simulation."""

    def __init__ (self, filename="dist_O3p_PA.dat", nwindows=11):
        lines = open (filename).readlines ()
        data  = []
        for line in lines:
            tokens = line.split ()
            data.append (float (tokens[1]))
        self.nwindows     = nwindows
        self.evbDistances = data


    def CalculateAverages (self):
        distPerWindow = len (self.evbDistances) / self.nwindows
        collect       = []
        averages      = []
        wcounter      = 0
        for counter, distance in enumerate (self.evbDistances):
            collect.append (distance)
            wcounter += 1
            if wcounter == distPerWindow:
                averageDistance = sum (collect) / distPerWindow
                wcounter        = 0
                collect         = []
                averages.append (averageDistance)
        self.evbAverages = averages


    def FitPolynomial (self, degree=4):
        steps           = range (1, self.nwindows + 1)
        self.polynomial = numpy.polyfit (numpy.array (steps), numpy.array (self.evbAverages), degree)


    def CalculateFitted (self, nuseWindows=0):
        if nuseWindows < 1:
            nuseWindows = self.nwindows
        lambd  = 1.
        delta  = float (self.nwindows - 1) / (nuseWindows - 1)
        fitted = []
        for nw in range (nuseWindows):
            fitted.append (numpy.polyval (self.polynomial, lambd))
            lambd += delta
        self.delta     = delta
        self.evbFitted = fitted


    def PrintAverages (self, filename="dist_from_evb.dat"):
        output = open (filename, "w")
        lambd  = 1.
        for average in self.evbAverages:
            output.write ("%f  %f\n" % (lambd, average))
            lambd += 1.
        output.close ()


    def PrintFitted (self, filename="dist_for_pmf.dat"):
        output = open (filename, "w")
        lambd  = 1.
        for fit in self.evbFitted:
            output.write ("%f  %f\n" % (lambd, fit))
            lambd += self.delta
        output.close ()


#===============================================================================
# . Main program
#===============================================================================
distances = Distances ()
distances.CalculateAverages ()
distances.FitPolynomial ()
distances.CalculateFitted (nuseWindows=31)

distances.PrintAverages ()
distances.PrintFitted ()
