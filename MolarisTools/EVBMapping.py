#-------------------------------------------------------------------------------
# . File      : EVBMapping.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2017)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
"""EVBMapping is a class that finds the extrema (RS, transition state, PS) on a free energy profile."""

import exceptions
import numpy
import sys
import os

DEFAULT_INPUT   =  "dG_dE.graph"
DEFAULT_OUTPUT  =  "check.dat"
DEFAULT_DEGREE  =  6
DEFAULT_REF     =  1
DEFAULT_X       =  79
DEFAULT_Y       =  30
IMAG_TOL        =  0.0001


class EVBMapping (object):
    """A class for the mapping of EVB."""

    def __init__ (self, inputfile=DEFAULT_INPUT):
        """Constructor."""
        if not os.path.exists (inputfile):
            raise exceptions.StandardError ("File %s not found" % inputfile)
        # . Load the mapping data
        lines = open (inputfile).readlines ()
        Gs    = []
        gaps  = []
        for line in lines[1:]:
            tokens = line.split ()
            if not tokens:
                break
            gap, G = map (float, tokens)
            Gs.append   (G)
            gaps.append (gap)
        # . Save relevant variables
        self.Gs            =  Gs
        self.gaps          =  gaps
        self.isLoaded      =  True
        self.isWritten     =  False
        self.isCalculated  =  False


    def CalculateExtrema (self, polynomialDegree=DEFAULT_DEGREE, reference=DEFAULT_REF):
        """Find minima and maxima on the energy profile."""
        if self.isLoaded:
            # . Fit a polynomial function to the data
            coeff = numpy.polyfit (numpy.array (self.gaps), numpy.array (self.Gs), polynomialDegree)
            # . Find a derivative function of the polynomial
            deriv = numpy.polyder (coeff)
            # . Find the roots of the polynomial function
            roots = numpy.roots (deriv)
            # . Select the roots between RS and PS
            start, end   = min (self.gaps), max (self.gaps)
            rootsOK      = filter (lambda r: r >= start and r <= end, roots)
            rootsOK.sort ()
            # . Calculate the values of G at the found extrema
            extrema = []
            for root in rootsOK:
                extrema.append (numpy.polyval (coeff, root))
            # . Select a reference extremum
            Gmin = extrema[reference - 1]
            # . Save relevant variables
            self.coeff        = coeff
            self.deriv        = deriv
            self.rootsOK      = rootsOK
            self.Gmin         = Gmin
            self.isCalculated = True


    def _CheckImag (self, value, label=""):
        """Make sure that the imaginary component is zero."""
        if abs (value.imag) > IMAG_TOL:
            raise exceptions.StandardError (("%s has a nonzero imaginary component: %f." % (label, value.imag)) if label else "Nonzero imaginary component: %f." % value.imag)
        return value.real


    def WriteGnuplotData (self, filename=DEFAULT_OUTPUT):
        """Write a file to view in Gnuplot."""
        if self.isCalculated:
            f     = open (filename, "w")
            f.write ("# %8s %14s %14s %14s %14s\n" % ("gap", "G", "Gfit", "Grel", "dG"))
    
            # . For each gap, write G, G from the fitted polynomial, G from the polynomial relative to Gmin, the value of the derivative function
            for gap, G in zip (self.gaps, self.Gs):
                Gfit   = numpy.polyval (self.coeff, gap)
                Grel   = Gfit - self.Gmin
                Gderiv = numpy.polyval (self.deriv, gap)
                f.write ("%14.3f %14.3f %14.3f %14.3f %14.3f\n" % (gap, G, self._CheckImag (Gfit, label="Gfit"), self._CheckImag (Grel, label="Grel"), self._CheckImag (Gderiv, label="Gderiv")))
            f.close ()
            self.filename  = filename
            self.isWritten = True


    def Plot (self, sizex=DEFAULT_X, sizey=DEFAULT_Y):
        """Print a plot with Gnuplot."""
        if self.isWritten:
            os.system ("gnuplot -e \"set term dumb %d %d ; set title 'dG (Egap)' ; plot '%s' u 1:2 w lp\"" % (sizex, sizey, self.filename))


    def Summary (self):
        """Print a summary."""
        if self.isCalculated:
            print ("---- Found extrema ----".center (57))
            for serial, root in enumerate (self.rootsOK, 1):
                G = numpy.polyval (self.coeff, root)
                print ("%-2d: gap = %8.2f      G = %8.2f      Grel = %8.2f" % (serial, self._CheckImag (root, label="root"), self._CheckImag (G, label="G"), self._CheckImag (G - self.Gmin, label="Grel")))


    def DoAll (self, polynomialDegree=DEFAULT_DEGREE, reference=DEFAULT_REF, filename=DEFAULT_OUTPUT, sizex=DEFAULT_X, sizey=DEFAULT_Y):
        """Execute all tasks."""
        # . Calculate the extrema
        self.CalculateExtrema (polynomialDegree=polynomialDegree, reference=reference)
        # . Write a file for Gnuplot
        self.WriteGnuplotData (filename=filename)
        # . Print a plot
        self.Plot (sizex=sizex, sizey=sizey)
        # . Print a summary
        self.Summary ()


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
