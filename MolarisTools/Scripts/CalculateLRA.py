#-------------------------------------------------------------------------------
# . File      : CalculateLRA.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2018)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import os, glob

from MolarisTools.Parser  import GapFile, GapFileEVB


def CalculateLRA (patha="lra_RS", pathb="lra_RS_qmmm", logging=True, verbose=False, skip=None, trim=None, returnTerms=False, gapFormat="QM"):
    """Calculate LRA for two endpoint simulations."""
    points = []
    for path in (patha, pathb):
        gapfiles = []
        logs     = glob.glob (os.path.join (path, "evb_equil_*out"))
        logs.sort ()
        for log in logs:
            (logfile, logext) = os.path.splitext (log)
            gapfile = os.path.join (logfile, "gap.out")
            gapfiles.append (gapfile)
        if logging:
            ngap = len (gapfiles)
            print ("# . Found %d gap files at location %s" % (ngap, path))
        points.append (gapfiles)
    (filesa, filesb) = points
    lena, lenb = map (len, points)
    ngap = lena
    if lena > lenb:
        ngap = lenb
    if logging:
        print ("# . Using %d gap files" % ngap)
    points = []
    for gapfiles in (filesa, filesb):
        if (gapFormat == "QM"):
            gap = GapFile (gapfiles[0], logging=(True if (logging and verbose) else False))
        else:
            gap = GapFileEVB (gapfiles[0], logging=(True if (logging and verbose) else False))
        for nextfile in gapfiles[1:ngap]:
            gap.Extend (nextfile)
        points.append (gap)
    (gapa, gapb) = points
    if logging:
        print ("# . Number of steps in each endpoint is (%d, %d)" % (gapa.nsteps, gapb.nsteps))
        if   isinstance (skip, int):
            print ("# . Skipping first %d configurations" % skip)
        elif isinstance (skip, float):
            percent  = int (skip * 100.)
            (na, nb) = int (gapa.nsteps * skip), int (gapb.nsteps * skip)
            print ("# . Skipping first %d%% (%d, %d) of configurations" % (percent, na, nb))
        if   isinstance (trim, int):
            print ("# . Skipping last %d configurations" % trim)
        elif isinstance (trim, float):
            percent  = int (trim * 100.)
            (na, nb) = int (gapa.nsteps * trim), int (gapb.nsteps * trim)
            print ("# . Skipping last %d%% (%d, %d) of configurations" % (percent, na, nb))
    a = gapa.CalculateLRATerm (skip=skip, trim=trim)
    b = gapb.CalculateLRATerm (skip=skip, trim=trim)
    lra = .5 * (a + b)
    if logging:
        print ("# . Calculated LRA = %f  (%f, %f)" % (lra, a, b))
    if returnTerms:
        return (lra, a, b)
    return lra


def CalculateOneSidedLRA (path="lra_RS_qmmm", logging=True, verbose=False, skip=None, trim=None, gapFormat="QM"):
    """Calculate LRA for one endpoint simulation."""
    gapfiles = []
    logs     = glob.glob (os.path.join (path, "evb_equil_*out"))
    logs.sort ()
    for log in logs:
        (logfile, logext) = os.path.splitext (log)
        gapfile = os.path.join (logfile, "gap.out")
        gapfiles.append (gapfile)
    ngap = len (gapfiles)
    if logging:
        print ("# . Found %d gap files at location %s" % (ngap, path))
    if logging:
        print ("# . Using %d gap files" % ngap)
    if (gapFormat == "QM"):
        gapa = GapFile (gapfiles[0], logging=(True if (logging and verbose) else False))
    else:
        gapa = GapFileEVB (gapfiles[0], logging=(True if (logging and verbose) else False))
    for nextfile in gapfiles[1:ngap]:
        gapa.Extend (nextfile)
    if logging:
        print ("# . Number of steps in endpoint is %d" % gapa.nsteps)
        if   isinstance (skip, int):
            print ("# . Skipping first %d configurations" % skip)
        elif isinstance (skip, float):
            percent  = int (skip * 100.)
            na       = int (gapa.nsteps * skip)
            print ("# . Skipping first %d%% (%d) of configurations" % (percent, na))
        if   isinstance (trim, int):
            print ("# . Skipping last %d configurations" % trim)
        elif isinstance (trim, float):
            percent  = int (trim * 100.)
            na       = int (gapa.nsteps * trim)
            print ("# . Skipping last %d%% (%d) of configurations" % (percent, na))
    lra = gapa.CalculateLRATerm (skip=skip, trim=trim)
    if logging:
        print ("# . Calculated LRA = %f" % lra)
    return lra


#===============================================================================
# . Main program
#===============================================================================
if (__name__ == "__main__"): pass
