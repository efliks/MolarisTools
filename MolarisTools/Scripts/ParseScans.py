#-------------------------------------------------------------------------------
# . File      : ParseScans.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2017)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import os, math, glob

from MolarisTools.Parser  import MolarisOutputFile


def ParsePESScan (pattern="evb_scan_", filenameTotal="total_e.dat", filenameTotalRelative="total_e_rel.dat", patternChanges="changes_", logging=True):
    """Parse a potential energy surface scan."""
    steps = []
    files = glob.glob ("%s*.out" % pattern)
    files.sort ()
    for fn in files:
        mof = MolarisOutputFile (filename=fn, logging=False)
        if logging:
            print ("# . Parsing file %s ..." % fn)
        collect = []
        for step in mof.qmmmComponentsI:
            Eqmmm = step.Eqmmm
            collect.append (Eqmmm)
        if logging:
            nsteps = len (collect)
            print ("# . Found %d steps" % nsteps)
        steps.append (collect)
    if logging:
        ntotal = 0
        for step in steps:
            ntotal += len (step)
        print ("Found %d total steps in %d files" % (ntotal, len (files)))
    base  = steps[0][-1]
    m, n  = (0, 0)
    baseRel, topRel = (base, base)
    for (i, step) in enumerate (steps, 1):
        if step[-1] < baseRel:
            (baseRel, m) = (step[-1], i)
        if step[-1] > topRel:
            (topRel,  n) = (step[-1], i)
    barrier = topRel - baseRel
    if logging:
        print ("Found a barrier between points (%d, %d) of %.1f kcal/mol" % (m, n, barrier))
    for (filename, subtract) in ((filenameTotal, base), (filenameTotalRelative, baseRel)):
        fo   = open (filename, "w")
        for (i, step) in enumerate (steps, 1):
            last = step[-1]
            fo.write ("%d   %f\n" % (i, last - subtract))
        fo.close ()
        if logging:
            print ("Wrote file %s" % filename)
    for (i, step) in enumerate (steps, 1):
        filename = "%s%03d.dat" % (patternChanges, i)
        fo       = open (filename, "w")
        previous = step[0]
        for iteration in step[1:]:
            change   = (iteration - previous)
            previous = iteration
            fo.write ("%f\n" % change)
        fo.close ()
        if logging:
            print ("Wrote file %s" % filename)


def ParsePESScan2D (pattern="evb_scan_", filenameTotal="total_e.dat", filenameTotalRelative="total_e_rel.dat", useDistances=False, zigzag=False, logging=True):
    """Parse a two-dimensional potential energy surface scan."""
    files  = glob.glob ("%s*.inp" % pattern)
    nfiles = len (files)
    size   = int (math.sqrt (nfiles))
    base   = 0.
    found  = False
    for i in range (1, size + 1):
        for j in range (1, size + 1):
            filename = "%s%02d_%02d.out" % (pattern, 1, 1)
            if os.path.exists (filename):
                if logging:
                    print ("# . Parsing file %s ..." % filename)
                mof   = MolarisOutputFile (filename=filename, logging=False)
                last  = mof.qmmmComponentsI[-1]
                base  = last.Eqmmm
                found = True
                break
        if found:
            break
    if logging:
        print ("# . Base energy is %f" % base)

    rows   = []
    pairs  = []
    nlogs  = 0
    for i in range (1, size + 1):
        columns = []
        gather  = []
        for j in range (1, size + 1):
            je = j
            if zigzag:
                # . In a zig-zag scan, the second coordinate 
                # . increases and decreases for odd and even iterations 
                # . of the first coordinate, respectively.
                if ((i % 2) == 0):
                    je = size - j + 1
            filename = "%s%02d_%02d.out" % (pattern, i, je)
            Eqmmm    = base
            if os.path.exists (filename):
                if logging:
                    print ("# . Parsing file %s ..." % filename)
                mof   = MolarisOutputFile (filename=filename, logging=False)
                last  = mof.qmmmComponentsI[-1]
                Eqmmm = last.Eqmmm
                if logging:
                    nsteps = len (mof.qmmmComponentsI)
                    print ("# . Found %d steps" % nsteps)
                nlogs += 1
            columns.append (Eqmmm)
            if useDistances:
                filename   = "%s%02d_%02d.inp" % (pattern, i, je)
                molaris    = MolarisInputFile (filename, logging=False)
                (coi, coj) = molaris.constrainedPairs[:2]
                (ri, rj)   = (coi.req, coj.req)
                pair = (ri, rj)
                gather.append (pair)
        rows.append (columns)
        if useDistances:
            pairs.append (gather)
    if logging:
        print ("Parsed %d log files out of %d calculations" % (nlogs, nfiles))

    low, high = (base, base)
    for (i, row) in enumerate (rows):
        for (j, column) in enumerate (row):
            if (column < low):
                low      = column
                (li, lj) = (i, j)
            if (column > high):
                high     = column
                (hi, hj) = (i, j)
    if logging:
        if useDistances:
            (ra, rb) = pairs[li][lj]
            print ("Lowest  point at (%.2f, %.2f) with energy of %.2f kcals" % (ra, rb, low ))
            (ra, rb) = pairs[hi][hj]
            print ("Highest point at (%.2f, %.2f) with energy of %.2f kcals" % (ra, rb, high))
        else:
            print ("Lowest  point at (%d, %d) with energy of %.2f kcals" % (li + 1, lj + 1, low ))
            print ("Highest point at (%d, %d) with energy of %.2f kcals" % (hi + 1, hj + 1, high))
        difference = high - low
        print ("Energy difference is %.1f kcals" % difference)

    for (filename, subtract) in ((filenameTotal, 0.), (filenameTotalRelative, low)):
        fo = open (filename, "w")
        for (i, row) in enumerate (rows):
            for (j, column) in enumerate (row):
                if useDistances:
                    (idist, jdist) = pairs[i][j]
                    line = "%6.2f   %6.2f   %f" % (idist, jdist, rows[i][j] - subtract)
                else:
                    line = "%2d   %2d   %f" % (j + 1, i + 1, rows[i][j] - subtract)
                fo.write ("%s\n" % line)
            fo.write ("\n")
        fo.close ()
        if logging:
            print ("Wrote file %s" % filename)

    nrows = len (rows)
    for irow in range (1, nrows - 2):
        ncolumns = len (rows[irow])
        for jcolumn in range (1, ncolumns - 2):
            up     = rows[irow - 1][jcolumn    ]
            down   = rows[irow + 1][jcolumn    ]
            left   = rows[irow    ][jcolumn - 1]
            right  = rows[irow    ][jcolumn + 1]
            center = rows[irow    ][jcolumn    ]
            checks = (center < up, center < down, center < left, center < right)
            if all (checks):
                if logging:
                    if useDistances:
                        (idist, jdist) = pairs[irow][jcolumn]
                        print ("Found minimum at (%.2f, %.2f) with energy of %.2f kcals" % (idist, jdist, center - low))
                    else:
                        print ("Found minimum at (%d, %d) with energy of %.2f kcals" % (irow, jcolumn, center - low))
            checks = (center > up, center > down, center > left, center > right)
            if all (checks):
                if logging:
                    if useDistances:
                        (idist, jdist) = pairs[irow][jcolumn]
                        print ("Found maximum at (%.2f, %.2f) with energy of %.2f kcals" % (idist, jdist, center - low))
                    else:
                        print ("Found maximum at (%d, %d) with energy of %.2f kcals" % (irow, jcolumn, center - low))


#===============================================================================
# . Main program
#===============================================================================
if (__name__ == "__main__"): pass
