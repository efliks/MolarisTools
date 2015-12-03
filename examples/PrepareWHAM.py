#!/usr/bin/env python

import sys, os.path
sys.path.append ("/home/mikolaj/devel/MolarisTools")

from MolarisTools import MolarisInputFile, DistanceFile


# . Collect equilibrium distances
nwindows = 21
equil    = []
pair     = (2, 6)

for wcount in range (1, nwindows + 1):
    inputfile = MolarisInputFile (os.path.join ("..", "w%d" % wcount, "window_%d.inp" % wcount), pair=pair)
    forceConst, equilDist = inputfile.equil
    equil.append (equilDist)

# . Collect distances from simulations
simul    = []
for wcount in range (1, nwindows + 1):
    distfile = DistanceFile (os.path.join ("..", "w%d" % wcount, "window_%d" % wcount, "dist.dat"))
    data = distfile.pairs[pair]
    simul.append (data)

# . Calculate average distances from simulations
aver     = []
for series in simul:
    average = sum (series) / len (series)
    aver.append (average)

# . Write meta file
metafile = open ("wham.meta", "w")
for wcount in range (1, nwindows + 1):
    metafile.write ("wham%02d.ts    %8.3f    200.000\n" % (wcount, aver[wcount - 1]))
metafile.close ()

# . Write timeseries files
for wcount in range (1, nwindows + 1):
    timeseries = open ("wham%02d.ts" % wcount, "w")
    series = simul[wcount - 1]
    for si, s in enumerate (series, 1):
        timeseries.write ("%4d    %8.3f\n" % (si, s))
    timeseries.close ()
