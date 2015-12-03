#!/usr/bin/env python

import sys, os.path
sys.path.append ("/home/mikolaj/devel/MolarisTools")

from MolarisTools import CHELPGCharges


argv       =  sys.argv
groups     =  []
extraCol   =  False
extraLine  =  False
if len (argv) < 2:
    print ("Usage: %s gaussian_or_orca.out [--average 1,2,3] [--merge 4,5] [--group 6,7,8,9,-1] [--fix 1,-0.5] [--column] [--line]" % os.path.basename (argv[0]))
    sys.exit ()


# . Load the log file
charges = CHELPGCharges (argv[1])

if len (argv) > 2:
    prev = ""
    for arg in argv[2:]:
        if   arg in ("-c", "--column"):
            extraCol    = True
        elif arg in ("-l", "--line"):
            extraLine   = True
        elif arg in ("--average", "-a", "--merge", "-m", "--group", "-g", "--fix", "-f"):
            prev = arg
        else:
            if not prev:
                raise exceptions.StandardError ("Unrecognized option: %s" % arg)
            tokens = arg.split (",")
            # . Convert tokens
            items = map (lambda token: int (token) if token.isdigit () else float (token), tokens)
            # items  = []
            # for token in tokens:
            #     if token.isdigit ():
            #         items.append (int (token))
            #     else:
            #         items.append (float (token))

            if   prev in ("--average", "-a"):
                # . Average charges
                charges.AverageCharges (items)
            elif prev in ("--merge", "-m"):
                # . Merge charges
                charges.MergeCharges (items)
            elif prev in ("--group", "-g"):
                # . Create a group of charges (currently, only one group is possible)
                atoms = items[:-1]
                charges.GroupCharges (atoms, items[-1], groups)
                groups.extend (atoms)
            elif prev in ("--fix", "-f"):
                # . Fix charges of selected atoms to predefined values
                charges.FixCharge (items[0], items[1])

# . Write out the final table
charges.WriteGeometryCharges (extraCol=extraCol, extraLine=extraLine, line="# %12s" % os.path.basename (sys.argv[1]))
