#!/usr/bin/env python

import sys, os.path, exceptions
sys.path.append ("/home/mikolaj/devel/MolarisTools")

from MolarisTools import EVBMapping


reference = 1
if len (sys.argv) >= 2:
    prev = ""
    for arg in sys.argv[1:]:
        if arg in ("-r", "--reference"):
            prev = arg
        else:
            if not prev:
                raise exceptions.StandardError ("Unrecognized option: %s" % arg)
            reference = int (arg)
mapping = EVBMapping ()
mapping.DoAll (reference=reference)
