#!/usr/bin/python

import os


lines = open ("equil_template.inp").readlines ()

# range (10, 310, 10)
for i, T in enumerate (range (1, 21), 1):
    output = []
    for line in lines:
        if line.count ("@PREVNUM@"):
            if i > 1:
                line = line.replace ("@PREVNUM@", "%03d" % (i - 1))
            else:
                if os.path.exists ("evb_step30.res"):
                    line = "%s  ./evb_step30.res\n" % line[:line.find ("rest_in") + len ("rest_in")]
                else:
                    line = ""
        elif line.count ("@NUM@"):
            line = line.replace ("@NUM@", "%03d" % i)
        elif line.count ("@TEMP@"):
            line = line.replace ("@TEMP@", "%.1f" % T)
        output.append (line)

    towrite = open ("evb_equil_%03d.inp" % i, "w")
    towrite.writelines (output)
    towrite.close ()
    print ("Wrote input file %d" % i)
