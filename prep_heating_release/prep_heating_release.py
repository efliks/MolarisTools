#!/usr/bin/python
lines = open ("heat_template.inp").readlines ()

for i, T in enumerate (range (10, 310, 10), 1):
    output = []
    for line in lines:
        if line.count ("@PREVNUM@"):
            if i > 1:
                line = line.replace ("@PREVNUM@", "%02d" % (i - 1))
            else:
                line = ""
        elif line.count ("@NUM@"):
            line = line.replace ("@NUM@", "%02d" % i)
        elif line.count ("@TEMP@"):
            line = line.replace ("@TEMP@", "%.1f" % T)
        output.append (line)

    towrite = open ("evb_heat_%02d.inp" % i, "w")
    towrite.writelines (output)
    towrite.close ()
    print ("Wrote input file %d (temperature %dK)" % (i, T))


# . Prepeare files for gradual release of some of the restraints
templ = open ("evb_heat_%02d.inp" % i).readlines ()

for count, (inputFile, prevRestart, currRestart, toUnfreeze) in enumerate ((
    ( "release1.inp" , "evb_heat_%02d/evb_step%02d.res" % (i, i) , "release1.res" , ("H5'1" ,              ) ),
    ( "release2.inp" , "release1/release1.res"                   , "release2.res" , ("H5'1" , "H4'" ,      ) ),
    ( "release3.inp" , "release2/release2.res"                   , "release3.res" , ("H5'1" , "H4'" , "H1'") ),
), 1):
    output = []
    for line in templ:
        if   line.count ("rest_in"):
            testStr = "rest_in"
            line    = "%s    %s\n" % (line[:line.find (testStr) + len (testStr)], prevRestart)

        elif line.count ("rest_out"):
            testStr = "rest_out"
            line    = "%s    %s\n" % (line[:line.find (testStr) + len (testStr)], currRestart)

        elif line.count ("constraint_post"):
            tokens = line.split ()
            for unfreeze in toUnfreeze:
                if tokens[-1] == unfreeze:
                    line = line.replace ("constraint_post", "# constraint_post")
        output.append (line)

    towrite = open (inputFile, "w")
    towrite.writelines (output)
    towrite.close ()
    print ("Wrote release file %d" % count)
