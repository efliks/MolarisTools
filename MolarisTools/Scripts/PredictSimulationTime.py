#-------------------------------------------------------------------------------
# . File      : PredictSimulationTime.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2018)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import os, glob, time, datetime


def PredictSimulationTime (pattern="evb_*out"):
    """Calculate remaining time for a simulation consisting of multiple files."""
    files = glob.glob (pattern)
    files.sort ()
    nfiles = len (files)
    print ("Found %d output files." % nfiles)
    pattern = "%a %b %d %H:%M:%S %Y"
    
    times = []
    for (i, log) in enumerate (files):
        mtime = int (os.path.getmtime (log))
        times.append (mtime)
    
    #    lines = open (log).readlines ()
    #    for line in lines:
    #        if line.startswith ("  Job finishing time:"):
    #            break
    #    string = " ".join (line.split ()[3:])
    #    epoch = int (time.mktime (time.strptime (string, pattern)))
    #    if (i < 1):
    #        print ("File %s: converting %s -> %d" % (log, string, epoch))
    #    else:
    #        seconds = epoch - times[i - 1]
    #        delta = str (datetime.timedelta (seconds=seconds))
    #        print ("File %s: converting %s -> %d  (%s)" % (log, string, epoch, delta))
    #    times.append (epoch)
    gradients = []
    for (i, t) in enumerate (times[1:], 1):
        grad = (t - times[i - 1])
        gradients.append (grad)
    
    gradient = sum (gradients) / len (gradients)
    intercept = times[0]
    print ("Gradient = %d, intercept = %d" % (gradient, intercept))
    average = str (datetime.timedelta (seconds=gradient))
    print ("Average time per file is %s" % average)
    
    inputs = glob.glob ("evb_*inp")
    inputs.sort ()
    nleft = len (inputs) - nfiles
    
    lines  = open (inputs[0]).readlines ()
    nsteps = 0
    for line in lines:
        # if line.count ("nsteps"):
        tokens = line.split ()
        if len (tokens) >= 2:
            if (tokens[0] == "nsteps"):
                nsteps = int (line.split ()[1])
                break
    if (nsteps > 0):
        totalSteps = float (nfiles * nsteps)
        delta = float (times[-1] - times[0])
        perTime = totalSteps / delta * 3600.
        print ("Simulation was running at the rate of %d steps per hour." % perTime)
        perStep = delta / totalSteps
        print ("One step takes %d seconds." % perStep)
    
    
    if (nleft > 0):
        print ("There are %d files left." % nleft)
        
        ends = times[-1] + gradient * nleft
        final = datetime.datetime.fromtimestamp (ends)
        fmt = final.strftime (pattern)
        
        now = int (time.time ())
        seconds = ends - now
        delta = str (datetime.timedelta (seconds=seconds))
        print ("Job will end on: %s (%s from now)" % (fmt, delta))


#===============================================================================
# . Main program
#===============================================================================
if (__name__ == "__main__"): pass
