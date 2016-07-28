#!/usr/bin/python

import sys, os, socket, collections, threading
sys.path.append (os.path.join (os.environ["HOME"], "devel", "MolarisTools"))

from MolarisTools import QMCallerGaussian, MolarisAtomsFile


charge      = -1
method      = "B3LYP/6-31G*"
parallelJob = True

states  = (
    ("_I",  "mol.in1",  "d.o1")  ,
    ("_II", "mol.in2",  "d.o2")  ,
    )
nstates = len (states)


#-------------------------------------------------------------------------------
class GaussianThread (threading.Thread):
    """Run a Gaussian job as a thread."""

    def __init__ (self, qmCaller):
        """Constructor."""
        threading.Thread.__init__ (self)
        if not isinstance (qmCaller, QMCallerGaussian):
            raise exceptions.StandardError ("Not a Gaussian caller.")
        self.qmCaller = qmCaller

    def run (self):
        self.qmCaller.Run ()
#-------------------------------------------------------------------------------


# . Determine which machine it is
Machine = collections.namedtuple ("Machine", "label  ncpu  memory  executable  scratch")

convert = {
            "eigen" : Machine (
                label       =   "eigen"  ,
                ncpu        =   1        ,
                memory      =   1        ,
                executable  =   os.path.join (os.environ["GAUSS_EXEDIR"], "g03") ,
                scratch     =   "."      ,
                ),

            "comet" : Machine (
                label       =   "comet"  ,
                ncpu        =   8        ,
                memory      =   8        ,
                executable  =   os.path.join (os.environ["GAUSS_EXEDIR"], "g09") ,
                scratch     =   os.environ["SCRDIR"] ,
                ),
}
hostname = socket.gethostname ()[:5]
machine  = convert[hostname]


# . Parallel job
if parallelJob:
    threads = []
    for (stateLabel, stateMolin, stateDo) in states:
        callerGaussian = QMCallerGaussian (
            ncpu      =   machine.ncpu     ,
            memory    =   machine.memory   ,
            archive   =   True             ,
            restart   =   True             ,
            qmmm      =   True             ,
            charge    =   charge           ,
            method    =   method           ,
            pathGaussian            =  machine.executable        ,
            fileAtoms               =  stateMolin                ,
            fileForces              =  stateDo                   ,
            fileGaussianError       =  "job%s.err" % stateLabel  ,
            fileGaussianInput       =  "job%s.inp" % stateLabel  ,
            fileGaussianOutput      =  "job%s.log" % stateLabel  ,
            fileGaussianCheckpoint  =  os.path.join (machine.scratch, "job%s.chk" % stateLabel) ,
            )
        thread = GaussianThread (callerGaussian)
        threads.append (thread)
    
    # . Run Gaussian threads
    for thread in threads: thread.start ()
    
    for thread in threads: thread.join ()


# . Sequential job
else:
        molin = MolarisAtomsFile ()
        (stateLabel, stateMolin, stateDo) = states[molin.stateID - 1]

        callerGaussian = QMCallerGaussian (
            ncpu      =   machine.ncpu     ,
            memory    =   machine.memory   ,
            archive   =   True             ,
            restart   =   True             ,
            qmmm      =   True             ,
            charge    =   charge           ,
            method    =   method           ,
            pathGaussian            =  machine.executable        ,
            fileGaussianError       =  "job%s.err" % stateLabel  ,
            fileGaussianInput       =  "job%s.inp" % stateLabel  ,
            fileGaussianOutput      =  "job%s.log" % stateLabel  ,
            fileGaussianCheckpoint  =  os.path.join (machine.scratch, "job%s.chk" % stateLabel) ,
            )
        callerGaussian.Run ()


# . All done, return to Molaris
