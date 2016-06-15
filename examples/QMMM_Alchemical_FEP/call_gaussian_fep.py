#!/usr/bin/python

import sys, os, socket, collections
sys.path.append (os.path.join (os.environ["HOME"], "devel", "MolarisTools"))

from MolarisTools import QMCallerGaussian, MolarisAtomsFile


# . Determine which state it is
STATE_I    =   26
STATE_II   =   24
convert    = {STATE_I : "_I" , STATE_II : "_II"}

molin      = MolarisAtomsFile (filename="mol.in")
nqatoms    = len (molin.qatoms)
fileGaussianError       =  "job%s.err" % convert[nqatoms]
fileGaussianInput       =  "job%s.inp" % convert[nqatoms]
fileGaussianOutput      =  "job%s.log" % convert[nqatoms]
fileGaussianCheckpoint  =  "job%s.chk" % convert[nqatoms]


# . Determine which machine it is
Machine = collections.namedtuple ("Machine", "label  ncpu  memory  executable  scratch")

convert = {
            "eigen" : Machine (
                label       =   "eigen"  ,
                ncpu        =   2        ,
                memory      =   2        ,
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


# . Set up a job
callerGaussian = QMCallerGaussian (
    ncpu      =   machine.ncpu     ,
    memory    =   machine.memory   ,
    archive   =   True             ,
    restart   =   True             ,
    qmmm      =   True             ,
    charge    =    -1              ,
    method    =   "B3LYP/6-31G*"   ,
    pathGaussian            =  machine.executable    ,
    fileGaussianError       =  fileGaussianError     ,
    fileGaussianInput       =  fileGaussianInput     ,
    fileGaussianOutput      =  fileGaussianOutput    ,
    fileGaussianCheckpoint  =  os.path.join (machine.scratch, fileGaussianCheckpoint) ,
    )

# . Run it
callerGaussian.Run ()
