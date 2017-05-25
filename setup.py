#!/usr/bin/python

from distutils.core import setup


setup (
    name          =  "MolarisTools" ,
    version       =  "0.1" ,
    description   =  "A library to automate routine tasks in Molaris-XG." ,
    author        =  "Mikolaj Feliks" ,
    author_email  =  "mikolaj.feliks@gmail.com" ,
    license       =  "GPLv3" ,
    url           =  "https://github.com/mfx9/MolarisTools" ,
    packages      =  ["MolarisTools", ] ,
    package_data  =  {"MolarisTools" : ["data/*", ], } ,
)
