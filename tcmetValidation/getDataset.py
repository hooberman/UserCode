#!/usr/bin/python

#-----------------------------------------------------------------------------------
#usage: python getDataset.py DATASETTYPE RELEASE
#This file returns a dataset name of the given type and in the given release
#-----------------------------------------------------------------------------------

import readline
import os
import commands
import re
import optparse
import sys

parser = optparse.OptionParser ("Usage: %prog [options] dataset globaltag key")
(options, args) = parser.parse_args()
neededArgs = 2

if len(args) < neededArgs:
    parser.error("Must provide 2 arguments")
    sys.exit(0)
                
datatype         = args[0]
release          = args[1]

print 'Find dataset of type       :  ' + datatype + ' for release ' + release

slashRE   = re.compile(r'^/')
fastsimRE = re.compile(r'FastSim')
first = True
nfiles = 0

datasetstring = '/RelVal' + datatype + '/' + release + '-*-RECO'

dataset = ''

for line in commands.getoutput('dbsql find dataset where dataset like %s' % datasetstring).split('\n'):
    if slashRE.search(line):
        line = line.strip()
        if not fastsimRE.search(line):
            nfiles = nfiles + 1
            if first:
                dataset = line
                first = False
                

print 'Found ' + str(nfiles) + ' files, choosing    :   ' + dataset                
print dataset

