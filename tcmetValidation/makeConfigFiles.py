#!/usr/bin/python

import readline
import os
import commands
import re
import optparse
import sys

parser = optparse.OptionParser ("Usage: %prog [options] dataset globaltag key metcollection")
(options, args) = parser.parse_args()
neededArgs = 4

if len(args) < neededArgs:
    parser.error("Must provide 4 arguments")
    sys.exit(0)
                
dataset         = args[0]
globaltag       = args[1]
mykey           = args[2]
metcoll         = args[3]

print '-------------------------------------------------------------------------------------'
print 'Making cfg file for dataset    :  ' + dataset
print 'Using global tag               :  ' + globaltag
print 'Writing cfg file               :  ' + mykey + '.py'
print 'MET collection                 :  ' + metcoll

slashRE = re.compile(r'^/')
first = True

listOfFiles = ''
file = ''

for line in commands.getoutput('dbsql find file where dataset like %s' % dataset).split('\n'):
    if slashRE.search(line):
        line = line.strip()
        if first:
            listOfFiles += "'" + line + "'"
            first = False
            file = line
        else:
            listOfFiles += ",\n    '" + line + "'"

print 'List of LFNs:'
print listOfFiles

#print 'SET NEVENTS = 100!!!'
vars = { 'MYFILENAMES' : listOfFiles, 'MYNEVENTS' : -1, 'MYGLOBALTAG' : globaltag, 'MYROOTFILE' : mykey + '.root', 'MYSAMPLE' : mykey , 'METCOLLECTION' : metcoll }

template = open('template.py','r')
templatestring = ''
for line in template:
    templatestring += line

template.close()

target = open(mykey + '.py','w')
target.write(templatestring % vars)
target.close()


