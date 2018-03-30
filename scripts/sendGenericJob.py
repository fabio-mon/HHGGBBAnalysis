#!/usr/bin/env python
import os, re
import commands
import math, time
import sys
import argparse
import subprocess

parser = argparse.ArgumentParser(description='This script splits CMSSW tasks in multiple jobs and sends them on HERCULES')

parser.add_argument("-l", "--label",          required=True,     type=str,  help="job label")
parser.add_argument("-b", "--baseFolder",     required=True,     type=str,  help="base folder")
parser.add_argument("-e", "--exeName",        required=True,     type=str,  help="absolute path of executable")
parser.add_argument("-I", "--inputFolder",    required=True,     type=str,  help="folder where input files are stored")
parser.add_argument("-i", "--inputFileName",  required=True,     type=str,  help="input file name(s)")
parser.add_argument("-t", "--inputTreeName",  required=True,     type=str,  help="input tree name")
parser.add_argument("-O", "--outputFolder",   required=True,     type=str,  help="folder where to store output files")
parser.add_argument("-o", "--outputFileName", required=True,     type=str,  help="output file name")
parser.add_argument("-c", "--configFile",     required=True,     type=str,  help="config file")
parser.add_argument("-q", "--queue",          default="1nh",     type=str,  help="lxbatch queue to use")
parser.add_argument("-n", "--nJobs",          default=1,         type=int,  help="number of jobs")
parser.add_argument("-s", "--submit",                                       help="submit jobs", action='store_true')
parser.add_argument("-v", "--verbose",                                      help="increase output verbosity", action='store_true')

args = parser.parse_args()


print 'START'

currDir = os.getcwd()

print

try:
   subprocess.check_output(['mkdir','jobs'])
except subprocess.CalledProcessError as e:
   print e.output
try:
   subprocess.check_output(['mkdir','jobs/'+args.label])
except subprocess.CalledProcessError as e:
   print e.output


##### creates directory and file list for job #######
jobDir = currDir+'/jobs/'+args.label+'/'
os.system('mkdir '+jobDir)
os.chdir(jobDir)

for x in range(1, args.nJobs+1):
   
   ##### creates config file #######
   with open(args.baseFolder+'/'+args.configFile) as fi:
      contents = fi.read()
      replaced_contents = contents.replace('INPUTFOLDER', args.inputFolder)
      replaced_contents = replaced_contents.replace('INPUTFILENAME', args.inputFileName)
      replaced_contents = replaced_contents.replace('INPUTTREENAME', args.inputTreeName)
      replaced_contents = replaced_contents.replace('INPUTTREENAME', args.inputTreeName)
      replaced_contents = replaced_contents.replace('OUTPUTFOLDER', args.outputFolder)
      replaced_contents = replaced_contents.replace('OUTPUTFILENAME', args.outputFileName)
      replaced_contents = replaced_contents.replace('NJOBS', str(args.nJobs))
      replaced_contents = replaced_contents.replace('JOBID', str(x))
      with open(jobDir+"/config_"+str(x)+".cfg", "w") as fo:
         fo.write(replaced_contents)
      
   ##### creates jobs #######
   jobFile = jobDir + '/job_' + args.label + '_' + str(x) + '.sh'
   with open(jobFile, 'w') as fout:
      fout.write("#!/bin/sh\n")
      fout.write("echo\n")
      fout.write("echo 'START---------------'\n")
      fout.write("echo 'current dir: ' ${PWD}\n")
      fout.write("cd /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/\n")
      fout.write("eval `scramv1 runtime -sh`\n")
      fout.write("cd "+str(args.baseFolder)+"\n")
      fout.write("echo 'current dir: ' ${PWD}\n")
      fout.write("source scripts/setup.sh\n")
      fout.write(args.baseFolder+'/'+args.exeName+" "+jobDir+"/config_"+str(x)+".cfg > "+jobDir+"/out_"+str(x)+".log\n")
      fout.write("echo 'STOP---------------'\n")
      fout.write("echo\n")
      fout.write("echo\n")
   os.system("chmod 755 "+jobFile)
   
   ###### sends bjobs ######
   if args.submit:
      command = "bsub -q "+args.queue+" -cwd "+jobDir+" "+jobFile
      print command
      os.system(command)
      #print "job " + jobFile + " submitted"
      command = "sleep 0.5"
      os.system(command)
      
os.chdir("../..")

print
print "your jobs:"
os.system("bjobs")
print
print 'END'
print
