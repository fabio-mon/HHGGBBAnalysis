#!/usr/bin/env python
import os, re
import commands
import math, time
import sys
import argparse
import subprocess

parser = argparse.ArgumentParser(description='This script splits Geant4 tasks in multiple jobs and sends them on LXBATCH')

parser.add_argument("-l", "--label",                 required=True,     type=str,    help="job label")
#parser.add_argument("-e", "--exe",                   required=True,     type=str,    help="executable")
parser.add_argument("-o", "--outputFolder",          required=True,     type=str,    help="folder where to store output files")
#parser.add_argument("-c", "--configFile",            required=True,     type=str,    help="config file to be run")
parser.add_argument("-kl_min",   "--klambda_min",    required=True,     type=float,  help="lower value of klambda to be considered")
parser.add_argument("-kl_max",   "--klambda_max",    required=True,     type=float,  help="upper value of klambda to be considered")
parser.add_argument("-N_kl",     "--N_klambda",      required=True,     type=int,    help="granularity of klambda scan(=nJobs to submit)")
parser.add_argument("-q", "--queue",                 default="1nd",     type=str,    help="hercules queue to use")
parser.add_argument("-s", "--submit",                                                help="submit jobs", action='store_true')
parser.add_argument("-v", "--verbose",                                               help="increase output verbosity", action='store_true')


args = parser.parse_args()

HHGGBBselector = "/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/bin/HHGGBB_selector.exe"
cfg_HHGGBBselector = "/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/cfg/fmonti_HHggbb_Delphes_klambdatemplate.cfg"
addMVA = "/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/bin/addMVA.exe"
MakeWorkspace = "/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/bin/MakeWorkspace.exe"
cfg_addMVA_HT = "/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/cfg/addMVA_HHggbb_HT_klambdatemplate.cfg"
cfg_addMVA_LT = "/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/cfg/addMVA_HHggbb_LT_klambdatemplate.cfg"


print 
print 'START'
print 

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
#try:
#   subprocess.check_output(['mkdir',args.outputFolder+"/"+args.label+"/"])
#except subprocess.CalledProcessError as e:
#   print e.output

for i_kl in range(0,args.N_klambda):

   #get the corresponding value of klambda
   kl = args.klambda_min + (i_kl+0.5)*(args.klambda_max-args.klambda_min)/args.N_klambda
   print(kl)
   
   ##### creates directory and file list for job #######
   jobDir = currDir+'/jobs/'+args.label+'/job_'+str(i_kl)
   os.system('mkdir '+jobDir)
   os.chdir(jobDir)
   
   ##### copy executable to the jobDir ######
   #os.system('cp '+HHGGBBselector+' '+jobDir+"/executable.exe")
   
   ##### creates config file for HHGGBBselector #######
   with open(cfg_HHGGBBselector) as fi:
      contents = fi.read()
      replaced_contents = contents.replace('KLAMBDA', str(kl)).replace('OUTFOLDER',args.outputFolder)
   with open(jobDir+"/HHggbb_config.cfg", "w") as fo:
      fo.write(replaced_contents)

   ##### creates config file for addMVA high mass #######
   with open(cfg_addMVA_HT) as fi2:
      contents = fi2.read()
      replaced_contents = contents.replace('KLAMBDA', str(kl)).replace('OUTFOLDER',args.outputFolder)
   with open(jobDir+"/addMVA_HT_config.cfg", "w") as fo:
      fo.write(replaced_contents)
  
 ##### creates config file for addMVA low mass #######
   with open(cfg_addMVA_LT) as fi3:
      contents = fi3.read()
      replaced_contents = contents.replace('KLAMBDA', str(kl)).replace('OUTFOLDER',args.outputFolder)
   with open(jobDir+"/addMVA_LT_config.cfg", "w") as fo:
      fo.write(replaced_contents)
 
   ##### creates jobs #######
   with open('job_'+str(i_kl)+'.sh', 'w') as fout:
      fout.write("#!/bin/sh\n")
      fout.write("echo\n")
      fout.write("echo 'START---------------'\n")
      fout.write("echo 'current dir: ' ${PWD}\n")
      fout.write("cd /afs/cern.ch/user/f/fmonti//work/myPhaseTwoAnalysis/CMSSW_9_3_2/src\n")
      fout.write("eval 'scram runtime -sh'\n")
      fout.write("export LD_LIBRARY_PATH=/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/lib:/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/DynamicTTree/lib/:/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/CfgManager/lib/:$LD_LIBRARY_PATH\n")
      fout.write("export DYLD_LIBRARY_PATH=/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/lib:/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/DynamicTTree/lib/:/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/CfgManager/lib/:$DYLD_LIBRARY_PATH\n")
      fout.write("export ROOT_INCLUDE_PATH=/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/interface:/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/DynamicTTree/interface/:/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/CfgManager/interface/:$ROOT_INCLUDE_PATH\n")
      #fout.write("source /afs/cern.ch/user/f/fmonti/bin/HHGGBBsetup.sh\n")
      fout.write("cd "+str(jobDir)+"\n")
      fout.write("echo 'current dir: ' ${PWD}\n")
      fout.write(HHGGBBselector+" HHggbb_config.cfg\n")
      fout.write(addMVA+" addMVA_HT_config.cfg\n")
      fout.write(addMVA+" addMVA_LT_config.cfg\n")
      fout.write(MakeWorkspace+" "+args.outputFolder+"/plotTree_HHggbb_HT_klambda_"+str(kl)+"_withMVA.root all_highMx intermediate\n")
      fout.write(MakeWorkspace+" "+args.outputFolder+"/plotTree_HHggbb_LT_klambda_"+str(kl)+"_withMVA.root all_lowMx intermediate\n")

      fout.write("echo 'STOP---------------'\n")
      fout.write("echo\n")
      fout.write("echo\n")
   os.system("chmod 755 job_"+str(i_kl)+".sh")
   
   ###### sends bjobs ######
   if args.submit:
      os.system("bsub -q "+args.queue+" job_"+str(i_kl)+".sh")
      print "job nr. " + str(i_kl) + " submitted"
   
   os.chdir("../..")
   
print
print "your jobs:"
os.system("bjobs")
print
print 'END'
print
