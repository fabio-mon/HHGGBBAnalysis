#!/usr/bin/env python
import os, re
import commands
import math, time
import sys
import argparse
import subprocess

parser = argparse.ArgumentParser(description='This script prepare workspace for klambda scan')

parser.add_argument("-l",              "--label",    required=True,     type=str,    help="job label")
parser.add_argument("-o",       "--outputFolder",    required=True,     type=str,    help="folder where to store output files")
parser.add_argument("-kl_min",   "--klambda_min",    required=True,     type=float,  help="lower value of klambda to be considered")
parser.add_argument("-kl_max",   "--klambda_max",    required=True,     type=float,  help="upper value of klambda to be considered")
parser.add_argument("-N_kl",       "--N_klambda",    required=True,     type=int,    help="granularity of klambda scan(=nJobs to submit)")
parser.add_argument("-bkg",       "--background",                                    help="reproduce background",      action='store_true')
parser.add_argument("-q",              "--queue",    default="1nd",     type=str,    help="hercules queue to use")
parser.add_argument("-s",             "--submit",                                    help="submit jobs",               action='store_true')
parser.add_argument("-v",            "--verbose",                                    help="increase output verbosity", action='store_true')


args = parser.parse_args()

HHGGBBselector = "/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/bin/HHGGBB_selector.exe"
addMVA = "/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/bin/addMVA.exe"
MakeWorkspace = "/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/bin/MakeWorkspace.exe"
#cfg_HHGGBBselector = "/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/cfg/fmonti_HHggbb_Delphes_klambdatemplate.cfg"
cfg_folder = "/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/cfg/"
cfg_addMVA_HT = "/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/cfg/addMVA_HHggbb_HT_klambdatemplate.cfg"
cfg_addMVA_LT = "/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/cfg/addMVA_HHggbb_LT_klambdatemplate.cfg"
normalization_file = "/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/data/klambda_rew.txt"

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
try:
   subprocess.check_output(['mkdir',args.outputFolder+"/"+args.label+"/"])
except subprocess.CalledProcessError as e:
   print e.output
outdir=args.outputFolder+"/"+args.label+"/"

if args.background:
   print("produce background ntuples")
   bkg_samples=["ggH", "gg", "qqH", "VH", "ttH", "bbH", "ttgg", "ttghad", "ttgsemilepfromt", "ttgsemilepfromtbar", "ttglep", "tt"]
   scenarios=["pessimistic", "intermediate", "optimistic"]
   for sample in bkg_samples:
      jobDir = currDir+"/jobs/"+args.label+"/job_"+sample
      os.system('mkdir '+jobDir)
      os.chdir(jobDir)

      ##### creates config file for HHGGBBselector #######
      with open(cfg_folder+"/template_common.cfg") as fi:
         contents = fi.read()
         replaced_contents = contents.replace('OUTFOLDER',outdir)
      with open(jobDir+"/common.cfg", "w") as fo:
         fo.write(replaced_contents)

      with open(cfg_folder+"fmonti_"+sample+"_Delphes.cfg") as fi:
         contents = fi.read()
         replaced_contents = contents.replace("importCfg /afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/cfg/common.cfg","importCfg "+jobDir+"/common.cfg")
      with open(jobDir+"/HHggbb_config.cfg", "w") as fo:
         fo.write(replaced_contents)

      ##### creates config file for addMVA high mass #######
      with open(cfg_folder+"/template_HT_addMVA_"+sample+".cfg") as fi2:
         contents = fi2.read()
         replaced_contents = contents.replace('INFOLDER', outdir).replace('OUTFOLDER',outdir)
      with open(jobDir+"/addMVA_HT_config.cfg", "w") as fo:
         fo.write(replaced_contents)
  
      ##### creates config file for addMVA low mass #######
      with open(cfg_folder+"/template_LT_addMVA_"+sample+".cfg") as fi3:
         contents = fi3.read()
         replaced_contents = contents.replace('INFOLDER', outdir).replace('OUTFOLDER',outdir)
      with open(jobDir+"/addMVA_LT_config.cfg", "w") as fo:
         fo.write(replaced_contents)
   
      ##### creates jobs #######
      with open("job_"+sample+".sh", 'w') as fout:
         fout.write("#!/bin/sh\n")
         fout.write("echo\n")
         fout.write("echo 'START---------------'\n")
         fout.write("workdir=${PWD}\n")
         fout.write("echo 'current dir: ' $workdir\n")
         fout.write("source /afs/cern.ch/user/f/fmonti/bin/HHGGBBsetup.sh\n")
         #fout.write("cd "+str(jobDir)+"\n")
         fout.write("cd $workdir\n")
         fout.write("echo 'current dir: ' ${PWD}\n")
         fout.write(HHGGBBselector+" "+jobDir+"/HHggbb_config.cfg\n")
         fout.write(addMVA+" "+jobDir+"/addMVA_HT_config.cfg\n")
         fout.write(addMVA+" "+jobDir+"/addMVA_LT_config.cfg\n")
         for scenario in scenarios:
            fout.write(MakeWorkspace+" "+outdir+"/plotTree_"+sample+"_HT_withMVA.root all_highMx "+scenario+" \n")
            fout.write(MakeWorkspace+" "+outdir+"/plotTree_"+sample+"_LT_withMVA.root all_lowMx  "+scenario+" \n")
         fout.write("echo 'STOP---------------'\n")
         fout.write("echo\n")
         fout.write("echo\n")
      os.system("chmod 755 job_"+sample+".sh")
      ###### sends bjobs ######
      if args.submit:
         os.system("bsub -q "+args.queue+" job_"+sample+".sh")
         print "job nr. " + sample + " submitted"


for i_kl in range(0,args.N_klambda):

   #get the corresponding value of klambda
   kl = args.klambda_min + i_kl*(args.klambda_max-args.klambda_min)/args.N_klambda
   print('%.3f' % kl)
   
   ##### creates directory and file list for job #######
   jobDir = currDir+'/jobs/'+args.label+'/job_'+str(i_kl)
   os.system('mkdir '+jobDir)
   os.chdir(jobDir)
   
   ##### copy executable to the jobDir ######
   #os.system('cp '+HHGGBBselector+' '+jobDir+"/executable.exe")
   
   ##### creates config file for HHGGBBselector #######
   BSM_CrossSection = 1.
   weight_sum = 1.
   with open(normalization_file) as f:
      for line in f:
         line_words = line.split('\t')
         if(line_words[0]!="KLAMBDA"):
            if(kl==float(line_words[0])):
               weight_sum=float(line_words[1])
               BSM_CrossSection=float(line_words[2])

   print ("weight sum = "+str(weight_sum)) 
   print ("BSM XS = "+str(BSM_CrossSection)) 
               
   with open(cfg_folder+"/fmonti_HHggbb_Delphes_klambdatemplate.cfg") as fi:
      contents = fi.read()
      replaced_contents = contents.replace('KLAMBDA', str('%.3f' % kl)).replace('OUTFOLDER',outdir).replace("BSMXS",str(BSM_CrossSection)).replace("WEIGHTSUM",str(weight_sum))
   with open(jobDir+"/HHggbb_config.cfg", "w") as fo:
      fo.write(replaced_contents)

   ##### creates config file for addMVA high mass #######
   with open(cfg_addMVA_HT) as fi2:
      contents = fi2.read()
      replaced_contents = contents.replace('KLAMBDA', str('%.3f' % kl)).replace('OUTFOLDER',outdir)
   with open(jobDir+"/addMVA_HT_config.cfg", "w") as fo:
      fo.write(replaced_contents)
  
 ##### creates config file for addMVA low mass #######
   with open(cfg_addMVA_LT) as fi3:
      contents = fi3.read()
      replaced_contents = contents.replace('KLAMBDA', str('%.3f' % kl)).replace('OUTFOLDER',outdir)
   with open(jobDir+"/addMVA_LT_config.cfg", "w") as fo:
      fo.write(replaced_contents)
 
   ##### creates jobs #######
   with open('job_'+str(i_kl)+'.sh', 'w') as fout:
      fout.write("#!/bin/sh\n")
      fout.write("echo\n")
      fout.write("echo 'START---------------'\n")
      fout.write("workdir=${PWD}\n")
      fout.write("echo 'current dir: ' $workdir\n")
      #fout.write("cd /afs/cern.ch/user/f/fmonti//work/myPhaseTwoAnalysis/CMSSW_9_3_2/src\n")
      #fout.write("eval 'scram runtime -sh'\n")
      #fout.write("export LD_LIBRARY_PATH=/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/lib:/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/DynamicTTree/lib/:/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/CfgManager/lib/:$LD_LIBRARY_PATH\n")
      #fout.write("export DYLD_LIBRARY_PATH=/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/lib:/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/DynamicTTree/lib/:/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/CfgManager/lib/:$DYLD_LIBRARY_PATH\n")
      #fout.write("export ROOT_INCLUDE_PATH=/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/interface:/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/DynamicTTree/interface/:/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/CfgManager/interface/:$ROOT_INCLUDE_PATH\n")
      fout.write("source /afs/cern.ch/user/f/fmonti/bin/HHGGBBsetup.sh\n")
      #fout.write("cd "+str(jobDir)+"\n")
      fout.write("cd $workdir\n")
      fout.write("echo 'current dir: ' ${PWD}\n")
      fout.write(HHGGBBselector+" "+jobDir+"/HHggbb_config.cfg\n")
      fout.write(addMVA+" "+jobDir+"/addMVA_HT_config.cfg\n")
      fout.write(addMVA+" "+jobDir+"/addMVA_LT_config.cfg\n")
      fout.write(MakeWorkspace+" "+outdir+"/plotTree_HHggbb_HT_klambda_"+str('%.3f' % kl)+"_withMVA.root all_highMx intermediate "+str(kl)+"\n")
      fout.write(MakeWorkspace+" "+outdir+"/plotTree_HHggbb_LT_klambda_"+str('%.3f' % kl)+"_withMVA.root all_lowMx intermediate "+str(kl)+"\n")
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
if args.background:
   print("to be done after bkg jobs termination:")
   print("cd "+outdir+" && mkdir -p bkg_LT/pessimistic bkg_LT/intermediate bkg_LT/optimistic bkg_HT/pessimistic bkg_HT/intermediate bkg_HT/optimistic")

   print("cd "+outdir+" && ls pessimistic_subset_plotTree_*_LT_withMVA.root | grep -v HHggbb | awk '{print \"mv \"$0\" bkg_LT/pessimistic/\"}' | /bin/sh")
   print("cd "+outdir+" && ls intermediate_subset_plotTree_*_LT_withMVA.root | grep -v HHggbb | awk '{print \"mv \"$0\" bkg_LT/intermediate/\"}' | /bin/sh")
   print("cd "+outdir+" && ls optimistic_subset_plotTree_*_LT_withMVA.root | grep -v HHggbb | awk '{print \"mv \"$0\" bkg_LT/optimistic/\"}' | /bin/sh")
   print("cd "+outdir+" && ls pessimistic_subset_plotTree_*_HT_withMVA.root | grep -v HHggbb | awk '{print \"mv \"$0\" bkg_HT/pessimistic/\"}' | /bin/sh")
   print("cd "+outdir+" && ls intermediate_subset_plotTree_*_HT_withMVA.root | grep -v HHggbb | awk '{print \"mv \"$0\" bkg_HT/intermediate/\"}' | /bin/sh")
   print("cd "+outdir+" && ls optimistic_subset_plotTree_*_HT_withMVA.root | grep -v HHggbb | awk '{print \"mv \"$0\" bkg_HT/optimistic/\"}' | /bin/sh")

   print("cd "+outdir+"/bkg_LT/pessimistic && hadd LT_DoubleEG.root *_gg_* *_tt_* *_ttgg_* *_ttghad_* *_ttglep_* *_ttgsemilepfromt_* *_ttgsemilepfromtbar_*")
   print("cd "+outdir+"/bkg_LT/intermediate && hadd LT_DoubleEG.root *_gg_* *_tt_* *_ttgg_* *_ttghad_* *_ttglep_* *_ttgsemilepfromt_* *_ttgsemilepfromtbar_*")
   print("cd "+outdir+"/bkg_LT/optimistic && hadd LT_DoubleEG.root *_gg_* *_tt_* *_ttgg_* *_ttghad_* *_ttglep_* *_ttgsemilepfromt_* *_ttgsemilepfromtbar_*")
   print("cd "+outdir+"/bkg_HT/pessimistic && hadd LT_DoubleEG.root *_gg_* *_tt_* *_ttgg_* *_ttghad_* *_ttglep_* *_ttgsemilepfromt_* *_ttgsemilepfromtbar_*")
   print("cd "+outdir+"/bkg_HT/intermediate && hadd LT_DoubleEG.root *_gg_* *_tt_* *_ttgg_* *_ttghad_* *_ttglep_* *_ttgsemilepfromt_* *_ttgsemilepfromtbar_*")
   print("cd "+outdir+"/bkg_HT/optimistic && hadd LT_DoubleEG.root *_gg_* *_tt_* *_ttgg_* *_ttghad_* *_ttglep_* *_ttgsemilepfromt_* *_ttgsemilepfromtbar_*")

   print("cd "+outdir+"/bkg_LT/pessimistic && mv *_bbH_* LT_output_bbHToGG_M-125_13TeV_amcatnlo.root && mv *_ttH_* LT_output_ttHToGG_M125_13TeV_powheg_pythia8_v2.root && mv *_VH_* LT_output_VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root && mv *_qqH_* LT_output_VBFHToGG_M-125_13TeV_powheg_pythia8.root && mv *_ggH_* LT_output_GluGluHToGG_M-125_13TeV_powheg_pythia8.root")
   print("cd "+outdir+"/bkg_LT/intermediate && mv *_bbH_* LT_output_bbHToGG_M-125_13TeV_amcatnlo.root && mv *_ttH_* LT_output_ttHToGG_M125_13TeV_powheg_pythia8_v2.root && mv *_VH_* LT_output_VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root && mv *_qqH_* LT_output_VBFHToGG_M-125_13TeV_powheg_pythia8.root && mv *_ggH_* LT_output_GluGluHToGG_M-125_13TeV_powheg_pythia8.root")
   print("cd "+outdir+"/bkg_LT/optimistic && mv *_bbH_* LT_output_bbHToGG_M-125_13TeV_amcatnlo.root && mv *_ttH_* LT_output_ttHToGG_M125_13TeV_powheg_pythia8_v2.root && mv *_VH_* LT_output_VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root && mv *_qqH_* LT_output_VBFHToGG_M-125_13TeV_powheg_pythia8.root && mv *_ggH_* LT_output_GluGluHToGG_M-125_13TeV_powheg_pythia8.root")
   print("cd "+outdir+"/bkg_HT/pessimistic && mv *_bbH_* LT_output_bbHToGG_M-125_13TeV_amcatnlo.root && mv *_ttH_* LT_output_ttHToGG_M125_13TeV_powheg_pythia8_v2.root && mv *_VH_* LT_output_VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root && mv *_qqH_* LT_output_VBFHToGG_M-125_13TeV_powheg_pythia8.root && mv *_ggH_* LT_output_GluGluHToGG_M-125_13TeV_powheg_pythia8.root")
   print("cd "+outdir+"/bkg_HT/intermediate && mv *_bbH_* LT_output_bbHToGG_M-125_13TeV_amcatnlo.root && mv *_ttH_* LT_output_ttHToGG_M125_13TeV_powheg_pythia8_v2.root && mv *_VH_* LT_output_VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root && mv *_qqH_* LT_output_VBFHToGG_M-125_13TeV_powheg_pythia8.root && mv *_ggH_* LT_output_GluGluHToGG_M-125_13TeV_powheg_pythia8.root")
   print("cd "+outdir+"/bkg_HT/optimistic && mv *_bbH_* LT_output_bbHToGG_M-125_13TeV_amcatnlo.root && mv *_ttH_* LT_output_ttHToGG_M125_13TeV_powheg_pythia8_v2.root && mv *_VH_* LT_output_VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root && mv *_qqH_* LT_output_VBFHToGG_M-125_13TeV_powheg_pythia8.root && mv *_ggH_* LT_output_GluGluHToGG_M-125_13TeV_powheg_pythia8.root")

print 'END'
print
