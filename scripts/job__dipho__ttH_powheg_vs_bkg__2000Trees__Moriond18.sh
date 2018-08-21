#!/bin/sh

cd /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/MVATraining/
PWD=$(pwd)
echo $PWD

export SCRAM_ARCH=slc6_amd64_gcc530
eval $(scram runtime -sh)

cd -

myCommand="/afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/MVATraining/classify.py -i /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ -T /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/MVATraining/dipho__ttH_powheg_vs_bkg__2000Trees__Moriond18.json"

echo "*** trying to execute "
echo $myCommand
echo "***"
eval "$myCommand"

ls -l
echo "cp tmva*.root /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/MVATraining/"
cp tmva*.root /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/MVATraining/
echo "cp weights/* /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/MVATraining/weights/"
cp weights/* /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/MVATraining/weights/
