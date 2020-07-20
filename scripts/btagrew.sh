#!/bin/sh

WORKDIR="/afs/cern.ch/work/f/fmonti/HHGGBBAnalysis" 
SAMPLEDIR="/eos/user/f/fmonti/HHGGBB/data/study_WITHMTD_v7/"
REWEIGHTSETUP="DelphesMTD_PabloMTD"
#REWEIGHTSETUP="DelphesMTD_DelphesNoMTD"

cd $WORKDIR
source scripts/original_setup.sh

#scenarios=("pessimistic" "intermediate" "optimistic")
scenarios=("intermediate")

for scenario in "${scenarios[@]}"
do
  echo $WORKDIR/bin/ReweightForBTagMTD_SECCA.exe $SAMPLEDIR/${scenario}//LT_DoubleEG.root                                              TCVARS $REWEIGHTSETUP
  $WORKDIR/bin/ReweightForBTagMTD_SECCA.exe $SAMPLEDIR/${scenario}//LT_DoubleEG.root                                              TCVARS $REWEIGHTSETUP
  echo $WORKDIR/bin/ReweightForBTagMTD_SECCA.exe $SAMPLEDIR/${scenario}//LT_output_GluGluHToGG_M-125_13TeV_powheg_pythia8.root         TCVARS $REWEIGHTSETUP
  $WORKDIR/bin/ReweightForBTagMTD_SECCA.exe $SAMPLEDIR/${scenario}//LT_output_GluGluHToGG_M-125_13TeV_powheg_pythia8.root         TCVARS $REWEIGHTSETUP
  echo $WORKDIR/bin/ReweightForBTagMTD_SECCA.exe $SAMPLEDIR/${scenario}//LT_output_GluGluToHHTo2B2G_node_SM_13TeV-madgraph.root        TCVARS $REWEIGHTSETUP
  $WORKDIR/bin/ReweightForBTagMTD_SECCA.exe $SAMPLEDIR/${scenario}//LT_output_GluGluToHHTo2B2G_node_SM_13TeV-madgraph.root        TCVARS $REWEIGHTSETUP
  echo $WORKDIR/bin/ReweightForBTagMTD_SECCA.exe $SAMPLEDIR/${scenario}//LT_output_VBFHToGG_M-125_13TeV_powheg_pythia8.root            TCVARS $REWEIGHTSETUP
  $WORKDIR/bin/ReweightForBTagMTD_SECCA.exe $SAMPLEDIR/${scenario}//LT_output_VBFHToGG_M-125_13TeV_powheg_pythia8.root            TCVARS $REWEIGHTSETUP
  echo $WORKDIR/bin/ReweightForBTagMTD_SECCA.exe $SAMPLEDIR/${scenario}//LT_output_VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root TCVARS $REWEIGHTSETUP
  $WORKDIR/bin/ReweightForBTagMTD_SECCA.exe $SAMPLEDIR/${scenario}//LT_output_VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root TCVARS $REWEIGHTSETUP
  echo $WORKDIR/bin/ReweightForBTagMTD_SECCA.exe $SAMPLEDIR/${scenario}//LT_output_bbHToGG_M-125_13TeV_amcatnlo.root                   TCVARS $REWEIGHTSETUP
  $WORKDIR/bin/ReweightForBTagMTD_SECCA.exe $SAMPLEDIR/${scenario}//LT_output_bbHToGG_M-125_13TeV_amcatnlo.root                   TCVARS $REWEIGHTSETUP
  echo $WORKDIR/bin/ReweightForBTagMTD_SECCA.exe $SAMPLEDIR/${scenario}//LT_output_ttHToGG_M125_13TeV_powheg_pythia8_v2.root           TCVARS $REWEIGHTSETUP
  $WORKDIR/bin/ReweightForBTagMTD_SECCA.exe $SAMPLEDIR/${scenario}//LT_output_ttHToGG_M125_13TeV_powheg_pythia8_v2.root           TCVARS $REWEIGHTSETUP

  echo mkdir $SAMPLEDIR/rewbtag_${scenario}
  mkdir $SAMPLEDIR/rewbtag_${scenario}
  echo mv $SAMPLEDIR/${scenario}/${REWEIGHTSETUP}_reweight_*.root $SAMPLEDIR/rewbtag_${scenario}/
  mv $SAMPLEDIR/${scenario}/${REWEIGHTSETUP}_reweight_*.root $SAMPLEDIR/rewbtag_${scenario}/
done


