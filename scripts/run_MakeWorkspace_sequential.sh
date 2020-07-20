#!/bin/sh

WORKDIR="/afs/cern.ch/work/f/fmonti/HHGGBBAnalysis" 
SAMPLEDIR="/eos/user/f/fmonti/HHGGBB/data/study_WITHMTD_v7/"

cd $WORKDIR
source scripts/original_setup.sh

labels=("HHggbb" "ggH" "gg" "qqH" "VH" "ttH" "bbH" "ttgg" "ttghad" "ttgsemilepfromt" "ttgsemilepfromtbar" "ttglep" "tt" "gjet" "qcd" "qcd30-40")

echo "HIGH MASS"
for var in "${labels[@]}"
do
  echo "${var}"
  printf "________________________\n"
  printf "${var}\n"
  $WORKDIR/bin/MakeWorkspace.exe $SAMPLEDIR/plotTree_${var}_LT_withMVA.root all_lowMx pessimistic  
  $WORKDIR/bin/MakeWorkspace.exe $SAMPLEDIR/plotTree_${var}_LT_withMVA.root all_lowMx intermediate 
  $WORKDIR/bin/MakeWorkspace.exe $SAMPLEDIR/plotTree_${var}_LT_withMVA.root all_lowMx optimistic   
  $WORKDIR/bin/MakeWorkspace.exe $SAMPLEDIR/plotTree_${var}_HT_withMVA.root all_highMx pessimistic
  $WORKDIR/bin/MakeWorkspace.exe $SAMPLEDIR/plotTree_${var}_HT_withMVA.root all_highMx intermediate
  $WORKDIR/bin/MakeWorkspace.exe $SAMPLEDIR/plotTree_${var}_HT_withMVA.root all_highMx optimistic
  printf "\n\n"
done

scenarios=("pessimistic" "intermediate" "optimistic")

#cd $WORKDIR
echo $PWD
for scenario in "${scenarios[@]}"
do
  mkdir -p $SAMPLEDIR/${scenario}/LT $SAMPLEDIR/${scenario}/HT
  hadd -f $SAMPLEDIR/${scenario}/LT/LT_DoubleEG.root \
       $SAMPLEDIR/${scenario}_subset_plotTree_tt_LT_withMVA.root \
       $SAMPLEDIR/${scenario}_subset_plotTree_ttglep_LT_withMVA.root \
       $SAMPLEDIR/${scenario}_subset_plotTree_ttgsemilepfromtbar_LT_withMVA.root \
       $SAMPLEDIR/${scenario}_subset_plotTree_ttgsemilepfromt_LT_withMVA.root \
       $SAMPLEDIR/${scenario}_subset_plotTree_ttghad_LT_withMVA.root \
       $SAMPLEDIR/${scenario}_subset_plotTree_ttgg_LT_withMVA.root \
       $SAMPLEDIR/${scenario}_subset_plotTree_gg_LT_withMVA.root \
       $SAMPLEDIR/${scenario}_subset_plotTree_ggH_LT_withMVA.root \
       $SAMPLEDIR/${scenario}_subset_plotTree_qqH_LT_withMVA.root \
       $SAMPLEDIR/${scenario}_subset_plotTree_VH_LT_withMVA.root \
       $SAMPLEDIR/${scenario}_subset_plotTree_bbH_LT_withMVA.root \
       $SAMPLEDIR/${scenario}_subset_plotTree_ttH_LT_withMVA.root

  hadd -f $SAMPLEDIR/${scenario}/HT/LT_DoubleEG.root \
       $SAMPLEDIR/${scenario}_subset_plotTree_tt_HT_withMVA.root \
       $SAMPLEDIR/${scenario}_subset_plotTree_ttglep_HT_withMVA.root \
       $SAMPLEDIR/${scenario}_subset_plotTree_ttgsemilepfromtbar_HT_withMVA.root \
       $SAMPLEDIR/${scenario}_subset_plotTree_ttgsemilepfromt_HT_withMVA.root \
       $SAMPLEDIR/${scenario}_subset_plotTree_ttghad_HT_withMVA.root \
       $SAMPLEDIR/${scenario}_subset_plotTree_ttgg_HT_withMVA.root \
       $SAMPLEDIR/${scenario}_subset_plotTree_gg_HT_withMVA.root \
       $SAMPLEDIR/${scenario}_subset_plotTree_ggH_HT_withMVA.root \
       $SAMPLEDIR/${scenario}_subset_plotTree_qqH_HT_withMVA.root \
       $SAMPLEDIR/${scenario}_subset_plotTree_VH_HT_withMVA.root \
       $SAMPLEDIR/${scenario}_subset_plotTree_bbH_HT_withMVA.root \
       $SAMPLEDIR/${scenario}_subset_plotTree_ttH_HT_withMVA.root

  cp $SAMPLEDIR/${scenario}_subset_plotTree_HHggbb_LT_withMVA.root   $SAMPLEDIR/${scenario}/LT/LT_output_GluGluToHHTo2B2G_node_SM_13TeV-madgraph.root
  cp $SAMPLEDIR/${scenario}_subset_plotTree_ggH_LT_withMVA.root      $SAMPLEDIR/${scenario}/LT/LT_output_GluGluHToGG_M-125_13TeV_powheg_pythia8.root
  cp $SAMPLEDIR/${scenario}_subset_plotTree_qqH_LT_withMVA.root      $SAMPLEDIR/${scenario}/LT/LT_output_VBFHToGG_M-125_13TeV_powheg_pythia8.root
  cp $SAMPLEDIR/${scenario}_subset_plotTree_VH_LT_withMVA.root       $SAMPLEDIR/${scenario}/LT/LT_output_VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root 
  cp $SAMPLEDIR/${scenario}_subset_plotTree_bbH_LT_withMVA.root      $SAMPLEDIR/${scenario}/LT/LT_output_bbHToGG_M-125_13TeV_amcatnlo.root
  cp $SAMPLEDIR/${scenario}_subset_plotTree_ttH_LT_withMVA.root      $SAMPLEDIR/${scenario}/LT/LT_output_ttHToGG_M125_13TeV_powheg_pythia8_v2.root

  cp $SAMPLEDIR/${scenario}_subset_plotTree_HHggbb_HT_withMVA.root   $SAMPLEDIR/${scenario}/HT/LT_output_GluGluToHHTo2B2G_node_SM_13TeV-madgraph.root
  cp $SAMPLEDIR/${scenario}_subset_plotTree_ggH_HT_withMVA.root      $SAMPLEDIR/${scenario}/HT/LT_output_GluGluHToGG_M-125_13TeV_powheg_pythia8.root
  cp $SAMPLEDIR/${scenario}_subset_plotTree_qqH_HT_withMVA.root      $SAMPLEDIR/${scenario}/HT/LT_output_VBFHToGG_M-125_13TeV_powheg_pythia8.root
  cp $SAMPLEDIR/${scenario}_subset_plotTree_VH_HT_withMVA.root       $SAMPLEDIR/${scenario}/HT/LT_output_VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root 
  cp $SAMPLEDIR/${scenario}_subset_plotTree_bbH_HT_withMVA.root      $SAMPLEDIR/${scenario}/HT/LT_output_bbHToGG_M-125_13TeV_amcatnlo.root
  cp $SAMPLEDIR/${scenario}_subset_plotTree_ttH_HT_withMVA.root      $SAMPLEDIR/${scenario}/HT/LT_output_ttHToGG_M125_13TeV_powheg_pythia8_v2.root
done
