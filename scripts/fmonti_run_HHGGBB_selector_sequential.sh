#!/bin/sh

WORKDIR="/afs/cern.ch/work/f/fmonti/HHGGBBAnalysis/" 
OUTDIR="/eos/user/f/fmonti/HHGGBB/data/study_WITHMTD_v7/"

cd $WORKDIR
source scripts/original_setup.sh

#labels=("gg" "HHggbb" "ggH" "qqH" "VH" "ttH" "bbH" "ttgg" "ttghad" "ttgsemilepfromt" "ttgsemilepfromtbar" "ttglep" "tt" "qcd" "qcd30-40")
#labels=("gjet" "gg" "HHggbb" "ggH" "qqH" "VH" "ttH" "bbH" "ttgg" "ttghad" "ttgsemilepfromt" "ttgsemilepfromtbar" "ttglep" "tt" "qcd" "qcd30-40")
labels=("ttglep" "tt" "qcd" "qcd30-40")

for var in "${labels[@]}"
do
  echo "-----------------------------------------------------------------------------------"
  echo "-----------------------------------------------------------------------------------"
  echo "${var}"
  #printf "________________________\n" >> $OUTDIR/HHGGBB_selector_report.txt 2>&1
  #printf "${var}\n" >> $OUTDIR/HHGGBB_selector_report.txt 2>&1
#  $WORKDIR/bin/HHGGBB_selector.exe $WORKDIR/cfg/fmonti_${var}_Delphes.cfg | grep "acceptance" >> $OUTDIR/HHGGBB_selector_report.txt 2>&1
  $WORKDIR/bin/HHGGBB_selector.exe $WORKDIR/cfg/fmonti_${var}_Delphes.cfg
  #printf "\n\n" >> $OUTDIR/HHGGBB_selector_report.txt 2>&1
done
