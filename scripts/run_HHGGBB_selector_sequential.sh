#!/bin/sh

WORKDIR="/afs/cern.ch/work/a/abenagli/HHGGBB/HHGGBBAnalysis" 
OUTDIR="/afs/cern.ch/user/a/abenagli/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200/"

cd $WORKDIR
#source scripts/setup.sh

labels=("gjet" "gg" "HHggbb" "ggH" "qqH" "VH" "ttH" "bbH" "ttgg" "ttghad" "ttgsemilepfromt" "ttgsemilepfromtbar" "ttglep" "tt" "qcd")

for var in "${labels[@]}"
do
  echo "${var}"
  printf "________________________\n" >> $OUTDIR/HHGGBB_selector_report.txt 2>&1
  printf "${var}\n" >> $OUTDIR/HHGGBB_selector_report.txt 2>&1
# $WORKDIR/bin/HHGGBB_selector.exe $WORKDIR/cfg/fmonti_${var}_Delphes.cfg | grep "acceptance" >> $OUTDIR/HHGGBB_selector_report.txt 2>&1
  $WORKDIR/bin/HHGGBB_selector.exe $WORKDIR/cfg/fmonti_${var}_Delphes.cfg
  printf "\n\n" >> $OUTDIR/HHGGBB_selector_report.txt 2>&1
done
