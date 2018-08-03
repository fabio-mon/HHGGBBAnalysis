#!/bin/sh

printf "________________________\n" >> HHGGBB_selector_report.txt
printf "HHGGBB\n" >> HHGGBB_selector_report.txt
/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/bin/HHGGBB_selector.exe /afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/cfg/HHGGBB_Delphes.cfg | grep "acceptance" >> HHGGBB_selector_report.txt
printf "\n\n" >> HHGGBB_selector_report.txt

printf "HGG\n" >> HHGGBB_selector_report.txt
/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/bin/HHGGBB_selector.exe /afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/cfg/HGG_Delphes.cfg | grep "acceptance" >> HHGGBB_selector_report.txt 
printf "\n\n" >> HHGGBB_selector_report.txt

printf "VBFHGG" >> HHGGBB_selector_report.txt
/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/bin/HHGGBB_selector.exe /afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/cfg/VBFHGG_Delphes.cfg | grep "acceptance" >> HHGGBB_selector_report.txt
printf "\n\n" >> HHGGBB_selector_report.txt

printf "ttH" >> HHGGBB_selector_report.txt
/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/bin/HHGGBB_selector.exe /afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/cfg/ttH_Delphes.cfg | grep "acceptance" >> HHGGBB_selector_report.txt
printf "\n\n" >> HHGGBB_selector_report.txt

printf "ZHqqgg" >> HHGGBB_selector_report.txt
/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/bin/HHGGBB_selector.exe /afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/cfg/ZHqqgg_Delphes.cfg | grep "acceptance" >> HHGGBB_selector_report.txt
printf "\n\n" >> HHGGBB_selector_report.txt

printf "ttgamma_hadr" >> HHGGBB_selector_report.txt
/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/bin/HHGGBB_selector.exe /afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/cfg/ttGammaHadr_Delphes.cfg | grep "acceptance" >> HHGGBB_selector_report.txt
printf "\n\n" >> HHGGBB_selector_report.txt

printf "GJet" >> HHGGBB_selector_report.txt
/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/bin/HHGGBB_selector.exe /afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/cfg/GJet_Delphes.cfg | grep "acceptance" >> HHGGBB_selector_report.txt
printf "\n\n" >> HHGGBB_selector_report.txt

printf "DiPhotonJetsBox" >> HHGGBB_selector_report.txt
/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/bin/HHGGBB_selector.exe /afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/cfg/DiPhotonJetsBox_Delphes.cfg | grep "acceptance" >> HHGGBB_selector_report.txt
printf "\n\n" >> HHGGBB_selector_report.txt

printf "QCD doubleEMenriched" >> HHGGBB_selector_report.txt
/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/bin/HHGGBB_selector.exe /afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/cfg/QCDdoubleEMenriched_Delphes.cfg | grep "acceptance" >> HHGGBB_selector_report.txt
printf "\n\n" >> HHGGBB_selector_report.txt



