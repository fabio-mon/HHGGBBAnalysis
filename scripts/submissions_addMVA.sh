./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l ttgg \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_bkg_v4/           -i TTGG\*.root   -t TTHGenericTagDumper/trees/ttGG_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_withMVAs_new_v4//   -o TTGG \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l ttgjets \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_bkg_v4/           -i TTGJets\*.root   -t TTHGenericTagDumper/trees/ttGJets_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_withMVAs_new_v4//   -o TTGJets \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l ttjets \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_bkg_v4/           -i TTJets\*.root   -t TTHGenericTagDumper/trees/ttJets_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_withMVAs_new_v4//   -o TTJets \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l gg \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_bkg_v4/           -i DiPhotonJetsBox\*.root   -t TTHGenericTagDumper/trees/diphoton_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_withMVAs_new_v4//   -o DiPhotonJetsBox \
    -n 20 \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l fake-fake \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_bkg_v4/           -i fake-fake.root   -t fake-fake \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_withMVAs_new_v4//   -o fake-fake \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l prompt-fake \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_bkg_v4/           -i prompt-fake.root   -t prompt-fake \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_withMVAs_new_v4//   -o prompt-fake \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l ttH \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_signal_v4/           -i ttHJetToGG_M125_13TeV_\*.root   -t TTHGenericTagDumper/trees/tth_125_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_withMVAs_new_v4//   -o ttHJetToGG_M125_13TeV \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l ggH \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_signal_v4/           -i GluGluHToGG_M125_13TeV_\*.root   -t TTHGenericTagDumper/trees/ggh_125_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_withMVAs_new_v4//   -o GluGluHToGG_M125_13TeV \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l VBF \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_signal_v4/           -i VBFHToGG_M125_13TeV\*.root   -t TTHGenericTagDumper/trees/vbf_125_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_withMVAs_new_v4//   -o VBFHToGG_M125_13TeV \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l VH \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_signal_v4/           -i VHToGG_M125_13TeV_\*.root   -t TTHGenericTagDumper/trees/vh_125_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_withMVAs_new_v4//   -o VHToGG_M125_13TeV \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l bbH \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_signal_v4/           -i bbHToGG_M125_\*.root   -t TTHGenericTagDumper/trees/bbh_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_withMVAs_new_v4//   -o bbHToGG_M125 \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l tHq \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_signal_v4/           -i THQ_M125_HToGG_13TeV-madgraph-pythia8_\*.root   -t TTHGenericTagDumper/trees/thq_125_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_withMVAs_new_v4//   -o THQ_M125_HToGG_13TeV-madgraph-pythia8 \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l tHW \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_signal_v4/           -i THW_M125_HToGG_13TeV-madgraph-pythia8_\*.root   -t TTHGenericTagDumper/trees/thw_125_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_withMVAs_new_v4//   -o THW_M125_HToGG_13TeV-madgraph-pythia8 \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l data \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_data_v4/           -i DoubleEG_Run2016\*.root   -t TTHGenericTagDumper/trees/Data_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_withMVAs_new_v4//   -o DoubleEG_Run2016 \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l CS_oneCategory \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_bkg_v4/           -i controlSample.root   -t ControlSample_oneCategory \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_withMVAs_new_v4//   -o controlSample_oneCategory \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l CS_diLepton \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_bkg_v4/           -i controlSample.root   -t ControlSample_diLepton \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_withMVAs_new_v4//   -o controlSample_diLepton \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l CS_singleLepton \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_bkg_v4/           -i controlSample.root   -t ControlSample_singleLepton \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_8_0_28/src/flashgg/Jobs/ntuples_HggPresel_2jets_withMVAs_new_v4//   -o controlSample_singleLepton \
    --submit


