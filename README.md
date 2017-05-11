# phase2L1TPatterns
Repository for Phase 2 L1 Trigger Pattern Creation


Checkout instructions:
```
scramv1 project CMSSW CMSSW_9_1_0_pre2
cd CMSSW_9_1_0_pre2/src
eval `scramv1 runtime -sh`
cmsenv
git cms-merge-topic isobelojalvo:CMSSW_9_0_X_L1PF_devel
git cms-merge-topic skinnari:Tracklet_91X
cd L1Trigger
git clone git@github.com:isobelojalvo/phase2L1TPatterns.git
cd ../
nohup scramv1 b -j 8
```