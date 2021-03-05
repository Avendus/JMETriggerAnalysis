### Tools for JME studies on the Run-3 HLT reconstruction

* [Tests on HLT Tracking for Run-3](#tests-on-hlt-tracking-for-run-3)

----------

### Tests on HLT Tracking for Run-3

```
cmsrel CMSSW_11_2_0_Patatrack
cd CMSSW_11_2_0_Patatrack/src
cmsenv
git cms-merge-topic missirol:devel_1120pa_kineParticleFilter -u
git cms-merge-topic missirol:devel_puppiPUProxy_1120patatrack -u
git cms-merge-topic mmasciov:tracking-allPVs -u
git clone https://github.com/missirol/JMETriggerAnalysis.git -o missirol -b run3

# Run-3 PFHC: copy preliminary HLT-PFHC for Run-3
mkdir -p ${CMSSW_BASE}/src/JMETriggerAnalysis/NTuplizers/data
cp /afs/cern.ch/user/c/chuh/public/PFCalibration/run3/PFCalibration.db ${CMSSW_BASE}/src/JMETriggerAnalysis/NTuplizers/data/PFHC_Run3Winter20_HLT_v01.db

git clone https://github.com/missirol/JetMETAnalysis.git -o missirol -b run3_jrantuples

scram b -j 12
```

The baseline HLT menu for Run-3 in 11_2_X can be found in
[Common/python/configs/HLT_dev_CMSSW_11_2_0_GRun_V19_configDump.py](https://github.com/missirol/JMETriggerAnalysis/blob/run3_devel_112X/Common/python/configs/HLT_dev_CMSSW_11_2_0_GRun_V19_configDump.py).

It was created with `hltGetConfiguration` via the commands listed in
[`Common/test/dumpHLTMenus_mcRun3.sh`](https://github.com/missirol/JMETriggerAnalysis/blob/run3_devel_112X/Common/test/dumpHLTMenus_mcRun3.sh).

An example of how to enable the different tracking customisations can be found in
[`NTuplizers/test/jmeTriggerNTuple_cfg.py`](https://github.com/missirol/JMETriggerAnalysis/blob/run3_devel_112X/NTuplizers/test/jmeTriggerNTuple_cfg.py)
(see option `reco`).
Test commands:
```
# (a) Run-3 tracking: standard
cmsRun jmeTriggerNTuple_cfg.py maxEvents=1 reco=HLT_Run3TRK output=out_HLT_Run3TRK.root

# (b) Run-3 tracking: all pixel vertices
cmsRun jmeTriggerNTuple_cfg.py maxEvents=1 reco=HLT_Run3TRKWithPU output=out_HLT_Run3TRKWithPU.root
```

----------
