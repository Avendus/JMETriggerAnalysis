----------
----------

**Note**: for instructions on HLT PF-Hadron calibrations and Jet Energy Scale Corrections,
please ignore this `readme`, and follow the instructions in the dedicated `readme` files:

 * [`readme` for HLT PF-Hadron Calibrations](https://github.com/missirol/JMETriggerAnalysis/blob/run3/PFHadronCalibrations/readme.md)
 * [`readme` for HLT Jet Energy Scale Corrections](https://github.com/missirol/JMETriggerAnalysis/blob/run3/JESCorrections/readme.md)

----------
----------

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

# external data
mkdir -p ${CMSSW_BASE}/src/JMETriggerAnalysis/NTuplizers/data

# PFHC: preliminary HLT-PFHC for Run-3
cp /afs/cern.ch/work/m/missirol/public/run3/PFHC/PFHC_Run3Winter20_HLT_v01.db \
   ${CMSSW_BASE}/src/JMETriggerAnalysis/NTuplizers/data/PFHC_Run3Winter20_HLT_v01.db

# JESC: preliminary HLT-JESCs for Run-3
cp /afs/cern.ch/work/m/missirol/public/run3/JESC/Run3Winter20_V1_MC/Run3Winter20_V1_MC.db \
   ${CMSSW_BASE}/src/JMETriggerAnalysis/NTuplizers/data/JESC_Run3Winter20_V1_MC.db

git clone https://github.com/missirol/JetMETAnalysis.git -o missirol -b run3_jrantuples

scram b -j 12
```

The baseline HLT menu for Run-3 in 11_2_X can be found in
[Common/python/configs/HLT_dev_CMSSW_11_2_0_GRun_V19_configDump.py](https://github.com/missirol/JMETriggerAnalysis/blob/run3_devel_112X/Common/python/configs/HLT_dev_CMSSW_11_2_0_GRun_V19_configDump.py).

It was created with `hltGetConfiguration` via the commands listed in
[`Common/test/dumpHLTMenus_mcRun3.sh`](https://github.com/missirol/JMETriggerAnalysis/blob/6a5807010c9934c25b0cb7bea22d7c90bbb87bc1/Common/test/dumpHLTMenus_mcRun3.sh).

An example of how to enable the different tracking customisations can be found in
[`NTuplizers/test/jmeTriggerNTuple_cfg.py`](https://github.com/missirol/JMETriggerAnalysis/blob/6a5807010c9934c25b0cb7bea22d7c90bbb87bc1/Common/test/dumpHLTMenus_mcRun3.sh)
(see option `reco`).
Test commands:
```
# (a) Run-3 tracking: standard
cmsRun jmeTriggerNTuple_cfg.py maxEvents=1 reco=HLT_Run3TRK output=out_HLT_Run3TRK.root

# (b) Run-3 tracking: all pixel vertices
cmsRun jmeTriggerNTuple_cfg.py maxEvents=1 reco=HLT_Run3TRKWithPU output=out_HLT_Run3TRKWithPU.root
```

----------
