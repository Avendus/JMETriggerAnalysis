This git repository was created to compare the online and offline of PFHC. It was forked from the run_13_0_X branch in cghuh's repository.

- Installation:
$ git clone https://github.com/Avendus/JMETriggerAnalysis.git -b run3_13_0_X

- EDAnalyzer
PATH: JMETriggerAnalysis/PFHadronCalibration/plugins/PFHadCalibNTuple.cc
Here, variables related to online and offline are received, and the online cut is applied to offline as well. However, the offline uses the offline's PFCandidates variable for the cut to be applied. Both online and offline must pass through the cut for variables to be filled in the ntuple.

- CMSSW Configure
PATH: JMETriggerAnalysis/PFHadronCalibration/test/pfHadCalibNTuple_cfg.py
Added online/offline Sim Particle and online/offline PF Candidates to the collection. (However, the online/offline Sim Particle uses the same input tag, so they are essentially the "same" collection.)

- Analyzer code
PATH: JMETriggerAnalysis/PFHadronCalibration/test/Analyzer_PFHC_Sample.cpp
Intends to compare the energy measured in the calorimeter using Online and Offline PF Candidates on an event-by-event basis.
