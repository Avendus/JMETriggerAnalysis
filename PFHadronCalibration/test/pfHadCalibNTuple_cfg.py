###
### command-line arguments
###
import FWCore.ParameterSet.VarParsing as vpo
opts = vpo.VarParsing('analysis')

opts.register('skipEvents', 0,
              vpo.VarParsing.multiplicity.singleton,
              vpo.VarParsing.varType.int,
              'number of events to be skipped')

opts.register('numThreads', 1,
              vpo.VarParsing.multiplicity.singleton,
              vpo.VarParsing.varType.int,
              'number of threads')

opts.register('numStreams', 0,
              vpo.VarParsing.multiplicity.singleton,
              vpo.VarParsing.varType.int,
              'number of streams')

opts.register('dumpPython', None,
              vpo.VarParsing.multiplicity.singleton,
              vpo.VarParsing.varType.string,
              'path to python file with content of cms.Process')

opts.register('output', 'pfHCNTuple.root',
              vpo.VarParsing.multiplicity.singleton,
              vpo.VarParsing.varType.string,
              'path to output ROOT file')

opts.register('verbosity', 0,
              vpo.VarParsing.multiplicity.singleton,
              vpo.VarParsing.varType.int,
              'verbosity level')

opts.parseArguments()

###
### HLT configuration
###

from HLT_dev_CMSSW_13_0_0_GRun_configDump import *
#from HLT_dev_CMSSW_13_0_0_GRun_configDump import cms, process

# remove cms.OutputModule objects from HLT config-dump
for _modname in process.outputModules_():
  _mod = getattr(process, _modname)
  if type(_mod) == cms.OutputModule:
    process.__delattr__(_modname)
    if opts.verbosity > 0:
      print('> removed cms.OutputModule:', _modname)

# remove cms.EndPath objects from HLT config-dump
for _modname in process.endpaths_():
    _mod = getattr(process, _modname)
    if type(_mod) == cms.EndPath:
       process.__delattr__(_modname)
       if opts.verbosity > 0:
          print('> removed cms.EndPath:', _modname)

# remove selected cms.Path objects from HLT config-dump
print('-'*108)
print('{:<99} | {:<4} |'.format('cms.Path', 'keep'))
print('-'*108)
for _modname in sorted(process.paths_()):
  _keepPath = _modname.startswith('MC_') and ('Jets' in _modname or 'MET' in _modname)
  _keepPath |= _modname.startswith('MC_ReducedIterativeTracking')
  if _keepPath:
    print('{:<99} | {:<4} |'.format(_modname, '+'))
    continue
  _mod = getattr(process, _modname)
  if type(_mod) == cms.Path:
    process.__delattr__(_modname)
    print('{:<99} | {:<4} |'.format(_modname, ''))
print('-'*108)

# remove FastTimerService
del process.FastTimerService

###
### PFHC-specific paths
###
# For PFSimparticle
process.particleFlowSimParticle = cms.EDProducer('PFSimParticleProducer',
  Fitter = cms.string('KFFittingSmoother'),
  MCTruthMatchingInfo = cms.untracked.bool(False),
  ParticleFilter = cms.PSet(
    EMin = cms.double(0),
    chargedPtMin = cms.double(0),
    etaMax = cms.double(5.3),
    invisibleParticles = cms.vint32(),
    protonEMin = cms.double(5000.0),
    rMax = cms.double(129.0),
    zMax = cms.double(317.0),
  ),
  Propagator = cms.string('PropagatorWithMaterial'),
  RecTracks = cms.InputTag('trackerDrivenElectronSeeds'),
  TTRHBuilder = cms.string('WithTrackAngle'),
  ecalRecHitsEB = cms.InputTag('caloRecHits','EcalRecHitsEB'),
  ecalRecHitsEE = cms.InputTag('caloRecHits','EcalRecHitsEE'),
  fastSimProducer = cms.untracked.InputTag('fastSimProducer','EcalHitsEB'),
  process_Particles = cms.untracked.bool(True),
  process_RecTracks = cms.untracked.bool(False),
  sim = cms.InputTag('g4SimHits'),
  verbose = cms.untracked.bool(False)
)

process.pfHadCalibNTuple = cms.EDAnalyzer('PFHadCalibNTuple',
  TTreeName = cms.string('Candidates'),

  genParticles = cms.InputTag('genParticles'), # Generater level info
  pfSimParticles = cms.InputTag('particleFlowSimParticle'), # Sim Particle (for Sim REAL INFO references)
  PFSimParticles = cms.InputTag('particleFlowSimParticle'), # Sim Particle (for Sim REAL INFO references)
  recoPFCandidates = cms.InputTag('hltParticleFlow'), # Online PF Candidates
  PFCandidates = cms.InputTag("particleFlow"), # Offline PF Candidates

  genParticleStatus = cms.int32(1), # status code of selected GEN particles
  genParticlePdgId = cms.int32(-211), # pdgID of selected GEN particles
  genParticleIsoMinDeltaR = cms.double(1.0), # min deltaR between "isolated" GEN particles and other GEN particles

  minPt = cms.double(1.00), # min pt
  minTrackP = cms.double(1.00), # min track momentum
  minTrackPt = cms.double(1.00), # min track transverse momentum
  minCaloEnergy = cms.double(0.5), # min ecal+hcal energy
  maxECalEnergy = cms.double(1e12), # max ecal energy

  # min nb of pixel and pixel+strip hits (per track-eta range)
  minPixelHits = cms.vuint32(2, 2, 2, 2, 2),
  #minTrackerHits = cms.vuint32(14, 17, 20, 17, 10),
  minTrackerHits = cms.vuint32(2, 2, 2, 2, 2),
  maxEtaForMinTrkHitsCuts = cms.vdouble(1.4, 1.6, 2.0, 2.4, 2.6),

  usePFBlockElements = cms.bool(True),
)

# for offline
process.load("RecoParticleFlow.PFProducer.particleFlowSimParticle_cff")
process.particleFlowSimParticle.ParticleFilter = cms.PSet(    
        # Allow *ALL* protons with energy > protonEMin                                                                               
        protonEMin = cms.double(5000.0),
        # Particles must have abs(eta) < etaMax (if close enough to 0,0,0)                                                           
        etaMax = cms.double(5.3),
        # Charged particles with pT < pTMin (GeV/c) are not simulated                                                                
        chargedPtMin = cms.double(0.0),
        # Particles must have energy greater than EMin [GeV]                                                                         
        EMin = cms.double(0.0),
        # half-length of the ECAL endcap inner surface                                                                                                                                                      
        rMax = cms.double(129.),
        zMax = cms.double(317.),
        # List of invisible particles (abs of pdgid)                                                                                                                                                        
        invisibleParticles = cms.vint32()
)

process.genReReco = cms.Sequence(process.particleFlowSimParticle) #for offline
process.genReRecoPath = cms.Path(process.genReReco) #for offline
process.schedule.append(process.genReRecoPath) #for offline

process.pfSimParticleSeq = cms.Sequence(process.particleFlowSimParticle)
process.pfSimParticlePath = cms.Path(process.pfSimParticleSeq)
process.schedule.append(process.pfSimParticlePath)

process.pfHadCalibNTupleSeq = cms.Sequence(process.pfHadCalibNTuple)
process.pfHadCalibNTupleEndPath = cms.EndPath(process.pfHadCalibNTupleSeq)
process.schedule.append(process.pfHadCalibNTupleEndPath)

###
### options
###

# number of events
process.maxEvents.input = opts.maxEvents

# number of events to be skipped
process.source.skipEvents = cms.untracked.uint32(opts.skipEvents)

# number of threads/streams
process.options.numberOfThreads = opts.numThreads
process.options.numberOfStreams = opts.numStreams

process.options.wantSummary = cms.untracked.bool(False)

# input files
if opts.inputFiles:
  process.source.fileNames = opts.inputFiles
else:
  #process.source.fileNames = [
  #              '/store/mc/Run3Winter23Reco/SinglePionGun_E0p2to200/GEN-SIM-RECO/EpsilonPU_126X_mcRun3_2023_forPU65_v1-v2/2550000/007001ee-f0d5-4161-99ea-f286ad7136e6.root', # offline sample
  #		'/store/mc/Run3Winter23Digi/SinglePionGun_E0p2to200/GEN-SIM-RAW/NoPUGTv3_126X_mcRun3_2023_forPU65_v3-v2/2560000/01723c2c-b1fa-456f-9102-c8781ee0e36d.root',
  #]
  process.source.fileNames = [
###                '/store/mc/Run3Winter23Reco/SinglePionGun_E0p2to200/GEN-SIM-RECO/NoPUGTv4_GTv4_126X_mcRun3_2023_forPU65_v4-v2/30000/470b7207-b55e-4e1f-93ee-2cb623899a0d.root', # offline sample
                '/store/mc/Run3Winter23Reco/SinglePionGun_E0p2to200/GEN-SIM-RECO/NoPUGTv4_GTv4_126X_mcRun3_2023_forPU65_v4-v2/2810000/40c98c3b-5f1c-46e4-ae43-1bbe231ad485.root',
  ]
# Real path of online sample current crab version
#/eos/cms/store/mc/Run3Winter23Digi/SinglePionGun_E0p2to200/GEN-SIM-RAW/NoPUGTv4_126X_mcRun3_2023_forPU65_v4-v2/30001

# input EDM files [secondary]
if not hasattr(process.source, 'secondaryFileNames'):
  process.source.secondaryFileNames = cms.untracked.vstring()

if opts.secondaryInputFiles:
  process.source.secondaryFileNames = opts.secondaryInputFiles
else:
  process.source.secondaryFileNames = [
###		'/store/mc/Run3Winter23Digi/SinglePionGun_E0p2to200/GEN-SIM-RAW/NoPUGTv4_126X_mcRun3_2023_forPU65_v4-v2/30000/012505de-55c8-489e-8fed-2d76f4617492.root',
                '/store/mc/Run3Winter23Digi/SinglePionGun_E0p2to200/GEN-SIM-RAW/NoPUGTv4_126X_mcRun3_2023_forPU65_v4-v2/2810000/2b9c157c-12e8-4344-b604-999366a9d91f.root',
  ]
# You have to check the secondary File is dauther of RECO Sample file at DAS


# output file
process.TFileService = cms.Service('TFileService', fileName = cms.string(opts.output))

# dump content of cms.Process to python file
if opts.dumpPython is not None:
  process.prune()
  open(opts.dumpPython, 'w').write(process.dumpPython())
