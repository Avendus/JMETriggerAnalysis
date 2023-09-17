#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFSimParticle.h"
#include "DataFormats/ParticleFlowReco/interface/PFSimParticleFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include <cassert>

#include <TFile.h>
#include <TTree.h>

class PFHadCalibNTuple : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit PFHadCalibNTuple(edm::ParameterSet const&);
  ~PFHadCalibNTuple() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  // name of module
  std::string const moduleLabel_;

  // name of TTree
  std::string const ttreeName_;

  // input EDM collections
  edm::EDGetTokenT<reco::GenParticleCollection> const genPartsToken_;// Gen Sim Particles
  edm::EDGetTokenT<reco::PFSimParticleCollection> const pfSimPartsToken_; // Online Sim Particles [ online & offline's Sim Particles using same tag, So Idealy they are same ]
  edm::EDGetTokenT<reco::PFSimParticleCollection> const offline_pfSimPartsToken_; // Offline Sim Particles [ online & offline's Sim Particles using same tag, So Idealy they are same ]
  edm::EDGetTokenT<reco::PFCandidateCollection> const recoPFCandsToken_; // Online PF Candidates
  edm::EDGetTokenT<reco::PFCandidateCollection> const offline_recoPFCandsToken_; // Offline PF Candidates

  // status code of selected GEN particles
  int const genParticleStatus_;

  // PDG-ID of selected GEN particles
  int const genParticlePdgId_;

  // min-DeltaR of isolated GEN particles from other GEN particles
  double const genParticleIsoMinDeltaR_;

  // min pt for charged hadrons
  double const minPt_;

  // min track-p for charged hadrons
  double const minTrackP_;

  // min track-pT for charged hadrons
  double const minTrackPt_;

  // min ecal+hcal raw energy for charged hadrons
  double const minCaloEnergy_;

  // max ecal raw energy to define a MIP
  double const maxECalEnergy_;

  // min number of pixel and pixel+strip hits for charged hadrons
  std::vector<uint> const minPixelHits_;
  std::vector<uint> const minTrackerHits_;
  std::vector<double> const maxEtaForMinTrkHitsCuts_;

  // use PFBlockElements to count tracks associated to reco::PFCandidate
  bool const usePFBlockElements_;

  // Number of tracks after cuts
  std::vector<uint> globalCounter_;

  TTree* ttree_ = nullptr;

  // Gen Sim Particles
  std::vector<float> gen_energy_; 

  // Online Sim Particles
  std::vector<float> true_energy_;
  std::vector<float> true_eta_;
  std::vector<float> true_phi_;
  std::vector<float> true_dr_;
  std::vector<int> true_charge_;

  // Online PF Candidates
  std::vector<float> pfc_ecal_;
  std::vector<float> pfc_hcal_;
  std::vector<float> pfc_eta_;
  std::vector<float> pfc_phi_;
  std::vector<int> pfc_charge_;
  std::vector<float> pfc_id_;
  std::vector<float> pfc_trackRef_p_;
  std::vector<float> pfc_trackRef_pt_;
  std::vector<float> pfc_trackRef_eta_;
  std::vector<float> pfc_trackRef_phi_;
  std::vector<uint> pfc_trackRef_nValidPixelHits_;
  std::vector<uint> pfc_trackRef_nValidTrackerHits_;
  
  std::vector<std::vector<int>> Is_Online_pfc_pass_cuts_;

  // Offline Sim Particles
  std::vector<float> true_energy_offline_;
  std::vector<float> true_eta_offline_;
  std::vector<float> true_phi_offline_;
  std::vector<float> true_dr_offline_;
  std::vector<int> true_charge_offline_;

  // Offline PF Candidates
  std::vector<float> pfc_ecal_offline_;
  std::vector<float> pfc_hcal_offline_;
  std::vector<float> pfc_eta_offline_;
  std::vector<float> pfc_phi_offline_;
  std::vector<float> pfc_charge_offline_;
  std::vector<float> pfc_id_offline_;
  std::vector<float> pfc_trackRef_p_offline_;
  std::vector<float> pfc_trackRef_pt_offline_;
  std::vector<float> pfc_trackRef_eta_offline_;
  std::vector<float> pfc_trackRef_phi_offline_;
  std::vector<float> pfc_trackRef_nValidPixelHits_offline_;
  std::vector<float> pfc_trackRef_nValidTrackerHits_offline_;

  ////std::vector<std::vector<int>> Is_Offline_pfc_pass_cuts_;

  int Online_pfc_size_;
  int Offline_pfc_size_;
  int diff_pfc_size_;

  std::vector<float> pfc_ecal_diff_;
  std::vector<float> pfc_hcal_diff_;
  std::vector<float> pfc_eta_diff_;
  std::vector<float> pfc_phi_diff_;
  std::vector<float> pfc_charge_diff_;
  std::vector<float> pfc_id_diff_;
  std::vector<float> pfc_trackRef_p_diff_;
  std::vector<float> pfc_trackRef_pt_diff_;
  std::vector<float> pfc_trackRef_eta_diff_;
  std::vector<float> pfc_trackRef_phi_diff_;
  std::vector<float> pfc_trackRef_nValidPixelHits_diff_;
  std::vector<float> pfc_trackRef_nValidTrackerHits_diff_;

  void reset_variables();

};

PFHadCalibNTuple::PFHadCalibNTuple(const edm::ParameterSet& iConfig)
    : moduleLabel_(iConfig.getParameter<std::string>("@module_label")),
      ttreeName_(iConfig.getParameter<std::string>("TTreeName")),
      genPartsToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))), // Token for Gen Particles collection
      pfSimPartsToken_(consumes<reco::PFSimParticleCollection>(iConfig.getParameter<edm::InputTag>("pfSimParticles"))), // Token for Online Sim Particles collection
      offline_pfSimPartsToken_(consumes<reco::PFSimParticleCollection>(iConfig.getParameter<edm::InputTag>("PFSimParticles"))), // Token for Offline Sim Particles collection
      recoPFCandsToken_(consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("recoPFCandidates"))), // Token for Online PF Candidates
      offline_recoPFCandsToken_(consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("PFCandidates"))), // Token for Offline PF Candidates
      genParticleStatus_(iConfig.getParameter<int>("genParticleStatus")),
      genParticlePdgId_(iConfig.getParameter<int>("genParticlePdgId")),
      genParticleIsoMinDeltaR_(iConfig.getParameter<double>("genParticleIsoMinDeltaR")),
      minPt_(iConfig.getParameter<double>("minPt")),
      minTrackP_(iConfig.getParameter<double>("minTrackP")),
      minTrackPt_(iConfig.getParameter<double>("minTrackPt")),
      minCaloEnergy_(iConfig.getParameter<double>("minCaloEnergy")),
      maxECalEnergy_(iConfig.getParameter<double>("maxECalEnergy")),
      minPixelHits_(iConfig.getParameter<std::vector<uint>>("minPixelHits")),
      minTrackerHits_(iConfig.getParameter<std::vector<uint>>("minTrackerHits")),
      maxEtaForMinTrkHitsCuts_(iConfig.getParameter<std::vector<double>>("maxEtaForMinTrkHitsCuts")),
      usePFBlockElements_(iConfig.getParameter<bool>("usePFBlockElements")) {
  assert(minPixelHits_.size() == minTrackerHits_.size());
  assert(minPixelHits_.size() == maxEtaForMinTrkHitsCuts_.size());
  for (uint idx = 0; 1 + idx < maxEtaForMinTrkHitsCuts_.size(); ++idx) {
    assert(maxEtaForMinTrkHitsCuts_[idx] < maxEtaForMinTrkHitsCuts_[idx + 1]);
  }

  // Variable for counting event number which is passed each cut
  globalCounter_ = std::vector<uint>(14, 0);

  usesResource(TFileService::kSharedResource);

  edm::Service<TFileService> fs;

  if (not fs) {
    throw edm::Exception(edm::errors::Configuration, "TFileService is not registered in cfg file");
  }

  ttree_ = fs->make<TTree>(ttreeName_.c_str(), ttreeName_.c_str());

  if (not ttree_) {
    throw edm::Exception(edm::errors::Configuration, "failed to create TTree via TFileService::make<TTree>");
  }

  // Branch for Gen Particles
  ttree_->Branch("gen_energy", &gen_energy_);

  // Branch for Online Sim Particles
  ttree_->Branch("true_energy", &true_energy_);
  ttree_->Branch("true_eta", &true_eta_);
  ttree_->Branch("true_phi", &true_phi_);
  ttree_->Branch("true_dr", &true_dr_);
  ttree_->Branch("true_charge", &true_charge_);

  // Branch for Online PF Candidates
  ttree_->Branch("pfc_ecal", &pfc_ecal_);
  ttree_->Branch("pfc_hcal", &pfc_hcal_);
  ttree_->Branch("pfc_eta", &pfc_eta_);
  ttree_->Branch("pfc_phi", &pfc_phi_);
  ttree_->Branch("pfc_charge", &pfc_charge_);
  ttree_->Branch("pfc_id", &pfc_id_);
  ttree_->Branch("pfc_trackRef_p", &pfc_trackRef_p_);
  ttree_->Branch("pfc_trackRef_pt", &pfc_trackRef_pt_);
  ttree_->Branch("pfc_trackRef_eta", &pfc_trackRef_eta_);
  ttree_->Branch("pfc_trackRef_phi", &pfc_trackRef_phi_);
  ttree_->Branch("pfc_trackRef_nValidPixelHits", &pfc_trackRef_nValidPixelHits_);
  ttree_->Branch("pfc_trackRef_nValidTrackerHits", &pfc_trackRef_nValidTrackerHits_);

  ttree_->Branch("Is_Online_pfc_pass_cuts", &Is_Online_pfc_pass_cuts_);

  // Branch for Offline Sim Particles
  ttree_->Branch("true_energy_offline", &true_energy_offline_);
  ttree_->Branch("true_dr_offline", &true_dr_offline_);
  ttree_->Branch("pfc_ecal_offline", &pfc_ecal_offline_);
  ttree_->Branch("pfc_hcal_offline", &pfc_hcal_offline_);
  ttree_->Branch("true_energy_offline", &true_energy_offline_);
  ttree_->Branch("true_eta_offline", &true_eta_offline_);
  ttree_->Branch("true_phi_offline", &true_phi_offline_);
  ttree_->Branch("true_dr_offline", &true_dr_offline_);
  ttree_->Branch("true_charge_offline", &true_charge_offline_);

  // Branch for Offline PF Candidates
  ttree_->Branch("pfc_ecal_offline", &pfc_ecal_offline_);
  ttree_->Branch("pfc_hcal_offline", &pfc_hcal_offline_);
  ttree_->Branch("pfc_eta_offline", &pfc_eta_offline_);
  ttree_->Branch("pfc_phi_offline", &pfc_phi_offline_);
  ttree_->Branch("pfc_charge_offline", &pfc_charge_offline_);
  ttree_->Branch("pfc_id_offline", &pfc_id_offline_);
  ttree_->Branch("pfc_trackRef_p_offline", &pfc_trackRef_p_offline_);
  ttree_->Branch("pfc_trackRef_pt_offline", &pfc_trackRef_pt_offline_);
  ttree_->Branch("pfc_trackRef_eta_offline", &pfc_trackRef_eta_offline_);
  ttree_->Branch("pfc_trackRef_phi_offline", &pfc_trackRef_phi_offline_);
  ttree_->Branch("pfc_trackRef_nValidPixelHits_offline", &pfc_trackRef_nValidPixelHits_offline_);
  ttree_->Branch("pfc_trackRef_nValidTrackerHits_offline", &pfc_trackRef_nValidTrackerHits_offline_);
  ttree_->Branch("Online_pfc_size", &Online_pfc_size_);
  ttree_->Branch("Offline_pfc_size", &Offline_pfc_size_);
  ttree_->Branch("diff_pfc_size", &diff_pfc_size_);

  ttree_->Branch("pfc_ecal_diff", &pfc_ecal_diff_);
  ttree_->Branch("pfc_hcal_diff", &pfc_hcal_diff_);
  ttree_->Branch("pfc_eta_diff", &pfc_eta_diff_);
  ttree_->Branch("pfc_phi_diff", &pfc_phi_diff_);
  ttree_->Branch("pfc_charge_diff", &pfc_charge_diff_);
  ttree_->Branch("pfc_id_diff", &pfc_id_diff_);
  ttree_->Branch("pfc_trackRef_p_diff", &pfc_trackRef_p_diff_);
  ttree_->Branch("pfc_trackRef_pt_diff", &pfc_trackRef_pt_diff_);
  ttree_->Branch("pfc_trackRef_eta_diff", &pfc_trackRef_eta_diff_);
  ttree_->Branch("pfc_trackRef_phi_diff", &pfc_trackRef_phi_diff_);
  ttree_->Branch("pfc_trackRef_nValidPixelHits_diff", &pfc_trackRef_nValidPixelHits_diff_);
  ttree_->Branch("pfc_trackRef_nValidTrackerHits_diff", &pfc_trackRef_nValidTrackerHits_diff_);

}  

PFHadCalibNTuple::~PFHadCalibNTuple() {
  edm::LogPrint("") << "----------------------------------------------------------";
  edm::LogPrint("") << moduleLabel_;
  edm::LogPrint("") << "----------------------------------------------------------";
  edm::LogPrint("") << "Total number of events: " << globalCounter_[0];
  edm::LogPrint("") << "Number of GEN isolated pions: " << globalCounter_[1];
  edm::LogPrint("") << "Number of pfSIM particles within dR = 0.01 of an isolated GEN pion: " << globalCounter_[2];
  edm::LogPrint("") << "Number of pfSIM particles with valid extrapolation point to ECAL surface: "
                    << globalCounter_[3];
  edm::LogPrint("") << "Number of PF candidates: " << globalCounter_[4];
  edm::LogPrint("") << "Number of PF Charged Hadrons: " << globalCounter_[5];
  edm::LogPrint("") << " - with pt > " << minPt_ << " GeV: " << globalCounter_[6];
  edm::LogPrint("") << " - with E_ECAL+E_HCAL > " << minCaloEnergy_ << " GeV: " << globalCounter_[7];
  if (usePFBlockElements_)
    edm::LogPrint("") << " - with only 1 track in the block: " << globalCounter_[8];
  edm::LogPrint("") << " - with track-p > " << minTrackP_ << " GeV and track-pT > " << minTrackPt_
                    << " GeV: " << globalCounter_[9];
  edm::LogPrint("") << " - with min nb of pixel hits: " << globalCounter_[10];
  edm::LogPrint("") << " - with min nb of pixel+strip hits: " << globalCounter_[11];
  edm::LogPrint("") << " - with E_ECAL < " << maxECalEnergy_ << " GeV: " << globalCounter_[12];
  edm::LogPrint("") << " Pass both online & offline cuts: " << globalCounter_[13];
}

void PFHadCalibNTuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  auto const& genParts = iEvent.get(genPartsToken_);
  auto const& pfSimParts = iEvent.get(pfSimPartsToken_);
  auto const& offline_pfSimParts = iEvent.get(offline_pfSimPartsToken_);
  auto const& recoPFCands = iEvent.get(recoPFCandsToken_);
  auto const& offline_recoPFCands = iEvent.get(offline_recoPFCandsToken_);

  reset_variables();
  const int NUM_CUTS = 13;
  std::vector<int> Is_Online_Single_pfc_pass_cuts(NUM_CUTS, false);
////  std::vector<int> Is_Offline_Single_pfc_pass_cuts(NUM_CUTS, false);

  ++globalCounter_[0];
  Is_Online_Single_pfc_pass_cuts[0] = true;
////  Is_Offline_Single_pfc_pass_cuts[0] = true;

  LogTrace("") << "----------------------------------------------------------";

  for (uint genpIdx_i = 0; genpIdx_i < genParts.size(); ++genpIdx_i) {
    auto const& genp_i = genParts.at(genpIdx_i);

    LogTrace("") << " genParticle[" << genpIdx_i << "]:"
                 << " pt=" << genp_i.pt() << " eta=" << genp_i.eta() << " phi=" << genp_i.phi()
                 << " pdgId=" << genp_i.pdgId() << " status=" << genp_i.status();

    // GEN and true(simPF) particle selection
    if (genp_i.status() == genParticleStatus_ and genp_i.pdgId() == genParticlePdgId_) {
      auto mindR = 999999.f;
      for (uint genpIdx_j = 0; genpIdx_j != genpIdx_i and genpIdx_j < genParts.size(); ++genpIdx_j) {
        auto const& genp_j = genParts.at(genpIdx_j);
        if (genp_j.status() == genParticleStatus_) {
          auto const dR = reco::deltaR(genp_i.eta(), genp_i.phi(), genp_j.eta(), genp_j.phi());
          if (dR < mindR)
            mindR = dR;
        }
      }
      // End of roof for calculating minimum dR

      // Isolated GEN particle
      if (mindR > genParticleIsoMinDeltaR_) {
        ++globalCounter_[1];
        Is_Online_Single_pfc_pass_cuts[1] = true;
////        Is_Offline_Single_pfc_pass_cuts[1] = true;

        for (auto const& ptc : pfSimParts) {
          // only consider negatively charged particles
          if (ptc.charge() >= 0)
            continue;

          // require true particle within deltaR = 0.01 of the gen particle
          auto const& gen = ptc.trajectoryPoint(reco::PFTrajectoryPoint::ClosestApproach);
          auto const dR = reco::deltaR(genp_i.eta(), genp_i.phi(), gen.momentum().Eta(), gen.momentum().Phi());
          if (dR > 0.01)
            continue;
          ++globalCounter_[2];
          Is_Online_Single_pfc_pass_cuts[2] = true;

          auto const& tpatecal = ptc.extrapolatedPoint(reco::PFTrajectoryPoint::ECALEntrance);
          if (not tpatecal.isValid())
            continue;

          auto const eta = tpatecal.positionREP().Eta();
          auto const phi = tpatecal.positionREP().Phi();
          auto const trueE = std::sqrt(tpatecal.momentum().Vect().Mag2());
          auto const true_dr = reco::deltaR(gen.momentum().Eta(), gen.momentum().Phi(), eta, phi);
          //auto const genE = std::sqrt(genp_i.momentum().Vect().Mag2());
          auto const genE = genp_i.energy();

          ++globalCounter_[3];
          Is_Online_Single_pfc_pass_cuts[3] = true;

          true_energy_.emplace_back(trueE);
          true_eta_.emplace_back(eta);
          true_phi_.emplace_back(phi);
          true_dr_.emplace_back(true_dr);
          true_charge_.emplace_back(ptc.charge());
          gen_energy_.emplace_back(genE);
        } 
	// End of [ pfSimParts ] roof : roof for Online PF Sim Particles 
        

        for (auto const& offline_ptc : offline_pfSimParts){
          ////if (offline_ptc.charge() >= 0)
          ////  continue;

          auto const& offline_gen = offline_ptc.trajectoryPoint(reco::PFTrajectoryPoint::ClosestApproach);
          auto const& offline_dR = reco::deltaR(genp_i.eta(), genp_i.phi(), offline_gen.momentum().Eta(), offline_gen.momentum().Phi());
          //if (offline_dR > 0.01)
            ////continue;
          //  Is_Offline_Single_pfc_pass_cuts[2] = false;
          //else Is_Offline_Single_pfc_pass_cuts[2] = true;

          auto const& offline_tpatecal = offline_ptc.extrapolatedPoint(reco::PFTrajectoryPoint::ECALEntrance);
          //if (not offline_tpatecal.isValid())
          ////  continue;
          //Is_Offline_Single_pfc_pass_cuts[3] = false;
          //else Is_Offline_Single_pfc_pass_cuts[3] = true;

          auto const offline_eta = offline_tpatecal.positionREP().Eta();
          auto const offline_phi = offline_tpatecal.positionREP().Phi();
          auto const offline_trueE = std::sqrt(offline_tpatecal.momentum().Vect().Mag2());
          auto const offline_true_dr = reco::deltaR(offline_gen.momentum().Eta(), offline_gen.momentum().Phi(), offline_eta, offline_phi);

          true_energy_offline_.emplace_back(offline_trueE);
          true_eta_offline_.emplace_back(offline_eta);
          true_phi_offline_.emplace_back(offline_phi);
          true_dr_offline_.emplace_back(offline_true_dr);
          true_charge_offline_.emplace_back(offline_ptc.charge());

        } 
	// End of [ offline_pfSimParts ] roof : roof for Offline PF Sim Particles
      }
      // End of [ mindR > genParticleIsoMinDeltaR_) ] if statement : if statement for selecting Isolated GEN particle
    }
    // End of [ genp_i.status() == genParticleStatus_ and genp_i.pdgId() == genParticlePdgId_ ] if statement
    // : if statement for selecting GEN and true(simPF) particle
  }
  // End of Gen Particles roof

  LogTrace("") << "----------------------------------------------------------";


  ////bool is_pass_cut_online = false;
  ////bool is_pass_cut_offline = false;


  // roof for Online PF Candidates: Reconstructed pion(pi-) selection
  
  Online_pfc_size_ = static_cast<int>(recoPFCands.size());
  Offline_pfc_size_ = static_cast<int>(offline_recoPFCands.size());
  diff_pfc_size_ = Online_pfc_size_ - Offline_pfc_size_;

  //for (auto const& pfc : recoPFCands) {
  for (size_t i = 0; i < recoPFCands.size(); i++) {

    if( static_cast<int>(i) >= Offline_pfc_size_ ) continue;

    const auto& pfc = recoPFCands[i];
    const auto& offline_pfc = offline_recoPFCands[i];

    ++globalCounter_[4];
    Is_Online_Single_pfc_pass_cuts[4] = true;

    if (pfc.particleId() != 1)
      continue;
    ++globalCounter_[5];
    Is_Online_Single_pfc_pass_cuts[5] = true;

    auto const ecalRaw = pfc.rawEcalEnergy();
    auto const hcalRaw = pfc.rawHcalEnergy();
    auto const ecalRaw_offline = offline_pfc.rawEcalEnergy();
    auto const hcalRaw_offline = offline_pfc.rawHcalEnergy();
    auto ecalRaw_diff = ecalRaw - ecalRaw_offline;
    auto hcalRaw_diff = hcalRaw - hcalRaw_offline;

    if ((ecalRaw + hcalRaw) < minCaloEnergy_)
      continue;
    ++globalCounter_[7];
    Is_Online_Single_pfc_pass_cuts[7] = true;

    auto nTracks = 0u;
    auto const& theElements = pfc.elementsInBlocks();
    if (theElements.empty()) {
      if (not usePFBlockElements_)
        nTracks = 1;  //!! hack for pfTICL (ref: https://github.com/cms-sw/cmssw/pull/32202)
    } else {
      auto const& elements = theElements[0].first->elements();
      for (unsigned iEle = 0; iEle < elements.size(); ++iEle) {
        if (elements[iEle].type() == reco::PFBlockElement::TRACK) {
          ++nTracks;
        }
      }
    }

    auto trackRef = pfc.trackRef();
    auto trackRef_offline = offline_pfc.trackRef();

    auto const track_p = trackRef->p();
    auto const track_pt = trackRef->pt();
    auto const track_eta = trackRef->eta();
    auto const track_phi = trackRef->phi();

    float track_p_offline = -50000.0;
    float track_pt_offline = -50000.0;
    float track_eta_offline = -50000.0;
    float track_phi_offline = -50000.0;
    if(trackRef_offline.isNonnull()){
        track_p_offline = trackRef_offline->p();
        track_pt_offline = trackRef_offline->pt();
        track_eta_offline = trackRef_offline->eta();
        track_phi_offline = trackRef_offline->phi();
    }
    auto track_p_diff = track_p - track_p_offline;
    auto track_pt_diff = track_pt - track_pt_offline;
    auto track_eta_diff = track_eta - track_eta_offline;
    auto track_phi_diff = track_phi - track_phi_offline;


    auto const& hp = trackRef->hitPattern();
    uint const track_nValidPixelHits = hp.numberOfValidPixelHits();
    uint const track_nValidTrackerHits = trackRef->numberOfValidHits();

    uint track_nValidPixelHits_offline = 1000;
    uint track_nValidTrackerHits_offline = 1000;
    if(trackRef_offline.isNonnull()){
        auto const& hp_offline = trackRef_offline->hitPattern();
        track_nValidPixelHits_offline = hp_offline.numberOfValidPixelHits();
        track_nValidTrackerHits_offline = trackRef_offline->numberOfValidHits();
    }
    uint track_nValidPixelHits_diff = track_nValidPixelHits - track_nValidPixelHits_offline;
    uint track_nValidTrackerHits_diff = track_nValidTrackerHits - track_nValidTrackerHits_offline;

    LogTrace("") << "----------------------------------------------------------";
    LogTrace("") << moduleLabel_ << " trackRef:";
    LogTrace("") << "     pt=" << trackRef->pt();
    LogTrace("") << "     eta=" << trackRef->eta();
    LogTrace("") << "     phi=" << trackRef->phi();
    LogTrace("") << "     algo=" << trackRef->algo();
    LogTrace("") << "     numberOfValidStripTOBHits=" << hp.numberOfValidStripTOBHits();
    LogTrace("") << "     numberOfValidStripTECHits=" << hp.numberOfValidStripTECHits();
    LogTrace("") << "     numberOfValidStripTIBHits=" << hp.numberOfValidStripTIBHits();
    LogTrace("") << "     numberOfValidStripTIDHits=" << hp.numberOfValidStripTIDHits();
    LogTrace("") << "     numberOfValidPixelBarrelHits=" << hp.numberOfValidPixelBarrelHits();
    LogTrace("") << "     numberOfValidPixelEndcapHits=" << hp.numberOfValidPixelEndcapHits();
    LogTrace("") << "     numberOfValidPixelHits=" << hp.numberOfValidPixelHits();
    LogTrace("") << "     numberOfValidStripHits=" << hp.numberOfValidStripHits();
    LogTrace("") << "     numberOfValidHits=" << trackRef->numberOfValidHits();
    LogTrace("") << "     numberOfLostHits=" << trackRef->numberOfLostHits();
    LogTrace("") << "----------------------------------------------------------";

    if (pfc.pt() < minPt_)
      continue;
    ++globalCounter_[6];
    Is_Online_Single_pfc_pass_cuts[6] = true;
    if (nTracks != 1)
      continue;
    ++globalCounter_[8];
    Is_Online_Single_pfc_pass_cuts[8] = true;
    if (track_p < minTrackP_ or track_pt < minTrackPt_)
      continue;
    ++globalCounter_[9];
    Is_Online_Single_pfc_pass_cuts[9] = true;


    auto hasMinPixelHits = false;
    auto hasMinTrackerHits = false;
    for (uint ieta = 0; ieta < maxEtaForMinTrkHitsCuts_.size(); ++ieta) {
      auto const etaMin = ieta ? maxEtaForMinTrkHitsCuts_[ieta - 1] : 0.;
      auto const etaMax = maxEtaForMinTrkHitsCuts_[ieta];

      if (std::abs(track_eta) >= etaMin and std::abs(track_eta) < etaMax) {
        hasMinPixelHits = track_nValidPixelHits >= minPixelHits_[ieta];
        hasMinTrackerHits = track_nValidTrackerHits >= minTrackerHits_[ieta];
        break;
      }
    }

    if (not hasMinPixelHits)
      continue;
    ++globalCounter_[10];
    Is_Online_Single_pfc_pass_cuts[10] = true;

    if (not hasMinTrackerHits)
      continue;
    ++globalCounter_[11];
    Is_Online_Single_pfc_pass_cuts[11] = true;

    if (ecalRaw > maxECalEnergy_)
      continue;
    ++globalCounter_[12];
    Is_Online_Single_pfc_pass_cuts[12] = true;

    // After Online cuts
    ////is_pass_cut_online = true;

    pfc_ecal_.emplace_back(ecalRaw);
    pfc_hcal_.emplace_back(hcalRaw);
    pfc_eta_.emplace_back(pfc.eta());
    pfc_phi_.emplace_back(pfc.phi());
    pfc_charge_.emplace_back(pfc.charge());
    pfc_id_.emplace_back(pfc.particleId());
    pfc_trackRef_p_.emplace_back(track_p);
    pfc_trackRef_pt_.emplace_back(track_pt);
    pfc_trackRef_eta_.emplace_back(track_eta);
    pfc_trackRef_phi_.emplace_back(track_phi);
    pfc_trackRef_nValidPixelHits_.emplace_back(track_nValidPixelHits);
    pfc_trackRef_nValidTrackerHits_.emplace_back(track_nValidTrackerHits);
    Is_Online_pfc_pass_cuts_.emplace_back(Is_Online_Single_pfc_pass_cuts);

    pfc_ecal_offline_.emplace_back(ecalRaw_offline);
    pfc_hcal_offline_.emplace_back(hcalRaw_offline);
    pfc_eta_offline_.emplace_back(offline_pfc.eta());
    pfc_phi_offline_.emplace_back(offline_pfc.phi());
    pfc_charge_offline_.emplace_back(offline_pfc.charge());
    pfc_id_offline_.emplace_back(offline_pfc.particleId());
    pfc_trackRef_p_offline_.emplace_back(track_p_offline);
    pfc_trackRef_pt_offline_.emplace_back(track_pt_offline);
    pfc_trackRef_eta_offline_.emplace_back(track_eta_offline);
    pfc_trackRef_phi_offline_.emplace_back(track_phi_offline);
    pfc_trackRef_nValidPixelHits_offline_.emplace_back(track_nValidPixelHits_offline);
    pfc_trackRef_nValidTrackerHits_offline_.emplace_back(track_nValidTrackerHits_offline);


    auto eta_diff = pfc.eta() - offline_pfc.eta();
    auto phi_diff = pfc.phi() - offline_pfc.phi();
    auto charge_diff = pfc.charge() - offline_pfc.charge();
    auto pfc_id_diff = pfc.particleId() - offline_pfc.particleId();
    pfc_ecal_diff_.emplace_back(ecalRaw_diff);
    pfc_hcal_diff_.emplace_back(hcalRaw_diff);
    pfc_eta_diff_.emplace_back(eta_diff);
    pfc_phi_diff_.emplace_back(phi_diff);
    pfc_charge_diff_.emplace_back(charge_diff);
    pfc_id_diff_.emplace_back(pfc_id_diff);
    pfc_trackRef_p_diff_.emplace_back(track_p_diff);
    pfc_trackRef_pt_diff_.emplace_back(track_pt_diff);
    pfc_trackRef_eta_diff_.emplace_back(track_eta_diff);
    pfc_trackRef_phi_diff_.emplace_back(track_phi_diff);
    pfc_trackRef_nValidPixelHits_diff_.emplace_back(track_nValidPixelHits_diff);
    pfc_trackRef_nValidTrackerHits_diff_.emplace_back(track_nValidTrackerHits_diff);

  } 
  // End of Online PF Candidates roof

//////////////////////////////////////////////////////////////////////////////
/*
  // roof for Offline PF Candidates: Reconstructed pion(pi-) selection
  for (auto const& pfc : offline_recoPFCands) {
    //Is_Offline_Single_pfc_pass_cuts[4] = true;

    ////if (pfc.particleId() != 1)
      ////continue;
    //Is_Offline_Single_pfc_pass_cuts[5] = false;
    //else Is_Offline_Single_pfc_pass_cuts[5] = true;

    auto const offline_ecalRaw = pfc.rawEcalEnergy();
    auto const offline_hcalRaw = pfc.rawHcalEnergy();
    ////if ((offline_ecalRaw + offline_hcalRaw) < minCaloEnergy_)
      ////continue;

    ////auto nTracks = 0u;
    ////auto const& theElements = pfc.elementsInBlocks();
    ////if (theElements.empty()) {
    ////  if (not usePFBlockElements_)
    ////    nTracks = 1;  //!! hack for pfTICL (ref: https://github.com/cms-sw/cmssw/pull/32202)
    ////} else {
    ////  auto const& elements = theElements[0].first->elements();
    ////  for (unsigned iEle = 0; iEle < elements.size(); ++iEle) {
    ////    if (elements[iEle].type() == reco::PFBlockElement::TRACK) {
    ////      ++nTracks;
    ////    }
    ////  }
    ////}

    auto trackRef = pfc.trackRef();

    float track_p = 34.0;
    float track_pt = 200.;
    float track_eta = 3.1;
    float track_phi = 3.1;

    if (trackRef.isNonnull() && trackRef.isAvailable()) {
     track_p = trackRef->p();
     track_pt = trackRef->pt();
     track_eta = trackRef->eta();
     track_phi = trackRef->phi();
    }

    //auto const& hp = trackRef->hitPattern();
    //uint const track_nValidPixelHits = hp.numberOfValidPixelHits();
    //uint const track_nValidTrackerHits = trackRef->numberOfValidHits();

    //if (pfc.pt() < minPt_)
      ////continue;
    //Is_Offline_Single_pfc_pass_cuts[6] = false;
    //else Is_Offline_Single_pfc_pass_cuts[6] = true;
    //if (nTracks != 1)
    //  ////continue;
    //Is_Offline_Single_pfc_pass_cuts[8] = false;
    //else Is_Offline_Single_pfc_pass_cuts[8] = true;
    //if (track_p < minTrackP_ or track_pt < minTrackPt_)
    //  ////continue;
    //Is_Offline_Single_pfc_pass_cuts[9] = false;
    //else Is_Offline_Single_pfc_pass_cuts[9] = true;


    ////auto hasMinPixelHits = false;
    ////auto hasMinTrackerHits = false;
    ////for (uint ieta = 0; ieta < maxEtaForMinTrkHitsCuts_.size(); ++ieta) {
    ////  auto const etaMin = ieta ? maxEtaForMinTrkHitsCuts_[ieta - 1] : 0.;
    ////  auto const etaMax = maxEtaForMinTrkHitsCuts_[ieta];

    ////  if (std::abs(track_eta) >= etaMin and std::abs(track_eta) < etaMax) {
    ////    hasMinPixelHits = track_nValidPixelHits >= minPixelHits_[ieta];
    ////    hasMinTrackerHits = track_nValidTrackerHits >= minTrackerHits_[ieta];
    ////    break;
    ////  }
    ////}

    ////if (not hasMinPixelHits)
    ////  ////continue;
    ////Is_Offline_Single_pfc_pass_cuts[10] = false;
    ////else Is_Offline_Single_pfc_pass_cuts[10] = true;

    ////if (not hasMinTrackerHits)
    ////  ////continue;
    ////Is_Offline_Single_pfc_pass_cuts[11] = false;
    ////else Is_Offline_Single_pfc_pass_cuts[11] = true;

    ////if (offline_ecalRaw > maxECalEnergy_)
    ////  ////continue;
    ////Is_Offline_Single_pfc_pass_cuts[12] = false;
    ////else Is_Offline_Single_pfc_pass_cuts[12] = true;

    // After Offline cuts (Offline cuts are same condition to Online. Only differences is offline just use Offline PF Candidates)
    ////is_pass_cut_offline = true;

    pfc_ecal_offline_.emplace_back(offline_ecalRaw);
    pfc_hcal_offline_.emplace_back(offline_hcalRaw);
    pfc_eta_offline_.emplace_back(pfc.eta());
    pfc_phi_offline_.emplace_back(pfc.phi());
    pfc_charge_offline_.emplace_back(pfc.charge());
    pfc_id_offline_.emplace_back(pfc.particleId());


    pfc_trackRef_p_offline_.emplace_back(track_p);
    pfc_trackRef_pt_offline_.emplace_back(track_pt);
    pfc_trackRef_eta_offline_.emplace_back(track_eta);
    pfc_trackRef_phi_offline_.emplace_back(track_phi);
    //pfc_trackRef_nValidPixelHits_offline_.emplace_back(track_nValidPixelHits);
    //pfc_trackRef_nValidTrackerHits_offline_.emplace_back(track_nValidTrackerHits);

    ////Is_Offline_pfc_pass_cuts_.emplace_back(Is_Offline_Single_pfc_pass_cuts);
  } 
  // End of Online PF Candidates roof

*/
/////////////////////////////////////////////////////////////////////////////

  ////if(!(is_pass_cut_online && is_pass_cut_offline)) return;
  ////// Event which pass both Online & Offine cuts,
  ////// will be filled at ntuple
  ////++globalCounter_[13];
 
  ttree_->Fill();
}

void PFHadCalibNTuple::reset_variables() {
  
  // Gen Particles
  gen_energy_.clear();

  // Online Sim Particles
  true_energy_.clear();
  true_eta_.clear();
  true_phi_.clear();
  true_dr_.clear();
  true_charge_.clear();

  // Online PF Candidates
  pfc_ecal_.clear();
  pfc_hcal_.clear();
  pfc_eta_.clear();
  pfc_phi_.clear();
  pfc_charge_.clear();
  pfc_id_.clear();
  pfc_trackRef_p_.clear();
  pfc_trackRef_pt_.clear();
  pfc_trackRef_eta_.clear();
  pfc_trackRef_phi_.clear();
  pfc_trackRef_nValidPixelHits_.clear();
  pfc_trackRef_nValidTrackerHits_.clear();

  Is_Online_pfc_pass_cuts_.clear();


  // Offline Sim Particles
  true_energy_offline_.clear();
  true_eta_offline_.clear();
  true_phi_offline_.clear();
  true_dr_offline_.clear();
  true_charge_offline_.clear();

  // Offline PF Candidates
  pfc_ecal_offline_.clear();
  pfc_hcal_offline_.clear();
  pfc_ecal_offline_.clear();
  pfc_hcal_offline_.clear();
  pfc_eta_offline_.clear();
  pfc_phi_offline_.clear();
  pfc_charge_offline_.clear();
  pfc_id_offline_.clear();
  pfc_trackRef_p_offline_.clear();
  pfc_trackRef_pt_offline_.clear();
  pfc_trackRef_eta_offline_.clear();
  pfc_trackRef_phi_offline_.clear();
  pfc_trackRef_nValidPixelHits_offline_.clear();
  pfc_trackRef_nValidTrackerHits_offline_.clear();

  ////Is_Offline_pfc_pass_cuts_.clear();
  Online_pfc_size_ = 0;
  Offline_pfc_size_ = 0;
  diff_pfc_size_ = 0;

  pfc_ecal_diff_.clear();
  pfc_hcal_diff_.clear();
  pfc_eta_diff_.clear();
  pfc_phi_diff_.clear();
  pfc_charge_diff_.clear();
  pfc_id_diff_.clear();
  pfc_trackRef_p_diff_.clear();
  pfc_trackRef_pt_diff_.clear();
  pfc_trackRef_eta_diff_.clear();
  pfc_trackRef_phi_diff_.clear();
  pfc_trackRef_nValidPixelHits_diff_.clear();
  pfc_trackRef_nValidTrackerHits_diff_.clear();

}

void PFHadCalibNTuple::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(PFHadCalibNTuple);
