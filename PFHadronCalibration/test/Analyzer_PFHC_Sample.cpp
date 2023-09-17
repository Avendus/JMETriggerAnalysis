#include <vector>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TF1.h>
#include <TF2.h>
#include "TH2F.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "Math/SMatrix.h"
#include "Math/SVector.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include <string>
#include <iostream>
#include <math.h>

using namespace std;

int main(){

  // User defined input variable
  TString SamplePath = "Sample/PFHadCalibration.root";
  TString TreePath_insideSample = "pfHadCalibNTuple/Candidates";
  TString OutputName = "PFHC_output.root";
  bool fullRun = true; 
  //bool fullRun = false;
  int SettedEntriesNum; if(!fullRun) SettedEntriesNum = 500000;  

//========================================================================================================


  // gen particle variable from [ genParticles ] tag
  vector<float>* gen_energy_= nullptr;

  // online variable
      // [ particleFlowSimParticle ] tag variable for true simulation information reference
      vector<float>* true_energy_= nullptr; 
      vector<float>* true_eta_= nullptr; 
      vector<float>* true_phi_= nullptr; 
      vector<float>* true_dr_= nullptr;  
      vector<float>* true_charge_= nullptr;
     
      // [ hltParticleFlow ] tag variable for reco PF Candidates (Online)
      vector<float>* pfc_ecal_= nullptr; // Reco clustered energy using PF algorithm at ECAL
      vector<float>* pfc_hcal_= nullptr; // Reco clustered energy using PF algorithm at HCAL
      vector<float>* pfc_eta_ = nullptr;
      vector<float>* pfc_phi_ = nullptr;
      vector<float>* pfc_charge_ = nullptr;
      vector<float>* pfc_id_ = nullptr;
      vector<float>* pfc_trackRef_p_= nullptr;
      vector<float>* pfc_trackRef_pt_ = nullptr;
      vector<float>* pfc_trackRef_eta_ = nullptr;
      vector<float>* pfc_trackRef_phi_ = nullptr;
      vector<float>* pfc_trackRef_nValidPixelHits_ = nullptr;
      vector<float>* pfc_trackRef_nValidTrackerHits_ = nullptr;


  // Offline variable
      // [ particleFlowSimParticle ] tag variable for true simulation information reference
      // same tag to online
      vector<float>* true_energy_offline_= nullptr; 
      vector<float>* true_eta_offline_= nullptr; 
      vector<float>* true_phi_offline_= nullptr; 
      vector<float>* true_dr_offline_= nullptr;  
      vector<float>* true_charge_offline_= nullptr;

      // [ particleFlow ] tag variable for reco PF Candidates (Offline)
      vector<float>* pfc_ecal_offline_= nullptr; // Reco clustered energy using PF algorithm at ECAL
      vector<float>* pfc_hcal_offline_= nullptr; // Reco clustered energy using PF algorithm at HCAL
      vector<float>* pfc_eta_offline_ = nullptr;
      vector<float>* pfc_phi_offline_ = nullptr;
      vector<float>* pfc_charge_offline_ = nullptr;
      vector<float>* pfc_id_offline_ = nullptr;
      vector<float>* pfc_trackRef_p_offline_= nullptr;
      vector<float>* pfc_trackRef_pt_offline_ = nullptr;
      vector<float>* pfc_trackRef_eta_offline_ = nullptr;
      vector<float>* pfc_trackRef_phi_offline_ = nullptr;
      vector<float>* pfc_trackRef_nValidPixelHits_offline_ = nullptr;
      vector<float>* pfc_trackRef_nValidTrackerHits_offline_ = nullptr;

//========================================================================================================

  TFile *fRead = new TFile( SamplePath, "READ" );
  TTree* sTree=(TTree*)fRead->Get( TreePath_insideSample );

  sTree->SetMakeClass(1);
 
  sTree->SetBranchAddress("gen_energy", &gen_energy_);

  sTree->SetBranchAddress("true_energy", &true_energy_);
  sTree->SetBranchAddress("true_eta", &true_eta_);
  sTree->SetBranchAddress("true_phi", &true_phi_);
  sTree->SetBranchAddress("true_dr", &true_dr_);
  sTree->SetBranchAddress("true_charge", &true_charge_);
  sTree->SetBranchAddress("pfc_ecal", &pfc_ecal_);
  sTree->SetBranchAddress("pfc_hcal", &pfc_hcal_);
  sTree->SetBranchAddress("pfc_eta", &pfc_eta_);
  sTree->SetBranchAddress("pfc_phi", &pfc_phi_);
  sTree->SetBranchAddress("pfc_charge", &pfc_charge_);
  sTree->SetBranchAddress("pfc_id", &pfc_id_);
  sTree->SetBranchAddress("pfc_trackRef_p", &pfc_trackRef_p_);
  //sTree->SetBranchAddress("pfc_trackRef_pt", &pfc_trackRef_pt_);
  //sTree->SetBranchAddress("pfc_trackRef_eta", &pfc_trackRef_eta_);
  //sTree->SetBranchAddress("pfc_trackRef_phi", &pfc_trackRef_phi_);
  //sTree->SetBranchAddress("pfc_trackRef_nValidPixelHits", &pfc_trackRef_nValidPixelHits_);
  //sTree->SetBranchAddress("pfc_trackRef_nValidTrackerHits", &pfc_trackRef_nValidTrackerHits_);

  sTree->SetBranchAddress("true_energy_offline", &true_energy_offline_);
  sTree->SetBranchAddress("true_eta_offline", &true_eta_offline_);
  sTree->SetBranchAddress("true_phi_offline", &true_phi_offline_);
  sTree->SetBranchAddress("true_dr_offline", &true_dr_offline_);
  sTree->SetBranchAddress("true_charge_offline", &true_charge_offline_);
  sTree->SetBranchAddress("pfc_ecal_offline", &pfc_ecal_offline_);
  sTree->SetBranchAddress("pfc_hcal_offline", &pfc_hcal_offline_);
  sTree->SetBranchAddress("pfc_eta_offline", &pfc_eta_offline_);
  sTree->SetBranchAddress("pfc_phi_offline", &pfc_phi_offline_);
  sTree->SetBranchAddress("pfc_charge_offline", &pfc_charge_offline_);
  sTree->SetBranchAddress("pfc_id_offline", &pfc_id_offline_);
  sTree->SetBranchAddress("pfc_trackRef_p_offline", &pfc_trackRef_p_offline_);
  //sTree->SetBranchAddress("pfc_trackRef_pt_offline", &pfc_trackRef_pt_offline_);
  //sTree->SetBranchAddress("pfc_trackRef_eta_offline", &pfc_trackRef_eta_offline_);
  //sTree->SetBranchAddress("pfc_trackRef_phi_offline", &pfc_trackRef_phi_offline_);
  //sTree->SetBranchAddress("pfc_trackRef_nValidPixelHits_offline", &pfc_trackRef_nValidPixelHits_offline_);
  //sTree->SetBranchAddress("pfc_trackRef_nValidTrackerHits_offline", &pfc_trackRef_nValidTrackerHits_offline_);


  TFile* output = new TFile(OutputName , "RECREATE");

  const int maxCadidatesNum = 10;
  TH1I* NumOfOnlinePFCandidates_eachEvt = new TH1I("NumOfOnlinePFCandidates_eachEvt", "Online PF candidates number at ecah event;Number of Online PF candidates;number of events", maxCadidatesNum, 1, maxCadidatesNum+1);
  TH1I* NumOfOfflinePFCandidates_eachEvt = new TH1I("NumOfOfflinePFCandidates_eachEvt", "Offline PF candidates number at ecah event;Number of Offline PF candidates;number of events", maxCadidatesNum, 1, maxCadidatesNum+1);

  TH2F* true_to_pfc_Dr_online_vs_offline = new TH2F("true_to_pfc_Dr_online_vs_offline", "Dr(true & pfc diff) distribution of online vs offline;offline dr;online dr",200, 0, 0.2, 200, 0, 0.2);

  TH1F* true_energy = new TH1F("true_energy", "True Energy;Energy [GeV];Events", 100, 0, 500);
  TH1F* ecal_energy = new TH1F("ecal_energy", "ECAL Energy[raw];Energy [GeV];Events", 100, 0, 500);
  TH1F* hcal_energy = new TH1F("hcal_energy", "HCAL Energy[raw];Energy [GeV];Events", 100, 0, 500);

  TH2F* Ratio_trueE_to_genE = new TH2F("Ratio_trueE_to_genE", "[ Ratio of true energy to gen energy ] vs genE;Energy [GeV];Ratio", 100, 0, 500, 100, 0, 1);

  TH2F* ecal_OnOffRatio_vs_Etrue_LogScale = new TH2F("ecal_OnOffRatio_vs_Etrue_LogScale", "ECAL Online Offline Logscale ratio vs Etrue;Energy [GeV];Online ECAL Energy / Offline ECAL Energy [LogScale]", 100, 0, 500, 105, -10.5, 10.5);
  TH2F* hcal_OnOffRatio_vs_Etrue_LogScale = new TH2F("hcal_OnOffRatio_vs_Etrue_LogScale", "HCAL Online Offline Logscale ratio vs Etrue;Energy [GeV];Online ECAL Energy / Offline ECAL Energy [LogScale]", 100, 0, 500, 105, -10.5, 10.5);
  TH2F* cal_OnOffRatio_vs_Etrue_LogScale = new TH2F("cal_OnOffRatio_vs_Etrue_LogScale", "CAL Online Offline Logscale ratio vs Etrue;Energy [GeV];Online ECAL Energy / Offline ECAL Energy [LogScale]", 100, 0, 500, 105, -10.5, 10.5);
  TH2F* ecal_OnOffRatio_vs_Etrue = new TH2F("ecal_OnOffRatio_vs_Etrue", "ECAL Online Offline ratio vs Etrue;Energy [GeV];Online ECAL Energy / Offline ECAL Energy", 100, 0, 500, 100, 0, 10);
  TH2F* hcal_OnOffRatio_vs_Etrue = new TH2F("hcal_OnOffRatio_vs_Etrue", "HCAL Online Offline ratio vs Etrue;Energy [GeV];Online HCAL Energy / Offline HCAL Energy", 100, 0, 500, 100, 0, 10);
  TH2F* cal_OnOffRatio_vs_Etrue = new TH2F("cal_OnOffRatio_vs_Etrue", "CAL Online Offline ratio vs Etrue;Energy [GeV];Online CAL Energy / Offline CAL Energy", 100, 0, 500, 100, 0, 10);

  // It's not mean REAL ERROR, this histogram is for compare on/offline measured energy
  // Used formula : e.g. (ecal - ecal_off) / ecal_off
  TH2F* ecal_OnOffRelativeError_vs_Etrue = new TH2F("ecal_OnOffRelativeError_Etrue", "ECAL Online Offline Relative Error vs Etrue; Energy [GeV];(online ecal - offline ecal) / offline ecal", 100, 0, 500, 105, -10.5, 10.5);
  TH2F* hcal_OnOffRelativeError_vs_Etrue = new TH2F("hcal_OnOffRelativeError_Etrue", "HCAL Online Offline Relative Error vs Etrue; Energy [GeV];(online hcal - offline hcal) / offline hcal", 100, 0, 500, 105, -10.5, 10.5);
  TH2F* cal_OnOffRelativeError_vs_Etrue = new TH2F("cal_OnOffRelativeError_Etrue", "CAL Online Offline Relative Error vs Etrue; Energy [GeV];(online cal - offline cal) / offline cal", 100, 0, 500, 105, -10.5, 10.5);

  TH2F* ecal_pfc_e_sub_true_e = new TH2F("ecal_pfc_e_sub_true_e", "ecal on/offline CAL Energy - True Energy;Online delta E [GeV];Offline delta E [GeV]", 100, 0, 50, 100, 0, 50);
  TH2F* hcal_pfc_e_sub_true_e = new TH2F("hcal_pfc_e_sub_true_e", "hcal on/offline CAL Energy - True Energy;Energy difference [GeV];Ratio", 100, 0, 50, 100, 0, 50);

  TH2F* cal_pfc_e_sub_true_e = new TH2F("cal_pfc_e_sub_true_e", "cal on/offline CAL Energy - True Energy;Online delta E [GeV];Offline delta E [GeV]", 100, 0, 50, 100, 0, 50);
  TH2F* cal_pfc_e_sub_true_e_above30 = new TH2F("cal_pfc_e_sub_true_e_above30", "cal on/offline CAL Energy above 30 - True Energy;Online delta E [GeV];Offline delta E [GeV]", 100, 0, 500, 100, 0, 500);


  bool flag[10] = {0,0,0,0,0,0,0,0,0,0};
  float e_gen, e_true;
  float minDr, Dr, true_index, minDr_offline, Dr_offline, true_index_offline;
  float ecal, hcal, ecal_off, hcal_off;
  float ecal_OnOffRatio, hcal_OnOffRatio, ecal_RelativeError, hcal_RelativeError, cal_OnOffRatio, cal_RelativeError;
  double sigmaEcalHcal = 1.0, sigmaEcalHcal_offline = 1.0;
  int NumberOfover1Candidates = 0;
  float online_ecal_e_diff_to_trueE = 0, offline_ecal_e_diff_to_trueE = 0;
  float online_hcal_e_diff_to_trueE = 0, offline_hcal_e_diff_to_trueE = 0;
  float ratio_ecal_e_diff_to_trueE = 0, ratio_hcal_e_diff_to_trueE = 0;
  float online_cal_e_diff_to_trueE = 0, offline_cal_e_diff_to_trueE = 0;

  vector<double> sigmas, sigmas_offline;

  cout << "Total Entries: "<< sTree->GetEntries() << endl;
  
  if(fullRun) SettedEntriesNum = sTree->GetEntries();
  else cout<< "You Setted Entries Num: "<< SettedEntriesNum << endl;

  for(int entry = 0; entry < SettedEntriesNum ; entry++) {
  // roof for event iteration
  
      sTree->GetEntry(entry);

      // Check for online PFCandidates & PFSimParticles existence
      if(pfc_id_->size() == 0 || true_energy_->size() == 0) continue;
      // Check for offline PFCandidates & PFSimParticles existence
      if(pfc_id_offline_->size() == 0 || true_energy_offline_->size() == 0) continue;

      // PF Candidates interation
      for(int i = 0; i < static_cast<int>(pfc_id_->size()); ++i){

          if(i > 0) {NumberOfover1Candidates++; continue;}
          
	  if(static_cast<int>(true_energy_->size()) != static_cast<int>(true_energy_offline_->size())){
	      cout<<"ERROR: Event num [ "<<entry+1<<" ]'s true energy online/offline size: "<<static_cast<int>(true_energy_->size())<<", "<<static_cast<int>(true_energy_offline_->size())<<endl;
              cout<<"ERROR: true variable of online/offline must be same, because they use same input tag"<<endl;
	      cout<<"ERROR: Check the tag & ntuple variable"<<endl;
	      exit(216);
	  }
	  NumOfOnlinePFCandidates_eachEvt->Fill(static_cast<int>(pfc_id_->size()));
	  NumOfOfflinePFCandidates_eachEvt->Fill(static_cast<int>(pfc_id_offline_->size()));
	  
	  // roof for finding online & offline minimum Dr & its index
          minDr = 99.;  Dr = 999.; true_index = -1;
          minDr_offline = 99.;  Dr_offline = 999.; true_index_offline = -1;
          for(int j = 0; j < static_cast<int>(true_energy_->size()); ++j){

              Dr = (pfc_eta_->at(i) - true_eta_->at(j))*(pfc_eta_->at(i) - true_eta_->at(j)) + (pfc_phi_->at(i) - true_phi_->at(j))*(pfc_phi_->at(i) - true_phi_->at(j));
              Dr_offline = (pfc_eta_offline_->at(i) - true_eta_offline_->at(j))*(pfc_eta_offline_->at(i) - true_eta_offline_->at(j)) + (pfc_phi_offline_->at(i) - true_phi_offline_->at(j))*(pfc_phi_offline_->at(i) - true_phi_offline_->at(j));

              if((minDr > Dr) && ((fabs(true_eta_->at(j)) <= 2.4 && pfc_trackRef_p_->at(i) > 0) || fabs(true_eta_->at(j)) > 2.4)) { minDr = Dr; true_index = j; }
              if((minDr_offline > Dr_offline) && ((fabs(true_eta_offline_->at(j)) <= 2.4 && pfc_trackRef_p_offline_->at(i) > 0) || fabs(true_eta_offline_->at(j)) > 2.4)) { minDr_offline = Dr_offline; true_index_offline = j; }

          }
	  // End of finding minimum Dr & its index roof


          if(true_index >= 0 && true_index_offline >= 0) {
              if(Dr < 1. && Dr_offline < 1.){
		  // Using true index, set the [ true information of candidates ]
		  e_gen = gen_energy_->at(true_index);
                  e_true = true_energy_->at(true_index);
                  ////etas.push_back(eta_->at(true_index));
                  ////phis.push_back(phi_->at(true_index));

                  if (e_true > 500) continue;
                  ecal = pfc_ecal_->at(i);
                  hcal = pfc_hcal_->at(i);
                  ////ecalEnergies.push_back(ecal);
                  ////hcalEnergies.push_back(hcal);

		  ecal_off = pfc_ecal_offline_->at(i);
		  hcal_off = pfc_hcal_offline_->at(i);

		  ecal_OnOffRatio = ecal / ecal_off;
		  hcal_OnOffRatio = hcal / hcal_off;
		  cal_OnOffRatio = (ecal + hcal) / (ecal_off + hcal_off);
		  ecal_RelativeError = (ecal - ecal_off) / ecal_off;
		  hcal_RelativeError = (hcal - hcal_off) / hcal_off;
                  cal_RelativeError = ((ecal+hcal) - (ecal_off+hcal_off)) / (ecal_off+hcal_off);

		  online_ecal_e_diff_to_trueE = fabs(ecal - e_true);
		  offline_ecal_e_diff_to_trueE = fabs(ecal_off - e_true);

		  online_hcal_e_diff_to_trueE = fabs(hcal - e_true);
		  offline_hcal_e_diff_to_trueE = fabs(hcal_off - e_true);

		  online_cal_e_diff_to_trueE = fabs((hcal+ecal) - e_true);
		  offline_cal_e_diff_to_trueE = fabs((hcal_off+ecal_off) - e_true);
		   

                  true_to_pfc_Dr_online_vs_offline->Fill(minDr_offline, minDr);
                  true_energy->Fill(e_true);
                  ecal_energy->Fill(ecal);
                  hcal_energy->Fill(hcal);           

		  ecal_OnOffRatio_vs_Etrue_LogScale->Fill(e_true, TMath::Log10(ecal_OnOffRatio));
                  hcal_OnOffRatio_vs_Etrue_LogScale->Fill(e_true, TMath::Log10(hcal_OnOffRatio));
                  cal_OnOffRatio_vs_Etrue_LogScale->Fill(e_true, TMath::Log10(cal_OnOffRatio));

		  ecal_OnOffRatio_vs_Etrue->Fill(e_true, ecal_OnOffRatio);
                  hcal_OnOffRatio_vs_Etrue->Fill(e_true, hcal_OnOffRatio);
                  cal_OnOffRatio_vs_Etrue->Fill(e_true, cal_OnOffRatio);

		  ////if(TMath::Log10(hcal_OnOffRatio) < -50){
		  //////if(TMath::Log10(hcal_OnOffRatio) > 50){
                  ////        cout<<"This might be the situation about particle deposit their energy only one calorimeter"<<endl;
		  ////        cout<<"ecal, hcal e: "<<ecal<<", "<<hcal<<endl;
		  ////        cout<<"off ecal, hcal e: "<<ecal_off<<", "<<hcal_off<<endl;
		  ////}

                  ecal_OnOffRelativeError_vs_Etrue->Fill(e_true, ecal_RelativeError);
                  hcal_OnOffRelativeError_vs_Etrue->Fill(e_true, hcal_RelativeError);
                  cal_OnOffRelativeError_vs_Etrue->Fill(e_true, cal_RelativeError);

		  Ratio_trueE_to_genE->Fill(e_gen, e_true/e_gen);

		  if(e_true < 30){
                      ecal_pfc_e_sub_true_e->Fill(online_ecal_e_diff_to_trueE, offline_ecal_e_diff_to_trueE);
                      hcal_pfc_e_sub_true_e->Fill(online_hcal_e_diff_to_trueE, offline_hcal_e_diff_to_trueE);
                      cal_pfc_e_sub_true_e->Fill(online_cal_e_diff_to_trueE, offline_cal_e_diff_to_trueE);
                  }
		  else{
                      cal_pfc_e_sub_true_e_above30->Fill(online_cal_e_diff_to_trueE, offline_cal_e_diff_to_trueE);
	          }

                  if(fabs(true_eta_->at(true_index)) < 1.5) sigmaEcalHcal = sqrt(0.08*0.08 + 1.04*1.04*(std::max((double)(ecal + hcal), 1.0)));
                  else sigmaEcalHcal = sqrt(0.04*0.04 + 1.80*1.80*(std::max((double)(ecal + hcal), 1.0)));
                  sigmas.push_back(sigmaEcalHcal);

		  if(fabs(true_eta_offline_->at(true_index)) < 1.5) sigmaEcalHcal_offline = sqrt(0.08*0.08 + 1.04*1.04*(std::max((double)(ecal_off + hcal_off), 1.0)));
                  else sigmaEcalHcal_offline = sqrt(0.04*0.04 + 1.80*1.80*(std::max((double)(ecal_off + hcal_off), 1.0)));

              }
	      // End of online/offline Dr cut
          }
	  // End of [ online/offline true index >= 0 ] cut
      }
      // End of PF Candidates roof

      unsigned N = sTree->GetEntriesFast();
      int frac = ((double)entry/N)*100;
      switch(frac) {
          case 10 : if (!flag[0]) { cout<<"Progress: 10%" <<endl; flag[0] = 1; } break;
	  case 20 : if (!flag[1]) { cout<<"Progress: 20%" <<endl; flag[1] = 1; } break;
	  case 30 : if (!flag[2]) { cout<<"Progress: 30%" <<endl; flag[2] = 1; } break;
	  case 40 : if (!flag[3]) { cout<<"Progress: 40%" <<endl; flag[3] = 1; } break;
	  case 50 : if (!flag[4]) { cout<<"Progress: 50%" <<endl; flag[4] = 1; } break;
	  case 60 : if (!flag[5]) { cout<<"Progress: 60%" <<endl; flag[5] = 1; } break;
	  case 70 : if (!flag[6]) { cout<<"Progress: 70%" <<endl; flag[6] = 1; } break;
	  case 80 : if (!flag[7]) { cout<<"Progress: 80%" <<endl; flag[7] = 1; } break;
	  case 90 : if (!flag[8]) { cout<<"Progress: 90%" <<endl; flag[8] = 1; } break;
	  case 99 : if (!flag[9]) { cout<<"Progress: 100%"<<endl; flag[9] = 1; } break;
          default : break;
      }

  }
  // End of [ entry ] for roof

  cout<<"Number of Evt of over 1 candidate : "<<NumberOfover1Candidates<<endl;

  output->Write();  

  return 0;

}
// End of main function
