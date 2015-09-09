#define anaGJet_cxx
#include "anaGJet.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include <TMath.h>

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include "SPRING15_50ns.cc"
#include <fstream>
using namespace std;
float pho_pt;
float pho_phi;
float pho_eta;
float pho_SCEta;
float pho_Sihih;
float pho_ChHadIso;
float gen_weight;
float evt_weight =0;
int pho_hasPixelSeed ;
int pho_EleVeto ;

void anaGJet::clearVariables(){
  pho_pt = 0;
  pho_phi = 0;
  pho_eta = 0;
  pho_SCEta = 0;
  pho_Sihih = 0;
  pho_ChHadIso = 0;
  pho_hasPixelSeed = -9;
  pho_EleVeto = -9;
  gen_weight = -999;
}
void anaGJet::Loop(TString name)
{
  //
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  TFile* tmp = TFile::Open(name, "RECREATE");
  TTree* miniTree = new TTree("miniTree", "miniTree");
  miniTree->Branch("event", &event);
  miniTree->Branch("run", &run);
  miniTree->Branch("evt_weight", &evt_weight, "evt_weight/F");
  miniTree->Branch("gen_weight",&gen_weight,"gen_weight/F");
  miniTree->Branch("pho_pt", &pho_pt, "pho_pt/F");
  miniTree->Branch("pho_eta", &pho_eta, "pho_eta/F");
  miniTree->Branch("pho_SCEta", &pho_SCEta, "pho_SCEta/F");
  miniTree->Branch("pho_Sihih", &pho_Sihih, "pho_Sihih/F");
  miniTree->Branch("pho_ChHadIso", &pho_ChHadIso, "pho_ChHadIso/F");
  miniTree->Branch("pho_hasPixelSeed",&pho_hasPixelSeed,"pho_hasPixelSeed/I");
  miniTree->Branch("pho_EleVeto",&pho_EleVeto,"pho_EleVeto/I");

  TH1F* hCounter=new TH1F("hCounter","hCounter",20,0.,20);
  TH1F* hCounter_wt=new TH1F("hCounter_wt","hCounter_wt",20,0.,20);
  TH1D*  hGenWeightedEvents_ = new TH1D("hGenWeightedEvents", "total weighted processed", 2, 0, 2);
  TH1D*  hEvents_ = new TH1D("hEvents", "total processed", 2, 0, 2);
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    hEvents_->Fill(1);
    hGenWeightedEvents_->Fill(1,gen_weight);
    if(jentry % 1000 == 0) cout << "Processed " << jentry
                                << " events out of " <<nentries<< endl;
    clearVariables();
    event = event;
    run = run;
    gen_weight = genWeight;
    hCounter->Fill(0);
    hCounter_wt->Fill(0.,genWeight);
    evt_weight = evt_weight + gen_weight;
    if(nPho<1)continue;
    hCounter->Fill(1);
    hCounter_wt->Fill(1.,genWeight);
    vector <int> iphotons;
    for (int ipho = 0; ipho < nPho; ++ipho){
      // PRE-PHOTON SELECTION
      if((*phoEt)[ipho] < 15) continue;
      if( fabs((*phoSCEta)[ipho])>2.5) continue;
      if( fabs((*phoSCEta)[ipho])<1.566 && fabs((*phoSCEta)[ipho])>1.4442) continue;
      if (!passPhotonID(ipho, 1)) continue;
      iphotons.push_back(ipho);
    }
    if(iphotons.size() < 1 ) continue;
    
    hCounter_wt->Fill(2.,genWeight);
    hCounter->Fill(2);
    
    TLorentzVector pho;
    int j = iphotons[0];
    pho.SetPtEtaPhiM((*phoEt)[j], (*phoEta)[j], (*phoPhi)[j], 0);
    
    pho_pt = pho.Pt();
    pho_eta = pho.Eta();
    pho_SCEta = (*phoSCEta)[j];
    pho_ChHadIso = (*phoPFChIso)[j]   - rho * phoEffArea03ChHad(pho_eta);
    pho_Sihih = (*phoSigmaIEtaIEtaFull5x5)[j];
    pho_hasPixelSeed = phohasPixelSeed->at(j);
    pho_EleVeto = (*phoEleVeto)[j];
    
    miniTree->Fill();
    
  }
  miniTree->Write();
  cout<<evt_weight<<endl;
  hGenWeightedEvents_->Write();
  hEvents_->Write();
  
  tmp->Write();
  tmp->Close();
}


int main(int argc, char* argv[]) {
  TString fName = argv[1];
  
  anaGJet t(fName);
  t.Loop("minitree_"+fName+".root");
  return 0;
}


