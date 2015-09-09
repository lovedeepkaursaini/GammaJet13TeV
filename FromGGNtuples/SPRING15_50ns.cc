const int    photonID_IsConv[2][3]                = { {0, 0, 0} ,             {0, 0, 0}             };
const double photonID_HoverE[2][3]                = { {0.05, 0.05, 0.05} ,    {0.05, 0.05, 0.05}    };
const double photonID_SigmaIEtaIEta[2][3]         = { {0.0103, 0.0100, 0.0100} , {0.0277,0.0267,0.0267} };
const double photonID_RhoCorrR03ChHadIso[2][3]    = { {2.44,1.31,0.91} ,       {1.84,1.25,0.65}       };
const double photonID_RhoCorrR03NeuHadIso_0[2][3] = { {2.57,0.60,0.33} ,       {4.00,1.65,0.93}       };
const double photonID_RhoCorrR03NeuHadIso_1[2][3] = { {0.0044, 0.0044, 0.0044} ,    {0.0040, 0.0040, 0.0040}    };
const double photonID_RhoCorrR03NeuHadIso_2[2][3] = { {0.5809,0.5809,0.5809} ,    {0.9402,0.9402,0.9402}    };

const double photonID_RhoCorrR03PhoIso_0[2][3]    = { {1.92,1.33,0.61} ,       {2.15, 1.02, 0.54}       };
const double photonID_RhoCorrR03PhoIso_1[2][3]    = { {0.0043, 0.0043, 0.0043} , {0.0041, 0.0041, 0.0041} };

double anaGJet::dR(double eta1, double phi1, double eta2, double phi2) {
  double dphi = acos(cos(phi1 - phi2));
  double deta = eta1 - eta2;
  return sqrt(dphi*dphi + deta*deta);
}

bool anaGJet::fidEtaPass(double Eta){

  double fabsEta = TMath::Abs(Eta);
  if( fabsEta > 2.5) return false;
  if( 1.4442 < fabsEta && fabsEta < 1.566) return false;
  return true;
}

int anaGJet::phoRegion(double absEta){
  int region = 0;
  if( absEta >= 1.0  ) region++;
  if( absEta >= 1.479) region++;
  if( absEta >= 2.0  ) region++;
  if( absEta >= 2.2  ) region++;
  if( absEta >= 2.3  ) region++;
  if( absEta >= 2.4  ) region++;
  return region;
}
//https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonIdentificationRun2#SPRING15_selections_bunch_crossi
//// Selection implementation details for SPRING15 
double anaGJet::phoEffArea03ChHad(double phoEta){
  double eta = TMath::Abs(phoEta);
  static const double area[7] = {0.0158,0.0143,0.0115,0.0094,0.0095,0.0068,0.0053};
  return area[phoRegion(eta)];
}

double anaGJet::phoEffArea03NeuHad(double phoEta){
  double eta = TMath::Abs(phoEta);
  static const double area[7] = {0.0143,0.0210,0.0147,0.0082,0.0124,0.0186,0.0320};
  return area[phoRegion(eta)];
}

double anaGJet::phoEffArea03Pho(double phoEta){
  double eta = TMath::Abs(phoEta);
  static const double area[7] = {0.0725,0.0604,0.0320,0.0512,0.0766,0.0949,0.1160};
  return area[phoRegion(eta)];
}

bool anaGJet::passPhotonID(int phoInd, int pho_ID_ind = 0) {
//cout<<phoInd<<endl;
  // phoInd - index of the photon in the tree
  // pho_ID_ind: 0 -- loose, 1 -- medium, 2 -- tight
  double eta = (*phoSCEta)[phoInd];
  double et = (*phoEt)[phoInd];
  
  // rho corrected isolations
  double Pho03ChHadIso =     (*phoPFChIso)[phoInd]   - rho * phoEffArea03ChHad(eta);
  double Pho03NeuHadIso =    (*phoPFNeuIso)[phoInd]  - rho * phoEffArea03NeuHad(eta);
  double Pho03PhoIso =       (*phoPFPhoIso)[phoInd]  - rho * phoEffArea03Pho(eta);
  
  
  int region = 1; //barrel
  if(TMath::Abs( eta )< 1.4442) region = 0; //endcap
  bool phoPresel = fidEtaPass( eta ) &&
     et > 15 && //  Et cut
     ((*phoHoverE)[phoInd] < photonID_HoverE[region][pho_ID_ind]) &&
//     ((*phoSigmaIEtaIEta)[phoInd]<photonID_SigmaIEtaIEta[region][pho_ID_ind]) &&
     (Pho03NeuHadIso < (photonID_RhoCorrR03NeuHadIso_0[region][pho_ID_ind] + exp(photonID_RhoCorrR03NeuHadIso_1[region][pho_ID_ind] * et + photonID_RhoCorrR03NeuHadIso_2[region][pho_ID_ind]))) &&
     (Pho03PhoIso < (photonID_RhoCorrR03PhoIso_0[region][pho_ID_ind] + et * photonID_RhoCorrR03PhoIso_1[region][pho_ID_ind])) ;
  //   (Pho03ChHadIso < photonID_RhoCorrR03ChHadIso[region][pho_ID_ind])   ;
 /*cout<<" phoSCEta: "<< eta
<<", phoEt: "<<et
<<"phoEleVeto: "<<(*phoEleVeto)[phoInd]  
<<", phoHoverE12: "<<(*phoHoverE12)[phoInd]
<<", Pho03NeuHadIso: "<<Pho03NeuHadIso
<<", Pho03PhoIso: "<<Pho03PhoIso
<<endl;
*/
  return phoPresel;
}
