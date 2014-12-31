#ifndef _HelperLeptonCore_h
#define _HelperLeptonCore_h

#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH3.h"
#include "TH2F.h"
#include "TF1.h"
#include "TF2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TAxis.h"
#include "TMath.h"
#include "TRandom3.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <exception>
#include <cmath>
#include <random>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
//Root includes
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include <TRandom3.h>

#if !defined(__CINT__) && !defined(__MAKECINT__)

#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"

#include "FWCore/Framework/interface/Event.h"

//Headers for the data items
#include "BEAN/Collections/interface/BNelectron.h"
#include "BEAN/Collections/interface/BNevent.h"
#include "BEAN/Collections/interface/BNjet.h"
#include "BEAN/Collections/interface/BNmcparticle.h"
#include "BEAN/Collections/interface/BNmet.h"
#include "BEAN/Collections/interface/BNmuon.h"
#include "BEAN/Collections/interface/BNphoton.h"
#include "BEAN/Collections/interface/BNsupercluster.h"
#include "BEAN/Collections/interface/BNtau.h"
#include "BEAN/Collections/interface/BNtrack.h"
#include "BEAN/Collections/interface/BNtrigger.h"
#include "BEAN/Collections/interface/BNskimbits.h"
#include "BEAN/Collections/interface/BNtrigobj.h"
#include "BEAN/Collections/interface/BNprimaryvertex.h"

#include "BEAN/BEANmaker/interface/BtagWeight.h"
#include "BEAN/BEANmaker/interface/BEANhelper.h"

#include "BEAN/BEANmaker/interface/AnglesUtil.h"
#include "PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h"

#endif

#include "ttHMultileptonAnalysis/TemplateMakers/interface/CERN/RochCor2012.h"
#include "ttHMultileptonAnalysis/TemplateMakers/interface/BEANFileInterface.h"
#include "ttHMultileptonAnalysis/TemplateMakers/interface/KinematicVariable.h"
#include "ttHMultileptonAnalysis/TemplateMakers/interface/TwoObjectKinematic.h"
#include "ttHMultileptonAnalysis/TemplateMakers/interface/GenericCollectionMember.h"
#include "ttHMultileptonAnalysis/TemplateMakers/interface/JobParameters.h"

class HelperLeptonCore  {
public:
  HelperLeptonCore();

  BEANhelper * setupAnalysisParameters(string year, string sampleName); // initialize
  void detectData(string sampleName);
  int convertSampleNameToNumber(string sampleName);
  void initializePUReweighting(); // setup some PU reweighting flags
  void initializeInputCollections(edm::EventBase&, bool, BEANFileInterface&);

  // Handle the gymnastics of tight and loose collection definitions
  void getTightLoosePreselectedElectrons(electronID::electronID tightID,
                                         electronID::electronID looseID,
                                         electronID::electronID preselectedID,
                                         BEANFileInterface* selectedCollections);

  void getTightLoosePreselectedMuons(muonID::muonID tightID,
                                     muonID::muonID looseID,
                                     muonID::muonID preselectedID,
                                     BEANFileInterface* selectedCollections);

  void getTightLoosePreselectedTaus(tauID::tauID tightID,
                                         tauID::tauID looseID,
                                         tauID::tauID preselectedID,
                                         BEANFileInterface* selectedCollections);

  void getTightCorrectedJets(double ptCut,
                             double etaCut,
                             jetID::jetID tightID,
                             BEANFileInterface* selectedCollections);

  BNjetCollection * getCorrectedSelectedJets(double ptCut,
                                             double etaCut,
                                             jetID::jetID jetID,
                                             const char csvWorkingPoint);

  void getCorrectedMet(BEANFileInterface * selectedCollections,
                       sysType::sysType jetSyst = sysType::NA);

  void fillLepCollectionWithSelectedLeptons(BEANFileInterface * selectedCollections);
  void fillZLepCollectionWithSelectedLeptons(BEANFileInterface * selectedCollections,
                                             TwoObjectKinematic<BNleptonCollection,BNleptonCollection> * myZLikeMassLepLepSFOS_tight,
                                             TwoObjectKinematic<BNleptonCollection,BNleptonCollection> * myZLikeMassLepLepSFOS_tightLoose,
                                             TwoObjectKinematic<BNleptonCollection,BNleptonCollection> * myZLikeMassLepLepSFOS_all);

  bool isFromB(BNmcparticle particle);

  double scaleIPVarsMC(double ipvar, int pdgId, double pt, double eta, int mcMatchId, int mcMatchAny);
  double scaleSIPMC(double& sip, int& genID, double& pt, int& mcMatchID, int& mcMatchAny, double& eta);
  double scaleDZMC(double dz, int genID, double pt, double eta, int mcMatchID, int mcMatchAny);
  double scaleDXYMC(double dxy, int genID, double pt, double eta, int mcMatchID, int mcMatchAny);
  double scaleLepJetPtRatioMC(double jetPtRatio, int genID, double pt, double eta, int mcMatchID, int mcMatchAny);
  double scaleLepJetDRMC(double jetDR, int genID, double pt, double eta, int mcMatchID, int mcMatchAny);

  template <typename collectionType> void fillMCMatchID(collectionType& collection, BNmcparticleCollection& mcParticles, const double& maxDR);
  template <typename collectionType> void fillMCMatchAny(collectionType& collection, BNmcparticleCollection& mcParticles, const double& maxDR);
  template <typename collectionType> void fillSIP(collectionType& collection, bool applySmearing);
  template <typename collectionType> void fillLepJetPtRatio(collectionType& collection, BNjetCollection& jetCollection, bool applySmearing);
  template <typename collectionType> void fillLepJetDeltaR(collectionType& collection, BNjetCollection& jetCollection, bool applySmearing);
  template <typename collectionType> void fillLepJetBTagCSV(collectionType& collection, BNjetCollection& jetCollection);
  template <typename collectionType> void scaleMCCollectionDZ(collectionType& collection);
  template <typename collectionType> void scaleMCCollectionDXY(collectionType& collection);
  template <typename particleType> bool isPlausible(particleType particle, const BNmcparticle& mcParticle);
  template <typename particleType> BNmcparticleCollection getPlausibleParticles(particleType& particle, BNmcparticleCollection mcParticles);
  template <typename collectionType> void applyRochesterCorrections(collectionType& collection);
  template <typename collectionType> void apply_SMP_12_011_correction_mu(collectionType& collection, bool doSmear);
  template <typename collectionType> void apply_SMP_12_011_correction_ele(collectionType& collection, bool doSmear);
  template <typename collectionType> void apply_gsf_Et(collectionType& collection);
  template <typename collectionType> void apply_sc_Et(collectionType& collection);
  template <typename collectionType> void shiftLeptonEnergy(collectionType& collection, double shift_scale_B, double shift_scale_E1, double shift_scale_E2);
  template <typename collectionType> void smearLeptonEnergy(collectionType& collection, double smear_size);
  template <typename particleType> BNjet GetClosestJet(const BNjetCollection& jets, const particleType& particle, const double maxDeltaR);

  TRandom *gSmearer;

  //-------------------- Variables
  std::string analysisYear;

  double weight_Xsec;
  int nGen;
  float Xsec;
  int sampleNumber;
  string sampleName;

  bool isData;
  string dataset;

  sysType::sysType jetEnergyShift;
  sysType::sysType csvShift;

  BEANhelper bHelp;
  BEANFileInterface rawCollections;

  JobParameters config;

  bool verbose;

  std::string listOfCollisionDatasets;
  std::string datasetForBEANHelper;

  edm::Handle<BNeventCollection> h_event;
  BNeventCollection events;

  edm::Handle<BNmuonCollection> h_muons;
  BNmuonCollection muonsRaw;
  BNmuonCollection muonsTight;
  BNmuonCollection muonsLoose;
  BNmuonCollection muonsPreselected;
  BNmuonCollection muonsTightLoose;
  BNmuonCollection muonsLoosePreselected;
  BNmuonCollection muonsTightLoosePreselected;

  edm::Handle<BNmcparticleCollection> h_mcparticles;
  BNmcparticleCollection mcparticles;

  edm::Handle<BNjetCollection> h_pfjets;
  BNjetCollection sortedCorrSelJets;
  BNjetCollection pfjets;
  BNjetCollection jetsTight;
  BNjetCollection jetsByCSVTight;
  BNjetCollection jetsLooseCSV;
  BNjetCollection jetsTightCSV;
  BNjetCollection jetsNotLooseCSV;
  BNjetCollection jetsNotMediumCSV;
  BNjetCollection jetsNotTightCSV;
  BNjetCollection jetsMediumCSV;

  edm::Handle<BNmetCollection> h_pfmets;
  BNmetCollection pfmets;
  BNmetCollection metCorrected;
  edm::Handle<BNmetCollection> h_pfType1CorrectedMetBN;
  BNmetCollection metpfType1CorrectedMetBN;

  edm::Handle<BNtriggerCollection> h_hlt;
  BNtriggerCollection hltInfo;

  edm::Handle<BNprimaryvertexCollection> h_pvs;
  BNprimaryvertexCollection pvs;

  edm::Handle<BNelectronCollection> h_electrons;
  BNelectronCollection electronsRaw;
  BNelectronCollection electronsTight;
  BNelectronCollection electronsLoose;
  BNelectronCollection electronsPreselected;
  BNelectronCollection electronsTightLoose;
  BNelectronCollection electronsLoosePreselected;
  BNelectronCollection electronsTightLoosePreselected;

  BNleptonCollection leptonsRaw;
  BNleptonCollection leptonsTight;
  BNleptonCollection leptonsLoose;
  BNleptonCollection leptonsPreselected;
  BNleptonCollection leptonsTightLoose;
  BNleptonCollection leptonsLoosePreselected;
  BNleptonCollection leptonsTightLoosePreselected;
  BNleptonCollection leptonsTightZ;
  BNleptonCollection leptonsTightNonZ;
  BNleptonCollection leptonsTightLooseZ;
  BNleptonCollection leptonsTightLooseNonZ;
  BNleptonCollection leptonsTightLoosePreselectedZ;
  BNleptonCollection leptonsTightLoosePreselectedNonZ;

  edm::Handle<BNtauCollection> h_taus;
  BNtauCollection tausRaw;
  BNtauCollection tausTight;
  BNtauCollection tausLoose;
  BNtauCollection tausPreselected;
  BNtauCollection tausTightLoose;
  BNtauCollection tausLoosePreselected;
  BNtauCollection tausTightLoosePreselected;

  edm::Handle<BNjetCollection> h_lepMvaJets;
  BNjetCollection lepMvaJets;

  BNmcparticleCollection higgsParticles;

//-------------------- Inline functions
  inline double smearMC(TRandom* gSmearer, double x, double mu, double sigma) {
    if (x == 0) return gSmearer->Gaus(mu,sigma);
    else return (x/abs(x))*(abs(x) + gSmearer->Gaus(mu,sigma));
  }
  inline double logSmearMC(TRandom* gSmearer, double x, double mu, double sigma) {
    if (x == 0) return std::exp(gSmearer->Gaus(mu,sigma));
    else return (x/abs(x))*std::exp(std::log(abs(x)) + gSmearer->Gaus(mu,sigma));
  }
  inline double shiftMC(double x, double delta) {
    if (x == 0) return delta;
    else return (x/abs(x))*(abs(x) + delta);
  }
  inline double scaleShiftMC(double x, double scale, double shift) {
    if (x == 0) return shift;
    else return (x/abs(x))*(abs(x)*scale + shift);
  }

};

//-------------------- Template functions
template <typename collectionType>
void HelperLeptonCore::fillMCMatchID(collectionType& collection, BNmcparticleCollection& mcParticles, const double& maxDR) {
  BNmcparticleCollection plausibleParticles;
  BNmcparticle matchedParticle;

  for (auto& object: collection) {
    plausibleParticles = getPlausibleParticles(object, mcParticles);
    matchedParticle = bHelp.GetMatchedMCparticle(plausibleParticles, object, maxDR);
    object.mcMatchID = matchedParticle.id;
  }
}

template <typename collectionType>
void HelperLeptonCore::fillMCMatchAny(collectionType& collection, BNmcparticleCollection& mcParticles, const double& maxDR) {
  BNmcparticleCollection plausibleParticles;
  BNmcparticle matchedParticle;

  for (auto& object: collection) {
    plausibleParticles = getPlausibleParticles(object, mcParticles);
    matchedParticle = bHelp.GetMatchedMCparticle(plausibleParticles, object, maxDR);
    object.mcMatchAny = 1 + int(isFromB(matchedParticle));
  }
}

template <typename particleType>
BNmcparticleCollection HelperLeptonCore::getPlausibleParticles(particleType& particle, BNmcparticleCollection mcParticles) {
  BNmcparticleCollection plausibleParticles;

  for (auto mcParticle: mcParticles) {
    if (isPlausible(particle, mcParticle)) plausibleParticles.push_back(mcParticle);
  }
  return plausibleParticles;
}

template <typename particleType>
bool HelperLeptonCore::isPlausible(particleType particle, const BNmcparticle& mcParticle) {
  if ((abs(particle.genId) == 11) && (abs(mcParticle.id != 11))) return false;
  if ((abs(particle.genId) == 13) && (abs(mcParticle.id != 13))) return false;
  double dR = deltaR(particle.eta, particle.phi, mcParticle.eta, mcParticle.phi);
  if (dR < 0.3) return true;
  if ((particle.pt < 10) && (abs(particle.genId) == 13) && (mcParticle.id != particle.genId)) return false;
  if (dR<0.7) return true;
  if (min(particle.pt, mcParticle.pt) / max(particle.pt, mcParticle.pt) < 0.3) return false;
  return true;
}

template <typename collectionType>
void HelperLeptonCore::fillSIP(collectionType& collection, bool applySmearing) {
  double sip = -99.0;
  for (auto& object: collection) {
    if (object.IP != -99 && object.IPError != -99) sip = object.IP / object.IPError;
    if (applySmearing) sip = scaleSIPMC(sip, object.genId, object.pt, object.mcMatchID, object.mcMatchAny, object.eta);
    object.SIP = sip;
  }
}

template <typename collectionType>
void HelperLeptonCore::fillLepJetBTagCSV(collectionType& collection, BNjetCollection& jetCollection) {
  BNjet matchedJet;

  for (auto& object: collection) {
    matchedJet = GetClosestJet(jetCollection, object, 0.5);
    object.jetBTagCSV = matchedJet.btagCombinedSecVertex;
  }
}

template <typename collectionType>
void HelperLeptonCore::fillLepJetPtRatio(collectionType& collection, BNjetCollection& jetCollection, bool applySmearing) {
  double jetPtRatio = 1.5;
  BNjet matchedJet;

  for (auto& object: collection) {
    matchedJet = GetClosestJet(jetCollection, object, 0.5);
    jetPtRatio = object.pt/matchedJet.pt;
    if (applySmearing) jetPtRatio = scaleLepJetPtRatioMC(jetPtRatio, object.genId, object.pt, object.eta, object.mcMatchID, object.mcMatchAny);

    object.jetPtRatio = jetPtRatio;
  }
}

template <typename collectionType>
void HelperLeptonCore::fillLepJetDeltaR(collectionType& collection, BNjetCollection& jetCollection, bool applySmearing) {
  double jetDeltaR = -99;
  BNjet matchedJet;

  for (auto& object: collection) {
    matchedJet = GetClosestJet(jetCollection, object, 0.5);
    jetDeltaR = deltaR(object.eta, object.phi, matchedJet.eta, matchedJet.phi);
    if (applySmearing) jetDeltaR = scaleLepJetDRMC(jetDeltaR, object.genId, object.pt, object.eta, object.mcMatchID, object.mcMatchAny);

    object.jetDeltaR = jetDeltaR;
  }
}

template <typename collectionType>
void HelperLeptonCore::scaleMCCollectionDZ(collectionType& collection) {
  double dz = -99.0;
  for (auto& object: collection) {
    dz = scaleDZMC(object.correctedDZ, object.genId, object.pt, object.eta, object.mcMatchID, object.mcMatchAny);
    object.correctedDZ = dz;
  }
}

template <typename collectionType>
void HelperLeptonCore::scaleMCCollectionDXY(collectionType& collection) {
  double dxy = -99.0;
  for (auto& object: collection) {
    dxy = scaleDXYMC(object.correctedD0Vertex, object.genId, object.pt, object.eta, object.mcMatchID, object.mcMatchAny);
    object.correctedD0Vertex = dxy;
  }
}

template <typename collectionType>
void HelperLeptonCore::applyRochesterCorrections(collectionType& collection) {
  RochCor2012 rochesterCorrections2012;
  TLorentzVector p4;

  for (auto& object: collection) {
    p4.SetPxPyPzE(object.px, object.py, object.pz, object.energy);
    if (isData) rochesterCorrections2012.momcor_data(p4, object.tkCharge, 0.0, 0);
    if (!isData) rochesterCorrections2012.momcor_mc(p4, object.tkCharge, 0.0, 0);

    //std::cout << "uncorrected pt: " << object.pt << "  corrected pt: " << p4.Pt() << std::endl;

    object.px = p4.Px();
    object.py = p4.Py();
    object.pz = p4.Pz();
    object.energy = p4.Energy();
    object.pt = p4.Pt();
  }
}

template <typename collectionType>
void HelperLeptonCore::apply_SMP_12_011_correction_mu(collectionType& collection, bool doSmear) {

  //Correction factors from SMP-12-011 (AN-12-067)
  //http://cms.cern.ch/iCMS/analysisadmin/viewanalysis?id=862&field=id&value=862&name=Inclusive%20W/Z%20cross%20section%20at%208%20TeV
  //http://cds.cern.ch/record/1646590?ln=en

  TLorentzVector p4;
  double pi = 3.14159265358979323846;
  double bin = 2*pi/11; //11 bins over a 2 pi range
  double phi = -99.0;
  double scale = 1.0;

  double smear_size = 0.5; //Quoted value
  double smear_val = 0.0;
  double smear_scale = 1.0;

  for (auto& object: collection) {
    p4.SetPxPyPzE(object.px, object.py, object.pz, object.energy);
    phi = p4.Phi();
    if (isData) continue;
    if (object.tkCharge == -1) {
      if      (phi >= 0*bin - pi && phi < 1*bin - pi) scale = 1.005;
      else if (phi >= 1*bin - pi && phi < 2*bin - pi) scale = 1.0075;
      else if (phi >= 2*bin - pi && phi < 3*bin - pi) scale = 1.006;
      else if (phi >= 3*bin - pi && phi < 4*bin - pi) scale = 1.001;
      else if (phi >= 4*bin - pi && phi < 5*bin - pi) scale = 0.999;
      else if (phi >= 5*bin - pi && phi < 6*bin - pi) scale = 0.9985;
      else if (phi >= 6*bin - pi && phi < 7*bin - pi) scale = 1.0015;
      else if (phi >= 7*bin - pi && phi < 8*bin - pi) scale = 0.9995;
      else if (phi >= 8*bin - pi && phi < 9*bin - pi) scale = 0.999;
      else if (phi >= 9*bin - pi && phi < 10*bin - pi) scale = 0.9965;
      else if (phi >= 10*bin - pi && phi <= 11*bin - pi) scale = 1.002;
    }
    else if (object.tkCharge == 1) {
      if      (phi >= 0*bin - pi && phi < 1*bin - pi) scale = 1.001;
      else if (phi >= 1*bin - pi && phi < 2*bin - pi) scale = 1.0025;
      else if (phi >= 2*bin - pi && phi < 3*bin - pi) scale = 0.998;
      else if (phi >= 3*bin - pi && phi < 4*bin - pi) scale = 0.999;
      else if (phi >= 4*bin - pi && phi < 5*bin - pi) scale = 1.000;
      else if (phi >= 5*bin - pi && phi < 6*bin - pi) scale = 1.001;
      else if (phi >= 6*bin - pi && phi < 7*bin - pi) scale = 1.0015;
      else if (phi >= 7*bin - pi && phi < 8*bin - pi) scale = 1.004;
      else if (phi >= 8*bin - pi && phi < 9*bin - pi) scale = 1.003;
      else if (phi >= 9*bin - pi && phi < 10*bin - pi) scale = 1.002;
      else if (phi >= 10*bin - pi && phi <= 11*bin - pi) scale = 0.999;
    }

    smear_val = gSmearer->Gaus(0.0, smear_size);    
    if (doSmear) smear_scale = (p4.Energy() + smear_val)/p4.Energy();
      
    // MC correction divide by scale
    object.px = p4.Px()*smear_scale/scale;
    object.py = p4.Py()*smear_scale/scale;
    object.pz = p4.Pz()*smear_scale/scale;
    object.energy = p4.Energy()*smear_scale/scale;
    object.pt = p4.Pt()*smear_scale/scale;
    
  } //End loop over objects

} //End apply_SMP_12_011_correction_mu

template <typename collectionType>
void HelperLeptonCore::apply_SMP_12_011_correction_ele(collectionType& collection, bool doSmear) {

  //Correction factors from SMP-12-011 (AN-12-067)
  //http://cms.cern.ch/iCMS/analysisadmin/viewanalysis?id=862&field=id&value=862&name=Inclusive%20W/Z%20cross%20section%20at%208%20TeV
  //http://cds.cern.ch/record/1646590?ln=en

  TLorentzVector p4;
  double absEta = -99.0;
  double scale = 1.0;
  
  double smear_size = 0.0;
  double smear_val = 0.0;
  double smear_scale = 1.0;

  for (auto& object: collection) {
    p4.SetPxPyPzE(object.px, object.py, object.pz, object.energy);
    absEta = abs(p4.Eta());
    if (isData) continue;
    
    if      (absEta >= 0.0 && absEta < 0.4) { scale = 1.0024; smear_size = 0.46; }
    else if (absEta >= 0.4 && absEta < 0.8) { scale = 1.0055; smear_size = 0.01; } 
    else if (absEta >= 0.8 && absEta < 1.2) { scale = 1.0063; smear_size = 0.82; }
    else if (absEta >= 1.2 && absEta < 1.4442) { scale = 1.0067; smear_size = 0.71; }
    else if (absEta >= 1.566 && absEta < 2.0) { scale = 0.9985; smear_size = 1.36; }
    else if (absEta >= 2.0 && absEta <= 2.5) { scale = 0.9835; smear_size = 1.28; } 

    smear_val = gSmearer->Gaus(0.0, smear_size);    
    if (doSmear) smear_scale = (p4.Energy() + smear_val)/p4.Energy();
      
    // MC correction divide by scale
    object.px = p4.Px()*smear_scale/scale;
    object.py = p4.Py()*smear_scale/scale;
    object.pz = p4.Pz()*smear_scale/scale;
    object.energy = p4.Energy()*smear_scale/scale;
    object.pt = p4.Pt()*smear_scale/scale;
    
  } //End loop over objects

} //End apply_SMP_12_011_correction_ele

template <typename collectionType>
void HelperLeptonCore::apply_gsf_Et(collectionType& collection) {

  TLorentzVector p4;
  double gsf_scale = 1.0;
  
  for (auto& object: collection) {
    p4.SetPxPyPzE(object.px, object.py, object.pz, object.energy);

    gsf_scale = object.gsfEt/object.et;
      
    // MC correction divide by scale
    object.px = p4.Px()*gsf_scale;
    object.py = p4.Py()*gsf_scale;
    object.pz = p4.Pz()*gsf_scale;
    object.energy = p4.Energy()*gsf_scale;
    object.pt = p4.Pt()*gsf_scale;
    
  } //End loop over objects

} //End apply_gsf_Et

template <typename collectionType>
void HelperLeptonCore::apply_sc_Et(collectionType& collection) {

  TLorentzVector p4;
  double sc_scale = 1.0;
  
  for (auto& object: collection) {
    p4.SetPxPyPzE(object.px, object.py, object.pz, object.energy);

    sc_scale = object.scEt/object.et;
      
    // MC correction divide by scale
    object.px = p4.Px()*sc_scale;
    object.py = p4.Py()*sc_scale;
    object.pz = p4.Pz()*sc_scale;
    object.energy = p4.Energy()*sc_scale;
    object.pt = p4.Pt()*sc_scale;
    
  } //End loop over objects

} //End apply_sc_Et

template <typename collectionType>
void HelperLeptonCore::shiftLeptonEnergy(collectionType& collection, double shift_scale_B, double shift_scale_E1, double shift_scale_E2) {

  TLorentzVector p4;
  double shift_scale = 1.0;
  
  for (auto& object: collection) {
    p4.SetPxPyPzE(object.px, object.py, object.pz, object.energy);
    if (isData) continue;

    if (abs(p4.Eta()) < 1.5) shift_scale = shift_scale_B;
    else if (abs(p4.Eta()) < 2.0) shift_scale = shift_scale_E1;
    else shift_scale = shift_scale_E2;

    object.px = p4.Px()*shift_scale;
    object.py = p4.Py()*shift_scale;
    object.pz = p4.Pz()*shift_scale;
    object.energy = p4.Energy()*shift_scale;
    object.pt = p4.Pt()*shift_scale;

  } //End loop over objects

} //End shiftLeptonEnergy
    
template <typename collectionType>
void HelperLeptonCore::smearLeptonEnergy(collectionType& collection, double smear_size) {

  TLorentzVector p4;
  double smear_val = 0.0;
  double smear_scale = 1.0;

  for (auto& object: collection) {
    p4.SetPxPyPzE(object.px, object.py, object.pz, object.energy);
    if (isData) continue;

    smear_val = gSmearer->Gaus(0.0, smear_size);
    smear_scale = (p4.Energy() + smear_val)/p4.Energy();

    object.px = p4.Px()*smear_scale;
    object.py = p4.Py()*smear_scale;
    object.pz = p4.Pz()*smear_scale;
    object.energy = p4.Energy()*smear_scale;
    object.pt = p4.Pt()*smear_scale;

  } //End loop over objects

} //End smearLeptonEnergy  

template <typename particleType>
BNjet HelperLeptonCore::GetClosestJet(const BNjetCollection& jets, const particleType& object, const double maxDeltaR) {
  BNjet result;
  double minDeltaR = 999;
  for (auto& jet: jets) {
    double dR = deltaR(jet.eta, jet.phi, object.eta, object.phi);
    if ((dR <= maxDeltaR) && (dR < minDeltaR)) {
      result = jet;
      minDeltaR = dR;
    }
  }

  if (minDeltaR == 999) {
    result.phi = object.phi;
    result.eta = object.eta;
    result.pt = object.pt;
  }

  return result;
}

#endif // _HelperLeptonCore_h
