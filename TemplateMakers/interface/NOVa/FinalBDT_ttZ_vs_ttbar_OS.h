#ifndef _FinalBDT_ttZ_vs_ttbar_OS_h
#define _FinalBDT_ttZ_vs_ttbar_OS_h

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include "ttHMultileptonAnalysis/TemplateMakers/interface/KinematicVariable.h"
#include "ttHMultileptonAnalysis/TemplateMakers/interface/BranchInfo.h"
#include <typeinfo>

class FinalBDT_ttZ_vs_ttbar_OS: public KinematicVariable<double> {

public:
  //Store branch values so they are accessible to other classes
  vector<BranchInfo<double>> myVars;
  
  //Input variables
  Float_t varnumJets_40;
  Float_t varmht;
  Float_t varmass_leplep;
  Float_t varvectPt_leplep;
  Float_t vardR_leplep;
  Float_t varMT_of_jets_mass_of_jets;
  Float_t varMatch_ttbar_jj_Bq_bqq;
  Float_t varMatch_ttbar_jj_Bqq_bq;
  Float_t varMatch_ttbar_ll_Bb;
  Float_t varttbar_ll_B_CSV;
  Float_t varttbar_ll_b_CSV;
  Float_t varMatch_ttbar_jj_Bqq_bqq;
  
  vector<TMVA::Reader *> reader;

  BNjetCollection **jets;
  BNjetCollection **jets_40;
  TwoObjectKinematic<BNleptonCollection, BNjetCollection> * myMHT;
  TwoObjectKinematic<BNleptonCollection,BNleptonCollection> * myMassLepLep;
  TwoObjectKinematic<BNleptonCollection,BNleptonCollection> * myVectPtLepLep;
  TwoObjectKinematic<BNleptonCollection,BNleptonCollection> * myDeltaRLepLep;
  TwoObjectKinematic<BNjetCollection,BNjetCollection> * myMTOfJets;
  TwoObjectKinematic<BNjetCollection,BNjetCollection> * myMassOfJets;
  MatchTester_ttbar_jj * myMatchTester_ttbar_jj;
  MatchTester_ttbar_ll * myMatchTester_ttbar_ll;

  FinalBDT_ttZ_vs_ttbar_OS(BNjetCollection **_jets,
                           BNjetCollection **_jets_40,
                           TwoObjectKinematic<BNleptonCollection, BNjetCollection> * _myMHT,
                           TwoObjectKinematic<BNleptonCollection,BNleptonCollection> * _myMassLepLep,
                           TwoObjectKinematic<BNleptonCollection,BNleptonCollection> * _myVectPtLepLep,
                           TwoObjectKinematic<BNleptonCollection,BNleptonCollection> * _myDeltaRLepLep,
                           TwoObjectKinematic<BNjetCollection,BNjetCollection> * _myMTOfJets,
                           TwoObjectKinematic<BNjetCollection,BNjetCollection> * _myMassOfJets,
                           MatchTester_ttbar_jj * _myMatchTester_ttbar_jj,
                           MatchTester_ttbar_ll * _myMatchTester_ttbar_ll);

  void evaluate();

};

FinalBDT_ttZ_vs_ttbar_OS::FinalBDT_ttZ_vs_ttbar_OS(BNjetCollection **_jets,
                                                   BNjetCollection **_jets_40,
                                                   TwoObjectKinematic<BNleptonCollection, BNjetCollection> * _myMHT,
                                                   TwoObjectKinematic<BNleptonCollection,BNleptonCollection> * _myMassLepLep,
                                                   TwoObjectKinematic<BNleptonCollection,BNleptonCollection> * _myVectPtLepLep,
                                                   TwoObjectKinematic<BNleptonCollection,BNleptonCollection> * _myDeltaRLepLep,
                                                   TwoObjectKinematic<BNjetCollection,BNjetCollection> * _myMTOfJets,
                                                   TwoObjectKinematic<BNjetCollection,BNjetCollection> * _myMassOfJets,
                                                   MatchTester_ttbar_jj * _myMatchTester_ttbar_jj,
                                                   MatchTester_ttbar_ll * _myMatchTester_ttbar_ll):
  jets(_jets), jets_40(_jets_40), myMHT(_myMHT), myMassLepLep(_myMassLepLep), myVectPtLepLep(_myVectPtLepLep),
  myDeltaRLepLep(_myDeltaRLepLep), myMTOfJets(_myMTOfJets), myMassOfJets(_myMassOfJets),
  myMatchTester_ttbar_jj(_myMatchTester_ttbar_jj), myMatchTester_ttbar_ll(_myMatchTester_ttbar_ll) {

  //std::cout << "Setting up FinalBDT_ttZ_vs_ttbar_OS" << std::endl;
  
  branches["FinalBDT_ttZ_vs_ttbar_OS"] = BranchInfo<double>("FinalBDT_ttZ_vs_ttbar_OS");
  
  std::vector< TString >catList;
  catList.push_back("eq5j_ttbar"); //0
  catList.push_back("ge6j_ttbar"); //1

  for( unsigned int jj = 0 ; jj < 2 ; ++jj ) { //2 readers for 2 BDTs
    
    //std::cout << "Setting up reader " << jj << std::endl;

    reader.push_back( new TMVA::Reader( "!Color:!Silent" ));

    reader[jj]->AddVariable( "numJets_40", &varnumJets_40 );
    reader[jj]->AddVariable( "mht", &varmht );
    reader[jj]->AddVariable( "mass_leplep", &varmass_leplep );
    reader[jj]->AddVariable( "vectPt_leplep", &varvectPt_leplep );
    reader[jj]->AddVariable( "dR_leplep", &vardR_leplep );
    reader[jj]->AddVariable( "MT_of_jets_mass_of_jets", &varMT_of_jets_mass_of_jets );
    reader[jj]->AddVariable( "Match_ttbar_jj_Bq_bqq", &varMatch_ttbar_jj_Bq_bqq );
    reader[jj]->AddVariable( "Match_ttbar_jj_Bqq_bq", &varMatch_ttbar_jj_Bqq_bq );
    reader[jj]->AddVariable( "Match_ttbar_ll_Bb", &varMatch_ttbar_ll_Bb );
    reader[jj]->AddVariable( "ttbar_ll_B_CSV", &varttbar_ll_B_CSV );
    reader[jj]->AddVariable( "ttbar_ll_b_CSV", &varttbar_ll_b_CSV );
    if (jj == 1) {
      reader[jj]->AddVariable( "Match_ttbar_jj_Bqq_bqq", &varMatch_ttbar_jj_Bqq_bqq );
    }

    TString dir = (string(getenv("CMSSW_BASE"))+"/src/ttHMultileptonAnalysis/TemplateMakers/data/NOVa/Nov14/ttZ_vs_zjets_ttbar_OS/").c_str();
    TString label = catList[jj];
    TString file_name = "TMVAClassification_BDTG.weights.xml";
    //TString file_name = "TMVAClassification_CFMlpANN.weights.xml";
    TString weight_file_name = dir + label + "/" + file_name;

    reader[jj]->BookMVA( "BDTG method", weight_file_name );
    //reader[jj]->BookMVA( "CFMlpANN method", weight_file_name );

    std::cout << "Loading weight file " << weight_file_name << std::endl;
  }

}

void FinalBDT_ttZ_vs_ttbar_OS::evaluate() {
  if (this->evaluatedThisEvent) return;
  if ((*jets)->size() < 5) return;
  evaluatedThisEvent = true;

  //std::cout << "Inside FinalBDT_ttZ_vs_ttbar_OS::evaluate()" << std::endl;

  myMHT->evaluate();
  myMassLepLep->evaluate();
  myVectPtLepLep->evaluate();
  myDeltaRLepLep->evaluate();
  myMTOfJets->evaluate();
  myMassOfJets->evaluate();
  myMatchTester_ttbar_jj->evaluate();
  myMatchTester_ttbar_ll->evaluate();


  varnumJets_40 = (*jets_40)->size()*1.0;
  varmht = (*myMHT).myVars[0].branchVal;
  varmass_leplep = (*myMassLepLep).myVars[0].branchVal;
  varvectPt_leplep = (*myVectPtLepLep).myVars[0].branchVal;
  vardR_leplep = (*myDeltaRLepLep).myVars[0].branchVal;
  varMT_of_jets_mass_of_jets = (*myMTOfJets).myVars[0].branchVal/(*myMassOfJets).myVars[0].branchVal;

  std::string branchName = "";
  for (unsigned int ii = 0; ii < (*myMatchTester_ttbar_jj).myVars.size(); ii++) {
    branchName = (*myMatchTester_ttbar_jj).myVars[ii].branchName;
    //std::cout << branchName << std::endl;
    if (branchName == "Match_ttbar_jj_Bq_bqq") varMatch_ttbar_jj_Bq_bqq = (*myMatchTester_ttbar_jj).myVars[ii].branchVal;
    if (branchName == "Match_ttbar_jj_Bqq_bq") varMatch_ttbar_jj_Bqq_bq = (*myMatchTester_ttbar_jj).myVars[ii].branchVal;
    if (branchName == "Match_ttbar_jj_Bqq_bqq") varMatch_ttbar_jj_Bqq_bqq = (*myMatchTester_ttbar_jj).myVars[ii].branchVal;
  }
  for (unsigned int ii = 0; ii < (*myMatchTester_ttbar_ll).myVars.size(); ii++) {
    branchName = (*myMatchTester_ttbar_ll).myVars[ii].branchName;
    //std::cout << branchName << std::endl;
    if (branchName == "Match_ttbar_ll_Bb") varMatch_ttbar_ll_Bb = (*myMatchTester_ttbar_ll).myVars[ii].branchVal;
    if (branchName == "ttbar_ll_B_CSV") varttbar_ll_B_CSV = (*myMatchTester_ttbar_ll).myVars[ii].branchVal;
    if (branchName == "ttbar_ll_b_CSV") varttbar_ll_b_CSV = (*myMatchTester_ttbar_ll).myVars[ii].branchVal;
  }
  
  //std::cout << "Here" << std::endl;

  //   std::cout << "varnumJets_40: " << varnumJets_40 << std::endl;

  for( unsigned int jj = 0 ; jj < 2 ; ++jj ) {
    
    TMVA::Reader  *tmpReader = reader[jj];
    TString mvaName = "BDTG";
    //TString mvaName = "CFMlpANN";

    TString methodName = mvaName + TString(" method");
    Float_t annOut  = tmpReader->EvaluateMVA( methodName );

    if (jj == 0 && (*jets)->size() == 5) branches["FinalBDT_ttZ_vs_ttbar_OS"].branchVal = annOut;
    if (jj == 1 && (*jets)->size() >= 6) branches["FinalBDT_ttZ_vs_ttbar_OS"].branchVal = annOut;

  }

  //Clean out values from last event
  myVars.clear();

  for (typename map<TString, BranchInfo<double>>::iterator iBranch = branches.begin();
       iBranch != branches.end(); iBranch++) {
    myVars.push_back(iBranch->second);
  }
  
}
  

#endif 
