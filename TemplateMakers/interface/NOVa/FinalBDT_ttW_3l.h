#ifndef _FinalBDT_ttW_3l_h
#define _FinalBDT_ttW_3l_h

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include "ttHMultileptonAnalysis/TemplateMakers/interface/KinematicVariable.h"
#include "ttHMultileptonAnalysis/TemplateMakers/interface/BranchInfo.h"
#include <typeinfo>

class FinalBDT_ttW_3l: public KinematicVariable<double> {

public:
  //Store branch values so they are accessible to other classes
  vector<BranchInfo<double>> myVars;

  Float_t varmet_pt;
  Float_t varMT_of_everything;
  Float_t varjets_by_pt_1_pt;
  Float_t varjets_by_CSV_2_btagCombinedSecVertex;
  Float_t varttbar_fake_3l_top_mass_lep_blep;
  Float_t varMatch_ttbar_fake_3l_B;
  Float_t varMatch_ttbar_fake_3l_b;
  Float_t varmax_max_Match_ttbar_fake_3l_B_Match_ttbar_fake_3l_b___6_;
  Float_t varMatch_ttW_3l_Bb;
  Float_t varall_SS_leptons_by_pt_1_pt;
  Float_t varall_SS_leptons_by_pt_2_pt;
  
  vector<TMVA::Reader *> reader;

  BNmetCollection **met;
  ThreeObjectKinematic<BNmetCollection,BNleptonCollection,BNjetCollection> * myMTOfEverything;
  BNjetCollection **jets;
  BNjetCollection **jetsByCSV;
  MatchTester_ttbar_fake_3l * myMatchTester_ttbar_fake_3l;
  MatchTester_ttW_3l * myMatchTester_ttW_3l;
  BNleptonCollection **leptonsSS;
  BNleptonCollection **tightLoosePreselectedLeptons;  

  FinalBDT_ttW_3l(BNmetCollection **_met,
                  ThreeObjectKinematic<BNmetCollection,BNleptonCollection,BNjetCollection> * _myMTOfEverything,
                  BNjetCollection **_jets,
                  BNjetCollection **_jetsByCSV,
                  MatchTester_ttbar_fake_3l * _myMatchTester_ttbar_fake_3l,
                  MatchTester_ttW_3l * _myMatchTester_ttW_3l,
                  BNleptonCollection **_leptonsSS,
                  BNleptonCollection **_tightLoosePreselectedLeptons);

  void evaluate();

};

FinalBDT_ttW_3l::FinalBDT_ttW_3l(BNmetCollection **_met,
                                 ThreeObjectKinematic<BNmetCollection,BNleptonCollection,BNjetCollection> * _myMTOfEverything,
                                 BNjetCollection **_jets,
                                 BNjetCollection **_jetsByCSV,
                                 MatchTester_ttbar_fake_3l * _myMatchTester_ttbar_fake_3l,
                                 MatchTester_ttW_3l * _myMatchTester_ttW_3l,
                                 BNleptonCollection **_leptonsSS,
                                 BNleptonCollection **_tightLoosePreselectedLeptons):
  met(_met), myMTOfEverything(_myMTOfEverything), jets(_jets), jetsByCSV(_jetsByCSV), myMatchTester_ttbar_fake_3l(_myMatchTester_ttbar_fake_3l),
  myMatchTester_ttW_3l(_myMatchTester_ttW_3l), leptonsSS(_leptonsSS), tightLoosePreselectedLeptons(_tightLoosePreselectedLeptons) {

  //std::cout << "Setting up FinalBDT_ttW_3l" << std::endl;
  
  branches["FinalBDT_ttW_3l"] = BranchInfo<double>("FinalBDT_ttW_3l");
  
  std::vector< TString >catList;
  catList.push_back("eq1j"); //0
  catList.push_back("ge2j"); //1

  for( unsigned int jj = 0 ; jj < 2 ; ++jj ) { //2 readers for 2 BDTs
    
    //std::cout << "Setting up reader " << jj << std::endl;

    reader.push_back( new TMVA::Reader( "!Color:!Silent" ));
    
    if (jj == 0) {
      reader[jj]->AddVariable( "met_pt", &varmet_pt ); 
      reader[jj]->AddVariable( "MT_of_everything", &varMT_of_everything );
      reader[jj]->AddVariable( "jets_by_pt_1_pt", &varjets_by_pt_1_pt );
      reader[jj]->AddVariable( "max_max_Match_ttbar_fake_3l_B_Match_ttbar_fake_3l_b___6_", &varmax_max_Match_ttbar_fake_3l_B_Match_ttbar_fake_3l_b___6_ );
    }
    if (jj == 1) {
      reader[jj]->AddVariable( "MT_of_everything", &varMT_of_everything );
      reader[jj]->AddVariable( "jets_by_CSV_2_btagCombinedSecVertex", &varjets_by_CSV_2_btagCombinedSecVertex );
      reader[jj]->AddVariable( "ttbar_fake_3l_top_mass_lep_blep", &varttbar_fake_3l_top_mass_lep_blep );
      reader[jj]->AddVariable( "Match_ttW_3l_Bb", &varMatch_ttW_3l_Bb );
    }
    reader[jj]->AddVariable( "all_SS_leptons_by_pt_1_pt", &varall_SS_leptons_by_pt_1_pt ); 
    reader[jj]->AddVariable( "all_SS_leptons_by_pt_2_pt", &varall_SS_leptons_by_pt_2_pt ); 
    
    TString dir = (string(getenv("CMSSW_BASE"))+"/src/ttHMultileptonAnalysis/TemplateMakers/data/NOVa/Nov14/ttW_vs_ttbar_3l/").c_str();
    TString label = catList[jj];
    TString file_name = "TMVAClassification_BDTG.weights.xml";
    //TString file_name = "TMVAClassification_CFMlpANN.weights.xml";
    TString weight_file_name = dir + label + "/" + file_name;
    
    reader[jj]->BookMVA( "BDTG method", weight_file_name );
    //reader[jj]->BookMVA( "CFMlpANN method", weight_file_name );
    
    std::cout << "Loading weight file " << weight_file_name << std::endl;
  }

}

void FinalBDT_ttW_3l::evaluate() {
  if (this->evaluatedThisEvent) return;
  if ((*jets)->size() < 1) return;
  if ((*leptonsSS)->size() < 2) return;
  if ((*tightLoosePreselectedLeptons)->size() != 3) return;
  if (abs((*tightLoosePreselectedLeptons)->at(0)->tkCharge + (*tightLoosePreselectedLeptons)->at(1)->tkCharge + (*tightLoosePreselectedLeptons)->at(2)->tkCharge) != 1) return;
  evaluatedThisEvent = true;

  //std::cout << "Inside FinalBDT_ttW_3l::evaluate()" << std::endl;

  varmet_pt = KinematicVariableConstants::FLOAT_INIT;
  varMT_of_everything = KinematicVariableConstants::FLOAT_INIT;
  varjets_by_pt_1_pt = KinematicVariableConstants::FLOAT_INIT;
  varjets_by_CSV_2_btagCombinedSecVertex = KinematicVariableConstants::FLOAT_INIT;
  varttbar_fake_3l_top_mass_lep_blep = KinematicVariableConstants::FLOAT_INIT;
  varMatch_ttbar_fake_3l_B = KinematicVariableConstants::FLOAT_INIT;
  varMatch_ttbar_fake_3l_b = KinematicVariableConstants::FLOAT_INIT;
  varmax_max_Match_ttbar_fake_3l_B_Match_ttbar_fake_3l_b___6_ = KinematicVariableConstants::FLOAT_INIT;
  varMatch_ttW_3l_Bb = KinematicVariableConstants::FLOAT_INIT;
  varall_SS_leptons_by_pt_1_pt = KinematicVariableConstants::FLOAT_INIT;
  varall_SS_leptons_by_pt_2_pt = KinematicVariableConstants::FLOAT_INIT;
  
  myMTOfEverything->evaluate();
  myMatchTester_ttbar_fake_3l->evaluate();
  myMatchTester_ttW_3l->evaluate();

  varmet_pt = (*met)->at(0).pt;
  varMT_of_everything = (*myMTOfEverything).myVars[0].branchVal;
  varjets_by_pt_1_pt = (*jets)->at(0).pt;
  if ((*jetsByCSV)->size() < 2) varjets_by_CSV_2_btagCombinedSecVertex = -1;
  else varjets_by_CSV_2_btagCombinedSecVertex = (*jetsByCSV)->at(1).btagCombinedSecVertex;
  varall_SS_leptons_by_pt_1_pt = (*leptonsSS)->at(0)->pt;
  varall_SS_leptons_by_pt_2_pt = (*leptonsSS)->at(1)->pt;

  std::string branchName = "";
  for (unsigned int ii = 0; ii < (*myMatchTester_ttbar_fake_3l).myVars.size(); ii++) {
    branchName = (*myMatchTester_ttbar_fake_3l).myVars[ii].branchName;
    if (branchName == "ttbar_fake_3l_top_mass_lep_blep") varttbar_fake_3l_top_mass_lep_blep = (*myMatchTester_ttbar_fake_3l).myVars[ii].branchVal;
    if (branchName == "Match_ttbar_fake_3l_B") varMatch_ttbar_fake_3l_B = max((*myMatchTester_ttbar_fake_3l).myVars[ii].branchVal,-6.0);
    if (branchName == "Match_ttbar_fake_3l_b") varMatch_ttbar_fake_3l_b = max((*myMatchTester_ttbar_fake_3l).myVars[ii].branchVal,-6.0);
  }
  varmax_max_Match_ttbar_fake_3l_B_Match_ttbar_fake_3l_b___6_ = max(varMatch_ttbar_fake_3l_B, varMatch_ttbar_fake_3l_b);
  for (unsigned int ii = 0; ii < (*myMatchTester_ttW_3l).myVars.size(); ii++) {
    branchName = (*myMatchTester_ttW_3l).myVars[ii].branchName;
    if (branchName == "Match_ttW_3l_Bb") varMatch_ttW_3l_Bb = (*myMatchTester_ttW_3l).myVars[ii].branchVal;
  }

  if (varmet_pt < -99) {
    std::cout << "Error! varmet_pt: " << varmet_pt << std::endl; }
  if (varMT_of_everything < -99) {
    std::cout << "Error! varMT_of_everything: " << varMT_of_everything << std::endl; }
  if (varjets_by_pt_1_pt < -99) {
    std::cout << "Error! varjets_by_pt_1_pt: " << varjets_by_pt_1_pt << std::endl; }
  if (varjets_by_CSV_2_btagCombinedSecVertex < -99 && (*jets)->size() >= 2) {
    std::cout << "Error! varjets_by_CSV_2_btagCombinedSecVertex: " << varjets_by_CSV_2_btagCombinedSecVertex << std::endl; }
  if (varttbar_fake_3l_top_mass_lep_blep < -99) {
    std::cout << "Error! varttbar_fake_3l_top_mass_lep_blep: " << varttbar_fake_3l_top_mass_lep_blep << std::endl; }
  if (varMatch_ttbar_fake_3l_B < -99) {
    std::cout << "Error! varMatch_ttbar_fake_3l_B: " << varMatch_ttbar_fake_3l_B << std::endl; }
  if (varMatch_ttbar_fake_3l_b < -99) {
    std::cout << "Error! varMatch_ttbar_fake_3l_b: " << varMatch_ttbar_fake_3l_b << std::endl; }
  if (varmax_max_Match_ttbar_fake_3l_B_Match_ttbar_fake_3l_b___6_ < -99) {
    std::cout << "Error! varmax_max_Match_ttbar_fake_3l_B_Match_ttbar_fake_3l_b___6_: " << varmax_max_Match_ttbar_fake_3l_B_Match_ttbar_fake_3l_b___6_ << std::endl; }
  if (varMatch_ttW_3l_Bb < -99 && (*jets)->size() >= 2) {
    std::cout << "Error! varMatch_ttW_3l_Bb: " << varMatch_ttW_3l_Bb << std::endl; }
  if (varall_SS_leptons_by_pt_1_pt < -99) {
    std::cout << "Error! varall_SS_leptons_by_pt_1_pt: " << varall_SS_leptons_by_pt_1_pt << std::endl; }
  if (varall_SS_leptons_by_pt_2_pt < -99) {
    std::cout << "Error! varall_SS_leptons_by_pt_2_pt: " << varall_SS_leptons_by_pt_2_pt << std::endl; }
  
  for( unsigned int jj = 0 ; jj < 2 ; ++jj ) {
    
    TMVA::Reader  *tmpReader = reader[jj];
    TString mvaName = "BDTG";
    //TString mvaName = "CFMlpANN";

    TString methodName = mvaName + TString(" method");
    Float_t annOut  = tmpReader->EvaluateMVA( methodName );

    if (jj == 0 && (*jets)->size() == 1) branches["FinalBDT_ttW_3l"].branchVal = annOut;
    if (jj == 1 && (*jets)->size() >= 2) branches["FinalBDT_ttW_3l"].branchVal = annOut;

//     std::cout << "annOut[jj]: " << annOut << "[" << jj << "]" << std::endl;
  }

  //Clean out values from last event
  myVars.clear();

  for (typename map<TString, BranchInfo<double>>::iterator iBranch = branches.begin();
       iBranch != branches.end(); iBranch++) {
    myVars.push_back(iBranch->second);
  }
  
}
  

#endif 
