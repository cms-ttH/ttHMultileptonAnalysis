#ifndef _FinalBDT_ttZ_vs_zjets_ttbar_OS_h
#define _FinalBDT_ttZ_vs_zjets_ttbar_OS_h

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include "ttHMultileptonAnalysis/TemplateMakers/interface/KinematicVariable.h"
#include "ttHMultileptonAnalysis/TemplateMakers/interface/BranchInfo.h"
#include <typeinfo>

class FinalBDT_ttZ_vs_zjets_ttbar_OS: public KinematicVariable<double> {

public:
  //Store branch values so they are accessible to other classes
  vector<BranchInfo<double>> myVars;
  
  //Input variables
  Float_t varnumJets_40;
  Float_t varMT_of_leps_jets_mass_of_leps_jets;
  Float_t varjets_by_pt_5_pt;
  Float_t varjets_by_CSV_1_btagCombinedSecVertex;
  Float_t varjets_by_CSV_2_btagCombinedSecVertex;
  Float_t varMatch_ttbar_jj_Bq_bqq;
  Float_t varMatch_ttbar_jj_Bqq_bq;
  Float_t varMatch_ttbar_jj_Bqq_bqq;
  Float_t varlog_chiSquared_;
  Float_t varFinalBDT_ttZ_vs_ttbar_OS;
  
  vector<TMVA::Reader *> reader;

  BNjetCollection **jets;
  BNjetCollection **jets_40;
  TwoObjectKinematic<BNleptonCollection,BNjetCollection> * myMTOfLepsJets;
  TwoObjectKinematic<BNleptonCollection,BNjetCollection> * myMassOfLepsJets;
  BNjetCollection **jetsByCSV;
  MatchTester_ttbar_jj * myMatchTester_ttbar_jj;
  KinFitterttHadHad * chiSquareTopMasses;
  FinalBDT_ttZ_vs_ttbar_OS * myFinalBDT_ttZ_vs_ttbar_OS;

  FinalBDT_ttZ_vs_zjets_ttbar_OS(BNjetCollection **_jets,
                           BNjetCollection **_jets_40,
                           TwoObjectKinematic<BNleptonCollection,BNjetCollection> * _myMTOfLepsJets,
                           TwoObjectKinematic<BNleptonCollection,BNjetCollection> * _myMassOfLepsJets,
                           BNjetCollection **_jetsByCSV,
                           MatchTester_ttbar_jj * _myMatchTester_ttbar_jj,
                           KinFitterttHadHad * _chiSquareTopMasses,
                           FinalBDT_ttZ_vs_ttbar_OS * _myFinalBDT_ttZ_vs_ttbar_OS);

  void evaluate();

};

FinalBDT_ttZ_vs_zjets_ttbar_OS::FinalBDT_ttZ_vs_zjets_ttbar_OS(BNjetCollection **_jets,
                                                   BNjetCollection **_jets_40,
                                                   TwoObjectKinematic<BNleptonCollection,BNjetCollection> * _myMTOfLepsJets,
                                                   TwoObjectKinematic<BNleptonCollection,BNjetCollection> * _myMassOfLepsJets,
                                                   BNjetCollection **_jetsByCSV,
                                                   MatchTester_ttbar_jj * _myMatchTester_ttbar_jj,
                                                   KinFitterttHadHad * _chiSquareTopMasses,
                                                   FinalBDT_ttZ_vs_ttbar_OS * _myFinalBDT_ttZ_vs_ttbar_OS):

  jets(_jets), jets_40(_jets_40), myMTOfLepsJets(_myMTOfLepsJets), myMassOfLepsJets(_myMassOfLepsJets),
  jetsByCSV(_jetsByCSV), myMatchTester_ttbar_jj(_myMatchTester_ttbar_jj), chiSquareTopMasses(_chiSquareTopMasses),
  myFinalBDT_ttZ_vs_ttbar_OS(_myFinalBDT_ttZ_vs_ttbar_OS) {

  //std::cout << "Setting up FinalBDT_ttZ_vs_zjets_ttbar_OS" << std::endl;
  
  branches["FinalBDT_ttZ_vs_zjets_ttbar_OS"] = BranchInfo<double>("FinalBDT_ttZ_vs_zjets_ttbar_OS");
  
  std::vector< TString >catList;
  catList.push_back("eq5j_zjets_ttbar"); //0
  catList.push_back("ge6j_zjets_ttbar"); //1

  for( unsigned int jj = 0 ; jj < 2 ; ++jj ) { //2 readers for 2 BDTs
    
    //std::cout << "Setting up reader " << jj << std::endl;

    reader.push_back( new TMVA::Reader( "!Color:!Silent" ));

    reader[jj]->AddVariable( "numJets_40", &varnumJets_40 );
    reader[jj]->AddVariable( "MT_of_leps_jets_mass_of_leps_jets", &varMT_of_leps_jets_mass_of_leps_jets );
    reader[jj]->AddVariable( "jets_by_pt_5_pt", &varjets_by_pt_5_pt );
    reader[jj]->AddVariable( "jets_by_CSV_1_btagCombinedSecVertex", &varjets_by_CSV_1_btagCombinedSecVertex );
    reader[jj]->AddVariable( "jets_by_CSV_2_btagCombinedSecVertex", &varjets_by_CSV_2_btagCombinedSecVertex );
    reader[jj]->AddVariable( "Match_ttbar_jj_Bq_bqq", &varMatch_ttbar_jj_Bq_bqq );
    reader[jj]->AddVariable( "Match_ttbar_jj_Bqq_bq", &varMatch_ttbar_jj_Bqq_bq );
    if (jj == 1) {
      reader[jj]->AddVariable( "Match_ttbar_jj_Bqq_bqq", &varMatch_ttbar_jj_Bqq_bqq );
      reader[jj]->AddVariable( "log_chiSquared_", &varlog_chiSquared_ );
    }
    reader[jj]->AddVariable( "FinalBDT_ttZ_vs_ttbar_OS", &varFinalBDT_ttZ_vs_ttbar_OS );

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

void FinalBDT_ttZ_vs_zjets_ttbar_OS::evaluate() {
  if (this->evaluatedThisEvent) return;
  if ((*jets)->size() < 5) return;
  evaluatedThisEvent = true;

  //std::cout << "Inside FinalBDT_ttZ_vs_zjets_ttbar_OS::evaluate()" << std::endl;

  varnumJets_40 = KinematicVariableConstants::FLOAT_INIT;
  varMT_of_leps_jets_mass_of_leps_jets = KinematicVariableConstants::FLOAT_INIT;
  varjets_by_pt_5_pt = KinematicVariableConstants::FLOAT_INIT;
  varjets_by_CSV_1_btagCombinedSecVertex = KinematicVariableConstants::FLOAT_INIT;
  varjets_by_CSV_2_btagCombinedSecVertex = KinematicVariableConstants::FLOAT_INIT;
  varMatch_ttbar_jj_Bq_bqq = KinematicVariableConstants::FLOAT_INIT;
  varMatch_ttbar_jj_Bqq_bq = KinematicVariableConstants::FLOAT_INIT;
  varMatch_ttbar_jj_Bqq_bqq = KinematicVariableConstants::FLOAT_INIT;
  varlog_chiSquared_ = KinematicVariableConstants::FLOAT_INIT;
  varFinalBDT_ttZ_vs_ttbar_OS = KinematicVariableConstants::FLOAT_INIT;  

  myMTOfLepsJets->evaluate();
  myMassOfLepsJets->evaluate();
  myMatchTester_ttbar_jj->evaluate();
  chiSquareTopMasses->evaluate();
  myFinalBDT_ttZ_vs_ttbar_OS->evaluate();

  
  varnumJets_40 = (*jets_40)->size()*1.0;
  varMT_of_leps_jets_mass_of_leps_jets = (*myMTOfLepsJets).myVars[0].branchVal/(*myMassOfLepsJets).myVars[0].branchVal;
  varjets_by_pt_5_pt = (*jets)->at(4).pt;
  varjets_by_CSV_1_btagCombinedSecVertex = (*jetsByCSV)->at(0).btagCombinedSecVertex;
  varjets_by_CSV_2_btagCombinedSecVertex = (*jetsByCSV)->at(1).btagCombinedSecVertex;

  std::string branchName = "";
  for (unsigned int ii = 0; ii < (*myMatchTester_ttbar_jj).myVars.size(); ii++) {
    branchName = (*myMatchTester_ttbar_jj).myVars[ii].branchName;
    //std::cout << branchName << std::endl;
    if (branchName == "Match_ttbar_jj_Bq_bqq") varMatch_ttbar_jj_Bq_bqq = (*myMatchTester_ttbar_jj).myVars[ii].branchVal;
    if (branchName == "Match_ttbar_jj_Bqq_bq") varMatch_ttbar_jj_Bqq_bq = (*myMatchTester_ttbar_jj).myVars[ii].branchVal;
    if (branchName == "Match_ttbar_jj_Bqq_bqq") varMatch_ttbar_jj_Bqq_bqq = (*myMatchTester_ttbar_jj).myVars[ii].branchVal;
  }
  for (unsigned int ii = 0; ii < (*chiSquareTopMasses).myVars.size(); ii++) {
    branchName = (*chiSquareTopMasses).myVars[ii].branchName;
    //std::cout << branchName << std::endl;
    if (branchName == "chiSquared") varlog_chiSquared_ = log((*chiSquareTopMasses).myVars[ii].branchVal);
  }
  
  varFinalBDT_ttZ_vs_ttbar_OS = (*myFinalBDT_ttZ_vs_ttbar_OS).myVars[0].branchVal;
  
  //std::cout << "Here" << std::endl;

  //   std::cout << "varnumJets_40: " << varnumJets_40 << std::endl;

  if (varnumJets_40 < -99) {
  	std::cout << "Error! varnumJets_40 = " << varnumJets_40 << std::endl; }
  if (varMT_of_leps_jets_mass_of_leps_jets < -99) {
    std::cout << "Error! varMT_of_leps_jets_mass_of_leps_jets = " << varMT_of_leps_jets_mass_of_leps_jets << std::endl; }
  if (varjets_by_pt_5_pt < -99) {
  	std::cout << "Error! varjets_by_pt_5_pt = " << varjets_by_pt_5_pt << std::endl; }
  if (varjets_by_CSV_1_btagCombinedSecVertex < -99) {
  	std::cout << "Error! varjets_by_CSV_1_btagCombinedSecVertex = " << varjets_by_CSV_1_btagCombinedSecVertex << std::endl; }
  if (varjets_by_CSV_2_btagCombinedSecVertex < -99) {
  	std::cout << "Error! varjets_by_CSV_2_btagCombinedSecVertex = " << varjets_by_CSV_2_btagCombinedSecVertex << std::endl; }
  if (varMatch_ttbar_jj_Bq_bqq < -99) {
  	std::cout << "Error! varMatch_ttbar_jj_Bq_bqq = " << varMatch_ttbar_jj_Bq_bqq << std::endl; }
  if (varMatch_ttbar_jj_Bqq_bq < -99) {
  	std::cout << "Error! varMatch_ttbar_jj_Bqq_bq = " << varMatch_ttbar_jj_Bqq_bq << std::endl; }
  if (varMatch_ttbar_jj_Bqq_bqq < -99 && (*jets)->size() >= 6) {
  	std::cout << "Error! varMatch_ttbar_jj_Bqq_bqq = " << varMatch_ttbar_jj_Bqq_bqq << std::endl; }
  if (varlog_chiSquared_ < -99 && (*jets)->size() >= 6) {
  	std::cout << "Error! varlog_chiSquared_ = " << varlog_chiSquared_ << std::endl; }
  if (varFinalBDT_ttZ_vs_ttbar_OS < -99) {  
  	std::cout << "Error! varFinalBDT_ttZ_vs_ttbar_OS = " << varFinalBDT_ttZ_vs_ttbar_OS << std::endl; }  

  for( unsigned int jj = 0 ; jj < 2 ; ++jj ) {
    
    TMVA::Reader  *tmpReader = reader[jj];
    TString mvaName = "BDTG";
    //TString mvaName = "CFMlpANN";

    TString methodName = mvaName + TString(" method");
    Float_t annOut  = tmpReader->EvaluateMVA( methodName );

    if (jj == 0 && (*jets)->size() == 5) branches["FinalBDT_ttZ_vs_zjets_ttbar_OS"].branchVal = annOut;
    if (jj == 1 && (*jets)->size() >= 6) branches["FinalBDT_ttZ_vs_zjets_ttbar_OS"].branchVal = annOut;

  }

  //Clean out values from last event
  myVars.clear();

  for (typename map<TString, BranchInfo<double>>::iterator iBranch = branches.begin();
       iBranch != branches.end(); iBranch++) {
    myVars.push_back(iBranch->second);
  }
  
}
  

#endif 
