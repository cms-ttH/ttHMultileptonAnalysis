#ifndef _FinalBDT_ttZ_3l_h
#define _FinalBDT_ttZ_3l_h

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include "ttHMultileptonAnalysis/TemplateMakers/interface/KinematicVariable.h"
#include "ttHMultileptonAnalysis/TemplateMakers/interface/BranchInfo.h"
#include <typeinfo>

class FinalBDT_ttZ_3l: public KinematicVariable<double> {

public:
  //Store branch values so they are accessible to other classes
  vector<BranchInfo<double>> myVars;
  
  //Input variables for SS dilepton
  Float_t varMT_of_everything;
  Float_t varnumMediumBJets;
  Float_t varMatch_ttZ_3l_Bb;
  Float_t varMatch_ttZ_3l_Bq;
  Float_t varMatch_ttZ_3l_bq;
  Float_t varMatch_ttZ_3l_Bbq;
  Float_t varMatch_ttZ_3l_Bqq;
  Float_t varMatch_ttZ_3l_bqq;
  Float_t varMatch_ttZ_3l_Bbqq;
  Float_t varZLike_mass_leplep_SFOS_all;
  
  vector<TMVA::Reader *> reader;

  BNjetCollection **jets;
  ThreeObjectKinematic<BNmetCollection,BNleptonCollection,BNjetCollection> * myMTOfEverything;
  BNjetCollection **mediumCSVJets;
  MatchTester_ttZ_3l * myMatchTester_ttZ_3l;
  TwoObjectKinematic<BNleptonCollection,BNleptonCollection> * myZLikeMassLepLepSFOSAll;

  FinalBDT_ttZ_3l(BNjetCollection **_jets,
                  ThreeObjectKinematic<BNmetCollection,BNleptonCollection,BNjetCollection> * _myMTOfEverything,
                  BNjetCollection **_mediumCSVJets,
                  MatchTester_ttZ_3l * _myMatchTester_ttZ_3l,
                  TwoObjectKinematic<BNleptonCollection,BNleptonCollection> * _myZLikeMassLepLepSFOSAll);

  void evaluate();

};

FinalBDT_ttZ_3l::FinalBDT_ttZ_3l(BNjetCollection **_jets,
                                 ThreeObjectKinematic<BNmetCollection,BNleptonCollection,BNjetCollection> * _myMTOfEverything,
                                 BNjetCollection **_mediumCSVJets,
                                 MatchTester_ttZ_3l * _myMatchTester_ttZ_3l,
                                 TwoObjectKinematic<BNleptonCollection,BNleptonCollection> * _myZLikeMassLepLepSFOSAll):
  jets(_jets), myMTOfEverything(_myMTOfEverything), mediumCSVJets(_mediumCSVJets), myMatchTester_ttZ_3l(_myMatchTester_ttZ_3l),
  myZLikeMassLepLepSFOSAll(_myZLikeMassLepLepSFOSAll) {

  //std::cout << "Setting up FinalBDT_ttZ_3l" << std::endl;
  
  branches["FinalBDT_ttZ_3l"] = BranchInfo<double>("FinalBDT_ttZ_3l");
  
  std::vector< TString >catList;
  catList.push_back("eq3j"); //0
  catList.push_back("ge4j"); //1

  for( unsigned int jj = 0 ; jj < 2 ; ++jj ) { //2 readers for 2 BDTs
    
    //std::cout << "Setting up reader " << jj << std::endl;

    reader.push_back( new TMVA::Reader( "!Color:!Silent" ));

    reader[jj]->AddVariable( "MT_of_everything", &varMT_of_everything );
    reader[jj]->AddVariable( "numMediumBJets", &varnumMediumBJets );
    if (jj == 0) {
      reader[jj]->AddVariable( "Match_ttZ_3l_Bb", &varMatch_ttZ_3l_Bb );
      reader[jj]->AddVariable( "Match_ttZ_3l_Bq", &varMatch_ttZ_3l_Bq );
      reader[jj]->AddVariable( "Match_ttZ_3l_bq", &varMatch_ttZ_3l_bq ); }
    reader[jj]->AddVariable( "Match_ttZ_3l_Bbq", &varMatch_ttZ_3l_Bbq );
    reader[jj]->AddVariable( "Match_ttZ_3l_Bqq", &varMatch_ttZ_3l_Bqq );
    reader[jj]->AddVariable( "Match_ttZ_3l_bqq", &varMatch_ttZ_3l_bqq );
    if (jj == 1) {
      reader[jj]->AddVariable( "Match_ttZ_3l_Bbqq", &varMatch_ttZ_3l_Bbqq ); }
    reader[jj]->AddVariable( "ZLike_mass_leplep_SFOS_all", &varZLike_mass_leplep_SFOS_all );

    TString dir = (string(getenv("CMSSW_BASE"))+"/src/ttHMultileptonAnalysis/TemplateMakers/data/NOVa/Nov14/ttZ_vs_WZ_and_ttbar_3l/").c_str();
    TString label = catList[jj];
    TString file_name = "TMVAClassification_BDTG.weights.xml";
    //TString file_name = "TMVAClassification_CFMlpANN.weights.xml";
    TString weight_file_name = dir + label + "/" + file_name;

    reader[jj]->BookMVA( "BDTG method", weight_file_name );
    //reader[jj]->BookMVA( "CFMlpANN method", weight_file_name );

    std::cout << "Loading weight file " << weight_file_name << std::endl;
  }

}

void FinalBDT_ttZ_3l::evaluate() {
  if (this->evaluatedThisEvent) return;
  if ((*jets)->size() < 3) return;
  myZLikeMassLepLepSFOSAll->evaluate();
  if ((*myZLikeMassLepLepSFOSAll).myVars[0].branchVal < 0) return;
  evaluatedThisEvent = true;

  //std::cout << "Inside FinalBDT_ttZ_3l::evaluate()" << std::endl;

  varMT_of_everything = KinematicVariableConstants::FLOAT_INIT;
  varnumMediumBJets = KinematicVariableConstants::FLOAT_INIT;
  varMatch_ttZ_3l_Bb = KinematicVariableConstants::FLOAT_INIT;
  varMatch_ttZ_3l_Bq = KinematicVariableConstants::FLOAT_INIT;
  varMatch_ttZ_3l_bq = KinematicVariableConstants::FLOAT_INIT;
  varMatch_ttZ_3l_Bbq = KinematicVariableConstants::FLOAT_INIT;
  varMatch_ttZ_3l_Bqq = KinematicVariableConstants::FLOAT_INIT;
  varMatch_ttZ_3l_bqq = KinematicVariableConstants::FLOAT_INIT;
  varMatch_ttZ_3l_Bbqq = KinematicVariableConstants::FLOAT_INIT;
  varZLike_mass_leplep_SFOS_all = KinematicVariableConstants::FLOAT_INIT;
  
  myMTOfEverything->evaluate();
  myMatchTester_ttZ_3l->evaluate();

  varMT_of_everything = (*myMTOfEverything).myVars[0].branchVal;
  varnumMediumBJets = (*mediumCSVJets)->size()*1.0;

  std::string branchName = "";
  for (unsigned int ii = 0; ii < (*myMatchTester_ttZ_3l).myVars.size(); ii++) {
    branchName = (*myMatchTester_ttZ_3l).myVars[ii].branchName;
    //std::cout << branchName << std::endl;
    if (branchName == "Match_ttZ_3l_Bb") varMatch_ttZ_3l_Bb = (*myMatchTester_ttZ_3l).myVars[ii].branchVal;
    if (branchName == "Match_ttZ_3l_Bq") varMatch_ttZ_3l_Bq = (*myMatchTester_ttZ_3l).myVars[ii].branchVal;
    if (branchName == "Match_ttZ_3l_bq") varMatch_ttZ_3l_bq = (*myMatchTester_ttZ_3l).myVars[ii].branchVal;
    if (branchName == "Match_ttZ_3l_Bbq") varMatch_ttZ_3l_Bbq = (*myMatchTester_ttZ_3l).myVars[ii].branchVal;
    if (branchName == "Match_ttZ_3l_Bqq") varMatch_ttZ_3l_Bqq = (*myMatchTester_ttZ_3l).myVars[ii].branchVal;
    if (branchName == "Match_ttZ_3l_bqq") varMatch_ttZ_3l_bqq = (*myMatchTester_ttZ_3l).myVars[ii].branchVal;
    if (branchName == "Match_ttZ_3l_Bbqq") varMatch_ttZ_3l_Bbqq = (*myMatchTester_ttZ_3l).myVars[ii].branchVal;
  }
  //std::cout << "Here" << std::endl;
  varZLike_mass_leplep_SFOS_all = (*myZLikeMassLepLepSFOSAll).myVars[0].branchVal;

  

  if (varMT_of_everything < -99) {
    std::cout << "Error! varMT_of_everything = " << varMT_of_everything << std::endl; }
  if (varnumMediumBJets < -99) {
    std::cout << "Error! varnumMediumBJets = " << varnumMediumBJets << std::endl; }
  if (varMatch_ttZ_3l_Bb < -99) {
    std::cout << "Error! varMatch_ttZ_3l_Bb = " << varMatch_ttZ_3l_Bb << std::endl; }
  if (varMatch_ttZ_3l_Bq < -99) {
    std::cout << "Error! varMatch_ttZ_3l_Bq = " << varMatch_ttZ_3l_Bq << std::endl; }
  if (varMatch_ttZ_3l_bq < -99) {
    std::cout << "Error! varMatch_ttZ_3l_bq = " << varMatch_ttZ_3l_bq << std::endl; }
  if (varMatch_ttZ_3l_Bbq < -99) {
    std::cout << "Error! varMatch_ttZ_3l_Bbq = " << varMatch_ttZ_3l_Bbq << std::endl; }
  if (varMatch_ttZ_3l_Bqq < -99) {
    std::cout << "Error! varMatch_ttZ_3l_Bqq = " << varMatch_ttZ_3l_Bqq << std::endl; }
  if (varMatch_ttZ_3l_bqq < -99) {
    std::cout << "Error! varMatch_ttZ_3l_bqq = " << varMatch_ttZ_3l_bqq << std::endl; }
  if (varMatch_ttZ_3l_Bbqq < -99 && (*jets)->size() >= 4) {
    std::cout << "Error! varMatch_ttZ_3l_Bbqq = " << varMatch_ttZ_3l_Bbqq << std::endl; }
  if (varZLike_mass_leplep_SFOS_all < -99) {
    std::cout << "Error! varZLike_mass_leplep_SFOS_all = " << varZLike_mass_leplep_SFOS_all << std::endl; }
  
  for( unsigned int jj = 0 ; jj < 2 ; ++jj ) {
    
    TMVA::Reader  *tmpReader = reader[jj];
    TString mvaName = "BDTG";
    //TString mvaName = "CFMlpANN";

    TString methodName = mvaName + TString(" method");
    Float_t annOut  = tmpReader->EvaluateMVA( methodName );

    if (jj == 0 && (*jets)->size() == 3) branches["FinalBDT_ttZ_3l"].branchVal = annOut;
    if (jj == 1 && (*jets)->size() >= 4) branches["FinalBDT_ttZ_3l"].branchVal = annOut;

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
