#ifndef _KinFitterttHadHad_h
#define _KinFitterttHadHad_h
#include "ttHMultileptonAnalysis/TemplateMakers/interface/KinematicVariable.h"
#include "ttHMultileptonAnalysis/TemplateMakers/interface/BranchInfo.h"
#include<vector>

class KinFitterttHadHad: public KinematicVariable<double> {

public:

  //Store branch values so they are accessible to other classes
  vector<BranchInfo<double>> myVars;
  
  BNjetCollection ** jets;
  double chiSquared;
  double W_massTTHadHad1;
  double W_massTTHadHad2;
  double top_massTTHadHad1;
  double top_massTTHadHad2;
  TLorentzVector B1;
  TLorentzVector B2;
  TLorentzVector Q1;
  TLorentzVector Q2;
  TLorentzVector Q3;
  TLorentzVector Q4;
  int B1_mother;
  int B2_mother;
  int Q1_mother;
  int Q2_mother;
  int Q3_mother;
  int Q4_mother;

  KinFitterttHadHad(BNjetCollection **_jets);
  TLorentzVector getFourVector(const int i);
  int getMother(const int i);
  void calculateChiSquared(const int firstBIndex, const int secondBIndex, const int *v, const int size);
  void permute(int *v, const int start, const int n);
  void evaluate();
};

KinFitterttHadHad::KinFitterttHadHad(BNjetCollection **_jets) : jets(_jets) {
  this->resetVal = -KinematicVariableConstants::DOUBLE_INIT;

  branches["chiSquared"] = BranchInfo<double>("chiSquared");
  branches["W_massTTHadHad1"] = BranchInfo<double>("W_massTTHadHad1");
  branches["W_massTTHadHad2"] = BranchInfo<double>("W_massTTHadHad2");
  branches["top_massTTHadHad1"] = BranchInfo<double>("top_massTTHadHad1");
  branches["top_massTTHadHad2"] = BranchInfo<double>("top_massTTHadHad2");
  branches["matched_W_TTHadHad1"] = BranchInfo<double>("matched_W_TTHadHad1");
  branches["matched_W_TTHadHad2"] = BranchInfo<double>("matched_W_TTHadHad2");
  branches["matched_top_TTHadHad1"] = BranchInfo<double>("matched_top_TTHadHad1");
  branches["matched_top_TTHadHad2"] = BranchInfo<double>("matched_top_TTHadHad2");
}

void KinFitterttHadHad::evaluate () {
  if (this->evaluatedThisEvent) return;
  evaluatedThisEvent = true;

  if((*jets)->size() <= 10) {
    int v[(*jets)->size() - 2];
    for(unsigned int i = 0; i < ((*jets)->size() - 2); i++) {
      v[i] = i+2;
    }
    if((*jets)->size() >= 6){
      permute(v, 0, sizeof(v)/sizeof(int));
    }
  }

  else {
    int v[8];
    for(unsigned int i = 0; i < 8; i++) {
      v[i] = i+2;
    }
    permute(v, 0, sizeof(v)/sizeof(int));
  }


}

void KinFitterttHadHad::permute(int *v, const int start, const int n) {
  if (start == n-1) {
    calculateChiSquared(0, 1, v, n);
    calculateChiSquared(1, 0, v, n);
  }
  else {
    for (int i = start; i < n; i++) {
      int tmp = v[i];
      v[i] = v[start];
      v[start] = tmp;
      permute(v, start+1, n);
      v[start] = v[i];
      v[i] = tmp;

    }
  }
}

void KinFitterttHadHad::calculateChiSquared(const int firstBIndex, const int secondBIndex, const int *v, const int size) {
  if (v != 0) {
    B1 = getFourVector(firstBIndex);
    B2 = getFourVector(secondBIndex);
    Q1 = getFourVector(v[0]);
    Q2 = getFourVector(v[1]);
    Q3 = getFourVector(v[2]);
    Q4 = getFourVector(v[3]);

    B1_mother = getMother(firstBIndex);
    B2_mother = getMother(secondBIndex);
    Q1_mother = getMother(v[0]);
    Q2_mother = getMother(v[1]);
    Q3_mother = getMother(v[2]);
    Q4_mother = getMother(v[3]);

//     for  (int i = 0; i < size; i++) {
//       printf("%4d", v[i]+1 );
//     }
//     printf("\n");

    
    double stdev_W = 13.8;
    double stdev_top = 22.7;
    chiSquared = pow((((Q1+Q2).M() - 88)/stdev_W), 2) + pow((((Q3+Q4).M() - 88)/stdev_W),2) + pow((((B1+Q1+Q2).M() - 180)/stdev_top),2) +  pow((((B2+Q3+Q4).M() - 180)/stdev_top),2);
    // chiSquared = pow((((Q1+Q2).M() - 81)/stdev_W), 2) + pow((((Q3+Q4).M() - 81)/stdev_W),2) + pow((((B1+Q1+Q2).M() - 173)/stdev_top),2) +  pow((((B2+Q3+Q4).M() - 173)/stdev_top),2);
    //    printf("chiSquared = %f\n",chiSquared);
    W_massTTHadHad1 = (Q1+Q2).M();
    W_massTTHadHad2 = (Q3+Q4).M();
    top_massTTHadHad1 = (B1+Q1+Q2).M();
    top_massTTHadHad2 = (B2+Q3+Q4).M();
    //    printf("mass 1 = %f\n", W_massTTHadHad1);
    // printf("mass 2 = %f\n", W_massTTHadHad2);
    //    printf("mass 1 = %f\n", top_massTTHadHad1);
    //    printf("mass 2 = %f\n", top_massTTHadHad2);
    
    int W1_match =  (abs(Q1_mother) == 24 && abs(Q2_mother) == 24 && Q1_mother == Q2_mother);
    int W2_match = (abs(Q3_mother) == 24 && abs(Q4_mother) == 24 && Q3_mother == Q4_mother);
    int top1_match = (W1_match && abs(Q1_mother + B1_mother) == 30);
    int top2_match = (W2_match && abs(Q3_mother + B2_mother) == 30);
    
    if (chiSquared < branches["chiSquared"].branchVal) {
      branches["chiSquared"].branchVal = chiSquared;
      branches["W_massTTHadHad1"].branchVal = W_massTTHadHad1;
      branches["W_massTTHadHad2"].branchVal = W_massTTHadHad2;
      branches["top_massTTHadHad1"].branchVal = top_massTTHadHad1;
      branches["top_massTTHadHad2"].branchVal = top_massTTHadHad2;
      branches["matched_W_TTHadHad1"].branchVal = W1_match;
      branches["matched_W_TTHadHad2"].branchVal = W2_match;
      branches["matched_top_TTHadHad1"].branchVal = top1_match;
      branches["matched_top_TTHadHad2"].branchVal = top2_match;

    }
  }

  myVars.clear();
  
  for (typename map<TString, BranchInfo<double>>::iterator iBranch = branches.begin();
       iBranch != branches.end(); iBranch++) {
    myVars.push_back(iBranch->second);
  }
  
}

TLorentzVector KinFitterttHadHad::getFourVector(const int i) {
  BNjet jet = (*jets)->at(i);
  TLorentzVector result(jet.px, jet.py, jet.pz, jet.energy);

  return result;
}

int KinFitterttHadHad::getMother(const int i) {
  BNjet jet = (*jets)->at(i);
  int result = jet.genPartonMotherId;

  return result;
}

#endif
