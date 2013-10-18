#ifndef _TwoObjectKinematic_h
#define _TwoObjectKinematic_h

#include  "ttHMultileptonAnalysis/TemplateMakers/interface/KinematicVariable.h"
#include "ttHMultileptonAnalysis/TemplateMakers/interface/BranchInfo.h"
#include <typeinfo>

// //Anna's fix to make BNleptonCollection work just like any other collection
// template<typename ObjectType_temp>
// ObjectType_temp * ptr(ObjectType_temp & obj) { return &obj; } //turn reference into pointer!

// template<typename ObjectType_temp>
// ObjectType_temp * ptr(ObjectType_temp * obj) { return obj; } //obj is already pointer, return it!

template <class collectionType1, class collectionType2>
class TwoObjectKinematic: public KinematicVariable<double> {

public:
  //An "abs" prefix on a variable indicates the absolute value of a
  //"pair" variable (e.g. |mass|, |deltaR|), or the absolute value
  //of each element of a non-pair variable (e.g. |eta1| + |eta2|)
  //A "vect" prefix indicates the sum of vectors (vect1 + vect2)
  string variable; //mass, MT, deltaR, deltaEta, deltaPhi, absDeltaEta, absDeltaPhi,
                   //pt, pz, energy, eta, phi, absPz, absEta, absPhi,
                   //vectPt, vectPz, vectEnergy, vectEta, vectPhi
                   //absVectPt, absVectPz, absVectEnergy, absVectEta, absVectPhi
  //An "abs" suffix indicates the absolute value of each variable as a whole,
  //e.g. sum ( |a1+a2| + |b1+b2| + |c1+c2| ...)
  //An "abs" prefix indicates the absolute value of the final output,
  //e.g. | sum ( (a1+a2) + (b1+b2) + (c1+c2) ...)
  string which_pair; //all_pairs, max, min, avg, sum, closest_to,
                     //all_pairs_abs, maxAbs, minAbs, avgAbs, sumAbs,
                     //absMax, absMin, absAvg, absSum, vector_sum
  string new_name; //Replacement name if desired
  double target_value; //only matters for closest_to; otherwise, set to -99
  collectionType1 **selCollection1; //first selected collection
  collectionType2 **selCollection2; //second selected collection
  string branch_name_1; //first selected collection name
  string branch_name_2; //second collected collection name
  unsigned int max1; //number of first collection objects 
  unsigned int max2; //number of second collection objects
  
  TwoObjectKinematic(string input_variable, string input_which_pair, double input_target_value, string input_new_name,
                collectionType1 **input_selCollection1, string input_branch_name_1, int input_max1,
                collectionType2 **input_selCollection2, string input_branch_name_2, int input_max2);
  
  void evaluate();
  
};

//Template for any collections other than BNlepton and BNmet
template <class collectionType1, class collectionType2>
TwoObjectKinematic<collectionType1,collectionType2>::TwoObjectKinematic(string input_variable, string input_which_pair, double input_target_value, string input_new_name,
                                                                        collectionType1 **input_selCollection1, string input_branch_name_1, int input_max1,
                                                                        collectionType2 **input_selCollection2, string input_branch_name_2, int input_max2):
  variable(input_variable), which_pair(input_which_pair), target_value(input_target_value), new_name(input_new_name),
  selCollection1(input_selCollection1), branch_name_1(input_branch_name_1),max1(input_max1),
  selCollection2(input_selCollection2), branch_name_2(input_branch_name_2),max2(input_max2)
{

  //Convert target value to integers for branch name
  int target_value_whole = floor(target_value);
  int target_value_fraction = static_cast<int>(target_value*10)%10;

  // --Create the branches--
  //Creates branches for each pair of objects
  if (which_pair == "all_pairs" || which_pair == "all_pairs_abs") { 
    for (unsigned int i=0; i<max1; i++) {
      for (unsigned int j=0; j<max2; j++) {
        //Eliminates pairs of the same object (e.g. jet1_jet1) and redundant pairs (jet1_jet2 and jet2_jet1) 
        if ( typeid(*selCollection1).name() != typeid(*selCollection2).name() || i < j ) {
          TString bName = Form("%s_%d_%s_%d_%s", branch_name_1.c_str(), i+1, branch_name_2.c_str(), j+1, variable.c_str());
          if (which_pair == "all_pairs_abs") bName = Form("%s_%d_%s_%d_abs_%s", branch_name_1.c_str(), i+1, branch_name_2.c_str(), j+1, variable.c_str());
          if (new_name.length() > 0) bName = Form("%s_%d_%d", new_name.c_str(), i+1, j+2);
          branches[bName] = BranchInfo<double>(bName);
        }
      }
    }
  }
  //If you gave it a new name, use that name
  else if (new_name.length() > 0) {
    TString bName = Form("%s", new_name.c_str());
    branches[bName] = BranchInfo<double>(bName);
  }    
  //Creates a branch for the pair of objects
  else if (which_pair == "closest_to") { 
    TString bName = Form("%s_%s_%s_%dp%d_%s", branch_name_1.c_str(), branch_name_2.c_str(), which_pair.c_str(), target_value_whole, target_value_fraction, variable.c_str()); 
    branches[bName] = BranchInfo<double>(bName);
  }
  //Creates a branch for the pair of objects
  else {
    TString bName = Form("%s_%s_%s_%s", branch_name_1.c_str(), branch_name_2.c_str(), which_pair.c_str(), variable.c_str()); 
    branches[bName] = BranchInfo<double>(bName);
  }
  //Initializes branch values to default
  this->resetVal = KinematicVariableConstants::DOUBLE_INIT;

}

//Template for any collections other than BNlepton and BNmet
template <class collectionType1, class collectionType2>
void TwoObjectKinematic<collectionType1,collectionType2>::evaluate() {
  if (this->evaluatedThisEvent) return;
  evaluatedThisEvent = true;

  double thisValue = KinematicVariableConstants::DOUBLE_INIT; //Values of each individual object
  double thisValueSum = 0.0; //Sum over each object value
  double thisValueSumAbs = 0.0; //Sum over absolute value of each object value
  double thisValueIterator = 0.0; //Counting the total number of individual objects

  double thisPairValue = KinematicVariableConstants::DOUBLE_INIT; //Values of each pair of objects
  double thisPairValueSum = 0.0; //Sum over value of each object pair
  double thisPairValueSumAbs = 0.0; //Sum over absolute value of each object pair
  double thisPairValueIterator = 0.0; //Conting the total number of distinct pairs

  TLorentzVector thisVectorSum(0.,0.,0.,0.); //Sum of vectors of individual objects
  TLorentzVector thisVectorTransverse(0.,0.,0.,0.); //Transverse vectors of individual objects
  TLorentzVector thisVectorTransverseSum(0.,0.,0.,0.); //Sum of transverse vectors of individual objects
  
  TLorentzVector vect1; //Vector for first object
  TLorentzVector vect2; //Vector for second object
  TLorentzVector vect12; //Combined vector of two objects

  TLorentzVector vect1_transverse; //Transverse vector for first object
  TLorentzVector vect2_transverse; //Transverse vector for second object
  TLorentzVector vect12_transverse; //Combined transverse vector of two
  
  //Convert target value to integers for branch name
  int target_value_whole = floor(target_value);
  int target_value_fraction = static_cast<int>(target_value*10)%10;

  bool useValue = 0; //Use individual values for sum and average calculations
  if ( variable == "pt" || variable == "pz" || variable == "energy" || variable == "eta" || variable == "phi" ||
       variable == "absPz" || variable == "absEta" || variable == "absPhi") {
    useValue = 1;
  }
  bool usePairValue = 0; //Use pair values for sum and average calculations
  if ( variable == "mass" || variable == "MT" || variable == "deltaR" || variable == "deltaEta" || variable == "deltaPhi" ||
       variable == "absDeltaEta" || variable == "absDeltaPhi" || variable == "vectPt" || variable == "vectPz" ||
       variable == "vectEnergy" || variable == "vectEta" || variable == "vectPhi" || variable == "absVectPhi" ||
       variable == "absVectPz" || variable == "absVectEnergy" || variable == "absVectEta" || variable == "absVectPhi" ) {
    usePairValue = 1;
  }


  //Use individual values
  if ( useValue || which_pair == "vector_sum" ) {
    //Loop over first set of objects alone
    //for (typename collectionType1::const_iterator object1 = (*selCollection1)->begin(); object1 != (*selCollection1)->end(); ++object1) {
    for (unsigned int iObj1 = 0; iObj1 < (*selCollection1)->size(); iObj1++) {

      //Sets max number of objects
      if ( iObj1 < max1 ) {
        thisValueIterator += 1.0;
        vect1.SetPtEtaPhiE(ptr((*selCollection1)->at(iObj1))->pt,ptr((*selCollection1)->at(iObj1))->eta,ptr((*selCollection1)->at(iObj1))->phi,
                           max(ptr((*selCollection1)->at(iObj1))->energy,ptr((*selCollection1)->at(iObj1))->pt)); //Hack for "energy" of MET
        vect1_transverse.SetPtEtaPhiE(ptr((*selCollection1)->at(iObj1))->pt,0.0,ptr((*selCollection1)->at(iObj1))->phi,ptr((*selCollection1)->at(iObj1))->pt);
        thisVectorSum += vect1; thisVectorTransverseSum += vect1_transverse;
        
        if ( variable == "pt" ) thisValueSum += vect1.Pt();
        else if ( variable == "pz" ) thisValueSum += vect1.Pt();
        else if ( variable == "energy" ) thisValueSum += vect1.Pt();
        else if ( variable == "eta" ) thisValueSum += vect1.Pt();
        else if ( variable == "phi" ) thisValueSum += vect1.Pt();
        else if ( variable == "absPz" ) thisValueSum += vect1.Pt();
        else if ( variable == "absEta" ) thisValueSum += vect1.Pt();
        else if ( variable == "absPhi" ) thisValueSum += vect1.Pt();
        else if ( which_pair != "vector_sum" ) { std::cerr << " No valid entry for variable: " << variable << std::endl; continue; }
        else continue;
      }      
    }
    
    //Loop over second set of objects alone
    //for (typename collectionType2::const_iterator object2 = (*selCollection2)->begin(); object2 != (*selCollection2)->end(); ++object2) {
    for (unsigned int iObj2 = 0; iObj2 < (*selCollection2)->size(); iObj2++) {

      //Sets max number of objects, eliminates use of the same object twice 
      if ( iObj2 < max2 && (typeid(*selCollection1).name() != typeid(*selCollection2).name()
                            || iObj2 >= (*selCollection1)->size() || iObj2 >= max1 ) ) {
        thisValueIterator += 1.0;
        vect2.SetPtEtaPhiE(ptr((*selCollection2)->at(iObj2))->pt,ptr((*selCollection2)->at(iObj2))->eta,ptr((*selCollection2)->at(iObj2))->phi,
                           max(ptr((*selCollection2)->at(iObj2))->energy,ptr((*selCollection2)->at(iObj2))->pt)); //Hack for "energy" of MET
        vect2_transverse.SetPtEtaPhiE(ptr((*selCollection2)->at(iObj2))->pt,0.0,ptr((*selCollection2)->at(iObj2))->phi,ptr((*selCollection2)->at(iObj2))->pt);
        thisVectorSum += vect2; thisVectorTransverseSum += vect2_transverse;
        
        if ( variable == "pt" ) thisValueSum += vect2.Pt();
        else if ( variable == "pz" ) thisValueSum += vect2.Pt();
        else if ( variable == "energy" ) thisValueSum += vect2.Pt();
        else if ( variable == "eta" ) thisValueSum += vect2.Pt();
        else if ( variable == "phi" ) thisValueSum += vect2.Pt();
        else if ( variable == "absPz" ) thisValueSum += vect2.Pt();
        else if ( variable == "absEta" ) thisValueSum += vect2.Pt();
        else if ( variable == "absPhi" ) thisValueSum += vect2.Pt();
        else if ( which_pair != "vector_sum" ) { std::cerr << " No valid entry for variable: " << variable << std::endl; continue; }
        else continue;
      }
    }  

    //Fill the branches
    TString bName = Form("%s_%s_%s_%s", branch_name_1.c_str(), branch_name_2.c_str(), which_pair.c_str(), variable.c_str());
    if (new_name.length() > 0) bName = Form("%s", new_name.c_str());

    if ( which_pair == "avg" ) branches[bName].branchVal = thisValueSum / thisValueIterator;
    else if ( which_pair == "sum" ) branches[bName].branchVal = thisValueSum;
    else if ( which_pair == "avgAbs" ) branches[bName].branchVal = thisValueSumAbs / thisValueIterator;
    else if ( which_pair == "sumAbs" ) branches[bName].branchVal = thisValueSumAbs;
    else if ( which_pair == "absAvg" ) branches[bName].branchVal = abs(thisValueSum) / thisValueIterator;
    else if ( which_pair == "absSum" ) branches[bName].branchVal = abs(thisValueSum);
    else if ( which_pair == "vector_sum" ) {
      if ( variable == "mass" ) branches[bName].branchVal = thisVectorSum.M();
      else if ( variable == "MT" ) branches[bName].branchVal = thisVectorTransverseSum.M();
      else if ( variable == "pt" ) branches[bName].branchVal = thisVectorSum.Pt();
      else if ( variable == "pz" ) branches[bName].branchVal = thisVectorSum.Pz();
      else if ( variable == "energy" ) branches[bName].branchVal = thisVectorSum.Energy();
      else if ( variable == "eta" ) branches[bName].branchVal = thisVectorSum.Eta();
      else if ( variable == "phi" ) branches[bName].branchVal = thisVectorSum.Phi();
      else if ( variable == "absPz" ) branches[bName].branchVal = abs(thisVectorSum.Pz());
      else if ( variable == "absEta" ) branches[bName].branchVal = abs(thisVectorSum.Eta());
      else if ( variable == "absPhi" ) branches[bName].branchVal = abs(thisVectorSum.Phi());
      else std::cerr << " No valid entry for variable in vector_sum" << std::endl;
    }
    else std::cerr << " No valid entry for which_pair: " << which_pair << std::endl; 
  } //End if ( useValue || which_pair == "vector_sum" )
  
  //Loop over pairs of objects
  //for (typename collectionType1::const_iterator object1 = (*selCollection1)->begin(); object1 != (*selCollection1)->end(); ++object1) {
  for (unsigned int iObj1 = 0; iObj1 < (*selCollection1)->size(); iObj1++) {

    //for (typename collectionType2::const_iterator object2 = (*selCollection2)->begin(); object2 != (*selCollection2)->end(); ++object2) {
    for (unsigned int iObj2 = 0; iObj2 < (*selCollection2)->size(); iObj2++) {

      //Sets max number of objects, eliminates pairs of the same object (e.g. jet1_jet1) and redundant pairs (jet1_jet2 and jet2_jet1) 
      if (iObj1<max1 && iObj2<max2 && (typeid(*selCollection1).name() != typeid(*selCollection2).name() || iObj1 < iObj2) ) {
        thisPairValueIterator += 1.0;

        vect1.SetPtEtaPhiE(ptr((*selCollection1)->at(iObj1))->pt,ptr((*selCollection1)->at(iObj1))->eta,ptr((*selCollection1)->at(iObj1))->phi,
                           max(ptr((*selCollection1)->at(iObj1))->energy,ptr((*selCollection1)->at(iObj1))->pt)); //Hack for "energy" of MET
        vect2.SetPtEtaPhiE(ptr((*selCollection2)->at(iObj2))->pt,ptr((*selCollection2)->at(iObj2))->eta,ptr((*selCollection2)->at(iObj2))->phi,
                           max(ptr((*selCollection2)->at(iObj2))->energy,ptr((*selCollection2)->at(iObj2))->pt)); //Hack for "energy" of MET
        vect12 = vect1 + vect2;
        vect1_transverse.SetPtEtaPhiE(ptr((*selCollection1)->at(iObj1))->pt,0.0,ptr((*selCollection1)->at(iObj1))->phi,ptr((*selCollection1)->at(iObj1))->pt);
        vect2_transverse.SetPtEtaPhiE(ptr((*selCollection2)->at(iObj2))->pt,0.0,ptr((*selCollection2)->at(iObj2))->phi,ptr((*selCollection2)->at(iObj2))->pt);
        vect12_transverse = vect1 + vect2;

        if ( variable == "mass" ) thisPairValue = vect12.M();
        else if ( variable == "MT" ) thisPairValue = vect12_transverse.M();
        else if ( variable == "deltaR" ) thisPairValue = vect1.DeltaR(vect2);
        else if ( variable == "deltaEta" ) thisPairValue = vect1.Eta() - vect2.Eta();
        else if ( variable == "deltaPhi" ) thisPairValue = vect1.Phi() - vect2.Phi();
        else if ( variable == "absDeltaEta" ) thisPairValue = abs(vect1.Eta() - vect2.Eta());
        else if ( variable == "absDeltaPhi" ) thisPairValue = abs(vect1.Phi() - vect2.Phi());
        else if ( variable == "pt" ) thisPairValue = vect1.Pt() + vect2.Pt();
        else if ( variable == "pz" ) thisPairValue = vect1.Pz() + vect2.Pz();
        else if ( variable == "energy" ) thisPairValue = vect1.Energy() + vect2.Energy();
        else if ( variable == "eta" ) thisPairValue = vect1.Eta() + vect2.Eta();
        else if ( variable == "phi" ) thisPairValue = vect1.Phi() + vect2.Phi();
        else if ( variable == "absPz" ) thisPairValue = abs(vect1.Pz()) + abs(vect2.Pz());
        else if ( variable == "absEta" ) thisPairValue = abs(vect1.Eta()) + abs(vect2.Eta());
        else if ( variable == "absPhi" ) thisPairValue = abs(vect1.Phi()) + abs(vect2.Phi());
        else if ( variable == "vectPt" ) thisPairValue = vect12.Pt(); 
        else if ( variable == "vectPz" ) thisPairValue = vect12.Pz(); 
        else if ( variable == "vectEnergy" ) thisPairValue = vect12.Energy(); 
        else if ( variable == "vectEta" ) thisPairValue = vect12.Eta();
        else if ( variable == "vectPhi" ) thisPairValue = vect12.Phi();
        else if ( variable == "absVectPz" ) thisPairValue = abs(vect12.Pz());
        else if ( variable == "absVectEta" ) thisPairValue = abs(vect12.Eta());
        else if ( variable == "absVectPhi" ) thisPairValue = abs(vect12.Phi());
        else { std::cerr << " No valid entry for variable: " << variable << std::endl; continue; }

        thisPairValueSum += thisPairValue;
        thisPairValueSumAbs += abs(thisPairValue);

        // --Fill the branches--
        //Fills branches for each pair of objects
        if ( which_pair == "all_pairs" || which_pair == "all_pairs_abs") { 
          if ( which_pair == "all_pairs_abs" ) {
            TString bName = Form("%s_%d_%s_%d_abs_%s", branch_name_1.c_str(), iObj1+1, branch_name_2.c_str(), iObj2+1, variable.c_str());
            if (new_name.length() > 0) bName = Form("%s_%d_%d", new_name.c_str(), iObj1+1, iObj2+1);
            branches[bName].branchVal = abs(thisPairValue);
          }
          else {
            TString bName = Form("%s_%d_%s_%d_%s", branch_name_1.c_str(), iObj1+1, branch_name_2.c_str(), iObj2+1, variable.c_str());
            if (new_name.length() > 0) bName = Form("%s_%d_%d", new_name.c_str(), iObj1+1, iObj2+1);
            branches[bName].branchVal = thisPairValue;
          }
        }
        //Fills a branch for the pair of objects
        else if ( which_pair == "closest_to" ) {
          TString bName = Form("%s_%s_%s_%dp%d_%s", branch_name_1.c_str(), branch_name_2.c_str(), which_pair.c_str(), target_value_whole, target_value_fraction, variable.c_str());
          if (new_name.length() > 0) bName = Form("%s", new_name.c_str());
          
          if ( abs(thisPairValue - target_value) < abs(branches[bName].branchVal - target_value) ) {
            branches[bName].branchVal = thisPairValue;
          }
        }
        //Fills a branch for the pair of objects
        else {
          TString bName = Form("%s_%s_%s_%s", branch_name_1.c_str(), branch_name_2.c_str(), which_pair.c_str(), variable.c_str());
          if (new_name.length() > 0) bName = Form("%s", new_name.c_str());
          
          if ( which_pair == "max" ) {
            if ( thisPairValue > branches[bName].branchVal || branches[bName].branchVal == KinematicVariableConstants::DOUBLE_INIT ) {
              branches[bName].branchVal = thisPairValue;
            }
          }
          else if ( which_pair == "min" ) {
            if ( thisPairValue < branches[bName].branchVal || branches[bName].branchVal == KinematicVariableConstants::DOUBLE_INIT ) {
              branches[bName].branchVal = thisPairValue;
            }
          }
          else if ( which_pair == "maxAbs" ) { 
            if ( abs(thisPairValue) > branches[bName].branchVal || branches[bName].branchVal == KinematicVariableConstants::DOUBLE_INIT ) {
              branches[bName].branchVal = abs(thisPairValue);
            }
          }
          else if ( which_pair == "minAbs" ) {
            if ( abs(thisPairValue) < branches[bName].branchVal || branches[bName].branchVal == KinematicVariableConstants::DOUBLE_INIT )  {
              branches[bName].branchVal = abs(thisPairValue);
            }
          }
          else if ( which_pair == "absMax" ) {
            if ( thisPairValue > branches[bName].branchVal || branches[bName].branchVal == KinematicVariableConstants::DOUBLE_INIT ) {
              branches[bName].branchVal = abs(thisPairValue);
            }
          }
          else if ( which_pair == "absMin" ) {
            if ( thisPairValue < branches[bName].branchVal || branches[bName].branchVal == KinematicVariableConstants::DOUBLE_INIT ) {
              branches[bName].branchVal = abs(thisPairValue);
            }
          }
          else if ( which_pair == "avg" && usePairValue ) branches[bName].branchVal = thisPairValueSum / thisPairValueIterator;                                                                                                                      
          else if ( which_pair == "sum" && usePairValue ) branches[bName].branchVal = thisPairValueSum;
          else if ( which_pair == "avgAbs" && usePairValue ) branches[bName].branchVal = thisPairValueSumAbs / thisPairValueIterator;
          else if ( which_pair == "sumAbs" && usePairValue ) branches[bName].branchVal = thisPairValueSumAbs;
          else if ( which_pair == "absAvg" && usePairValue ) branches[bName].branchVal = abs(thisPairValueSum) / thisPairValueIterator;
          else if ( which_pair == "absSum" && usePairValue ) branches[bName].branchVal = abs(thisPairValueSum);
          else if ( which_pair != "vector_sum" && usePairValue ) { std::cerr << " No valid entry for variable " << variable << " with which_pair " << which_pair << std::endl; continue; }
          else continue;
        }
      } // end if (iObj1<max1 && iObj2<max2)
    } //End loop over object2
  } //End loop over object1

}

#endif 