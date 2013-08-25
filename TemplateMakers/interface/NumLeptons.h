#ifndef _NumLeptons_h
#define _NumLeptons_h


// just store the number of leptons

class NumLeptons: public KinematicVariable<int> {

public:
  
  NumLeptons();

  void evaluate ();


};

NumLeptons::NumLeptons  () {

  branches["numTightMuons"] = BranchInfo<int>("numTightMuons");
  branches["numLooseMuons"] = BranchInfo<int>("numLooseMuons");
  
  branches["numTightEle"] = BranchInfo<int>("numLooseMuons");
  branches["numLooseEle"] = BranchInfo<int>("numLooseMuons");
  
}

void NumLeptons::evaluate () {

  if (this->evaluatedThisEvent) return;
  evaluatedThisEvent = true;
  
  branches["numTightMuons"].branchVal = int(this->blocks->muonCollection->size());
  branches["numLooseMuons"].branchVal = int(this->blocks->looseMuonCollection->size());

  branches["numTightEle"].branchVal = int(this->blocks->eleCollection->size());
  branches["numLooseEle"].branchVal = int(this->blocks->looseEleCollection->size());
  
}



#endif 
