//--------------------------------------------------------------------------------------------------
// Reduce the ntuples and compute all the necessary weights for further studies and plots
//
// Authors: L. Di Matteo                                                                  (Dec 2014)
//--------------------------------------------------------------------------------------------------
#include <iostream>
#include <TSystem.h>
#include <TROOT.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPlots/Input/interface/TaskSamples.h"
#include "MitMonoJet/Core/interface/MitDMSTree.h"
#include "MitMonoJet/Core/interface/MitSynchTree.h"

#include "TLorentzVector.h"

using namespace std;
using namespace mithep;

//---
TString getEnv(const char* name);
//---
void fillOutNtuples(MitSynchTree &outtree, MitDMSTree &intree, double baseWeight);
//---
bool eventPassSelection(MitDMSTree &intree);
//==================================================================================================
void makeSynchTree(double lumi = 19700.0) 
{
  // Define output files mode
  TString outFileMode = "RECREATE";

  // Read all environment variables
  TString anaDir = getEnv("MIT_MONOJET_DIR");
  TString hstDir = getEnv("MIT_ANA_HIST");
  TString anaCfg = "boostedv-ana";
  TString prdCfg = getEnv("MIT_PROD_CFG");
  
  // Fix data list for photons
  anaCfg = "boostedv-synch";

  // Define samples
  TaskSamples* samples = new TaskSamples(prdCfg.Data(),hstDir.Data());
  samples->SetNameTxt(anaCfg.Data());
  samples->ReadFile((anaDir + TString("/config")).Data());
  vector<const Sample*> listOfSamples;
  for (UInt_t iSample=0; iSample < samples->NDataSamples(); iSample++) listOfSamples.push_back(samples->GetDataSample(iSample));
  for (UInt_t iSample=0; iSample < samples->NSamples(); iSample++) listOfSamples.push_back(samples->GetSample(iSample));  

  // Prepare pointer to outfile
  TFile *fin;
  TFile *fout;
  TString outFileName = "synch.root";
  fout = new TFile(outFileName,outFileMode);
  
  // Prepare object to store outtree
  MitSynchTree outTree;
  
  // Generate reduced trees
  // loop through the samples and produce the reduced trees
  for (UInt_t iSample=0; iSample < listOfSamples.size(); iSample++) {
    TString inFilePath = hstDir+"/"+*listOfSamples.at(iSample)->File();
    fin = new TFile(inFilePath,"READ");
    cout << "INFO -- Reading sample file: " << *listOfSamples.at(iSample)->File() << endl;
    // Prepare event weight
    double thisXsec = *listOfSamples.at(iSample)->Xsec();
    double nGenEvts = ((TH1D*)fin->FindObjectAny("hDAllEvents"))->GetEntries();
    double baseWeight = lumi*thisXsec/nGenEvts;
    // Read input tree
    MitDMSTree inTree;
    inTree.LoadTree(inFilePath,1);
    inTree.InitTree(1);
    // Start a new sample group according to cfg file legend
    if (*listOfSamples.at(iSample)->Legend() != " ") {
      // Close previous group
      if (iSample > 0) {
        fout->cd();
        outTree.tree_->Write();
        fout->Close();
        fout = TFile::Open(outFileName,"UPDATE");     
      } 
      fout->cd();
      TString outTreeName = *listOfSamples.at(iSample)->Legend();
      outTree.CreateTree(outTreeName.Data());
      outTree.InitVariables();
    }
    // Scan on input and fill output ntuple
    cout << "INFO ---> Number of events passing the selection is: ";
    fillOutNtuples(outTree,inTree,baseWeight);
    
    // Close last group
    if (iSample == (listOfSamples.size()-1)) {
      fout->cd();
      outTree.tree_->Write(); 
      fout->Close();
      break;
    }      

    fin->Close();
  }  
  
  return;
}


//==================================================================================================
TString getEnv(const char* name)
{
  if (! gSystem->Getenv(name)) {
    printf(" Environment variable: %s  not defined. EXIT!\n",name);
    return TString("");
  } 
  return TString(gSystem->Getenv(name));  
}

//==================================================================================================
void fillOutNtuples(MitSynchTree &outtree, MitDMSTree &intree, double baseWeight)
{
  // Loop over tree entries
  Int_t nEntries = intree.tree_->GetEntries();
  for ( Int_t iEntry = 0; iEntry < nEntries; iEntry++ ) {
    
    // Get this tree entry
    intree.tree_-> GetEntry(iEntry);

    // Determine if event passes selection
    if (!eventPassSelection(intree))
      continue;
          
    // Determine the outtree variables for this event
    outtree.event_= intree.event_;
    outtree.run_  = intree.run_;
    outtree.lumi_ = intree.lumi_;
    outtree.met_  = intree.metRaw_;
    outtree.phopt_ = intree.pho1_.Pt();
    outtree.phoeta_ = intree.pho1_.Eta();
    outtree.jetpt_ = intree.jet1_.Pt();
    outtree.jeteta_ = intree.jet1_.Eta();
    // Fill output tree
    outtree.tree_->Fill();
  }
  
  return;
}

//==================================================================================================
bool eventPassSelection(MitDMSTree &intree)
{
  // Trigger
  bool triggerBit = (intree.trigger_ & MitDMSTree::HLTPhoton);
                      
  // Narrow jets
  bool jetBit = (intree.jet1_.Pt() > 110 && abs(intree.jet1_.eta()) < 2.5);
    
  // Extra
  bool extraBit = true;
  extraBit = extraBit && 
            (intree.pho1_.Pt() > 160 && abs(intree.pho1_.Eta()) < 2.5);
        
  bool theDecision;
  theDecision = triggerBit && jetBit && extraBit;
  
  return theDecision;
}
