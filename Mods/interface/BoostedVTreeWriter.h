//--------------------------------------------------------------------------------------------------
// BoostedVTreeWriter
//
// This module writes a tree about jets corresponding to Vector bosons. See documentation here:
//
//
// Authors: C.Paus
//--------------------------------------------------------------------------------------------------
#ifndef MITMONOJET_MODS_BOOSTEDVTREEWRITER_H
#define MITMONOJET_MODS_BOOSTEDVTREEWRITER_H

#include <TFile.h>

#include "fastjet/JetDefinition.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/Pruner.hh"

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/TriggerObjectCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/PFTauCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/PhysicsUtils/interface/QGTagger.h"
#include "MitMonoJet/Core/MitGPBoostedVTree.h"

namespace mithep
{

  class BoostedVTreeWriter : public BaseMod
  {
  public:
    typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

    BoostedVTreeWriter(const char *name  = "BoostedVTreeWriter",
		       const char *title = "Vector Boson Tagging module");

    ~BoostedVTreeWriter();

    void                          SetIsData(Bool_t b)               { fIsData = b; }
    void                          SetMcPartsName(const char *n)     { fMcPartsName = n; }
    void                          SetTriggerObjsName(const char *n) { fTriggerObjsName = n; }
    void                          SetJetsName(const char *n)        { fJetsName = n; }
    void                          SetJetsFromBranch(bool b)         { fJetsFromBranch = b; }
    void                          SetPFCandidatesName(const char *n){ fPfCandidatesName = n; }
    void                          SetPFCandidatesFromBranch(bool b) { fPfCandidatesFromBranch = b; }
    void                          SetPhotonsName(const char *n)     { fPhotonsName= n; }
    void                          SetPhotonsFromBranch(bool b)      { fPhotonsFromBranch = b; }
    void                          SetPFTausName(const char *n)      { fPfTausName = n; }
    void                          SetPFTausFromBranch(bool b)       { fPfTausFromBranch = b; }
    void                          SetLeptonsName(const char *n)     { fLeptonsName  = n; } 
    void                          SetPFNoPileUpName(const char *n)  { fPfNoPileUpName  = n; } 
    void                          SetPFPileUpName(const char *n)    { fPfPileUpName  = n; }
    void                          SetQgTaggerCHS(bool b)            { fQgTaggerCHS = b; }
    void                          SetOutputName(const char *n)      { fOutputName = n; }
    void                          SetPruning(Int_t n)               { fPrune = n; }

  protected:
    void                          Process();
    void                          SlaveBegin();
    void                          SlaveTerminate();

  private:
    void                          ProcessMc();             // Deals with the MC separately
    void                          GetJetTriggerObjs();     // Fills the jet trigger objs
    Double_t                      MinTriggerDeltaR(LorentzVector jet);
    float                         GetTau(fastjet::PseudoJet &iJet,int iN, float iKappa);
    bool                          IsTightMuon(const Muon *muon);

    Bool_t                        fIsData;                 //is this data or MC?
    TString                       fMcPartsName;            //(i) name of MC particles
    const MCParticleCol          *fMcParts;	               //MC particle coll
    TString                       fTriggerObjsName;        //(i) name of trigger objects
    const TriggerObjectCol       *fTrigObjs;               //trigger objects coll handle
    TString                       fJetsName;               //(i) name of jets used to make trigger
    Bool_t                        fJetsFromBranch;         //are jets from Branch?
    const JetCol                 *fJets;                   //jets used to make the trigger
    TString                       fPfCandidatesName;       //(i) name of PF candidates coll
    Bool_t                        fPfCandidatesFromBranch; //are PF candidates from Branch?
    const PFCandidateCol         *fPfCandidates;           //particle flow candidates coll handle
    TString                       fPhotonsName;            //(i) name of photon coll
    Bool_t                        fPhotonsFromBranch;      //are photons from Branch?
    const PhotonCol              *fPhotons;                //photon coll handle
    TString                       fPfTausName;             //(i) name of PF tau coll
    Bool_t                        fPfTausFromBranch;       //are PF taus from Branch?
    const PFTauCol               *fPfTaus;                 //PF tau coll handle

    TString                       fLeptonsName;            //(i) name of the merged leptons
    TString                       fPfNoPileUpName;         //(i) name for PF no pileup candidates
    TString                       fPfPileUpName;           //(i) name for PF pileup candidates

    Bool_t                        fQgTaggerCHS;            //are the input jets CHS?
    QGTagger                     *fQgTagger;               //compute jets Q/G discrimination

    TString                       fPileUpDenName;          //(i) name for PF pileup energy density
    const PileupEnergyDensityCol *fPileUpDen;              //PU energy density handle             
    TString                       fVertexesName;           //(i) name for the good vertexes

    std::vector<const TriggerObject *>
                                  fJetTriggerObjs;         //jet trigger objects

    // Objects from fastjet we want to use
    double                        fConeSize;
    int                           fPrune;                  //apply pruning: 0-no, 1-standard CMS
    fastjet::Pruner              *fPruner;
    fastjet::JetDefinition       *fCAJetDef;
    fastjet::GhostedAreaSpec     *fActiveArea;
    fastjet::AreaDefinition      *fAreaDefinition;

    // Counters
    Int_t                         fNAnalyzed;              //number of events analyzed

    // Output tree
    TString                       fOutputName;             //(o) name of ntuple output
    TFile	                     *fOutputFile;
    MitGPBoostedVTree             fMitGPTree;

    ClassDef(BoostedVTreeWriter, 1) // Boosted Vector boson tree writer
  };
}
#endif
