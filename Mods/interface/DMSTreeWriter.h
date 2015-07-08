//--------------------------------------------------------------------------------------------------
// $Id: DMSTreeWriter.h,v 1.10 2013/10/21 19:34:02 dimatteo Exp $
//
// DMSTreeWriter
//
// Authors: L. Di Matteo        (original class from G.G.Ceballos)
//--------------------------------------------------------------------------------------------------
#ifndef MITMONOJET_MODS_DMSTreeWriter_H
#define MITMONOJET_MODS_DMSTreeWriter_H

#include <TH1D.h>
#include <TFile.h>

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/MCEventInfo.h"
#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/TriggerObjectCol.h"
#include "MitPhysics/Utils/interface/VertexTools.h"
#include "MitAna/DataTree/interface/XlJetCol.h"
#include "MitAna/DataTree/interface/XlFatJetCol.h"
#include "MitAna/DataTree/interface/XlSubJetCol.h"
#include "MitAna/DataTree/interface/XlEvtSelData.h"
#include "MitAna/DataTree/interface/XsIsoParticleCol.h"
#include "MitMonoJet/Utils/interface/DiJetMVA.h"

#include "MitMonoJet/Core/MitDMSTree.h"


namespace mithep 
{

  class DMSTreeWriter : public BaseMod
  {
  public:
    typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

    DMSTreeWriter(const char *name  = "DMSTreeWriter", 
                  const char *title = "DMS Tree writer");
    
    ~DMSTreeWriter();

    // setting all the input Names
    void    SetRawMetName(const char *n)            { fRawMetName= n;              }
    void    SetPhotonsName(const char *n)           { fPhotonsName= n;             }
    void    SetPhotonsFromBranch(bool b)            { fPhotonsFromBranch = b;      }
    void    SetElectronsName(const char *n)         { fElectronsName = n;          }
    void    SetElectronsFromBranch(bool b)          { fElectronsFromBranch = b;    }
    void    SetMuonsName(const char *n)             { fMuonsName = n;              }
    void    SetMuonsFromBranch(bool b)              { fMuonsFromBranch = b;        }
    void    SetTausName(const char *n)              { fTausName = n;               }
    void    SetTausFromBranch(bool b)               { fTausFromBranch = b;         }
    void    SetJetsName(const char *n)              { fJetsName = n;               }
    void    SetJetsFromBranch(bool b)               { fJetsFromBranch = b;         }
    void    SetFatJetsName(const char *n)           { fFatJetsName = n;            }
    void    SetFatJetsFromBranch(bool b)            { fFatJetsFromBranch = b;      }
    void    SetSubJetsName(const char *n)           { fSubJetsName = n;            }
    void    SetSubJetsFromBranch(bool b)            { fSubJetsFromBranch = b;      }
                                                    
    void    SetPVName(const char *n)                { fPVName = n;                 }
    void    SetPVFromBranch(bool b)                 { fPVFromBranch = b;           }
    void    SetPUInfoName(const char *n)            { fPileUpName = n;             }

    void    SetTriggerObjectsName(const char *n)    { fTriggerObjectsName = n;     }
    void    SetIsData (Bool_t b)                    { fIsData = b;                 }

    void    SetInPUHistoFileName(const char *n)     { fPUInputFileName = n;        }
    void    SetTargetPUHistoFileName(const char *n) { fPUTargetFileName = n;       }

  protected:
    void    Process();
    void    SlaveBegin();
    void    SlaveTerminate();
            
  private:              
    // Private auxiliary methods... Names of the input Collections

    TString                        fEvtSelDataName;
    TString                        fRawMetName;
    TString                        fPhotonsName;
    TString                        fElectronsName;
    TString                        fMuonsName;
    TString                        fTausName;
    TString                        fJetsName;
    TString                        fFatJetsName;
    TString                        fSubJetsName;
    TString                        fPVName;
    TString                        fPileUpDenName;    
    TString                        fPileUpName;
    TString                        fMCEventInfoName;
    TString                        fMCParticlesName;
    TString                        fTriggerObjectsName;
			           
    Bool_t                         fIsData;
    Bool_t                         fMetMVAFromBranch;
    Bool_t                         fPhotonsFromBranch;
    Bool_t                         fElectronsFromBranch;
    Bool_t                         fMuonsFromBranch;
    Bool_t                         fTausFromBranch;
    Bool_t                         fJetsFromBranch;
    Bool_t                         fFatJetsFromBranch;
    Bool_t                         fSubJetsFromBranch;
    Bool_t                         fPVFromBranch;
			         
    const PFMetCol                *fRawMet;
    const XsIsoParticleCol        *fPhotons;
    const XsIsoParticleCol        *fElectrons;
    const XsIsoParticleCol        *fMuons;
    const XsIsoParticleCol        *fPFTaus;
    const XlJetCol                *fJets;
    const XlFatJetCol             *fFatJets;
    const XlSubJetCol             *fSubJets;
    const VertexCol               *fPV;
    const PileupInfoCol           *fPileUp;    
    const PileupEnergyDensityCol  *fPileUpDen;
    const MCEventInfo             *fMCEventInfo;
    const MCParticleCol           *fMCParticles;
    const XlEvtSelData            *fEvtSelData;
    const TriggerObjectCol        *fTrigObj;

    // Auxiliary
    TString                        fPUInputFileName;     // input PU histo file
    TString                        fPUTargetFileName;    // target PU histo file
    TH1D                          *fPUInput;             // input PU histo
    TH1D                          *fPUTarget;            // target PU histo
    const TH1D                    *fPUWeight;            // target PU histo
    Float_t                        PUWeight(Float_t npu);// PU reweighting function
    // Met helpers  
    void                           CorrectMet(const float met, const float metPhi, 
                                              const LorentzVector& l1, const LorentzVector& l2,
                                              float& newMet, float& newMetPhi);
    void                           CorrectMet(const float met, const float metPhi, 
                                              const LorentzVector& l1,
                                              float& newMet, float& newMetPhi,
                                              Bool_t applyFprintCorrection = kFALSE);
    float                          GetMt(const LorentzVector& l1, const float met, const float metPhi);
    float                          GetFootprint(const float pt, const float eta);
    // Gen Level Info Analysis   
    void                           getGenLevelInfo(MitDMSTree& tree);
    // Hlt Objects Matcher   
    Bool_t                         IsHLTMatched(LorentzVector& v,
                                                TriggerObject::ETriggerObject,
                                                Float_t deltaR = 0.5,
                                                Float_t minPt = 0.,
                                                Float_t maxEta = 99.);
    // Jet cleaner   
    Bool_t                         fjetIsCleaned(LorentzVector& v,
                                                 Float_t deltaR = 0.5);
    // Jet-Parton Id Matcher   
    Int_t                          JetPartonMatch(LorentzVector& v,
                                                  Float_t deltaR = 0.5);
    // Jet-Jets delta phi calculation   
    Float_t                        GetJetJetsDphi(LorentzVector& v);

    // Get B-tag via FatJet-Standard Jet matching   
    Float_t                        GetFatJetBtag(LorentzVector& v,
                                                 Float_t deltaR = 0.3);

    DiJetMVA                      *fDiJetMVA; //MVA helper for resolved category
    
    // Get jet uncertainties   
    JetCorrectionUncertainty      *fJetUncertainties;
    JetCorrectionUncertainty      *fFatJetUncertainties;
                                                
    TFile                         *fOutputFile;
    MitDMSTree                     fMitDMSTree;

    ClassDef(DMSTreeWriter,3)
  };
}
#endif
