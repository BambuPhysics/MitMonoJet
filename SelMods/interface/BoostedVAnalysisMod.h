//--------------------------------------------------------------------------------------------------
// BoostedVAnalysisMod
//
// An analysis module for selecting jet+MET events.
//
//
// Authors: LDM
//
//--------------------------------------------------------------------------------------------------
#ifndef MITMONOJET_SELMODS_BOOSTEDVANALYSISMOD_H
#define MITMONOJET_SELMODS_BOOSTEDVANALYSISMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/EvtSelData.h"
#include "MitAna/DataTree/interface/TriggerObjectCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/PFMetCol.h"

#include "MitAna/DataTree/interface/XlEvtSelData.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class BoostedVAnalysisMod : public BaseMod
  {
  public:
    BoostedVAnalysisMod(const char *name  = "BoostedVAnalysisMod", 
                        const char *title = "Boosted V selection Module");
    ~BoostedVAnalysisMod();
    
    // setting all the input Names
    void                  SetInputMetName          (const char *n) { fMetBranchName= n;           }
    void                  SetFatJetsName           (const char *n) { fFatJetsName = n;            } 
    void                  SetJetsName              (const char *n) { fJetsName = n;               } 
    void                  SetElectronsName         (const char *n) { fElectronsName = n;          }
    void                  SetMuonsName             (const char *n) { fMuonsName = n;              }
    void                  SetPhotonsName           (const char *n) { fPhotonsName = n;            }
    void                  SetTriggerObjectsName    (const char *n) { fTriggerObjectsName = n;     }

    // decide whether to read from branch
    void                  SetMetFromBranch         (Bool_t b)      { fMetFromBranch = b;          }
    void                  SetFatJetsFromBranch     (Bool_t b)      { fFatJetsFromBranch = b;      }
    void                  SetJetsFromBranch        (Bool_t b)      { fJetsFromBranch = b;         }
    void                  SetElectronsFromBranch   (Bool_t b)      { fElectronsFromBranch = b;    }
    void                  SetMuonsFromBranch       (Bool_t b)      { fMuonsFromBranch = b;        }
    void                  SetPhotonsFromBranch     (Bool_t b)      { fPhotonsFromBranch = b;      }

    // decide whether to apply given preselections or not
    void                  ApplyResolvedPresel      (Bool_t b)      { fApplyResolvedPresel = b;    }
    void                  ApplyTopPresel           (Bool_t b)      { fApplyTopPresel = b;         }
    void                  ApplyWlepPresel          (Bool_t b)      { fApplyWlepPresel = b;        }
    void                  ApplyZmmPresel           (Bool_t b)      { fApplyZmmPresel = b;         }
    void                  ApplyZeePresel           (Bool_t b)      { fApplyZeePresel = b;         }
    void                  ApplyMetPresel           (Bool_t b)      { fApplyMetPresel = b;         }
    void                  ApplyVbfPresel           (Bool_t b)      { fApplyVbfPresel = b;         }
    void                  ApplyGjetPresel          (Bool_t b)      { fApplyGjetPresel = b;        }
    void                  ApplyFatJetPresel        (Bool_t b)      { fApplyFatJetPresel = b;      }

    // decide whether to to fill and publish the preseleciton word or not
    void                  FillAndPublishPresel     (Bool_t b)      { fFillAndPublishPresel = b;   }

    // decide whether to skip events or not if they fail preselection criteria
    void                  SkipEvents               (Bool_t b)      { fSkipEvents = b;             }

    // Setting cut values
    void                  SetMinResolvedMass       (Double_t x)    { fMinResolvedMass = x;        }
    void                  SetMaxResolvedMass       (Double_t x)    { fMaxResolvedMass = x;        }
    void                  SetMinFatJetPt           (Double_t x)    { fMinFatJetPt = x;            }
    void                  SetMinTagJetPt           (Double_t x)    { fMinTagJetPt = x;            }
    void                  SetMinVbfJetPt           (Double_t x)    { fMinVbfJetPt = x;            }
    void                  SetMinMet                (Double_t x)    { fMinMet = x;                 }
    void                  SetMinVbfMass            (Double_t x)    { fMinVbfMass = x;             }
    void                  SetMinPhotonPt           (Double_t x)    { fMinPhotonPt = x;            }
    
  protected:
    // Standard module methods
    void                  Begin();
    void                  Process();
    void                  SlaveBegin();
    void                  SlaveTerminate();
    void                  Terminate();      
    // Auxiliary methods
    void                  CorrectMet(const float met, const float metPhi,
                                     const Particle *l1, const Particle *l2,
                                     float &newMet, float &newMetPhi);
    int                   GetHLTWord( 
                          bool passSingleMuHLT,
                          bool passMonoJetHLT,
                          bool passVbfHLT,
                          bool passGjetHLT,
                          bool passDiEleHLT);
    int                   GetPreselWord( 
                          bool passTopPresel, 
                          bool passWlepPresel,
                          bool passZmmPresel,
                          bool passZeePresel,
                          bool passMetPresel, 
                          bool passVbfPresel, 
                          bool passGjetPresel,
                          bool passResolvedPresel);
                                     
    // names of the collections
    TString               fMetBranchName;
    TString               fFatJetsName;
    TString               fJetsName;
    TString               fElectronsName;
    TString               fMuonsName;
    TString               fPhotonsName;
    TString               fTriggerObjectsName;
    // logical whether to read from branch
    Bool_t                fMetFromBranch;
    Bool_t                fFatJetsFromBranch;
    Bool_t                fJetsFromBranch;
    Bool_t                fElectronsFromBranch;
    Bool_t                fMuonsFromBranch;
    Bool_t                fPhotonsFromBranch;
    // logical whether to apply given preselections or not
    Bool_t                fApplyResolvedPresel; 
    Bool_t                fApplyTopPresel; 
    Bool_t                fApplyWlepPresel;
    Bool_t                fApplyZmmPresel;
    Bool_t                fApplyZeePresel;
    Bool_t                fApplyMetPresel; 
    Bool_t                fApplyVbfPresel; 
    Bool_t                fApplyGjetPresel; 
    Bool_t                fApplyFatJetPresel;
    // logical whether to fill and publish the preseleciton word or not
    Bool_t                fFillAndPublishPresel; 
    // logical whether to skip events if they fail preselection criteria
    Bool_t                fSkipEvents;
    // hooks to the collections
    const PFMetCol       *fMet;
    const JetCol         *fFatJets; 
    const JetCol         *fJets; 
    const ElectronCol    *fElectrons;
    const MuonCol        *fMuons;
    const PhotonCol      *fPhotons;
    const TriggerObjectCol *fTrigObj;
    const EvtSelData     *fEvtSelData;
    XlEvtSelData         *fXlEvtSelData;  //extended event selection data object
    // Cuts
    Double_t              fMinResolvedMass;
    Double_t              fMaxResolvedMass;
    Double_t              fMinFatJetPt;
    Double_t              fMinTagJetPt;
    Double_t              fMinVbfJetPt;
    Double_t              fMinMet;
    Double_t              fMinVbfMass;
    Double_t              fMinVbfMet;
    Double_t              fMinPhotonPt;
    // Counters
    Long64_t              fAll;
    Long64_t              fPass;
        
    ClassDef(BoostedVAnalysisMod,5) // MonJet Selection Module
  };
}
#endif
