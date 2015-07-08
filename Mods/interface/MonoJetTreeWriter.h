//--------------------------------------------------------------------------------------------------
// $Id: MonoJetTreeWriter.h,v 1.10 2013/10/21 19:34:02 dimatteo Exp $
//
// MonoJetTreeWriter
//
// Authors: L. Di Matteo
//--------------------------------------------------------------------------------------------------
#ifndef MITMONOJET_MODS_MonoJetTreeWriter_H
#define MITMONOJET_MODS_MonoJetTreeWriter_H

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/MCEventInfo.h"
#include "MitAna/DataTree/interface/MCParticleFwd.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/GenJetCol.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitAna/DataTree/interface/EvtSelData.h"
#include "MitAna/DataTree/interface/TriggerObjectCol.h"
#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/PhotonFwd.h"
#include "MitAna/DataTree/interface/TrackCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitAna/DataTree/interface/PFTauCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/SuperClusterCol.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitPhysics/Utils/interface/MVAMet.h"
#include "MitPhysics/Utils/interface/MVAVBF.h"
#include "MitAna/PhysicsUtils/interface/QGTagger.h"
#include "MitPhysics/Utils/interface/RhoUtilities.h"
#include "MitPhysics/Utils/interface/VertexTools.h"
#include "MitPhysics/Utils/interface/ElectronIDMVA.h"

#include "MitMonoJet/Core/MitGPTree.h"


namespace mithep 
{

  class MonoJetTreeWriter : public BaseMod
  {
  public:
    MonoJetTreeWriter(const char *name  = "MonoJetTreeWriter", 
		      const char *title = "Selecting MonoJets");
    
    ~MonoJetTreeWriter();

    // setting all the input Names
    void                SetMetName(const char *n)         { fMetName= n;                 }
    void                SetRawMetName(const char *n)      { fRawMetName= n;              }
    void                SetMetFromBranch(bool b)          { fMetFromBranch = b;          }
    void                SetPhotonsName(const char *n)     { fPhotonsName= n;             }
    void                SetPhotonsFromBranch(bool b)      { fPhotonsFromBranch = b;      }
    void                SetElectronsName(const char *n)   { fElectronsName = n;          }
    void                SetElectronsFromBranch(bool b)    { fElectronsFromBranch = b;    }
    void                SetMuonsName(const char *n)       { fMuonsName = n;              }
    void                SetMuonsFromBranch(bool b)        { fMuonsFromBranch = b;        }
    void                SetTausName(const char *n)        { fTausName = n;               }
    void                SetTausFromBranch(bool b)         { fTausFromBranch = b;         }

    void                SetJetsName(const char *n)        { fJetsName = n;               }
    void                SetRawJetsName(const char *n)     { fRawJetsName = n;            }
    void                SetJetsFromBranch(bool b)         { fJetsFromBranch = b;         }
    void                SetQGTaggerCHS(bool b)            { fQGTaggerCHS = b;            }
    void                SetLeptonsName(const char *n)     { fLeptonsName = n;            }
    void                SetPFCandidatesName(const char *n){ fPFCandidatesName = n;       }
    void                SetPFCandidatesFromBranch(bool b) { fPFCandidatesFromBranch = b; }

    void                SetSuperClustersName(const char *n){ fSuperClustersName = n;     }
    void                SetTracksName(const char *n)      { fTracksName = n;             }
    void                SetPVName(const char *n)          { fPVName = n;                 }
    void                SetPVFromBranch(bool b)           { fPVFromBranch = b;           }
    void                SetPUInfoName(const char *n)      { fPileUpName = n;             }
    void                SetBeamspotName(const char *n)    { fBeamspotName = n;           }

    void                SetIsData (Bool_t b)              { fIsData = b;                 }
    void                SetTriggerObjectsName(const char *n) { fTriggerObjectsName = n;  }

    void                SetProcessID(Int_t n)             { fDecay = n;                  }
    void                SetFillNtupleType(Int_t d)        { fFillNtupleType= d;          }
    
    void               SetPFNoPileUpName(const char *n)   { fPFNoPileUpName  = n;        } 
    void               SetPFPileUpName(const char *n)     { fPFPileUpName  = n;          }


  protected:
    void                Process();
    void                SlaveBegin();
    void                SlaveTerminate();

  private:
    bool                IsTightMuon(const Muon *muon);
    void                CorrectMet(const float met, const float metPhi,
				   const Particle *l1, const Particle *l2,
                                   float &newMet, float &newMetPhi);

    // Private auxiliary methods... Names of the input Collections

    TString                        fEvtSelDataName;
    TString                        fRawMetName;
    TString                        fMetName;
    TString                        fPhotonsName;
    TString                        fElectronsName;
    TString                        fMuonsName;
    TString                        fTausName;
    TString                        fJetsName;
    TString                        fRawJetsName;
    TString                        fLeptonsName;
    TString                        fPFCandidatesName;
    TString                        fVertexName;
    TString                        fSuperClustersName;
    TString                        fTracksName;
    TString                        fPVName;
    TString                        fPileUpDenName;    
    TString                        fPileUpName;
    TString                        fBeamspotName;
    TString                        fMCEvInfoName;
    TString                        fMCPartName;
    TString                        fTriggerObjectsName;
    TString                        fPFNoPileUpName;
    TString                        fPFPileUpName;
			           
    Bool_t                         fIsData;
    Bool_t                         fMetFromBranch;
    Bool_t                         fPhotonsFromBranch;
    Bool_t                         fElectronsFromBranch;
    Bool_t                         fMuonsFromBranch;
    Bool_t                         fTausFromBranch;
    Bool_t                         fJetsFromBranch;
    Bool_t                         fPFCandidatesFromBranch;
    Bool_t                         fPVFromBranch;
    Bool_t                         fQGTaggerCHS;
			         
    QGTagger                      *qgTagger;

    std::vector<std::string>       fCorrectionFiles;   // list of jet correction files
    FactorizedJetCorrector        *fJetCorrector;
    JetCorrectionUncertainty      *fJetUncertainties;

    const PFMetCol                *fRawMet;
    const MetCol                  *fMet;
    //    MVAMet                        *fMVAMet;
    const PhotonCol               *fPhotons;
    const ElectronCol             *fElectrons;
    const MuonCol                 *fMuons;
    const PFTauCol                *fPFTaus;
    const JetCol                  *fJets;
    const JetCol                  *fRawJets;
    const TriggerObjectCol        *fTrigObj;
    const PFCandidateCol          *fPFCandidates;
    const TrackCol                *fTracks;
    const VertexCol               *fPV;
    const BeamSpotCol             *fBeamspot;
    const MCEventInfo             *fMCEventInfo;
    const PileupInfoCol           *fPileUp;    
    const PileupEnergyDensityCol  *fPileUpDen;
    const SuperClusterCol         *fSuperClusters; 
    const MCParticleCol           *fParticles;	        
    const EvtSelData              *fEvtSelData;
    const Vertex                  *fVertex;

    Int_t                          fDecay;
    Int_t                          fFillNtupleType;
    Int_t                          fNEventsSelected;

    TFile	                  *fOutputFile;
    MitGPTree                      fMitGPTree;

    ClassDef(MonoJetTreeWriter,1)
  };
}
#endif
