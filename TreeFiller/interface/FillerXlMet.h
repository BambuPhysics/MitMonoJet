//--------------------------------------------------------------------------------------------------
// $Id: FillerXlMet.h,v 1.9 2011/03/01 17:27:22 mzanetti Exp $
//
// FillerXlMet
//
// This module fill an extended version of the MET, based on the multivariate
// MET,with the associated  4-momentum and relative covariance matrix
//
// Authors: L.DiMatteo
//--------------------------------------------------------------------------------------------------

#ifndef MITMONOJET_TREEFILLER_FILLERXLMET_H
#define MITMONOJET_TREEFILLER_FILLERXLMET_H

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/PFTauCol.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitPhysics/Utils/interface/MVAMet.h"

#include "MitAna/DataTree/interface/XlMet.h"
#include "MitAna/DataTree/interface/XlMetFwd.h"
#include "MitAna/DataTree/interface/XlEvtSelData.h"

namespace mithep
{
  class FillerXlMet : public BaseMod
  {
    public:
      FillerXlMet(const char *name = "FillerXlMet",
                  const char *title = "XlMet Filler module");
      ~FillerXlMet();
 
      // set the data flag
      void SetIsData (Bool_t b)               { fIsData = b;                }
 
      // setting all the input Names
      void SetJetsName(const char *n)         { fJetsName = n;              }
      void SetJetsFromBranch(bool b)          { fJetsFromBranch = b;        }
      void SetMuonsName(const char *n)        { fMuonsName = n;             }
      void SetMuonsFromBranch(bool b)         { fMuonsFromBranch = b;       }
      void SetElectronsName(const char *n)    { fElectronsName = n;         }
      void SetElectronsFromBranch(bool b)     { fElectronsFromBranch = b;   }
      void SetTausName(const char *n)         { fTausName = n;              }
      void SetTausFromBranch(bool b)          { fTausFromBranch = b;        }      
      void SetPhotonsName(const char *n)      { fPhotonsName = n;           }
      void SetPhotonsFromBranch(bool b)       { fPhotonsFromBranch = b;     }      
      void SetPFCandidatesName(const char *n) { fPFCandidatesName = n;      }
      void SetPFCandidatesFromBranch(bool b)  { fPFCandidatesFromBranch = b;}
      void SetPVName(const char *n)           { fPVName = n;                }
      void SetPVFromBranch(bool b)            { fPVFromBranch = b;          }
      void SetRawMetName(const char *n)       { fRawMetName = n;            }

      // setting all the output collection name
      void SetXlMetName(const char *n)        { fXlMetName = n;             }
      void PublishOutput(Bool_t b)            { fPublishOutput = b;         }
      
    protected:
      void Process();
      void SlaveBegin();
      void SlaveTerminate();
  
    private:
      Bool_t fIsData;                      //is this data or MC?

      TString fJetsName;
      TString fMuonsName;
      TString fElectronsName;
      TString fTausName;
      TString fPhotonsName;
      TString fPFCandidatesName;
      TString fPVName;
      TString fPileUpDenName;
      TString fRawMetName;
      TString fXlMetName;

      Bool_t  fJetsFromBranch;
      Bool_t  fMuonsFromBranch;
      Bool_t  fElectronsFromBranch;
      Bool_t  fTausFromBranch;
      Bool_t  fPhotonsFromBranch;
      Bool_t  fPFCandidatesFromBranch;
      Bool_t  fPVFromBranch;
      Bool_t  fPublishOutput;
        
      const JetCol                 *fJets;
      const MuonCol                *fMuons;
      const ElectronCol            *fElectrons;
      const PFTauCol               *fPFTaus;
      const PhotonCol              *fPhotons;
      const PFCandidateCol         *fPFCandidates;
      const VertexCol              *fPV;
      const PileupEnergyDensityCol *fPileUpDen;
      const PFMetCol               *fRawMet;
      XlMetArr                     *fXlMet;      //array of XlMet
      const XlEvtSelData           *fEvtSelData; //array of XlMet

      MVAMet *fMVAMet;                           
      ClassDef(FillerXlMet, 0)                   //XlMet filler      
  };
}
#endif
