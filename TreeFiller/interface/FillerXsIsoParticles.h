//--------------------------------------------------------------------------------------------------
// $Id: FillerXsIsoParticles.h,v 1.9 2011/03/01 17:27:22 dimatteo Exp $
//
// FillerXsIsoParticles
//
// This module process a collection of input particles, extracts the minimal
// information and fill output collections of XsIsoParticles
//
// Authors: L.DiMatteo
//--------------------------------------------------------------------------------------------------

#ifndef MITMONOJET_TREEFILLER_FILLERXSISOPARTICLES_H
#define MITMONOJET_TREEFILLER_FILLERXSISOPARTICLES_H

#include "MitAna/DataTree/interface/XsIsoParticleFwd.h"
#include "MitAna/DataTree/interface/XsIsoParticle.h"

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/PFTauCol.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/Particle.h"

namespace mithep
{
  class FillerXsIsoParticles : public BaseMod
  {
    public:
      FillerXsIsoParticles(const char *name = "FillerXsIsoParticles",
                           const char *title = "XsIsoParticles Filler module");
      ~FillerXsIsoParticles();

      void IsData(Bool_t b)                { fIsData = b;             }
      void FillXsMuons (Bool_t b)          { fFillXsMuons = b;        }
      void FillXsElectrons (Bool_t b)      { fFillXsElectrons = b;    }
      void FillXsTaus (Bool_t b)           { fFillXsTaus = b;         }
      void FillXsPhotons (Bool_t b)        { fFillXsPhotons = b;      }
      void PublishOutput(Bool_t b)         { fPublishOutput = b;      }
                                                                                                                                          
      void SetMuonsName(const char *n)     { fMuonsName = n;          }
      void SetMuonsFromBranch(Bool_t b)    { fMuonsFromBranch = b;    }
      void SetIsoMuonsName(const char *n)  { fIsoMuonsName = n;       }
      void SetIsoMuonsFromBranch(Bool_t b) { fIsoMuonsFromBranch = b; }
      void SetElectronsName(const char *n) { fElectronsName = n;      }
      void SetElectronsFromBranch(Bool_t b){ fElectronsFromBranch=b;  }
      void SetTausName(const char *n)      { fTausName = n;           }
      void SetTausFromBranch(Bool_t b)     { fTausFromBranch = b;     }
      void SetPhotonsName(const char *n)   { fPhotonsName = n;        }
      void SetPhotonsFromBranch(Bool_t b)  { fPhotonsFromBranch = b;  }
                                                                      
      void SetXsMuonsName(const char *n)   { fXsMuonsName = n;        }
      void SetXsElectronsName(const char *n){ fXsElectronsName = n;   }
      void SetXsTausName(const char *n)    { fXsTausName = n;         }
      void SetXsPhotonsName(const char *n) { fXsPhotonsName = n;      }
                                                                     
    protected:
      void Process();
      void SlaveBegin();
      void SlaveTerminate();
 
      void FillXsIsoParticle(XsIsoParticleArr *pXsArr, const Particle *pParticle);
      // For XsMuons
      void FillXsIsoParticle(XsIsoParticleArr *pXsArr, const Particle *pParticle,
                             Bool_t isTight, Bool_t isIso);

      // Muon collection helpers
      Bool_t IsTightMuon(const Muon *muon);
      Bool_t IsIsoMuon(const Muon *muon);

    private:
      Bool_t fIsData;                      //is this data or MC?
      Bool_t fFillXsMuons;                 //=true if xs muons are stored
      Bool_t fFillXsElectrons;             //=true if xs electrons are stored
      Bool_t fFillXsTaus;                  //=true if xs taus are stored
      Bool_t fFillXsPhotons;               //==true if xs photons are stored
      Bool_t fPublishOutput;               //=true if output collection are published

      TString fMuonsName;                  //(i) name of input muons
      Bool_t fMuonsFromBranch;             //are input muons from Branch?
      const MuonCol *fMuons;               //input muons

      TString fIsoMuonsName;               //(i) name of input isolated muons
      Bool_t fIsoMuonsFromBranch;          //are input isolated muons from Branch?
      const MuonCol *fIsoMuons;            //input isolated muons

      TString fElectronsName;              //(i) name of input electrons
      Bool_t fElectronsFromBranch;         //are input electrons from Branch?
      const ElectronCol *fElectrons;       //input electrons

      TString fTausName;                   //(i) name of input taus
      Bool_t fTausFromBranch;              //are input taus from Branch?
      const PFTauCol *fTaus;               //input taus

      TString fPhotonsName;                //(i) name of input photons
      Bool_t fPhotonsFromBranch;           //are input photons from Branch?
      const PhotonCol *fPhotons;           //input photons
 
      TString fXsMuonsName;                //name of output fXsMuons collection
      XsIsoParticleArr *fXsMuons;          //array of fXsMuons
      TString fXsElectronsName;            //name of output fXsElectrons collection
      XsIsoParticleArr *fXsElectrons;      //array of fXsElectrons
      TString fXsTausName;                 //name of output fXsTaus collection
      XsIsoParticleArr *fXsTaus;           //array of fXsTaus
      TString fXsPhotonsName;              //name of output fXsPhotons collection
      XsIsoParticleArr *fXsPhotons;        //array of fXsPhotons
      
      ClassDef(FillerXsIsoParticles, 0)    //XsIsoParticles filler      
  };
}
#endif
