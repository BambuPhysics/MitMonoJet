//--------------------------------------------------------------------------------------------------
// $Id: FillerXlJets.h,v 1.9 2011/03/01 17:27:22 mzanetti Exp $
//
// FillerXlJets
//
// This module process a collection of input jet, and add extra info
// relevant for analysis (QGtagging and color pull)
//
// Authors: L.DiMatteo
//--------------------------------------------------------------------------------------------------

#ifndef MITMONOJET_TREEFILLER_FILLERXLJETS_H
#define MITMONOJET_TREEFILLER_FILLERXLJETS_H

#include <TVector2.h>

#include "MitAna/DataTree/interface/XlJetFwd.h"
#include "MitAna/DataTree/interface/XlJet.h"

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJet.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/PhysicsUtils/interface/QGTagger.h"

namespace mithep
{
  class FillerXlJets : public BaseMod
  {
    public:
      FillerXlJets(const char *name = "FillerXlJets",
                   const char *title = "XlJets Filler module");
      ~FillerXlJets();

      void IsData(Bool_t b)                { fIsData = b;           }
      void SetfQGTaggingOn(Bool_t b)       { fQGTaggingActive = b;  }
      void SetQGTaggerCHS(bool b)          { fQGTaggerCHS = b;      }
      void PublishOutput(Bool_t b)         { fPublishOutput = b;    }
                                                                                                                                        
      void SetJetsName(const char *n)      { fJetsName = n;         }
      void SetJetsFromBranch(Bool_t b)     { fJetsFromBranch = b;   }
                                                                    
      void SetXlJetsName(const char *n)    { fXlJetsName = n;       }
                                                                   
    protected:
      void Process();
      void SlaveBegin();
      void SlaveTerminate();
 
      void FillXlJet    (const PFJet *pPFJet);
      // Color pull helpers
      TVector2 GetPull  (const PFJet *inPFJet);

    private:
      Bool_t fIsData;                      //is this data or MC?
      Bool_t fQGTaggingActive;             //=true if QGTagging info is filled
      Bool_t fQGTaggerCHS;                 //=true if QGTagging weights are taken from CHS
      Bool_t fPublishOutput;               //=true if output collection are published

      TString fJetsName;                   //(i) name of input jets
      Bool_t fJetsFromBranch;              //are input jets from Branch?
      const JetCol *fJets;                 //input jets

      TString fPileUpDenName;              //(i) name of PU energy density coll
      Bool_t fPileUpDenFromBranch;         //is PU energy density from Branch?
      const PileupEnergyDensityCol *fPileUpDen; //PU energy density coll handle

      TString fVertexesName;               //(i) name of vertex coll
      Bool_t fVertexesFromBranch;          //are vertexex from Branch?
      const VertexCol *fVertexes;          //vertex coll handle
 
      TString fXlJetsName;                 //name of output fXlJets collection
      XlJetArr *fXlJets;                   //array of fXlJets
            
      // QG tagger 
      QGTagger *fQGTagger;                 //QGTagger calculator
      
      ClassDef(FillerXlJets, 0)            //XlJets, Fat and Sub, filler      
  };
}
#endif
