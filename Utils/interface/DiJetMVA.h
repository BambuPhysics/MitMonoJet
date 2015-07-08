//--------------------------------------------------------------------------------------------------
// DiJetMVA
//
// Helper Class for DiJet V-tagging MVA
//
// Authors: L. Di Matteo (porting from K.Sung code)
//--------------------------------------------------------------------------------------------------
#ifndef MITMONOJET_UTILS_DiJetMVA_H
#define MITMONOJET_UTILS_DiJetMVA_H

#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/VertexFwd.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/XlJetFwd.h"
#include "MitAna/DataTree/interface/XlJet.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

// for QG corrections
#include "MitMonoJet/Utils/interface/QGSyst.h"

#include "TMVA/Reader.h"

class TRandom3;

namespace mithep {
  class DiJetMVA final {
    public:
      DiJetMVA();
      ~DiJetMVA(); 

      void     SetLowPtMethodName(const char *n)  { fLowPtMethodName = n;    }
      void     SetHighPtMethodName(const char *n) { fHighPtMethodName = n;   }
      void     SetApplyQGCorrection(Bool_t b)     { fApplyQGCorrection = b;  }
      void     SetDiJetPtThr(double d)            { fDiJetPtThr = d;         }

      void     Initialize( TString lowPtWeights,
			                     TString highPtWeights, 
			                     std::string QGCorrDatabase );

      Double_t MVAValue( const XlJet *jet1, const XlJet *jet2,
                         const Float_t rho, 
                         const MCParticleCol *genParticles = 0, 
                         Bool_t printDebug = kFALSE, float *vars = 0);

    protected:      
      // Helper for pull angles
      Float_t       computePullAngle(const XlJet *basisJet, const XlJet *testJet);
      Int_t         JetPartonMatch(const XlJet *jet,
                                   const MCParticleCol *genParticles,
                                   Float_t deltaR);

      TString       fLowPtMethodName;      //name of MVA training method for low pt
      TString       fHighPtMethodName;     //name of MVA training method for low pt
      Bool_t        fIsInitialized;        //is reader properly initialized?
      Bool_t        fApplyQGCorrection;    //apply QG correction?

      double        fDiJetPtThr;           //threshold dividing low pt/high pt dijets
    
      TMVA::Reader *fLowPtReader;          //low pt MVA reader
      TMVA::Reader *fHighPtReader;         //high pt MVA reader
      
      QGSyst       *fQGSyst;               //pointer to QG corrector

      //MVA variables list
      Float_t                   fMVAVar_isSig;
      Float_t                   fMVAVar_ptmjj;
      Float_t                   fMVAVar_pull1ang;
      Float_t                   fMVAVar_pull2ang;
      Float_t                   fMVAVar_qgid1;
      Float_t                   fMVAVar_qgid2; 
      Float_t                   fMVAVar_mdrop;
  };
}

#endif
