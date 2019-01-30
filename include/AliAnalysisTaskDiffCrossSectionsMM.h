// ROOT Documentation about streamers:
// https://root.cern.ch/root/htmldoc/guides/users-guide/ROOTUsersGuide.html#streamers
//
//
// mikael.mieskolainen@cern.ch, 2018
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.
//
// Skeleton of this class is based on ALICE collaboration common grid code and macros

#ifndef ALIANALYSISTASKDIFFCROSSSECTIONSMM_H
#define ALIANALYSISTASKDIFFCROSSSECTIONSMM_H

class TH1;
class TTree;
class TList;

class AliESDHeader;
class AliESDAD;
class AliESDVZERO;
class AliESDZDC;
class AliMCEvent;
class AliESDEvent;

// C++
#include <vector>

// ROOT
#include <TObject.h>
#include <TString.h>
#include <TBits.h>
#include <TClonesArray.h>
#include <TParticle.h>

// AliROOT
#include "AliStack.h"
#include "AliESDVertex.h"
#include "AliAnalysisTaskSE.h"
#include "AliTriggerAnalysis.h"
#include "AliAnalysisUtils.h"


// Grid analysis task for diffractive cross sections
// 

class AliAnalysisTaskDiffCrossSectionsMM : public AliAnalysisTaskSE {

public:


  // Constructor and destructor
  AliAnalysisTaskDiffCrossSectionsMM(const char* name="AliAnalysisTaskDiffCrossSectionsMM");
  virtual ~AliAnalysisTaskDiffCrossSectionsMM();
  
  // -----------------------------------------------------------------------
  // AliROOT 
  virtual void NotifyRun();
  virtual void UserCreateOutputObjects();   // called once at the beginning of runtime
  virtual void UserExec(Option_t* option);  // called for each event
  virtual void Terminate(Option_t*);        // called at end of analysis
  // -----------------------------------------------------------------------


  // Setters
  void SetIsMC(Bool_t b = kTRUE) { fIsMC = b; }
  void SetMCType(TString s) { fMCType = s; }
  void SetTriggerSelection(TString ts) { fTriggerSelection = ts; }

  // ROOT tree name creator
  TString GetTreeName() const {
    TString s = "TE";
    if (!fIsMC && fTriggerSelection != "") {
      s += fTriggerSelection;
      s.ReplaceAll("|", "_");
    }
    return s;
  }
  
  TString GetResultsFileName() const { return "results.root"; }

  // -----------------------------------------------------------------------
  // Event information structure
  struct EventInfo {
    EventInfo()
      : fClassMask(0)
      , fClassMaskNext50(0)
      , fBCID(0)
      , fPeriod(0)
      , fTimeStamp(0)
      , fL0Inputs(0)
      , fL1Inputs(0)
      , fL2Inputs(0)
      , fRunNumber(0)
      , fnTrklet(0)
      , fOrbitID(0) {
      fnSPDClusters[0] = fnSPDClusters[1] = 0;
    }

    void Fill(const AliESDEvent*);

    ULong64_t fClassMask;
    ULong64_t fClassMaskNext50;
    UInt_t    fBCID;
    UInt_t    fPeriod;
    UInt_t    fTimeStamp;
    UInt_t    fL0Inputs;
    UInt_t    fL1Inputs;
    UInt_t    fnSPDClusters[2]; // 0 -> inner layer, 1 -> outer layer
    Int_t     fRunNumber;
    UShort_t  fnTrklet;
    UShort_t  fL2Inputs;
    UShort_t  fOrbitID;
  };

  // -----------------------------------------------------------------------
  // AD and V0 data structure
  struct ADV0 {
    enum {
      kCside = 0,
      kAside = 1
    };

    // constructor
    ADV0() {
      for (Int_t i = 0; i < 2; ++i) {
        fTime[i] = -10240.0f;
        fCharge[i] = 0.0;
        fBB[i] = fBG[i] = -1;
        fDecisionOnline[i] = fDecisionOffline[i] = -1;
      }
    }
    void FillAD(const AliESDEvent*, AliTriggerAnalysis&);
    void FillV0(const AliESDEvent*, AliTriggerAnalysis&);

    void FillInvalid();

    Float_t    fTime[2];            //
    Float_t    fCharge[2];          //
    Char_t     fBB[2];              // 
    Char_t     fBG[2];              //
    Double32_t fDecisionOnline[2];  //[-1,3,2]
    Double32_t fDecisionOffline[2]; //[-1,3,2]
  };
  
  // -----------------------------------------------------------------------
  // ZDC data structure
  struct ZDC {
    enum {
      kCside = 0,
      kAside = 1
    };

    // Constructor
    ZDC() {
      for (Int_t i = 0; i < 2; ++i) {
        fZNHit[i] = kFALSE;
        fZPHit[i] = kFALSE;
        fZNEnergy[i] = 0.0;
        fZPEnergy[i] = 0.0;
      }
    }
   void FillZDC(const AliESDEvent*, AliTriggerAnalysis&);

   void FillInvalid();

   Bool_t     fZNHit[2];              // Hit in Neutron part
   Bool_t     fZPHit[2];              // Hit in Proton part
   Float_t    fZNEnergy[2];           // Reconstructed energy in the neutron ZDC
   Float_t    fZPEnergy[2];           // Reconstructed energy in the proton ZDC
  };
  
  // -----------------------------------------------------------------------
  // Vertex structure
  struct VtxInfo {

    // Constructor
    VtxInfo()
      : fZ(0)
      , fNcontr(-4) {}

    void Fill(const AliESDVertex*);

    Double32_t fZ;      // [-32,32,7]
    Char_t     fNcontr; //
  };

  // -----------------------------------------------------------------------
  class TreeData : public TObject {
  public:

    // Constructor
    TreeData() 
      : TObject()
      , fEventInfo()
      , fVtxInfo()
      , fV0Info()
      , fADInfo()
      , fZDCInfo()
      , fPhysSelBits(0)
      , fIsIncompleteDAQ(kFALSE)
      , fIsSPDClusterVsTrackletBG(kFALSE) {}

    EventInfo fEventInfo;
    VtxInfo   fVtxInfo;
    ADV0      fV0Info;
    ADV0      fADInfo;
    ZDC       fZDCInfo;
    UInt_t    fPhysSelBits;
    Bool_t    fIsIncompleteDAQ;
    Bool_t    fIsSPDClusterVsTrackletBG;

    ClassDef(TreeData, 1);
  };

  // -----------------------------------------------------------------------
  // MC truth class
  class MCInfo : public TObject {
  public:
    enum {
      kInvalid = -1,
      kSDL,
      kSDR,
      kDD,
      kCD,
      kND,
      kElastic
    };

    // Structure for fiducial window generator level MC data
    struct DetGenLevel {

      // Constructor
      DetGenLevel()
        : NCharged(0)
        , NNeutral(0)
        , MeanPtCharged(-1024.0f)
        , MeanPtNeutral(-1024.0f)
        , MinPtCharged(1e9)      // <- Important init values
        , MinPtNeutral(1e9)      // <- Important init values
        , MaxPtCharged(-1e9)     // <- Important init values
        , MaxPtNeutral(-1e9) {}  // <- Important init values 
      
      // Number of particles
      UShort_t NCharged;
      UShort_t NNeutral;

      // Mean pt of particles
      Float_t MeanPtCharged;
      Float_t MeanPtNeutral;

      // Minimum pt of particles
      Float_t MinPtCharged;
      Float_t MinPtNeutral;

      // Maximum pt of particles
      Float_t MaxPtCharged;
      Float_t MaxPtNeutral;
    };

    // Structure for diffractive system(s) MC data
    struct DiffSystem {

      // Constructor (these init values do not matter for the algorithms)
      DiffSystem()
        : t(0) {

          for (UInt_t i = 0; i < 2; ++i) {

            Mass[i]     = -1024.0f;
            Mass_PDG[i] = -1024.0f;

            NCharged[i] = 0;
            NNeutral[i] = 0;

            MeanPtCharged[i] = 0;
            MeanPtNeutral[i] = 0;

            MinEtaCharged[i] = 0;
            MinEtaNeutral[i] = 0;

            MaxEtaCharged[i] = 0;
            MaxEtaNeutral[i] = 0;
          }
        }

      // 2->2 Kinematics, important especially in terms of Regge phenomenology
      Float_t  t;                 // Mandelstam, physical domain is negative values in GeV^2.
      Float_t  Mass[2];           // Diffractive system mass [0] -> L, [1] -> R (calculated from the final states 4-momenta)
      Float_t  Mass_PDG[2];       // (calculated from Pythia Diffractive Pseudoparticle 4-momentum with PDG code 9902210)

      // Particle production data. N.B. one can construct e.g. the mean |pt| of all particle by
      // multiplicity weighted average of charged and neutral mean |pt| values, by linearity of the expectation (sum) operator.
      UShort_t NCharged[2];       // Charged multiplicity
      UShort_t NNeutral[2];       // Neutral

      Float_t  MeanPtCharged[2];  // Mean |pt| of charged
      Float_t  MeanPtNeutral[2];  // Neutral

      Float_t  MinEtaCharged[2];  // Minimum pseudorapidity of charged particles in the system
      Float_t  MinEtaNeutral[2];  // Neutral

      Float_t  MaxEtaCharged[2];  // Maximum pseudorapidity of charged particles in the system
      Float_t  MaxEtaNeutral[2];  // Neutral
    };

    // Constructor
    MCInfo()
      : TObject()
      , fEventType(-1)
      , EM_check_SC(0)
      , EM_check_PP(0)
      , fDiffSys()
      
      , ZDNC()
      , ADC()
      , V0C()
      , SPDC()
      , SPDA()
      , V0A()
      , ADA()
      , ZDNA()
      
      , MinEtaVisibleCharged( 1024)
      , MaxEtaVisibleCharged(-1024)
      , MaxGapVisibleCharged(-1024)

      , MinEtaVisibleNeutral( 1024)
      , MaxEtaVisibleNeutral(-1024)
      , MaxGapVisibleNeutral(-1024)
      
      , MaxGapVisibleAll(-1024) {
    }
    //virtual ~MCInfo() {}
    
    void Fill(const AliMCEvent*, TString);
    void FillComb(const AliMCEvent*);
    void CheckFiducial(const TParticle* p, DetGenLevel& det, const Double_t eta_range[], const Double_t minPt);

    void ReCalc4Vector(TLorentzVector& v, const TParticle* p);
    void CalcInitialState(TLorentzVector& p_plus, TLorentzVector& p_minus, AliStack* stack);
    void CalcPseudoSystem(AliStack* stack);
    void EMCheck(AliStack* stack);

    void FindEtaGaps(AliStack* stack);
    void GapCalculus(const std::vector<Double_t>& eta, Float_t& MinEta, Float_t& MaxEta, Float_t& MaxGap);

    void CalcSD(const TLorentzVector& p_plus, const TLorentzVector& p_minus, AliStack* stack);
    void CalcSDslashDD(const TLorentzVector& p_plus, const TLorentzVector& p_minus, AliStack* stack);
    void CalcAssignmentSDslashDD(const TLorentzVector& p_plus, const TLorentzVector& p_minus, AliStack* stack);


    void vsort(std::vector<Double_t>& v, const std::vector<UInt_t>& sort_ind);
    std::vector<UInt_t> vsortind(const std::vector<Double_t>& v);
    Double_t VectorMean(const std::vector<Double_t>& vec);
    Double_t VectorMin(const std::vector<Double_t>& vec);
    Double_t VectorMax(const std::vector<Double_t>& vec);

    // Generic event level information
    Float_t fEventType;       // [-3,5,3]

    // Diffractive system kinematics & multiplicity data
    DiffSystem fDiffSys;

    // Energy-Momentum conservation check
    Float_t EM_check_SC;     
    Float_t EM_check_PP;

    // Fiducial generator level information
    DetGenLevel ZDNC;
    DetGenLevel ADC;
    DetGenLevel V0C;
    DetGenLevel SPDC;
    DetGenLevel SPDA;
    DetGenLevel V0A;
    DetGenLevel ADA;
    DetGenLevel ZDNA;

    // Pseudorapidity gap variables with no pt cutoffs (i.e. pt > 0)
    Float_t MinEtaVisibleCharged;
    Float_t MaxEtaVisibleCharged;
    Float_t MaxGapVisibleCharged;

    Float_t MinEtaVisibleNeutral;
    Float_t MaxEtaVisibleNeutral;
    Float_t MaxGapVisibleNeutral;

    Float_t MaxGapVisibleAll;

    ClassDef(MCInfo, 1);
  };


protected:
  void SetBranches(TTree* t);


private:  
  AliAnalysisTaskDiffCrossSectionsMM(const AliAnalysisTaskDiffCrossSectionsMM&);            // not implemented
  AliAnalysisTaskDiffCrossSectionsMM& operator=(const AliAnalysisTaskDiffCrossSectionsMM&); // not implemented

  Bool_t             fIsMC;                //
  TString            fMCType;              //
  TString            fTriggerSelection;    //

  AliTriggerAnalysis fTriggerAnalysis;     //!
  AliAnalysisUtils   fAnalysisUtils;       //!

  TTree*             fTE;                  //!
  TBits              fFastOrMap;           //!
  TBits              fFiredChipMap;        //!
  AliESDVertex       fVertexSPD;           //!
  TreeData           fTreeData;            //!
  MCInfo             fMCInfo;              //!
  
  ClassDef(AliAnalysisTaskDiffCrossSectionsMM, 1);
};

#endif // ALIANALYSISTASKDIFFCROSSSECTIONSMM_H
