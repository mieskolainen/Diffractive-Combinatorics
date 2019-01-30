// Analysis class for Combinatorial Diffractive Cross Sections
//
// mikael.mieskolainen@cern.ch, 2018
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.


#ifndef COMBINATORICS_H
#define COMBINATORICS_H


// C++ headers
#include <iostream>
#include <vector>

// ROOT headers
#include "TH1D.h"
#include "TProfile.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TArrayD.h"
#include <TObject.h>
#include <TString.h>
#include <TBits.h>
#include <TClonesArray.h>


// RooUnfold
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldTUnfold.h"
#include "RooUnfoldBinByBin.h"


// Own
#include "AliAnalysisTaskDiffCrossSectionsMM.h"
#include "CombEvent.h"
#include "VecOper.h"


#define EPS 1e-12



// DEFINE VECTOR SPACE DIMENSION
#define DIM 6
#define NCOMB 64 // 2^DIMENSION

//#define DIM 8
//#define NCOMB 256 // 2^DIMENSION

// Number of bootstrap samples in EM-cross section fits for non-linear statistical uncertainty
// Larger number gives better statistical coverage (put as high as you can)
extern UInt_t N_BOOTSTRAP;

const Double_t SQRTS = 13000; // CMS energy
const Double_t MP    = 0.938; // Proton mass

// -------------------------------------------------------------------
// Detector cuts (FIXED HERE, SLIDING BELOW)

// NOTE HERE: The cuts should match the cuts made in the van der Meer
// scans, otherwise one needs to take that into account!


const UInt_t  SPD_cut[2] = {2, 2};                // C-side and A-side > number of fastOR
const Double_t ZN_cut[2] = {500000.0, 500000.0};  // C-side and A-side Analog-to-digital units (A.U.)


// TIME WINDOW CUTS
const Double_t AD_cut_T[2][2] = {{63.0, 69.0},    // C-side [min,max] (nsec)
                         	     {54.0, 60.0}};   // A-side [min,max] (nsec)
const Double_t V0_cut_T[2][2] = {{0.0,  6.0},     // C-side [min,max] (nsec)
                         	     {7.0,  14.0}};   // A-side [min,max] (nsec)


// FROM debug output: Cutscale = 0.14, [ADC,V0C,SPDC,SPDA,V0A,ADA] = [164.3, 8.8, 2, 2, 5.7, 38.3]

// CHARGE SUM (Q) CUTS
const Double_t AD_cut_Q[2][2] = {{0.1, 1e12},       // C-side [min,max] Analog-to-digital units (A.U.)
                         	     {0.1, 1e12}};      // A-side [min,max] Analog-to-digital units (A.U.) 
const Double_t V0_cut_Q[2][2] = {{0.1, 1e12},       // C-side [min,max] Analog-to-digital units (A.U.)
                         	     {0.1, 1e12}};      // A-side [min,max] Analog-to-digital units (A.U.)

// SPD noisy fired chip list (masked away)
const std::vector<UInt_t> SPD_noisy = 
{50,69,151,188,208,230,232,233,320,329,331,421,460,541,719,827,1027,1059,1070,1090,1117,1190,
39,151,235,331,402,422,460,479,529,599,650,711,742,747,780,919,1027,1190,
79,389,518,552,572,590,751,757,767,781,810,827,828,840,890,900,930,1027,1033,1060,1072,1090,1127,1169,1176,1190};

//const std::vector<UInt_t> SPD_noisy = {230};


// Analysis mode control parameters
extern Bool_t FASTSTAT;
extern Int_t FOLDING_MODE;
extern Double_t PT_MIN;
extern Bool_t GENERATOR_LEVEL;
extern Bool_t HISTOGRAMS_ON;
extern Double_t MAXEVENTS_DATA;
extern Double_t MAXEVENTS_MC;
extern Bool_t GAPFLOW_ON;
extern UInt_t GAPFLOW_N;
extern Double_t GAPFLOW_MAXSCALE; 
extern Bool_t SKIP_CD;
extern Int_t DD_XIMAX_MODE;


// Nano event object
struct evec {

	UInt_t c     = 9999;   // Observed  ID
	UInt_t c_gen = 9999;   // Generator ID
	UInt_t proc  = 9999;   // Process class
	Double_t weight = 1.0; // By default 1.0 event weight

	Double_t  M2_L = -1.0; // Diffractive mass (left system)
	Double_t  M2_R = -1.0; // Diffractive mass (right system)
};



class Combinatorics {

public:

	// Constructor and destructor
	Combinatorics(const char* filename, const char* MCname, Bool_t fIsMC, Double_t vdm_scale, Double_t vdm_scale_error, UInt_t RunNumber, const TriggerData& trdata);
	~Combinatorics();

	// Operators
	void EmptyData();
	void ReadTree();
	void Printx();
	void Printx_unfolded(int iterations, int modelnumber);
	
	void PrintCombinatorics(const std::vector<Double_t>& input_x, const std::string& nametag, Double_t scale);
	void CorrectPileup();

	void VdmScale();
	void PrintData(const CombEvent* event);
	void CalculateF();
	void Plot1D(std::vector<TCanvas*>& can,  UInt_t color, UInt_t marker);
	void Plot1DA(std::vector<TCanvas*>& can, UInt_t color, UInt_t marker);
	void Plot1DB(std::vector<TCanvas*>& can, UInt_t color, UInt_t marker);
	
	void Plot();
	void PlotCodingScheme();
	void CreateDirectories();

	// Trigger functions
	Bool_t OnlineTrigger(const CombEvent* event);
	Bool_t BGVeto(const CombEvent* event);
	Bool_t EarlyInteractionVeto(const CombEvent* event, const Double_t V0_cut[2][2], const Double_t AD_cut[2][2]);
	Bool_t LateInteractionVeto(const CombEvent* event, const Double_t V0_cut[2][2], const Double_t AD_cut[2][2]);
	
	void SPDBitsSplit(const CombEvent* event, UInt_t& nbits_inner_minus, UInt_t& nbits_inner_plus, 
											  UInt_t& nbits_outer_minus, UInt_t& nbits_outer_plus);
	
	// Signal construction
	Double_t CalculateMCReWeight(const Double_t delta, const Double_t M2_L, const Double_t M2_R, const UInt_t fEventType);
	std::vector<Bool_t> ConstructVectorDet(const CombEvent* event, 
		const Double_t V0_cut_T[2][2], const Double_t AD_cut_T[2][2], const Double_t V0_cut_Q[2][2], 
		const Double_t AD_cut_Q[2][2], const UInt_t SPD_cut[2], const Double_t ZN_cut[2]);
	std::vector<Bool_t> ConstructVectorGen(const CombEvent* event);
	

	// Generate the bootstrapped sample
	void GenerateBootStrap();


	// SETTERS >>

	// Set unfolded rates in 						     [units of event counts]
	void SetXUnfolded(const std::vector<Double_t>& x_in) { x_unf_ = x_in; }
	// Set generator level rates   						 [units of event counts]
	void SetXGen(const std::vector<Double_t>& x_in)      { x_gen_ = x_in; }

	// Set unfolded bootstrap sample in
	void SetBOOTXUnfolded(const std::vector<std::vector<Double_t> >& X_in) { BOOTX_unf_ = X_in; }


	// Set Beam-Gas input
	void SetBeamGasX(const std::vector<Double_t>& xA, 
					 const std::vector<Double_t>& xC,
					 const std::vector<Double_t>& xE, Double_t scaleA, Double_t scaleC, Double_t scaleE) {
		xA_  = xA;
		xC_  = xC;
		xE_  = xE;
		BG_A = scaleA;
		BG_C = scaleC;
		BG_E = scaleE;

		// Set the Beam-Gas corrected event rates
		std::vector<Double_t> xcorrected = VecOper::BGSubstract(x_, xA_, xC_, xE_, scaleA, scaleC, scaleE);
		SetXCorrected(xcorrected);

		// Print out statistics
		VecOper::TriggerMaskStatistics(x_, xA_, xC_, xE_, BG_A, BG_C, BG_E);
	}


	// GETTERS >>
	const std::vector<Double_t>& 				 GetX()          { return x_; }
	const std::vector<Double_t>& 				 GetXCorrected() { return x_cor_; }
	const std::vector<Double_t>& 				 GetXUnfolded()  { return x_unf_; }
	const std::vector<Double_t>& 				 GetXGen()       { return x_gen_; }


	// Get bootstrapped samples
	const std::vector<std::vector<Double_t> >&   GetBOOTXCorrected() { return BOOTX_cor_; }
	const std::vector<std::vector<Double_t> >&   GetBOOTXUnfolded()  { return BOOTX_unf_; }


	// Get Monte Carlo process likelihood matrices
	const std::vector<std::vector<Double_t> >&   GetF()          { return F_; }
	const std::vector<std::vector<Double_t> >&   GetFGen()       { return F_gen_; }
	

	// -------------------------------------------------------------
	// Cross section variables (total != fiducial != visible)
	
	Double_t GetTotSigmaInel() { return sigma_inel_tot_; }  // MC and Unfolded+Extrapoled DATA
	Double_t GetFidSigmaInel() { return sigma_inel_fid_; }  // MC and Unfolded DATA

	Double_t GetVisSigmaInel() { return sigma_inel_vis_; }  		  // MC and DATA
	Double_t GetVisSigmaInelError() { return sigma_inel_vis_error_; } // Data (vdM scan error, Luminosity)


	// Unfolding and extrapolation
	void SetTotSigmaInelUnfolded(Double_t xsec) { sigma_inel_tot_unfolded_ = xsec; }
	void SetFidSigmaInelUnfolded(Double_t xsec) { sigma_inel_fid_unfolded_ = xsec; }
	
	Double_t GetTotSigmaInelUnfolded() { return sigma_inel_tot_unfolded_; } 
	Double_t GetFidSigmaInelUnfolded() { return sigma_inel_fid_unfolded_; } 


	// MC process class rates/fractions
	std::vector<Double_t> GetTotalProcessCount(Bool_t normalize);
	std::vector<Double_t> GetVisibleProcessCount(Bool_t normalize);


	TString  GetName();
	UInt_t   GetRunNumber() { return fRunNumber_; }


	// Trigger data
	TriggerData trdata_;


	// NOTE THAT THESE ARE SET WITH PREPROCESSOR MACRO ABOVE

	const Int_t d_ = DIM; 	              // Vector space dimension such as [ADC,V0C,SPD,V0A,ADA] = 5
	const Int_t N_ = NCOMB;               // The number of combinations 2^d

	// -------------------------------------------------------------------
	// Histograms made public

	// Unfolding
	TH2D* M;       // Combinatorial Folding matrix
	TH1D* hxDet;   // Detector level
	TH1D* hxGen;   // Generator level
	TH1D* hxEmpty; // Empty for template reasons

	// Coding matrices
	TH2F* h2a;   // Coding matrix <charge>
	TH2F* h2b;   // Coding matrix std(time)

	// Data and MC Histograms
	TH1F* hSPDbit[NCOMB];
	TH2F* h2SPDFOTR[NCOMB][2];

	TH1F* hSPDFO[NCOMB][2];
	TH1F* hSPDTR[NCOMB];

	TH1F* hADCharge[NCOMB][2];
	TH1F* hADTime[NCOMB][2];

	TH1F* hV0Charge[NCOMB][2];
	TH1F* hV0Time[NCOMB][2];

	TH1F* hZDN[NCOMB][2];
	TH1F* hZDP[NCOMB][2];

	TH2F* h2ADCT[NCOMB][2];
	TH2F* h2V0CT[NCOMB][2];


	// MC only histograms
	TH1F* hDiffDeltaYSD_gene[2];

	TH1F* hDiffDeltaYSD_seen_DET[DIM][2]; // 2 ~ SDL and SDR
	TH1F* hDiffDeltaYSD_seen_GEN[DIM][2]; // 2 ~ SDL and SDR


	TH1F* h1DiffDeltaYDD_gene;
	
	TH1F* h1DiffDeltaYDD_seen_DET[DIM];
	TH1F* h1DiffDeltaYDD_seen_GEN[DIM];

	
	TH1F* hDiffMassSDLowM_gene[2];
	TH1F* hDiffMassSD_gene[2];
	
	TH1F* hDiffMassSD_seen[DIM][2];   	  // 2 ~ SDL and SDR
	TH1F* hDiffMassSD_seen_DET[DIM][2];   // 2 ~ SDL and SDR
	TH1F* hDiffMassSD_seen_GEN[DIM][2];   // 2 ~ SDL and SDR
	
	
	TH2F* h2DiffMassDD_gene;
	
	TH2F* h2DiffMassDD_seen_DET[DIM];	  // 2 ~ DD
	TH2F* h2DiffMassDD_seen_GEN[DIM];	  // 2 ~ DD
	
	
	//  Detector level <observable> as a function of invariant mass
	TProfile* hAD_M_avgCharge[2][2];
	TProfile* hV0_M_avgCharge[2][2];
	TProfile* hSPD_M_avgFO[2][2];
	TProfile* hZDN_M_avgCharge[2][2];
	TProfile* hZDP_M_avgCharge[2][2];

	TH2F* h2AD_M_Charge[2][2];
	TH2F* h2V0_M_Charge[2][2];
	TH2F* h2SPD_M_FO[2][2];
	TH2F* h2ZDN_M_Charge[2][2];
	TH2F* h2ZDP_M_Charge[2][2];

	// Detector level <observable> as a function of fiducial generator level
	TProfile* hAD_N_avgCharge[2][2];
	TProfile* hAD_Nch_avgCharge[2][2];
	TProfile* hAD_Nn_avgCharge[2][2];
	
	TProfile* hV0_N_avgCharge[2][2];
	TProfile* hV0_Nch_avgCharge[2][2];
	TProfile* hV0_Nn_avgCharge[2][2];
	
	TProfile* hSPD_N_avgFO[2][2];
	TProfile* hSPD_Nch_avgFO[2][2];
	TProfile* hSPD_Nn_avgFO[2][2];

	TProfile* hZDN_Nn_avgCharge[2][2];


	// -------------------------------------------------------------
	// Cross-Correlations
	TH2F* h2XC[64][6][6];


	// RooUnfold object (~folding matrix)
	RooUnfoldResponse* response_;

	UInt_t C_; 				     		   // Number of process classes [SDL,SDR,DD,ND,CD etc.]

	Bool_t fIsMC_;			       		   // True for MC


	// Generate new re-weighted model
	void GenerateModel(Double_t POMERON_DELTA, Double_t XI_MAX);


private:

	// FIXED INDEXING FOR C-side (negative pseudorapidity) and A-side (positive pseudorapidity)
	const int C = 0;
	const int A = 1;


	// Nano event objects
	std::vector<evec> events;


	Bool_t   MassCutOff(const evec& ev, Double_t XI_MAX);
	Double_t CalculateMCReWeight(const evec& ev, Double_t POMERON_DELTA);


	// Output
	void FillSPDAscii(const CombEvent* event, Int_t eventnr, FILE* fp);
	void FillHistograms(const CombEvent* event, const evec& ev);
	void FillGapFlow(const CombEvent* event, const evec& ev);
	void FillCSVOutput(const CombEvent* event, const evec& ev, FILE* fp);


	// Observables
	Double_t GetObservable(const CombEvent* event, UInt_t index);


	// Unfolding construction functions
	Bool_t UnfoldInMemory = kFALSE;
	void NewUnfoldResponse();
	void ConstructUnfolding();


	// This is called after residual beam gas correction [units of event counts]
	void SetXCorrected(const std::vector<Double_t>& x_in) {
		x_cor_ = x_in;

		// For unfolding >>
		// Important, clear first to remove old data.
		hxDet->Reset();

		// Over 2^N combinations
		for (UInt_t c = 0; c < x_cor_.size(); ++c) {
			for (UInt_t events = 0; events < std::round(x_cor_.at(c)); ++events) {
				hxDet->Fill(c);
			}
		}

		// Call scaling to physical dimensions now, after beam gas correction!
		VdmScale(); 
	}


	// SPD functions
	Bool_t IsNoisySPD(UInt_t chipnumber);
	void   PrintHotSPD(UInt_t c);


	std::vector<Double_t> SPD_FO_z_pos_;   // SPD FO chip z-coordinates

	UInt_t fRunNumber_;                    // Run number
	TString fFilename_;                    // Filename of TTree
	
	std::vector<std::vector<Int_t> > B_;   // Binary detector combinations matrix [checker board matrix]
	std::vector<TString> det_labels_;      // Names (strings for plots) of the detectors

	// Detector indices in the vector format, to be initialized in constructor
	Int_t ZDN_ind_[2];
	Int_t AD_ind_[2];
	Int_t V0_ind_[2];
	Int_t SPD_ind_[2];

	TString fMCName_;		  		  	   		   // Name of the data/MC source (e.g. Pythia-6.x)
	std::vector<Double_t> sigma_dim_; 	   		   // Visible cross section per vector space dimension


	// -------------------------------------------------------------
	// Signal indices
	std::vector<UInt_t>   s_ind_;      			   // Signal indices


	// -------------------------------------------------------------
	// Event rates

	std::vector<Double_t> x_;          			   // Measured Beam-Beam           [counts]
	std::vector<Double_t> xA_;                     // Measured A-side beam - Empty [counts]
	std::vector<Double_t> xC_;                     // Measured C-side beam - Empty [counts]
	std::vector<Double_t> xE_;					   // Measured Empty - Empty       [counts]


	Double_t BG_A;								   // Beam-Gas Correction factors
	Double_t BG_C;								   // *
	double_t BG_E;								   // *


	std::vector<Double_t> x_cor_;                  // Beam-Gas Corrected           [counts]
	std::vector<Double_t> x_unf_; 			   	   // Unfolded event rates vector  [counts]
	std::vector<Double_t> x_gen_;      			   // Generator level rates vector [counts] (MC only)


	// -------------------------------------------------------------
	// Boostrapped samples for statistical uncertainty

	// Bootstrapped, Beam-Gas corrected vectors [2^N x Number of Bootstrap]
	// -> We take into account also the full statistical uncertainty of each beam-empty triggers
	std::vector<std::vector<Double_t> > BOOTX_cor_;

	// Bootstrapped, Unfolded vectors [2^N x Number of Bootstrap]
	std::vector<std::vector<Double_t> > BOOTX_unf_;


	// -------------------------------------------------------------
	// Cross section variables

	Double_t sigma_vdM_scan_ = 0.0; 		   	   // van der Meer scan (V0-AND) [mb]  (DATA only)
	Double_t sigma_vdM_scan_error_ = 0.0;          // luminosity scan error

	Double_t sigma_inel_tot_ = 0.0;                // Total inelastic x-section [mb]   (MC only)
	Double_t sigma_inel_tot_error_ = 0.0;          // luminosity scan error

	Double_t sigma_inel_fid_ = 0.0;	               // Total fiducial inelastic (MC only)
	Double_t sigma_inel_vis_ = 0.0;                // Total visible inelastic  (MC and DATA)
	Double_t sigma_inel_vis_error_ = 0.0;          // Total visible inelastic error (luminosity scan) (DATA)

	Double_t sigma_inel_tot_unfolded_ = 0.0;       // Unfolded and extrapolated total inelastic [mb] (MC and DATA)
	Double_t sigma_inel_fid_unfolded_ = 0.0;       // Unfolded total fiducial inelastic x-section [mb] (MC and DATA)


	// -------------------------------------------------------------
	// Class likelihood densities and default MC fractions
	std::vector<Double_t> fMCTotalProcessCount_;   // Total event counts per process

	std::vector<std::vector<Double_t> > F_;        // Class densities (probability densities as columns for each process)
	std::vector<std::vector<Double_t> > F_gen_;    // Class densities (at generator level)


	// Running mean
	Double_t N_ev     = 1e-12;
	Double_t SPDC_sum = 0.0;
	Double_t SPDA_sum = 0.0;
	Double_t V0C_sum  = 0.0;
	Double_t V0A_sum  = 0.0;
	Double_t ADC_sum  = 0.0;
	Double_t ADA_sum  = 0.0;
	Double_t ZDC_sum  = 0.0;
	Double_t ZDA_sum  = 0.0;

	Double_t SPDC_sum2 = 0.0;
	Double_t SPDA_sum2 = 0.0;
	Double_t V0C_sum2  = 0.0;
	Double_t V0A_sum2  = 0.0;
	Double_t ADC_sum2  = 0.0;
	Double_t ADA_sum2  = 0.0;
	Double_t ZDC_sum2  = 0.0;
	Double_t ZDA_sum2  = 0.0;

	// Rapidity-Combinatorics-Flow matrix
	std::vector<std::vector<Double_t> > gapflow1;
	std::vector<std::vector<Double_t> > gapflow2;
	std::vector<std::vector<Double_t> > gapflow3;
	

	// SPD bits (these are CPU costly to calculate from TBits object, so do only once per event)
	UInt_t nbits_inner_[2] = {0};				 // Inner layer C and A side
	UInt_t nbits_outer_[2] = {0};				 // Outer layer C and A side


	void PrintGapFlowMatrix(const std::vector<std::vector<Double_t>>& matrix, std::string filename);

	//ClassDef(Combinatorics,1);        	     // ROOT system integration
};


#endif
