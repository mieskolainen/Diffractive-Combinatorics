// Combinatorics analysis class for Combinatorial Diffractive Cross Sections
//
//
// mikael.mieskolainen@cern.ch, 2018
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.


// AliROOT's headers
#include "AliCDBManager.h"
#include "AliSPDGeometryUtils.h"

// ROOT headers
#include "TStopwatch.h"
#include "TLatex.h"
#include "TArrayD.h"
#include "TH2F.h"
#include "THStack.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TMatrixD.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "TVector3.h"
#include "TPrincipal.h"
#include "TF1.h"
#include "TLegend.h"


// Own headers
#include "AliAnalysisTaskDiffCrossSectionsMM.h"
#include "CombinatoricsSuper.h"
#include "Combinatorics.h"
#include "CombEvent.h"
#include "CutFlow.h"
#include "VecOper.h"
#include "PileUp.h"

// Own helper etc. functions
using namespace VecOper;

// -----------------------------------------------------------------------
// Number of bootstrap samples in EM-cross section fits for non-linear statistical uncertainty
// Larger number gives better statistical coverage (larger the better, default 300)
UInt_t N_BOOTSTRAP = 300;
// Exact multinomial (default = kTRUE) versus approximate Poisson (kFALSE) statistics
Bool_t FASTSTAT = kFALSE;
// -----------------------------------------------------------------------


// -----------------------------------------------------------------------
// Folding matrix construction mode (default 1, PT_MIN = 0.05)
Int_t FOLDING_MODE = 1;
Double_t PT_MIN = 0.05;
// -----------------------------------------------------------------------


// Analysis at generator (particle) level (default KFALSE) THIS IS ONLY FOR DEBUG!!
Bool_t GENERATOR_LEVEL = kFALSE;


// -----------------------------------------------------------------------
// Max number of events (fraction, 0...1, default 1.0)
Double_t MAXEVENTS_DATA = 1.0;
Double_t MAXEVENTS_MC   = 1.0;
// Histogramming/plotting (default kFALSE)
Bool_t HISTOGRAMS_ON = kFALSE;
// -----------------------------------------------------------------------


// -----------------------------------------------------------------------
// GAP-FLOW (default KFALSE)
Bool_t GAPFLOW_ON = kFALSE;
// Gapflow discretization (default = 100)
UInt_t GAPFLOW_N = 100;
// MAXSCALE times the mean of the observable in the GapFlow (fefault = 10)
Double_t GAPFLOW_MAXSCALE = 10.0; 
// -----------------------------------------------------------------------


// Skip central diffraction (default kTRUE)
Bool_t SKIP_CD = kTRUE;
// Double Diffraction kinematic cutoff mode (default = 0, 1 separate kinematics)
Int_t DD_XIMAX_MODE = 0;



// Powerlaw fit function: param[0]*1/(x)^(1+param[1])
static Double_t Reggefit(Double_t* x, Double_t* par) {

	Double_t val = par[0]*1/(EPS + pow(x[0], 1 + par[1]));

   	//		TMath::Max( 1.0e-10, 2*par[0]*par[0]) );
   	// To keep it positive definite when x < par[1]
   	// Due to shifting, the new domain is [shift,inf)
  	if (val < 0) {
  		val = 1.0e-10;
  	}
   	return val;
}


// Constructors
Combinatorics::Combinatorics(
	const char* filename, const char* MCname, Bool_t fIsMC, 
	Double_t vdm_scale, Double_t vdm_scale_error, UInt_t RunNumber, const TriggerData& trdata) {

	printf("\nCombinatorics::Combinatorics:: \n");

	// -------------------------------------------------------------------
	// Create a mapping of SPD FO chips
/*
	AliCDBManager::Instance()->SetDefaultStorage("raw://");
	//AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
	//AliCDBManager::Instance()->SetDefaultStorage("local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB");
	AliSPDGeometryUtils::LoadGeom(RunNumber);
	TVector3 p[4]; // 4 elements = 4 corners of the FO chips

	for (UInt_t i = 0; i < 1200; ++i) {
		AliSPDGeometryUtils::GetChipXYZ(i, p);
		SPD_FO_z_pos_.push_back( p[0].Z() );

		printf("FO Chip %d/%d : z = %0.3f \n", i, 1199, SPD_FO_z_pos_.at(i) );
	}
	// -------------------------------------------------------------------
*/
	
	// Temporary solution
	Double_t pos_array[20] = {13.160, 11.760, 10.360, 8.960, 7.570 ,6.069 ,4.669 ,3.269 ,1.869 ,0.479 ,
		                     -1.083 ,-2.483 ,-3.883 ,-5.283 ,-6.673 ,-8.174 ,-9.574 ,-10.974 ,-12.374 ,-13.764};

	Int_t k = 0;
	for (UInt_t i = 0; i < 1200; ++i) {
		SPD_FO_z_pos_.push_back( pos_array[k]);
		++k;
		if (k == 20) { k = 0; }
		//printf("FO Chip %d/%d : z = %0.3f \n", i, 1199, SPD_FO_z_pos_.at(i) );
	}

	if (fIsMC == kTRUE) {
		sigma_inel_tot_ 	  = vdm_scale; // Total inelastic in MC
		sigma_inel_tot_error_ = vdm_scale_error;
	} else {
		sigma_vdM_scan_ 	  = vdm_scale; // Visible inelastic in DATA
		sigma_vdM_scan_error_ = vdm_scale_error;
	}

	fFilename_  = filename;
	fMCName_    = MCname;
	fIsMC_      = fIsMC;
	trdata_     = trdata;
	fRunNumber_ = RunNumber;


	// Number of process classes
	C_ = 5; // SDL, SDR, DD, CD, ND

	// -------------------------------------------------------------------
	// Init matrix of size [2^d x d]
	std::vector<std::vector<Int_t> > temp(std::pow(2,d_), std::vector<Int_t>(d_));
	B_ = temp;


	// Init vector of size [2^d x 1]
	//TMatrixD temp2(pow(2, d_), 1);
	std::vector<Double_t> temp2(std::pow(2,d_), 0.0);
	x_ 			= temp2;
	x_cor_      = temp2;
	x_unf_      = temp2;
	x_gen_ 		= temp2;

	// --------------------------------------------------------------------
	// SET VECTOR SPACE PARAMETERS
/*
	// D = 8 Schema
	if (d_ == 8) {

		// Init detector <-> to combinatorics vector x element mapping
		ZDN_ind_[0] = 0;
		AD_ind_[0]  = 1;
		V0_ind_[0]  = 2;
		SPD_ind_[0] = 3;
		SPD_ind_[1] = 4;
		V0_ind_[1]  = 5;
		AD_ind_[1]  = 6;
		ZDN_ind_[1] = 7;

		// Detector labels
		std::vector<TString> det_labels = {"ZDNC","ADC","V0C","SPDC","SPDA","VOA","ADA","ZDNA"};
		det_labels_ = det_labels;
	}
*/
	// D = 6 Schema
	if (d_ == 6) {

		// Init detector <-> to combinatorics vector x element mapping
		AD_ind_[C]  = 0;
		V0_ind_[C]  = 1;
		SPD_ind_[C] = 2;
		SPD_ind_[A] = 3;
		V0_ind_[A]  = 4;
		AD_ind_[A]  = 5;

		// Detector labels
		std::vector<TString> det_labels = {"ADC","V0C","SPDC","SPDA","V0A","ADA"};
		det_labels_ = det_labels;
	}

/*
	// D=5 Schema
	if (d_ == 5) {
		
		// Init detector <-> to combinatorics vector x element mapping
		AD_ind_[0]  = 0;
		V0_ind_[0]  = 1;
		SPD_ind_    = 2;
		V0_ind_[1]  = 3;
		AD_ind_[1]  = 4;

		// Detector labels
		std::vector<TString> det_labels = {"ADC","V0C","SPD","VOA","ADA"};
		det_labels_ = det_labels;
	}
*/
	
	std::vector<Double_t> temp3(C_, 0);
	fMCTotalProcessCount_   = temp3;

	std::vector<Double_t> temp4(d_, 0);
	sigma_dim_ = temp4;

	// Init class density matrix F with zeros
	std::vector<std::vector<Double_t> > F_temp(std::pow(2,d_), std::vector<Double_t>(C_, 0.0));
	F_     = F_temp;
	F_gen_ = F_temp;
	
	
	// Create binary Matrix
	B_ = VecOper::ConstructB(d_);

	// Init Gap Flow Matrix
	std::vector<std::vector<Double_t> > gapflow_temp(GAPFLOW_N, std::vector<Double_t>(std::pow(2,d_), 0.0));
	gapflow1 = gapflow_temp;
	gapflow2 = gapflow_temp;
	gapflow3 = gapflow_temp;

	// -------------------------------------------------------------------
	// Initialize histograms

	// Unfolding histograms always on
	hxDet   = new TH1D(Form("hxDet_%s_%d", fMCName_.Data(), fRunNumber_),   "Detector level",  N_, -0.1, N_-0.1);
	hxGen   = new TH1D(Form("hxGen_%s_%d", fMCName_.Data(), fRunNumber_),   "Generator level", N_, -0.1, N_-0.1);
	hxEmpty = new TH1D(Form("hxEmpty_%s_%d", fMCName_.Data(), fRunNumber_), "Empty", N_, -0.1, N_-0.1);


	if (HISTOGRAMS_ON) {
	
		printf("Combinatorics::Combinatorics:: Initialize histograms..\n");

		// Loop over different 2^d combinations
		for (UInt_t c = 0; c < pow(2,d_); ++c) {

			hSPDbit[c] = new TH1F(Form("hSPDbit_%d_%s", c, fMCName_.Data()), Form("ID%d : SPD; SPD FiredChip chip number; Counts", c), 1200, 0, 1200);
			hSPDTR[c] = new TH1F(Form("hSPDTR_%d_%s", c, fMCName_.Data()), Form("ID%d : SPD; Tracklets [#]; Counts", c), 90, 0, 90);

			// Loop over C and A side
			for (UInt_t k = 0; k < 2; ++k) {

				TString side = "";
				side = k == 0 ? "C" : "A";

				hSPDFO[c][k] 	= new TH1F(Form("hSPDFO%d_%d_%s", k, c, fMCName_.Data()), Form("ID%d : SPD%s; Fired chips [#]", c, side.Data()), 300, 0, 600);

				hZDN[c][k] 		= new TH1F(Form("hZDN%d_%d_%s", k, c, fMCName_.Data()), Form("ID%d: ZDN%s;Signal [A/D-C units];", c, side.Data()), 500, 0, 3500);
				hZDP[c][k]   	= new TH1F(Form("hZDP%d_%d_%s", k, c, fMCName_.Data()), Form("ID%d: ZDP%s;Signal [A/D-C units];", c, side.Data()), 500, 0, 3500);

				hADCharge[c][k] = new TH1F(Form("hADCharge%d_%d_%s", k, c, fMCName_.Data()), Form("ID%d: AD%s;Charge [A/D-C units];", c, side.Data()), 1536, 0, 6144);
				hADTime[c][k]   = new TH1F(Form("hADTime%d_%d_%s", k, c, fMCName_.Data()), Form("ID%d: AD%s;Time [ns];", c, side.Data()), 500, 50, 75);

				hV0Charge[c][k] = new TH1F(Form("hV0Charge%d_%d_%s", k, c, fMCName_.Data()), Form("ID%d: V0%s;Charge [A/D-C units];", c, side.Data()), 512, 0, 1024);
				hV0Time[c][k]   = new TH1F(Form("hV0Time%d_%d_%s", k, c, fMCName_.Data()), Form("ID%d: V0%s;Time [ns];", c, side.Data()), 500, -5, 20);

				h2ADCT[c][k] 	= new TH2F(Form("h2ADCT%d_%d_%s", k, c, fMCName_.Data()), Form("ID%d: AD%s;Time [ns];Charge [A/D-C units]", c, side.Data()), 200, 50, 75, 200, 0, 6144);
				h2V0CT[c][k] 	= new TH2F(Form("h2V0CT%d_%d_%s", k, c, fMCName_.Data()), Form("ID%d: V0%s;Time [ns];Charge [A/D-C units]", c, side.Data()), 200, -5, 20, 200, 0, 1024);
			
				h2SPDFOTR[c][k] = new TH2F(Form("h2SPDFOTR%d_%d_%s", k, c, fMCName_.Data()), Form("ID%d : SPD%s; Tracklets (full SPD) [#]; Fired chips [#]", c, side.Data()), 100, 0, 100, 100, 0, 100);
			}
		}

		// 2D-Cross Correlations
		std::vector<std::string> det_XC = {"ADC","V0C","SPDC","SPDA","V0A","ADA"};
		std::vector<std::string> obs_XC = {"Charge","Charge","FiredChip","FiredChip","Charge","Charge"};

		std::vector<Double_t> min_XC = {0, 0, 0, 0, 0, 0};
		std::vector<Double_t> max_XC = {6144/10.0, 1024/10.0, 600/10.0, 600/10.0, 1024/10.0, 6144/10.0};
		std::vector<UInt_t> bin_XC   = {120, 120, 120, 120, 120, 120};


		// Loop over different 2^d combinations
		for (UInt_t c = 0; c < pow(2,d_); ++c) {

			// First observable
			for (Int_t i = 0; i < d_; ++i) {

				// Second observable
				for (Int_t j = i; j < d_; ++j) {

					h2XC[c][i][j] = new TH2F(Form("h2XC_%d_%d_%d_%s", c, i, j, fMCName_.Data()),
						Form("ID%d; %s %s; %s %s", c, 
							det_XC.at(i).c_str(), obs_XC.at(i).c_str(), 
							det_XC.at(j).c_str(), obs_XC.at(j).c_str()), 
							bin_XC.at(i), min_XC.at(i), max_XC.at(i), bin_XC.at(j), min_XC.at(j), max_XC.at(j));
				}
			}
		}
	

		// MC Histograms
		if (fIsMC_) {

			// Combinatorial Folding Matrix
			M = new TH2D(Form("FoldingMatrix_%s_%d", fMCName_.Data(), fRunNumber_), Form("%s;#Omega_{k} (Detector level);#Omega_{k} (Generator level)", fMCName_.Data()), N_, -0.5, N_-1+0.5, N_, -0.5, N_-1+0.5);

			// Left and right diffractive systems
			for (UInt_t k = 0; k < 2; ++k) {
				TString side = (k == 0 ? "SDL" : "SDR");
				TString edge = (k == 0 ? "right" : "left"); // The opposite edge than system dissociation direction
				
				hDiffMassSDLowM_gene[k] = new TH1F(Form("hSDDiffLowMGene%d_%s",   k, fMCName_.Data()), Form("%s %s;M^{2}_{X} (GeV^{2});dN/dM^{2}_{X} (events)", fMCName_.Data(), side.Data()), 100, 1.0, 25);
				hDiffMassSD_gene[k]     = new TH1F(Form("hSDDiffMassGene%d_%s",   k, fMCName_.Data()), Form("%s %s;M^{2}_{X} (GeV^{2});dN/dM^{2}_{X} (events)", fMCName_.Data(), side.Data()), 1000, 1.0, 10000);
				hDiffDeltaYSD_gene[k]   = new TH1F(Form("hSDDiffDeltaYGene%d_%s", k, fMCName_.Data()), Form("%s %s;System %s #LTY#GT edge;dN/d#LT#DeltaY#GT (events)", fMCName_.Data(), side.Data(), edge.Data()), 50, 0-10, 10);
			
				// <observable> as a function of mass
				// Loop over C and A side
				for (UInt_t kk = 0; kk < 2; ++kk) {

					TString CA = "";
					CA = kk == 0 ? "C" : "A";
					hAD_M_avgCharge[k][kk]  = new TProfile(Form("hAD_M_avgCharge%d_%d_%s",  k, kk, fMCName_.Data()), Form("%d %s;log10(M);AD%s  #LT #SigmaQ #GT (A-D/C units)",  fRunNumber_, side.Data(), CA.Data()), 100, std::log10(1.0), std::log10(SQRTS), 0, 6144);
					hV0_M_avgCharge[k][kk]  = new TProfile(Form("hV0_M_avgCharge%d_%d_%s",  k, kk, fMCName_.Data()), Form("%d %s;log10(M);V0%s  #LT #SigmaQ #GT (A-D/C units)",  fRunNumber_, side.Data(), CA.Data()), 100, std::log10(1.0), std::log10(SQRTS), 0, 1024);
					hSPD_M_avgFO[k][kk]     = new TProfile(Form("hSPD_M_avgFO%d_%d_%s", 	k, kk, fMCName_.Data()), Form("%d %s;log10(M);SPD%s  #LT FiredChip #GT (#)",            fRunNumber_, side.Data(), CA.Data()), 100, std::log10(1.0), std::log10(SQRTS), 0, 600);
					hZDN_M_avgCharge[k][kk] = new TProfile(Form("hZDN_M_avgCharge%d_%d_%s", k, kk, fMCName_.Data()), Form("%d %s;log10(M);ZDN%s  #LT #SigmaQ #GT (A-D/C units)", fRunNumber_, side.Data(), CA.Data()), 100, std::log10(1.0), std::log10(SQRTS), 0, 3500);
					hZDP_M_avgCharge[k][kk] = new TProfile(Form("hZDP_M_avgCharge%d_%d_%s", k, kk, fMCName_.Data()), Form("%d %s;log10(M);ZDP%s  #LT #SigmaQ #GT (A-D/C units)", fRunNumber_, side.Data(), CA.Data()), 100, std::log10(1.0), std::log10(SQRTS), 0, 3500);
					
					h2AD_M_Charge[k][kk]    = new TH2F(Form("h2AD_M_Charge%d_%d_%s",  k, kk, fMCName_.Data()), Form("%s %s;M_{X} (GeV);AD%s #SigmaQ (A-D/C units)",  fMCName_.Data(), side.Data(), CA.Data()), 100, 1.0, 200, 100, 0, 4000);
					h2V0_M_Charge[k][kk]    = new TH2F(Form("h2V0_M_Charge%d_%d_%s",  k, kk, fMCName_.Data()), Form("%s %s;M_{X} (GeV);V0%s #SigmaQ (A-D/C units)",  fMCName_.Data(), side.Data(), CA.Data()), 100, 1.0, 200, 100, 0, 400);
					h2SPD_M_FO[k][kk]       = new TH2F(Form("h2SPD_M_FO%d_%d_%s", 	  k, kk, fMCName_.Data()), Form("%s %s;M_{X} (GeV);SPD%s FiredChip (#)",  fMCName_.Data(), side.Data(), CA.Data()), 		   100, 1.0, 200, 100, 0, 100);
					h2ZDN_M_Charge[k][kk]   = new TH2F(Form("h2ZDN_M_Charge%d_%d_%s", k, kk, fMCName_.Data()), Form("%s %s;M_{X} (GeV);ZDN%s #SigmaQ (A-D/C units)", fMCName_.Data(), side.Data(), CA.Data()), 100, 1.0, 200, 100, 0, 3500);
					h2ZDP_M_Charge[k][kk]   = new TH2F(Form("h2ZDP_M_Charge%d_%d_%s", k, kk, fMCName_.Data()), Form("%s %s;M_{X} (GeV);ZDP%s #SigmaQ (A-D/C units)", fMCName_.Data(), side.Data(), CA.Data()), 100, 1.0, 200, 100, 0, 3500);
				}

				// <observable> as a function of all/charged/neutral particles
				// Loop over C and A side
				for (UInt_t kk = 0; kk < 2; ++kk) {

					TString CA = "";
					CA = kk == 0 ? "C" : "A";
					hAD_N_avgCharge[k][kk]    = new TProfile(Form("hAD_N_avgCharge%d_%d_%s",  k, kk, fMCName_.Data()), Form("%d %s (charged + neutral);N_{all};AD%s  #LT #SigmaQ #GT (A-D/C units)",  fRunNumber_, side.Data(), CA.Data()), 21, 0, 20, 0, 6144);
					hAD_Nch_avgCharge[k][kk]  = new TProfile(Form("hAD_Nch_avgCharge%d_%d_%s",  k, kk, fMCName_.Data()), Form("%d %s (charged exclusive);N_{ch};AD%s  #LT #SigmaQ #GT (A-D/C units)",  fRunNumber_, side.Data(), CA.Data()), 21, 0, 20, 0, 6144);
					hAD_Nn_avgCharge[k][kk]   = new TProfile(Form("hAD_Nn_avgCharge%d_%d_%s",  k, kk, fMCName_.Data()), Form("%d %s (neutral exclusive);N_{n} ;AD%s  #LT #SigmaQ #GT (A-D/C units)",  fRunNumber_, side.Data(), CA.Data()),   21, 0, 20, 0, 6144);
					
					hV0_N_avgCharge[k][kk]    = new TProfile(Form("hV0_N_avgCharge%d_%d_%s",  k, kk, fMCName_.Data()), Form("%d %s (charged + neutral);N_{};V0%s  #LT #SigmaQ #GT (A-D/C units)",  fRunNumber_, side.Data(), CA.Data()), 21, 0, 20, 0, 1024);
					hV0_Nch_avgCharge[k][kk]  = new TProfile(Form("hV0_Nch_avgCharge%d_%d_%s",  k, kk, fMCName_.Data()), Form("%d %s (charged exclusive);N_{ch};V0%s  #LT #SigmaQ #GT (A-D/C units)",  fRunNumber_, side.Data(), CA.Data()), 21, 0, 20, 0, 1024);
					hV0_Nn_avgCharge[k][kk]   = new TProfile(Form("hV0_Nn_avgCharge%d_%d_%s",  k, kk, fMCName_.Data()), Form("%d %s (neutral exclusive);N_{n};V0%s  #LT #SigmaQ #GT (A-D/C units)",  fRunNumber_, side.Data(), CA.Data()),   21, 0, 20, 0, 1024);
					
					hSPD_N_avgFO[k][kk]       = new TProfile(Form("hSPD_N_avgFO%d_%d_%s", 	k, kk, fMCName_.Data()), Form("%d %s (charged + neutral);N_{all};SPD%s  #LT FiredChip #GT (#)",            fRunNumber_, side.Data(), CA.Data()),   21, 0, 20, 0, 600);
					hSPD_Nch_avgFO[k][kk]     = new TProfile(Form("hSPD_Nch_avgFO%d_%d_%s", 	k, kk, fMCName_.Data()), Form("%d %s (charged exclusive);N_{ch};SPD%s  #LT FiredChip #GT (#)",            fRunNumber_, side.Data(), CA.Data()), 21, 0, 20, 0, 600);
					hSPD_Nn_avgFO[k][kk]      = new TProfile(Form("hSPD_Nn_avgFO%d_%d_%s", 	k, kk, fMCName_.Data()), Form("%d %s (neutral exclusive);N_{n};SPD%s  #LT FiredChip #GT (#)",            fRunNumber_, side.Data(), CA.Data()),      21, 0, 20, 0, 600);
					
					hZDN_Nn_avgCharge[k][kk]  = new TProfile(Form("hZDN_Nn_avgCharge%d_%d_%s",  k, kk, fMCName_.Data()), Form("%d %s;N_{n};ZDN%s  #LT #SigmaQ #GT (A-D/C units)",  fRunNumber_, side.Data(), CA.Data()), 21, 0, 20, 0, 3500);
				}
			}


			// Generated
			h2DiffMassDD_gene   = new TH2F(Form("h2DDDiffMassGene_%s",  fMCName_.Data()), Form("%s DD;M^{2}_{X} (GeV^{2});M^{2}_{Y} (GeV^{2})", fMCName_.Data()), 50, 1.0, 10000, 50, 1.0, 10000);
			h1DiffDeltaYDD_gene = new TH1F(Form("h1DDDifDeltaYGene_%s", fMCName_.Data()), Form("%s DD;#LT#DeltaY#GT;dN/d#LT#DeltaY#GT (events)", fMCName_.Data()), 50, 0, 19);


			// Different detectors
			for (Int_t i = 0; i < d_; ++i) {

					// Left and Right systems
					for (Int_t k = 0; k < 2; ++k) {
						TString side = (k == 0 ? "SDL" : "SDR");
						TString edge = (k == 0 ? "right" : "left"); // The opposite edge than system dissociation direction
				
						hDiffMassSD_seen_GEN[i][k]   = new TH1F(Form("hDiffMassSeen_GEN%d_%d_%s", i, k, fMCName_.Data()), Form("%s : %s Process %s;M^{2}_{X} (GeV^{2});Fiducial Acceptance (red), Detector level Eff x Acc (blue)", det_labels_[i].Data(), fMCName_.Data(), side.Data()), 1000, 1.0, 10000);
						hDiffDeltaYSD_seen_GEN[i][k] = new TH1F(Form("hDiffDeltaYSeen_GEN%d_%d_%s", i, k, fMCName_.Data()), Form("%s : %s Process %s;System %s #LTY#GT edge;Fiducial Acceptance (red), Detector level Eff x Acc (blue)", det_labels_[i].Data(), fMCName_.Data(), side.Data(), edge.Data()), 50, -10, 10);

						hDiffMassSD_seen_DET[i][k]   = new TH1F(Form("hDiffMassSeen_DET%d_%d_%s", i, k, fMCName_.Data()), Form("%s : %s Process %s;M^{2}_{X} (GeV^{2});Fiducial Acceptance (red), Detector level Eff x Acc (blue)", det_labels_[i].Data(), fMCName_.Data(), side.Data()), 1000, 1.0, 10000);
						hDiffDeltaYSD_seen_DET[i][k] = new TH1F(Form("hDiffDeltaYSeen_DET%d_%d_%s", i, k, fMCName_.Data()), Form("%s : %s Process %s;System %s #LTY#GT edge;Fiducial Acceptance (red), Detector level Eff x Acc (blue)", det_labels_[i].Data(), fMCName_.Data(), side.Data(), edge.Data()), 50, -10, 10);

						hDiffMassSD_seen_GEN[i][k]->SetLineColor(46); // Red
						hDiffDeltaYSD_seen_GEN[i][k]->SetLineColor(46); // Red
					}

					h2DiffMassDD_seen_GEN[i]   = new TH2F(Form("h2DiffMassSeen_GEN%d_%s", i, fMCName_.Data()), Form("%s DD - %s Efficiency;M^{2}_{X} (GeV^{2});M^{2}_{Y} (GeV^{2})", fMCName_.Data(), det_labels_[i].Data()), 50, 1.0, 10000, 50, 1.0, 10000);
					h1DiffDeltaYDD_seen_GEN[i] = new TH1F(Form("h1DiffDeltaYSeen_GEN%d_%s", i, fMCName_.Data()), Form("%s DD - %s;#LT#DeltaY#GT; Efficiency", fMCName_.Data(), det_labels_[i].Data()), 50, 0, 19);

					h2DiffMassDD_seen_DET[i]   = new TH2F(Form("h2DiffMassSeen_DET%d_%s", i, fMCName_.Data()), Form("%s DD - %s Efficiency;M^{2}_{X} (GeV^{2});M^{2}_{Y} (GeV^{2})", fMCName_.Data(), det_labels_[i].Data()), 50, 1.0, 10000, 50, 1.0, 10000);
					h1DiffDeltaYDD_seen_DET[i] = new TH1F(Form("h1DiffDeltaYSeen_DET%d_%s", i, fMCName_.Data()), Form("%s DD - %s;#LT#DeltaY#GT; Efficiency", fMCName_.Data(), det_labels_[i].Data()), 50, 0, 19);

					// Error saving
					/*
					hDiffMass_gen[i]->Sumw2();
					hDiffMass_seen[i]->Sumw2();

					hDiffMass_seen_ADC[i]->Sumw2();
					hDiffMass_seen_ADA[i]->Sumw2();

					hDiffMass_seen_V0C[i]->Sumw2();
					hDiffMass_seen_V0A[i]->Sumw2();
					*/
			}
		}
	} // Histograms on

}

// New unfold response function
void Combinatorics::NewUnfoldResponse() {

	if (UnfoldInMemory)
		delete response_;
	
	UnfoldInMemory = kTRUE;
	//RooUnfoldResponse response (nbins_measured, x_lo_measured, x_hi_measured,
	//                            nbins_true,     x_lo_true,     x_hi_true);
	response_ = new RooUnfoldResponse(N_, -0.1, N_-0.1); // These offsets are just to make binning right
	hxDet->Reset();
	hxGen->Reset();
}


// Construct unfolding matrix/response functions
void Combinatorics::ConstructUnfolding() {

	// UNWEIGHTED = Directly loop over MC events
	// All variables here unweighted
	for (UInt_t i = 0; i < events.size(); ++i) {

		Int_t c     = events.at(i).c;     // Detector level ID
		Int_t c_gen = events.at(i).c_gen; // Generator level ID
		//Int_t proc  = events.at(i).proc;  // Process

		// Visible event
		if (c != 0) {
			response_->Fill(c, c_gen);    // Measured signal
			hxDet->Fill(c);               // Detector level
			hxGen->Fill(c_gen);           // Generator level

		// Not measured due to efficiency/acceptance lost
		} else {
			response_->Miss(c_gen);    
			// hxDet->Fill(c);            // DO NOT FILL, AS IN DATA
			hxGen->Fill(c_gen);
		}
	}

	// RE-WEIGHTED DESCRIPTION
	// ** implement here **
}


// Destructor
Combinatorics::~Combinatorics() {

	// Always on
	delete hxDet;
	delete hxGen;
	delete hxEmpty;

	if (UnfoldInMemory)
		delete response_;

	// Delete histograms
	if (HISTOGRAMS_ON) {

		// Loop over different 2^d combinations
		for (UInt_t c = 0; c < std::pow(2,d_); ++c) {

			delete hSPDbit[c];
			delete hSPDTR[c];

			// Loop over C and A side
			for (UInt_t k = 0; k < 2; ++k) {

				delete hSPDFO[c][k];
				delete hZDN[c][k];
				delete hZDP[c][k];

				delete hADCharge[c][k];
				delete hADTime[c][k];

				delete hV0Charge[c][k];
				delete hV0Time[c][k];

				delete h2ADCT[c][k];
				delete h2V0CT[c][k];

				delete h2SPDFOTR[c][k];
			}
		}

		// Loop over different 2^d combinations
		for (UInt_t c = 0; c < pow(2,d_); ++c) {

			// First observable
			for (Int_t i = 0; i < d_; ++i) {

				// Second observable
				for (Int_t j = i; j < d_; ++j) {

					delete h2XC[c][i][j];
				}
			}
		}
	

		// MC Histograms
		if (fIsMC_) {

			// Combinatorial Folding Matrix
			delete M;

			// Left and right diffractive systems
			for (UInt_t k = 0; k < 2; ++k) {

				delete hDiffMassSDLowM_gene[k];
				delete hDiffMassSD_gene[k];
				delete hDiffDeltaYSD_gene[k];

				// <observable> as a function of mass
				// Loop over C and A side
				for (UInt_t kk = 0; kk < 2; ++kk) {

					delete hAD_M_avgCharge[k][kk];
					delete hV0_M_avgCharge[k][kk];
					delete hSPD_M_avgFO[k][kk];
					delete hZDN_M_avgCharge[k][kk];
					delete hZDP_M_avgCharge[k][kk];

					delete h2AD_M_Charge[k][kk];
					delete h2V0_M_Charge[k][kk];
					delete h2SPD_M_FO[k][kk];
					delete h2ZDN_M_Charge[k][kk];
					delete h2ZDP_M_Charge[k][kk];
				}

				// <observable> as a function of all/charged/neutral particles
				// Loop over C and A side
				for (UInt_t kk = 0; kk < 2; ++kk) {

					delete hAD_N_avgCharge[k][kk];
					delete hAD_Nch_avgCharge[k][kk];
					delete hAD_Nn_avgCharge[k][kk];

					delete hV0_N_avgCharge[k][kk];
					delete hV0_Nch_avgCharge[k][kk];
					delete hV0_Nn_avgCharge[k][kk];

					delete hSPD_N_avgFO[k][kk];
					delete hSPD_Nch_avgFO[k][kk];
					delete hSPD_Nn_avgFO[k][kk];

					delete hZDN_Nn_avgCharge[k][kk];
				}
			}

			// Generated
			delete h2DiffMassDD_gene;
			delete h1DiffDeltaYDD_gene;

			// Different detectors
			for (Int_t i = 0; i < d_; ++i) {

					// Left and Right systems
					for (Int_t k = 0; k < 2; ++k) {

						delete hDiffMassSD_seen_GEN[i][k];
						delete hDiffDeltaYSD_seen_GEN[i][k];

						delete hDiffMassSD_seen_DET[i][k];
						delete hDiffDeltaYSD_seen_DET[i][k];

						delete hDiffMassSD_seen_GEN[i][k];
						delete hDiffDeltaYSD_seen_GEN[i][k];
					}

					delete h2DiffMassDD_seen_GEN[i];
					delete h1DiffDeltaYDD_seen_GEN[i];

					delete h2DiffMassDD_seen_DET[i];
					delete h1DiffDeltaYDD_seen_DET[i];
			}
		}
	}
}


void Combinatorics::PlotCodingScheme() {

	if (HISTOGRAMS_ON) {

		printf("Combinatorics::PlotCodingScheme::\n");

	//const char* ID[n] = {"7","6","5","4","3","2","1","0"};

	std::vector<UInt_t> ID(std::pow(2,d_),0);
	for (UInt_t i = 0; i < ID.size()-1; ++i) {
		ID.at(i) = N_-i-1;
	}

	TCanvas* c1 = new TCanvas("c1","c1",700,5000);
	//c1->SetGrid();
	c1->SetLeftMargin(0.17);
	c1->SetTopMargin(0.04);
	c1->SetBottomMargin(0.04);

	h2a = new TH2F(Form("h2a_%s", fMCName_.Data()), Form("Vector space B^{%d}, %s, #sigma_{inel}^{vis} = %0.2f mb", d_, fMCName_.Data(), sigma_inel_vis_), d_,0,d_, N_,0,N_);
	h2b = new TH2F(Form("h2b_%s", fMCName_.Data()), Form("Vector space B^{%d}, %s, #sigma_{inel}^{vis} = %0.2f mb", d_, fMCName_.Data(), sigma_inel_vis_), d_,0,d_, N_,0,N_);

	Int_t k = 0;
	for (Int_t i = N_-1; i >= 0; --i) {

		// Get binary expansion (vector) for this k = 0...2^d-1
		std::vector<Bool_t> binvec = Ind2Vec(k, d_);

		for (Int_t j = 0; j < d_; ++j) {

			Double_t value1 = 0;
			Double_t value2 = 0;


			// Coding Scheme D = 8
			if (d_ == 8) {
				if (j == 0) {
				value1 = (!binvec.at(j) == kTRUE ? hZDN[N_-k-1][0]->GetMean() : 0);
				value2 = 0;
				}
				if (j == 1) {
				value1 = (!binvec.at(j) == kTRUE ? hADCharge[N_-k-1][0]->GetMean() : 0);
				value2 = (!binvec.at(j) == kTRUE ? hADTime[N_-k-1][0]->GetStdDev() : 0);
				}
				if (j == 2) {
				value1 = (!binvec.at(j) == kTRUE ? hV0Charge[N_-k-1][0]->GetMean() : 0);
				value2 = (!binvec.at(j) == kTRUE ? hV0Time[N_-k-1][0]->GetStdDev() : 0);
				}
				if (j == 3) {
				value1 = (!binvec.at(j) == kTRUE ? hSPDFO[N_-k-1][0]->GetMean() : 0);
				value2 = 0;
				}
				if (j == 4) {
				value1 = (!binvec.at(j) == kTRUE ? hSPDFO[N_-k-1][1]->GetMean() : 0);
				value2 = 0;
				}
				if (j == 5) {
				value1 = (!binvec.at(j) == kTRUE ? hV0Charge[N_-k-1][1]->GetMean() : 0);
				value2 = (!binvec.at(j) == kTRUE ? hV0Time[N_-k-1][1]->GetStdDev() : 0);
				}
				if (j == 6) {
				value1 = (!binvec.at(j) == kTRUE ? hADCharge[N_-k-1][1]->GetMean() : 0);
				value2 = (!binvec.at(j) == kTRUE ? hADTime[N_-k-1][1]->GetStdDev() : 0);
				}
				if (j == 7) {
				value1 = (!binvec.at(j) == kTRUE ? hZDN[N_-k-1][1]->GetMean() : 0);
				value2 = 0;
				}
			}

			// Coding Scheme 6
			if (d_ == 6) {
				if (j == 0) {
				value1 = (!binvec.at(j) == kTRUE ? hADCharge[N_-k-1][0]->GetMean() : 0);
				value2 = (!binvec.at(j) == kTRUE ? hADTime[N_-k-1][0]->GetStdDev() : 0);
				}
				if (j == 1) {
				value1 = (!binvec.at(j) == kTRUE ? hV0Charge[N_-k-1][0]->GetMean() : 0);
				value2 = (!binvec.at(j) == kTRUE ? hV0Time[N_-k-1][0]->GetStdDev() : 0);
				}
				if (j == 2) {
				value1 = (!binvec.at(j) == kTRUE ? hSPDFO[N_-k-1][0]->GetMean() : 0);
				value2 = 0;
				}
				if (j == 3) {
				value1 = (!binvec.at(j) == kTRUE ? hSPDFO[N_-k-1][1]->GetMean() : 0);
				value2 = 0;
				}
				if (j == 4) {
				value1 = (!binvec.at(j) == kTRUE ? hV0Charge[N_-k-1][1]->GetMean() : 0);
				value2 = (!binvec.at(j) == kTRUE ? hV0Time[N_-k-1][1]->GetStdDev() : 0);
				}
				if (j == 5) {
				value1 = (!binvec.at(j) == kTRUE ? hADCharge[N_-k-1][1]->GetMean() : 0);
				value2 = (!binvec.at(j) == kTRUE ? hADTime[N_-k-1][1]->GetStdDev() : 0);
				}
			}

			h2a->Fill(Form("%s", det_labels_[j].Data()), Form("#sigma_{%02d}", ID[k]), value1);
			h2b->Fill(Form("%s", det_labels_[j].Data()), Form("#sigma_{%02d}", ID[k]), value2);
		}	// x_sect_.at(n-k-1))
		++k;
	}

	h2a->SetMarkerSize(0.8);
	h2a->LabelsDeflate("X");
	h2a->LabelsDeflate("Y");
	h2a->LabelsOption("v");

	// Change y-axis text labels to include millibarns
	for (Int_t i = 1; i < N_+1; ++i) {
		h2a->GetYaxis()->SetBinLabel(i, Form("%02d = %0.2f mb", ID[i-1], x_.at(ID[i-1])) );
	}
	// Change x-axis text labels to include millibarns
	for (Int_t i = 1; i < d_+1; ++i) {
		h2a->GetXaxis()->SetBinLabel(i, Form("%s (%0.2f mb)", det_labels_[i-1].Data(), sigma_dim_.at(i-1)) );
	}

   	h2b->SetMarkerSize(0.8);
	h2b->LabelsDeflate("X");
	h2b->LabelsDeflate("Y");
	h2b->LabelsOption("v");


	h2a->Draw("TEXT SAME");

	// Add upper text box
	//TPaveText* pt = new TPaveText(0.1,0.8,0.9,0.9);
   	//pt->AddText("#LTcharge#GT    #LTcharge#GT    #LTcharge#GT    #FO#GT    #FO#GT    #LTcharge#GT    #LTcharge#GT    #LTcharge#GT");
   	//pt->Draw("SAME");


   	// Add Lower text box
	//TPaveText* pt2 = new TPaveText(0.01070689,0.04968919,1.569517,14.70807);
   	//pt2->AddText(Form("#sigma_{inel}^{vis} = %0.2f mb", sigma_inel_vis_));
   	//pt2->Draw("SAME");

	c1->SaveAs(Form("./figures_xsec/%d/Matrix/CodingScheme_%s.pdf", fRunNumber_, fMCName_.Data()) );

	// Change Y-axis text labels back to without millibarns
	for (Int_t i = 1; i < N_+1; ++i) {
		h2a->GetYaxis()->SetBinLabel(i, Form("#sigma_{%02d}", ID[i-1]) );
	}
	// Change X-axis text labels back to without millibarns
	for (Int_t i = 1; i < d_+1; ++i) {
		h2a->GetXaxis()->SetBinLabel(i, Form("%s", det_labels_[i-1].Data()) );
	}

	delete c1;
	}
}


// Generation of bootstrapped samples, taking into account the beam-gas combinations
void Combinatorics::GenerateBootStrap() {

	printf("Combinatorics::GenerateBootStrap:: \n");

	// Make sure it is clear
	std::vector<std::vector<Double_t> > temp(x_.size(), std::vector<Double_t>(N_BOOTSTRAP, 0));
	BOOTX_cor_ = temp;

	// Beam-Beam trigger mask
	std::vector<std::vector<Double_t> > BOOT_B = CreateBootStrapSample(nvec(x_),  std::round(vsum(x_)),  N_BOOTSTRAP, FASTSTAT);
	
	// For MC, only Beam-Beam -> no other
	if (fIsMC_) {
		BOOTX_cor_ = BOOT_B;
		return;
	}

	// A-Beam mask
	std::vector<std::vector<Double_t> > BOOT_A = CreateBootStrapSample(nvec(xA_), std::round(vsum(xA_)), N_BOOTSTRAP, FASTSTAT);
	// C-Beam mask
	std::vector<std::vector<Double_t> > BOOT_C = CreateBootStrapSample(nvec(xC_), std::round(vsum(xC_)), N_BOOTSTRAP, FASTSTAT);
	// Empty-Empty mask
	std::vector<std::vector<Double_t> > BOOT_E = CreateBootStrapSample(nvec(xE_), std::round(vsum(xE_)), N_BOOTSTRAP, FASTSTAT);

	printf("Beam-Gas Corrected sample generation... \n");

	Double_t counter = 0.1; // For printing out the progression
	printf("0 "); fflush(stdout);

	// Loop over boostrap samples
	for (UInt_t j = 0; j < N_BOOTSTRAP; ++j) {

		if ( (j /(double)N_BOOTSTRAP) > counter) {
			printf("%0.0f ", counter * 100); fflush(stdout);
			counter += 0.1;
		}

		// Get new realization
		std::vector<Double_t> Bx = VecOper::getcolvec(BOOT_B, j);
		std::vector<Double_t> Ax = VecOper::getcolvec(BOOT_A, j);
		std::vector<Double_t> Cx = VecOper::getcolvec(BOOT_C, j);
		std::vector<Double_t> Ex = VecOper::getcolvec(BOOT_E, j);

		// Beam-Gas substraction
		std::vector<Double_t> y  = VecOper::BGSubstract(Bx, Ax, Cx, Ex, BG_A, BG_C, BG_E);

		// Fill it as a column in the matrix
		VecOper::setcolvec(BOOTX_cor_, y, j);
	}
	printf("100 %%\n"); // at the end of event loop

}

// Correct pileup
void Combinatorics::CorrectPileup() {

	// Global OR rate
	const double R = trdata_.L0bLMb;

	// Pileup correction object
	PileUp pu;

	// Loop over columns
	for (uint j = 0; j < BOOTX_cor_.size(); ++j) {

		if (j > 0)
			pu.SetVerbose(false);

		// Get column vector out
		std::vector<double> y = VecOper::getcolvec(BOOTX_cor_, j);

		// Pileup inversion, returns 2^N count vector
		std::vector<double> y_corr = pu.MapCounts("Inverse", y, R);
		
		// Fill it as a column vector back
		VecOper::setcolvec(BOOTX_cor_, y_corr, j);
	}

}



// Print data
void Combinatorics::Printx() {

	printf("Combinatorics::Printx:: \n");

   	// OUTPUT .csv for external analysis
	TString output_file( Form("./figures_xsec/%d/Ascii/%s_x_rates.csv", fRunNumber_, fMCName_.Data()));
		
	FILE* fp;
	fp = fopen(output_file, "w");
	if (!fp) {
		printf("Combinatorics::Printx:: Error, could not open %s output! \n", output_file.Data());
		return;
	}
	fprintf(fp, "#ID (int, 0...2^N-1), event rate, sigma_inel_vis = %0.3f \n", sigma_inel_vis_);

	// Combination rates out
	for (Int_t c = 0; c < std::pow(2,d_); ++c) {
		fprintf(fp, "%d,%0.0f\n", c, x_.at(c) );
	}
	
	fclose(fp);
}


// Print unfolded rates
void Combinatorics::Printx_unfolded(int iterations, int modelnumber) {

	// printf("Combinatorics::Printx_unfolded:: \n");

   	// OUTPUT .csv for external analysis
	TString output_file( Form("./figures_xsec/%d/Ascii/%s_x_unfolded_rates_iter_%d_model_%d.csv", 
		fRunNumber_, fMCName_.Data(), iterations, modelnumber) );
	
	FILE* fp;
	fp = fopen(output_file, "w");
	if (!fp) {
		printf("Combinatorics::Printx_unfolded:: Error, could not open %s output! \n", output_file.Data());
		return;
	}
	fprintf(fp, "ID (int, 0...2^N-1), Events, sigma_inel_fid_unfolded_ = %0.2f, sigma_inel_tot_unfolded_ = %0.2f \n", sigma_inel_fid_unfolded_, sigma_inel_tot_unfolded_ );

	// Combination rates out
	for (Int_t c = 0; c < std::pow(2,d_); ++c) {
		fprintf(fp, "%d,%0.0f\n", c, x_unf_.at(c) );
	}
	
	fclose(fp);
}


// Calculate visible cross sections
void Combinatorics::VdmScale() {

	// Nullify
	sigma_inel_vis_ 		= 0.0;
	sigma_inel_vis_error_   = 0.0;

	// Choose the combinations
	Double_t sum         	   = EPS;
	Double_t sigma_scale 	   = 0.0;
	Double_t sigma_scale_error = 0.0;

	if (fIsMC_) { // MC
		sigma_scale       = sigma_inel_tot_;
		printf("Combinatorics::VdmScale:: MC mode \n");

		// Loop over all combinations
		for (UInt_t r = 0; r < x_cor_.size(); ++r) {
			sum += x_cor_.at(r);
		}

	} else { 	 // DATA
		sigma_scale        = sigma_vdM_scan_;
		sigma_scale_error  = sigma_vdM_scan_error_;
		printf("Combinatorics::VdmScale:: Data mode \n");
		printf("V0-AND vdM scan = %0.2f +- %0.2f mb \n",
						sigma_vdM_scan_, sigma_vdM_scan_error_);

		// Loop over all combinations and 
		for (UInt_t r = 0; r < x_cor_.size(); ++r) {

			// Choose V0C = 1 & V0A = 1 (V0_AND) combinations
			if (B_.at(r).at(V0_ind_[C]) == 1 && B_.at(r).at(V0_ind_[A]) == 1) {
				sum += x_cor_.at(r);
			}
		}
	}

	// Scale each combination
	for (UInt_t i = 0; i < x_cor_.size(); ++i) {

		// Visible
		if (i > 0) {
			sigma_inel_vis_     += x_cor_.at(i) * sigma_scale / sum;
		}
	}

	// Error on visible cross section
	sigma_inel_vis_error_ = sigma_inel_vis_ * (sigma_scale_error / sigma_scale);

/*
	printf("\nLeft-Right comparison sequence: \n");
	std::vector<UInt_t> lrseq = LRsequence(d_);

	for (UInt_t i = 0; i < x_.size(); ++i) {
		printf("[%d %d %d %d %d] <-> [%d %d %d %d %d] = ", B_[i][0], B_[i][1], B_[i][2], B_[i][3], B_[i][4],
			B_[lrseq.at(i)][0], B_[lrseq.at(i)][1], B_[lrseq.at(i)][2], B_[lrseq.at(i)][3], B_[lrseq.at(i)][4]);

		printf("[%02d <-> %02d] \t [%0.3f %0.3f] [mb] \t (%02d/%02d = %0.3f) \n", i, lrseq.at(i), x_sect_.at(i), x_sect_.at(lrseq.at(i)), 
			i, lrseq.at(i), x_sect_.at(i) / (x_sect_.at(lrseq.at(i)) + EPS) );
	}
	printf("Total = %0.2f [mb] \n\n", total);


	printf("CA ratio control: CA-parity reflected and the larger to sum:\n");
	total = 0;
	for (UInt_t i = 0; i < x_.size(); ++i) {
		// Choose the larger one
		UInt_t index = x_sect_.at(i) > x_sect_.at(lrseq.at(i)) ? i : lrseq.at(i);
		total += x_sect_.at(index);
	}
	printf("Total = %0.2f [mb] \n\n", total);
*/

	printf("\n");
	printf("\n");

	if (fIsMC_) {
	printf("%s :: Total inelastic:                        %0.3f mb \n",   fMCName_.Data(), sigma_inel_tot_);
	}
	
	printf("%s :: Visible inelastic (minbias):            %0.3f mb \n\n", fMCName_.Data(), sigma_inel_vis_);
	printf("\n");

	std::string nametag = "visible minbias";
	PrintCombinatorics(x_, nametag, sigma_inel_vis_);

}


void Combinatorics::PrintCombinatorics(const std::vector<Double_t>& input_x, const std::string& nametag, Double_t scale) {

  std::vector<UInt_t> xx(input_x.size(), 0);
  for (UInt_t i = 0; i < input_x.size(); ++i) {
  	xx.at(i) = (UInt_t) input_x.at(i);
  }

  // Set 0-bin to zero (it has values in MC (non-visible), don't want to have them here for comparison purposes)
  xx.at(0) = 0;
  printf("\n\nSystematic Subspace Combinatorics (mb): \n\n");


  printf("Full %s: \n", nametag.c_str());
  {
    std::vector<int> subind = {AD_ind_[C], V0_ind_[C], SPD_ind_[C], SPD_ind_[A], V0_ind_[A], AD_ind_[A]};
    printf("[ADC V0C SPDC SPDA V0A ADA]: \n");
    Subspace(xx, subind, scale);
  }
  printf("Subspace of %s: \n", nametag.c_str());
  {
    std::vector<int> subind = {AD_ind_[C], V0_ind_[C], V0_ind_[A], AD_ind_[A]};
    printf("[ADC V0C V0A ADA]: \n");
    Subspace(xx, subind, scale);
  }
  printf("Subspace of %s: \n", nametag.c_str());
  {
    std::vector<int> subind = {AD_ind_[C], V0_ind_[A]};
    printf("[ADC V0C]: \n");
    Subspace(xx, subind, scale);
  }
  printf("Subspace of %s: \n", nametag.c_str());
  {
    std::vector<int> subind = {V0_ind_[C], AD_ind_[A]};
    printf("[V0A ADA]: \n");
    Subspace(xx, subind, scale);
  }
  printf("Subspace of %s: \n", nametag.c_str());
  {
    std::vector<int> subind = {V0_ind_[C], V0_ind_[A]};
    printf("[V0C V0A]: \n");
    Subspace(xx, subind, scale);
  }
  printf("Subspace of %s: \n", nametag.c_str());
  {
    std::vector<int> subind = {AD_ind_[C], AD_ind_[A]};
    printf("[ADC ADA]: \n");
    Subspace(xx, subind, scale);
  }
  printf("Subspace of %s: \n", nametag.c_str());
  {
    std::vector<int> subind = {SPD_ind_[C], SPD_ind_[A]};
    printf("[SPDC SPDA]: \n");
    Subspace(xx, subind, scale);
  }  
  printf("Subspace of %s: \n", nametag.c_str());
  {
    std::vector<int> subind = {SPD_ind_[C]};
    printf("[SPDC]: \n");
    Subspace(xx, subind, scale);
  }
  printf("Subspace of %s: \n", nametag.c_str());
  {
    std::vector<int> subind = {SPD_ind_[A]};
    printf("[SPDA]: \n");
    Subspace(xx, subind, scale);
  }
  printf("Subspace of %s: \n", nametag.c_str());
  {
    std::vector<int> subind = {AD_ind_[C]};
    printf("[ADC]: \n");
    Subspace(xx, subind, scale);
  }
  printf("Subspace of %s: \n", nametag.c_str());
  {
    std::vector<int> subind = {V0_ind_[C]};
    printf("[V0C]: \n");
    Subspace(xx, subind, scale);
  }
  printf("Subspace of %s: \n", nametag.c_str());
  {
    std::vector<int> subind = {V0_ind_[A]};
    printf("[V0A]: \n");
    Subspace(xx, subind, scale);
  }
  printf("Subspace of %s: \n", nametag.c_str());
  {
    std::vector<int> subind = {AD_ind_[A]};
    printf("[ADA]: \n");
    Subspace(xx, subind, scale);
  }    

  printf("\n");


}



// DEBUG output printing
void Combinatorics::PrintData(const CombEvent* event) {

	printf("Combinatorics::PrintData:: \n");

	printf("\n---------------------------------------------------------\n");

	ULong64_t 	   c = event->fTreeData()->fEventInfo.fClassMask;
	ULong64_t      d = event->fTreeData()->fEventInfo.fClassMaskNext50;
	UInt_t         e = event->fTreeData()->fEventInfo.fBCID;
	UInt_t         f = event->fTreeData()->fEventInfo.fPeriod;
	UInt_t         g = event->fTreeData()->fEventInfo.fTimeStamp;
	UInt_t         h = event->fTreeData()->fEventInfo.fL0Inputs;
	UInt_t         i2 = event->fTreeData()->fEventInfo.fL1Inputs;

	Int_t          j = event->fTreeData()->fEventInfo.fRunNumber;
	UShort_t       k2 = event->fTreeData()->fEventInfo.fnTrklet;
	UShort_t       l = event->fTreeData()->fEventInfo.fL2Inputs;
	UShort_t       m = event->fTreeData()->fEventInfo.fOrbitID;

	printf("Binary codes:\n");

	printf(".fL0inputs:\n");
	PrintBits(sizeof(h), &h);
	printf(".fL1inputs:\n");
	PrintBits(sizeof(i2), &i2);
	printf(".fClassMask:\n");
	PrintBits(sizeof(c), &c);
	printf(".fClassMaskNext50:\n");
	PrintBits(sizeof(d), &d);
	printf("\n");

	// --------------------------------
	printf("Generic:\n");
	printf(".fClassMask: %llu \n.fClassMaskNext50: %llu \n.fBCID: %u \n.fPeriod: %u \n.fTimeStamp: %u \n.fL0Inputs: %u \n.fL1Inputs: %u \n.fRunNumber: %d \n.fnTrklet: %u \n.fL2Inputs: %u \n.fOrbitID: %u\n", 
		   c,d,e,f,g,h,i2,j,k2,l,m);

	UInt_t           n  = event->fTreeData()->fPhysSelBits;
	Bool_t           o  = event->fTreeData()->fIsIncompleteDAQ;
	Bool_t           p  = event->fTreeData()->fIsSPDClusterVsTrackletBG;

	printf(".fPhysSelBits: %u \n.fIsIncompleteDAQ: %d \n.fIsSPDClusterVsTrackletBG: %d\n\n", n, o, p);

	printf("\n");

	// --------------------------------
	for (UInt_t k = 0; k < 2; ++k) {

		Float_t 	 zz = event->fTreeData()->fZDCInfo.fZPEnergy[k];
		Float_t 	 yy = event->fTreeData()->fZDCInfo.fZNEnergy[k];

		Float_t      bb = event->fTreeData()->fV0Info.fTime[k];
		Float_t      cc = event->fTreeData()->fV0Info.fCharge[k];

		Char_t       dd = event->fTreeData()->fV0Info.fBB[k];
		Char_t       ee = event->fTreeData()->fV0Info.fBG[k];

		Int_t        ff = TMath::Nint(event->fTreeData()->fV0Info.fDecisionOnline[k]);
		Int_t     	 gg = TMath::Nint(event->fTreeData()->fV0Info.fDecisionOffline[k]);

		Float_t      hh = event->fTreeData()->fADInfo.fTime[k];
		Float_t      ii = event->fTreeData()->fADInfo.fCharge[k];

		Char_t       jj = event->fTreeData()->fADInfo.fBB[k];
		Char_t       kk = event->fTreeData()->fADInfo.fBG[k];

		Int_t     	 ll = TMath::Nint(event->fTreeData()->fADInfo.fDecisionOnline[k]);
		Int_t    	 mm = TMath::Nint(event->fTreeData()->fADInfo.fDecisionOffline[k]);

		if (k == 0) {
			printf("C-side: \n");
		}
		if (k == 1) {
			printf("A-side: \n");
		}
		printf(".fZDCInfo.fZPEnergy: %0.2f \n.fZDCInfo.fZNEnergy: %0.2f \n.fV0Info.fTime: %0.2f \n.fV0Info.fCharge: %0.2f \n.fV0Info.fBB: %d \n.fV0Info.fBG: %d \n.fV0Info.fDecisionOnline: %d \n.fV0Info.fDecisionOffline: %d \n.fADInfo.fTime: %0.2f \n.fADInfo.fCharge: %0.2f \n.fADInfo.fBB: %d \n.fADInfo.fBG: %d \n.fADInfo.fDecisionOnline: %d \n.fADInfo.fDecisionOffline: %d \n\n", 
			   zz, yy, bb,cc,dd,ee,ff,gg,hh,ii,jj,kk,ll,mm);
	}
	printf("\n");

	// --------------------------------

	for (UInt_t k = 0; k < 2; ++k) {

		UInt_t aa = event->fTreeData()->fEventInfo.fnSPDClusters[k];

		if (k == 0) {
			printf("SPD Layer 0: \n");
		}
		if (k == 1) {
			printf("SPD Layer 1: \n");
		}
		printf(".fnSPDClusters: %d \n", aa);
	}

	UInt_t bFastOR = event->fFastOrMap()->CountBits();

	printf("\nSPD:\n");
	printf(".fFastOrMap.CountBits(): %d \n\n", bFastOR);


	UInt_t bChip   = event->fFiredChipMap()->CountBits();

	printf("\nSPD:\n");
	printf(".fFiredChipMap.CountBits(): %d \n\n", bChip);

	// MC
	if (fIsMC_) {
		printf("MC variables: \n");

		Int_t ee       = event->fMCInfo()->fEventType;
		Float_t mass_0 = event->fMCInfo()->fDiffSys.Mass[0];
		Float_t mass_1 = event->fMCInfo()->fDiffSys.Mass[1];

		printf(".fMCInfo.ev.proc: %d \n.fMCInfo.fDiffSys.Mass[0]: %0.2f [GeV] \n.fMCInfo.fDiffSys.Mass[1]: %0.2f [GeV] \n\n", ee, mass_0, mass_1);
	}
}


// Calculate F-matrix (class density matrix) (ONLY MC (other than EPOS) USE THIS ROUTINE)
void Combinatorics::CalculateF() {

	if (fIsMC_ && !fMCName_.Contains("EPOS-LHC")) {

		printf("Combinatorics::CalculateF:: \n");

		printf("\n%s :: \n", fMCName_.Data());
		printf(" Generator Default Process Fractions and Cross Sections (XI_MAX cut always included): \n");
		Double_t norm = vsum(fMCTotalProcessCount_);
		for (UInt_t j = 0; j < C_; ++j) {

			printf("  %d | %0.3f (%0.2f mb)\n", j, 
				fMCTotalProcessCount_.at(j) / norm, fMCTotalProcessCount_.at(j) / norm * sigma_inel_tot_);
		}
		printf("TOT | 1.000 (%0.2f mb) \n", sigma_inel_tot_);
		printf("\n");

		// ---------------------------------------------------------------
		// ** GENERATOR level **

		// Normalize the columns to sum to one (0-bin included, because that keeps the ACCEPTANCES)
		// THIS TAKES INTO ACCOUNT THE MASS DISTRIBUTION RE-WEIGHTING
		for (UInt_t j = 0; j < C_; ++j) {
			Double_t sum = EPS;
			for (UInt_t i = 0; i < std::pow(2,d_); ++i) {
				sum += F_gen_.at(i).at(j);
			}
			for (UInt_t i = 0; i < std::pow(2,d_); ++i) {
				F_gen_.at(i).at(j) /= sum;
			}
		}
		// Generate synthetic x-vector (pure counts) - for comparisons with data
		for (UInt_t i = 0; i < x_.size(); ++i) {
			Double_t sum = 0;

			for (UInt_t j = 0; j < C_; ++j) {
				sum += fMCTotalProcessCount_.at(j) * F_gen_.at(i).at(j); // TotalProcessCount is the normalization
			}
			x_gen_.at(i) = TMath::Nint(sum); // Note, must be integers
		}

		// ---------------------------------------------------------------
		// ** DETECTOR level **

		// Normalize the columns to sum to one (0-bin included, because that keeps the EFFICIENCY x ACCEPTANCES)
		// THIS TAKES INTO ACCOUNT THE MASS DISTRIBUTION RE-WEIGHTING
		for (UInt_t j = 0; j < C_; ++j) {
			Double_t sum = EPS;
			for (UInt_t i = 0; i < std::pow(2,d_); ++i) {
				sum += F_.at(i).at(j);
			}
			for (UInt_t i = 0; i < std::pow(2,d_); ++i) {
				F_.at(i).at(j) /= sum;
			}
		}
		// Generate synthetic x-vector (counts) - for comparisons with data
		for (UInt_t i = 0; i < x_.size(); ++i) {
			Double_t sum = 0;

			for (UInt_t j = 0; j < C_; ++j) {
				sum += fMCTotalProcessCount_.at(j) * F_.at(i).at(j);  // TotalProcessCount is the normalization
			}
			x_.at(i) = TMath::Nint(sum); // Note, must be integers
		}

	} else {
		printf("Combinatorics::CalculateF() method is not for for DATA or EPOS! \n");
	}

}


// Plot 1D histograms
void Combinatorics::Plot1D(std::vector<TCanvas*>& can, UInt_t color, UInt_t marker) {

	if (HISTOGRAMS_ON) {

		// Schema d_ = 6 and d_ = 8 compatible
		Int_t ZDN_plot_ind[2] = {1,8};
		Int_t ZDP_plot_ind[2] = {9,16};

		Int_t AD_charge_plot_ind[2] = {2,7};
		Int_t V0_charge_plot_ind[2] = {3,6};

		Int_t AD_time_plot_ind[2] = {10,15};
		Int_t V0_time_plot_ind[2] = {11,14};

		Int_t SPD_FO_plot_ind[2]  = {4,5};


		// Loop over combinations
		for (UInt_t c = 0; c < std::pow(2,d_); ++c) {

			// Get corresponding vector
			std::vector<Bool_t> boolvec = Ind2Vec(c, d_);

			// Loop over C and A-side
			for (UInt_t k = 0; k < 2; ++k) {

				// Upper row
				//can[c]->cd(ZDN_plot_ind[k])->SetLogx();
				can[c]->cd(ZDN_plot_ind[k])->SetLogy();

				hZDN[c][k]->SetLineColor(color);
				NormHist(hZDN[c][k]);
				hZDN[c][k]->Draw("SAME");

				can[c]->cd(AD_charge_plot_ind[k])->SetLogx();
				can[c]->cd(AD_charge_plot_ind[k])->SetLogy();

				hADCharge[c][k]->SetLineColor(color);
				NormHist(hADCharge[c][k]);
				hADCharge[c][k]->Draw("SAME");

				can[c]->cd(V0_charge_plot_ind[k])->SetLogx();
				can[c]->cd(V0_charge_plot_ind[k])->SetLogy();

				hV0Charge[c][k]->SetLineColor(color);
				NormHist(hV0Charge[c][k]);
				hV0Charge[c][k]->Draw("SAME");

				can[c]->cd(SPD_FO_plot_ind[k])->SetLogx();
				can[c]->cd(SPD_FO_plot_ind[k])->SetLogy();

				hSPDFO[c][k]->SetLineColor(color);
				NormHist(hSPDFO[c][k]);
				hSPDFO[c][k]->Draw("SAME");


				// Lower row
				//can[c]->cd(ZDN_plot_ind[k])->SetLogx();
				can[c]->cd(ZDP_plot_ind[k])->SetLogy();

				hZDP[c][k]->SetLineColor(color);
				NormHist(hZDP[c][k]);
				hZDP[c][k]->Draw("SAME");

				can[c]->cd(AD_time_plot_ind[k])->SetLogy();

				hADTime[c][k]->SetLineColor(color);
				NormHist(hADTime[c][k]);
				hADTime[c][k]->Draw("SAME");

				can[c]->cd(V0_time_plot_ind[k])->SetLogy();

				hV0Time[c][k]->SetLineColor(color);
				NormHist(hV0Time[c][k]);
				hV0Time[c][k]->Draw("SAME");
			}


			//can[c]->cd(SPD_ind_+1)->SetLogx();
			//can[c]->cd(SPD_ind_+1)->SetLogy();

			//hSPDFO[c]->SetLineColor(color);
			//NormHist(hSPDFO[c]);
			//hSPDFO[c]->Draw("SAME");

			// Lower row SPD
			/*can[c]->cd(SPD_plot_ind2)->SetLogx();
			can[c]->cd(SPD_plot_ind2)->SetLogy();

			hSPDTR[c]->SetLineColor(color);
			NormHist(hSPDTR[c]);
			hSPDTR[c]->Draw("SAME");
			*/

			// Fill backgrounds for visualization
			for (Int_t k = 1; k < d_+1; ++k) {

				Double_t bg_color = (boolvec.at(k-1) == 1) ? kWhite : kGray;

				// Schema d = 6 compatible + ZDNs on sides
				if (d_ == 6) {
					can[c]->cd(k+1);      gPad->SetFillColor(bg_color); // Upper row
					can[c]->cd(k+1+d_+2); gPad->SetFillColor(bg_color); // Lower row
				}
				// Schema d = 8
				if (d_ == 8) {
					can[c]->cd(k);        gPad->SetFillColor(bg_color); // Upper row
					can[c]->cd(k+d_); 	  gPad->SetFillColor(bg_color); // Lower row
				}
			}

			// ZDC background color
			can[c]->cd(ZDN_plot_ind[0]); gPad->SetFillColor(kRed-10);
			can[c]->cd(ZDN_plot_ind[1]); gPad->SetFillColor(kRed-10);

			can[c]->cd(ZDP_plot_ind[0]); gPad->SetFillColor(kRed-10);
			can[c]->cd(ZDP_plot_ind[1]); gPad->SetFillColor(kRed-10);
		}
	}

}

// MC <observable> versus mass
void Combinatorics::Plot1DA(std::vector<TCanvas*>& can, UInt_t color, UInt_t marker) {

	if (fIsMC_ && HISTOGRAMS_ON) {

		Int_t kk = 1;
		for (Int_t i = 0; i < 2; ++i) {
			for (Int_t j = 0; j < 2; ++j) {
						
				can[0]->cd(kk);  // can[0]->cd(kk)->SetLogx();
				hAD_M_avgCharge[i][j]->SetLineColor(color);
				hAD_M_avgCharge[i][j]->Draw("SAME");  hAD_M_avgCharge[i][j]->GetYaxis()->SetRangeUser(0.0, 1000.0);
				
				can[1]->cd(kk); // can[1]->cd(kk)->SetLogx();
				hV0_M_avgCharge[i][j]->SetLineColor(color);
				hV0_M_avgCharge[i][j]->Draw("SAME");  hV0_M_avgCharge[i][j]->GetYaxis()->SetRangeUser(0.0, 50.0);
				
				can[2]->cd(kk); // can[2]->cd(kk)->SetLogx();
				hSPD_M_avgFO[i][j]->SetLineColor(color);
				hSPD_M_avgFO[i][j]->Draw("SAME");     hSPD_M_avgFO[i][j]->GetYaxis()->SetRangeUser(0.0, 20.0);

				can[3]->cd(kk); // can[3]->cd(kk)->SetLogx();
				hZDN_M_avgCharge[i][j]->SetLineColor(color);
				hZDN_M_avgCharge[i][j]->Draw("SAME"); hZDN_M_avgCharge[i][j]->GetYaxis()->SetRangeUser(0.0, 1400.0);
				
				can[4]->cd(kk); // can[4]->cd(kk)->SetLogx();
				hZDP_M_avgCharge[i][j]->SetLineColor(color);
				hZDP_M_avgCharge[i][j]->Draw("SAME"); hZDP_M_avgCharge[i][j]->GetYaxis()->SetRangeUser(0.0, 400.0);

				++kk;
			}
		}
	}
}

// MC <observable> versus fiducial tracks
void Combinatorics::Plot1DB(std::vector<TCanvas*>& can, UInt_t color, UInt_t marker) {

	if (fIsMC_ && HISTOGRAMS_ON) {

		Int_t kk = 1;
		for (Int_t i = 0; i < 2; ++i) {
			for (Int_t j = 0; j < 2; ++j) {
						
				// AD
				can[0]->cd(kk);  // can[0]->cd(kk)->SetLogx();
				hAD_N_avgCharge[i][j]->SetLineColor(color);
				hAD_N_avgCharge[i][j]->Draw("SAME");   hAD_N_avgCharge[i][j]->GetYaxis()->SetRangeUser(0.0, 1000.0);
				can[0]->cd(kk+1);  // can[0]->cd(kk)->SetLogx();
				hAD_Nch_avgCharge[i][j]->SetLineColor(color);
				hAD_Nch_avgCharge[i][j]->Draw("SAME"); hAD_Nch_avgCharge[i][j]->GetYaxis()->SetRangeUser(0.0, 1000.0);
				can[0]->cd(kk+2);  // can[0]->cd(kk)->SetLogx();
				hAD_Nn_avgCharge[i][j]->SetLineColor(color);
				hAD_Nn_avgCharge[i][j]->Draw("SAME");  hAD_Nn_avgCharge[i][j]->GetYaxis()->SetRangeUser(0.0, 1000.0);
				

				// V0
				can[1]->cd(kk); // can[1]->cd(kk)->SetLogx();
				hV0_N_avgCharge[i][j]->SetLineColor(color);
				hV0_N_avgCharge[i][j]->Draw("SAME");   hV0_N_avgCharge[i][j]->GetYaxis()->SetRangeUser(0.0, 50.0);
				can[1]->cd(kk+1); // can[1]->cd(kk)->SetLogx();
				hV0_Nch_avgCharge[i][j]->SetLineColor(color);
				hV0_Nch_avgCharge[i][j]->Draw("SAME"); hV0_Nch_avgCharge[i][j]->GetYaxis()->SetRangeUser(0.0, 50.0);
				can[1]->cd(kk+2); // can[1]->cd(kk)->SetLogx();
				hV0_Nn_avgCharge[i][j]->SetLineColor(color);
				hV0_Nn_avgCharge[i][j]->Draw("SAME");  hV0_Nn_avgCharge[i][j]->GetYaxis()->SetRangeUser(0.0, 50.0);
					

				// SPD
				can[2]->cd(kk); // can[2]->cd(kk)->SetLogx();
				hSPD_N_avgFO[i][j]->SetLineColor(color);
				hSPD_N_avgFO[i][j]->Draw("SAME");      hSPD_N_avgFO[i][j]->GetYaxis()->SetRangeUser(0.0, 20.0);
				can[2]->cd(kk+1); // can[2]->cd(kk)->SetLogx();
				hSPD_Nch_avgFO[i][j]->SetLineColor(color);
				hSPD_Nch_avgFO[i][j]->Draw("SAME");    hSPD_Nch_avgFO[i][j]->GetYaxis()->SetRangeUser(0.0, 20.0);
				can[2]->cd(kk+2); // can[2]->cd(kk)->SetLogx();
				hSPD_Nn_avgFO[i][j]->SetLineColor(color);
				hSPD_Nn_avgFO[i][j]->Draw("SAME");     hSPD_Nn_avgFO[i][j]->GetYaxis()->SetRangeUser(0.0, 20.0);
				

				// ZDN
				can[3]->cd(kk); // can[3]->cd(kk)->SetLogx();
				hZDN_Nn_avgCharge[i][j]->SetLineColor(color);
				hZDN_Nn_avgCharge[i][j]->Draw("SAME");  hZDN_Nn_avgCharge[i][j]->GetYaxis()->SetRangeUser(0.0, 1400.0);
				
				kk += 3;
			}
		}
	}
}



// Plotting
void Combinatorics::Plot() {

	// THIS IS ALWAYS PLOTTED
	// Folding matrix
	if (fIsMC_) {

		TCanvas* ccc = new TCanvas(Form("cM%s",fMCName_.Data()), "cM", 400, 400);
		ccc->cd();
		gPad->SetLogz(1);

		// Set X-axis tick divisions
		const Int_t n1 = 14;  // Primary divisions
		const Int_t n2 = 5;   // Secondary divisions
		const Int_t n = n1 + 100*n2;

		M->GetXaxis()->SetNdivisions(n);
		M->GetYaxis()->SetNdivisions(n);

		M->Draw("COLZ");
		ccc->Print(Form("./figures_xsec/%d/Unfolding/FoldingMatrix_%s.pdf", fRunNumber_, fMCName_.Data()));
		delete ccc;
	}
	


	if (HISTOGRAMS_ON) {

		printf("Combinatorics::Plot:: \n");

	   	// The number of elements in the vector space
	   	const UInt_t ncomb = (UInt_t) std::pow(2,d_);

	   	for (UInt_t c = 0; c < ncomb; ++c) {
			TCanvas* cSPD = new TCanvas(Form("cSPD_ID%d_%s",c,fMCName_.Data()), Form("hSPD_%d",c), 400, 300);
			cSPD->cd();
			cSPD->SetLogy();
			hSPDbit[c]->Draw();

			cSPD->Print(Form("./figures_xsec/%d/SPD/%s/SPD_ID%03u.pdf", fRunNumber_, fMCName_.Data(), c));
			delete cSPD;
		}


	   	// 2D-cross correlations

	   	// -------------------------------------------------------------------
	   	// Create canvases

	   	std::vector<TCanvas*> c2DXC(ncomb, NULL);

	   	for (UInt_t c = 0; c < ncomb; ++c) {
			c2DXC[c] = new TCanvas(Form("c2DXC_ID%d_%s",c,fMCName_.Data()), Form("h2Det_%d",c), 1200, 1200);
			c2DXC[c]->Divide(d_,d_); 			  // (horizontal, vertical) boxes
		}

		// 2D histograms
		// Loop over combinations
		for (UInt_t c = 0; c < ncomb; ++c) {
			UInt_t k = 1;

			// Get corresponding vector
			std::vector<Bool_t> boolvec = Ind2Vec(c, d_);

			for (Int_t i = 0; i < d_; ++i) {
				for (Int_t j = 0; j < d_; ++j) {

					if (j >= i) { // Plot only diagonal and upper right diagonal

						c2DXC[c]->cd(k);//->SetLogx();
						
						// Coloring
							Double_t color = (boolvec.at(i) == 1 && boolvec.at(j) == 1) ? kWhite : kGray;
							gPad->SetFillColor(color);

						// Pearson Cross Correlation
						Double_t xcorr = h2XC[c][i][j]->GetCorrelationFactor();
						h2XC[c][i][j]->SetTitle(Form("ID: %d, #rho = %0.2f", c, xcorr));					
	   					h2XC[c][i][j]->Draw("COLZ");
					}
					++k;
				}
			}

			c2DXC[c]->Print(Form("./figures_xsec/%d/XC/%s/Vector_ID%03u.pdf", fRunNumber_, fMCName_.Data(), c));
		}

	   	// -------------------------------------------------------------------
	   	// Create canvases

	   	std::vector<TCanvas*> c2D(ncomb, NULL);

	   	for (UInt_t c = 0; c < ncomb; ++c) {
			c2D[c] = new TCanvas(Form("c2D_ID%d_%s",c,fMCName_.Data()), Form("h2Det_%d",c), 1200, 200);
			c2D[c]->Divide(d_,1); 			  // (horizontal, vertical) boxes
		}

		// 2D histograms
		// Loop over combinations
		for (UInt_t c = 0; c < ncomb; ++c) {

			// Get corresponding vector
			std::vector<Bool_t> boolvec = Ind2Vec(c, d_);

			// Loop over C and A-side
			for (UInt_t k = 0; k < 2; ++k) {

				c2D[c]->cd(AD_ind_[k]+1)->SetLogy();
				h2ADCT[c][k]->Draw("COLZ");

				c2D[c]->cd(V0_ind_[k]+1)->SetLogy();
				h2V0CT[c][k]->Draw("COLZ");
				
				c2D[c]->cd(SPD_ind_[k]+1);//->SetLogx();
				h2SPDFOTR[c][k]->Draw("COLZ");
			}


			// Fill backgrounds for visualization
			for (Int_t k = 0; k < d_; ++k) {

				Double_t color = boolvec.at(k) == 1 ? kWhite : kGray;
				c2D[c]->cd(k+1); gPad->SetFillColor(color);

			}
			// Save PDFs
			c2D[c]->Print(Form("./figures_xsec/%d/Detector_2D/%s/Vector_ID%03u.pdf", fRunNumber_, fMCName_.Data(), c));
		}

/*
		// 2D-MC <observable> versus mass
		if (fIsMC_) {
			TCanvas* cxAD  = new TCanvas(Form("cxAD_%s",  fMCName_.Data()), "h2-AD",  500, 400); cxAD->Divide(2,2);
			TCanvas* cxV0  = new TCanvas(Form("cxV0_%s",  fMCName_.Data()), "h2-V0",  500, 400); cxV0->Divide(2,2);
			TCanvas* cxSPD = new TCanvas(Form("cxSPD_%s", fMCName_.Data()), "h2-SPD", 500, 400); cxSPD->Divide(2,2);
			TCanvas* cxZDN = new TCanvas(Form("cxZDN_%s", fMCName_.Data()), "h2-ZDN", 500, 400); cxZDN->Divide(2,2);
			TCanvas* cxZDP = new TCanvas(Form("cxZDP_%s", fMCName_.Data()), "h2-ZDP", 500, 400); cxZDP->Divide(2,2);

			Int_t kk = 1;
			for (Int_t i = 0; i < 2; ++i) {
				for (Int_t j = 0; j < 2; ++j) {
					
					cxAD->cd(kk);  cxAD->cd(kk)->SetLogx();  cxAD->cd(kk)->SetLogy();
					h2AD_M_Charge[i][j]->Draw();
					
					cxV0->cd(kk);  cxV0->cd(kk)->SetLogx();  cxV0->cd(kk)->SetLogy();
					h2V0_M_Charge[i][j]->Draw();
					
					cxSPD->cd(kk); cxSPD->cd(kk)->SetLogx(); cxSPD->cd(kk)->SetLogy();
					h2SPD_M_FO[i][j]->Draw();
					
					cxZDN->cd(kk); cxZDN->cd(kk)->SetLogx(); cxZDN->cd(kk)->SetLogy();
					h2ZDN_M_Charge[i][j]->Draw();

					cxZDP->cd(kk); cxZDP->cd(kk)->SetLogx(); cxZDP->cd(kk)->SetLogy();
					h2ZDP_M_Charge[i][j]->Draw();

					++kk;
				}
			}
			cxAD->Print(Form("./figures_xsec/%d/Detector/cx2D_AD_%s.pdf",   fRunNumber_, fMCName_.Data() ));
			cxV0->Print(Form("./figures_xsec/%d/Detector/cx2D_V0_%s.pdf",   fRunNumber_, fMCName_.Data() ));
			cxSPD->Print(Form("./figures_xsec/%d/Detector/cx2D_SPD_%s.pdf", fRunNumber_, fMCName_.Data() ));
			cxZDN->Print(Form("./figures_xsec/%d/Detector/cx2D_ZDN_%s.pdf", fRunNumber_, fMCName_.Data() ));
			cxZDP->Print(Form("./figures_xsec/%d/Detector/cx2D_ZDP_%s.pdf", fRunNumber_, fMCName_.Data() ));


			delete cxAD;
			delete cxV0;
			delete cxSPD;
			delete cxZDN;
			delete cxZDP;	
		}
*/

		TCanvas* c5 = new TCanvas(Form("c5_%s",fMCName_.Data()), "DiffMassDD", 230, 500);
		c5->Divide(1,2);


		TCanvas* c7 = new TCanvas(Form("c7_%s",fMCName_.Data()), "DiffDeltaYDD", 230, 500);
		c7->Divide(1,2);


		// FIT object for fitting the MC invariant mass squared histograms
		// FIT dsigma/dm^2 ~ 1/(m^2)^(1+POMERON_DELTA), here m^2 = x

		// HIGH_MASS fit
		TF1* M2fit = new TF1("M2fit",Reggefit, 10, 10000, 2); // 2 parameters
		M2fit->SetParameters(1, 1); // Initial values
		M2fit->SetParName(0,"C0");
		M2fit->SetParName(1,"POMERON_DELTA");

		// LOW-MASS fit
		TF1* M2Lowfit = new TF1("M2Lowfit",Reggefit, 2.5, 25, 2); // 2 parameters
		M2Lowfit->SetParameters(1, 1); // Initial values
		M2Lowfit->SetParName(0,"C0");
		M2Lowfit->SetParName(1,"POMERON_DELTA");

		// MC only, Mass histograms
		if (fIsMC_ && !fMCName_.Contains("EPOS-LHC")) {

			TCanvas* ccx = new TCanvas(Form("ccx_%s",fMCName_.Data()), "SDMass", 500, 280);
			ccx->Divide(2,1);
			TLegend legW2[2] = TLegend(0.35,0.73, 0.82,0.87);

			for (UInt_t kk = 0; kk < 2; ++kk) {
				ccx->cd(kk+1);

				hDiffMassSDLowM_gene[kk]->Fit("M2Lowfit","ER"); // "E" minos errors, R" restricts the fit
				hDiffMassSDLowM_gene[kk]->Draw();

				// legend (x1,y1, x2, y2) with FIT values
				legW2[kk].AddEntry(M2Lowfit, Form("C x 1/(M^{2})^{1+#Delta}, #Delta = %0.2f #pm %0.2f",
					M2Lowfit->GetParameter(1), M2Lowfit->GetParError(1)), "l");
				legW2[kk].SetBorderSize(0);
				legW2[kk].SetTextSize(0.03);
				legW2[kk].Draw();
			}

			ccx->SaveAs(Form("./figures_xsec/%d/Acceptance/SDLowM_%s.pdf", fRunNumber_, fMCName_.Data()) );
			delete ccx;


			TCanvas* ccx2 = new TCanvas(Form("ccx2_%s",fMCName_.Data()), "SDMass", 500, 280);
			ccx2->Divide(2,1);
			TLegend legW22[2] = TLegend(0.45,0.73, 0.82,0.87);

			for (UInt_t kk = 0; kk < 2; ++kk) {
				ccx2->cd(kk+1);
				ccx2->cd(kk+1)->SetLogy();
				ccx2->cd(kk+1)->SetLogx();

				hDiffMassSD_gene[kk]->Fit("M2fit","ER"); // "E" minos errors, R" restricts the fit
				hDiffMassSD_gene[kk]->Draw();

				// legend (x1,y1, x2, y2) with FIT values
				legW22[kk].AddEntry(M2fit, Form("C x 1/(M^{2})^{1+#Delta}, #Delta = %0.2f #pm %0.2f",
					M2fit->GetParameter(1), M2fit->GetParError(1)), "l");
				legW22[kk].SetBorderSize(0);
				legW22[kk].SetTextSize(0.03);
				legW22[kk].Draw();
			}

			ccx2->SaveAs(Form("./figures_xsec/%d/Acceptance/SDHighM_%s.pdf", fRunNumber_, fMCName_.Data()) );
			delete ccx2;


			// ----------------------------------------------------------	
			// Efficiency simulations >>

			for (Int_t side = 0; side < 2; ++side) {
				// Fit
   				hDiffMassSD_gene[side]->Fit("M2fit","ER"); // "E" minos errors, R" restricts the fit
			}


			// Loop over detectors >>

			// Calculate efficiencies
			for (Int_t i = 0; i < d_; ++i) {

				// Loop over C and A sides
				for (Int_t side = 0; side < 2; ++side) {

					// Calculate the efficiency by dividing with generated
					hDiffMassSD_seen_GEN[i][side]->Divide(hDiffMassSD_gene[side]);
					hDiffMassSD_seen_DET[i][side]->Divide(hDiffMassSD_gene[side]);

					hDiffMassSD_seen_GEN[i][side]->SetMaximum(1.0);
					hDiffMassSD_seen_DET[i][side]->SetMaximum(1.0);
					hDiffMassSD_seen_GEN[i][side]->SetMinimum(0.0);
					hDiffMassSD_seen_DET[i][side]->SetMinimum(0.0);

					// Calculate the efficiency by dividing with generated
					hDiffDeltaYSD_seen_GEN[i][side]->Divide(hDiffDeltaYSD_gene[side]);
					hDiffDeltaYSD_seen_DET[i][side]->Divide(hDiffDeltaYSD_gene[side]);

					hDiffDeltaYSD_seen_GEN[i][side]->SetMaximum(1.0);
					hDiffDeltaYSD_seen_DET[i][side]->SetMaximum(1.0);
					hDiffDeltaYSD_seen_GEN[i][side]->SetMinimum(0.0);
					hDiffDeltaYSD_seen_DET[i][side]->SetMinimum(0.0);
				}
			}


			// Loop over detectors
			TCanvas* c4 = new TCanvas(Form("c4_%s",fMCName_.Data()), "DiffMassSD", 300, 800);
			c4->Divide(2,d_);
			Int_t kk = 1;
			for (Int_t i = 0; i < d_; ++i) {

				// Loop over SDL and SDR
				for (Int_t side = 0; side < 2; ++side) {

					c4->cd(kk)->SetLogx(); 
					hDiffMassSD_seen_GEN[i][side]->Draw("SAME");
					hDiffMassSD_seen_DET[i][side]->Draw("SAME");
					++kk;
				}
			}
			c4->SaveAs(Form("./figures_xsec/%d/Acceptance/SD_%s_Acc_Mass_Vertical.pdf", fRunNumber_, fMCName_.Data() ));
			delete c4;

			// Loop over SDL and SDR
			TCanvas* c44 = new TCanvas(Form("c44_%s",fMCName_.Data()), "DiffMassSD", 1000, 300);
			c44->Divide(d_,2);
			kk = 1;
			for (Int_t side = 0; side < 2; ++side) {

				// Loop over detectors
				for (Int_t i = 0; i < d_; ++i) {

					c44->cd(kk)->SetLogx(); 
					hDiffMassSD_seen_GEN[i][side]->Draw("SAME");
					hDiffMassSD_seen_DET[i][side]->Draw("SAME");
					++kk;
				}
			}
			c44->SaveAs(Form("./figures_xsec/%d/Acceptance/SD_%s_Acc_Mass_Horizontal.pdf", fRunNumber_, fMCName_.Data() ));
			delete c44;

			/*
			// Regge-Mass distributions
			c4->cd(side+3)->SetLogy();
			c4->cd(side+3)->SetLogx();
			hDiffMassSD_gene[side]->Draw();	
	
			// legend (x1,y1, x2, y2) with FIT values
			TLegend legW = TLegend(0.5,0.73, 0.82,0.87);
			legW.AddEntry(M2fit, Form("C x 1/(M^{2})^{1+#epsilon}, #epsilon = %0.2f #pm %0.2f",
				M2fit->GetParameter(1), M2fit->GetParError(1)), "l");
			legW.SetBorderSize(0);
			legW.SetTextSize(0.03);
			legW.Draw();
			*/

			// Calculate the efficiency by dividing with the generated
			/*
			h2DiffMassDD_seen_GEN[i]->Divide(h2DiffMassDD_gene);
			h2DiffMassDD_seen_DET[i]->Divide(h2DiffMassDD_gene);

			for (Int_t side = 0; side < 2; ++side) {
				c5->cd(1)->SetLogy();
				c5->cd(1)->SetLogx();
				h2DiffMassDD_seen_GEN[i]->Draw("COLZ");
			}
			
			c5->cd(2)->SetLogy(); c5->cd(2)->SetLogx();
			h2DiffMassDD_gene->Draw("COLZ");

			c5->SaveAs(Form("./figures_xsec/%d/Acceptance/DD_%s_%s_Acc_Mass.pdf", fRunNumber_, fMCName_.Data(), det_labels_[i].Data() ));
			*/

			// ***********************************************************

			// Loop over detectors >>
			TCanvas* c6 = new TCanvas(Form("c6_%s",fMCName_.Data()), "DiffDeltaYSD", 300, 800);
			c6->Divide(2,d_);
			kk = 1;
			for (Int_t i = 0; i < d_; ++i) {

				// Loop over SDL and SDR
				for (Int_t side = 0; side < 2; ++side) {

					c6->cd(kk);
					hDiffDeltaYSD_seen_GEN[i][side]->Draw("SAME");
					hDiffDeltaYSD_seen_DET[i][side]->Draw("SAME");
					++kk;
				}
			}
			c6->SaveAs(Form("./figures_xsec/%d/Acceptance/SD_%s_Acc_DeltaY_Vertical.pdf", fRunNumber_, fMCName_.Data()) );
			delete c6;

			// Loop over SDL and SDR
			TCanvas* c66 = new TCanvas(Form("c66_%s",fMCName_.Data()), "DiffDeltaYSD", 1000, 300);
			c66->Divide(d_,2);
			kk = 1;
			for (Int_t side = 0; side < 2; ++side) {

				// Loop over detectors
				for (Int_t i = 0; i < d_; ++i) {

					c66->cd(kk);
					hDiffDeltaYSD_seen_GEN[i][side]->Draw("SAME");
					hDiffDeltaYSD_seen_DET[i][side]->Draw("SAME");
					++kk;
				}
			}
			c66->SaveAs(Form("./figures_xsec/%d/Acceptance/SD_%s_Acc_DeltaY_Horizontal.pdf", fRunNumber_, fMCName_.Data()) );
			delete c66;

			/*
			// Generated
			c6->cd(side+3);
			hDiffDeltaYSD_gene[side]->Draw();
			*/

			/*
			h1DiffDeltaYDD_seen_GEN[i]->Divide(h1DiffDeltaYDD_gene);

			c7->cd(1);
			h1DiffDeltaYDD_seen_GEN[i]->Draw();

			c7->cd(2);
			h1DiffDeltaYDD_gene->Draw();

			c7->SaveAs(Form("./figures_xsec/%d/Acceptance/DD_%s_%s_Acc_DeltaY.pdf", fRunNumber_, fMCName_.Data(), det_labels_[i].Data() ));
			*/
		}


		// Delete fit object
		delete M2fit;
		delete M2Lowfit;

		// Delete canvases
	   	for (UInt_t c = 0; c < ncomb; ++c) {
	   		delete c2DXC[c];
	   		delete c2D[c];
	   	}
	   	delete c5;
	   	delete c7;
	}

}


// Online-Trigger (re-emulated)
// (SPD || V0A || V0C || ADA || ADC)
Bool_t Combinatorics::OnlineTrigger(const CombEvent* event) {

	Bool_t online_trigger = kFALSE;

	UInt_t shift = 0; // Beam-Beam etc. collision masks
	if (trdata_.BCMask.compare("B") == 0) { shift = 0; }
	if (trdata_.BCMask.compare("A") == 0) { shift = 1; }
	if (trdata_.BCMask.compare("C") == 0) { shift = 2; }
	if (trdata_.BCMask.compare("E") == 0) { shift = 3; }

	ULong64_t class_mask = event->fTreeData()->fEventInfo.fClassMask;


	// ---------------------------------------------------------------
	// Here we check that the online trigger bits are valid

	// The magic trigger bits below are from the ALICE e-logbook
	if (fRunNumber_ == 226062) {
		// C0SMB-(B/A/C/E)-NOPF-ALLNOTRD || CINT10-(B/A/C/E)-NOPF-ALLNOTRD
		online_trigger = ( ((class_mask & (1ULL<<(12+shift))) !=0) ||
			             (  (class_mask & (1ULL<<(20+shift))) !=0) );
	} 
	else if (fRunNumber_ == 234039) {
		// CINT11-(B/A/C/E)-NOPF-ALLNOTRD
		online_trigger = (  (class_mask & (1ULL<<(28+shift))) !=0 );
	}
	else if (fRunNumber_ == 274593) {
		// CINT11-(B/A/C/E)-NOPF-CENTNOTRD
		const int offset = shift == 0 ? 0 : 1; // Due to I-class between B and A,C,E
		online_trigger = (  (class_mask & (1ULL<<(34+shift+offset))) !=0 );
	}
	else if (fRunNumber_ == 274594) {
		// CINT11-(B/A/C/E)-NOPF-CENTNOTRD
		const int offset = shift == 0 ? 0 : 1; // Due to I-class
		online_trigger = (  (class_mask & (1ULL<<(34+shift+offset))) !=0 );
	}
	else if (fRunNumber_ == 274595) {
		// CINT11-(B/A/C/E)-NOPF-CENTNOTRD
		const int offset = shift == 0 ? 0 : 1; // Due to I-class
		online_trigger = (  (class_mask & (1ULL<<(34+shift+offset))) !=0 );
	} else {

		printf("Combinatorics::OnlineTrigger:: Unknown RUN = %d \n", fRunNumber_);
	}


	// ---------------------------------------------------------------
	// For MC, simulated AD || V0 || SPD online "CINT11" trigger is here
	if (fIsMC_) {

		// Signal
		online_trigger = ( ( event->fFastOrMap()->CountBits() >= 1 ) ||
		  TMath::Nint(event->fTreeData()->fADInfo.fDecisionOnline[C]) == AliTriggerAnalysis::kADBB ||
		  TMath::Nint(event->fTreeData()->fV0Info.fDecisionOnline[C]) == AliTriggerAnalysis::kV0BB ||
		  TMath::Nint(event->fTreeData()->fV0Info.fDecisionOnline[A]) == AliTriggerAnalysis::kV0BB ||
		  TMath::Nint(event->fTreeData()->fADInfo.fDecisionOnline[A]) == AliTriggerAnalysis::kADBB );
	}

	return online_trigger;
}


// Beam-Gas interaction VETO
Bool_t Combinatorics::BGVeto(const CombEvent* event) {

	Bool_t filter_ok = kTRUE;

	if ( (TMath::Nint(event->fTreeData()->fADInfo.fDecisionOffline[C]) == AliTriggerAnalysis::kADBG) ||
		 (TMath::Nint(event->fTreeData()->fV0Info.fDecisionOffline[C]) == AliTriggerAnalysis::kV0BG) ||
	     (TMath::Nint(event->fTreeData()->fV0Info.fDecisionOffline[A]) == AliTriggerAnalysis::kV0BG) ||
		 (TMath::Nint(event->fTreeData()->fADInfo.fDecisionOffline[A]) == AliTriggerAnalysis::kADBG) ) {
		filter_ok = kFALSE;
	}

	return filter_ok;
}


// 1-dimensional Early interaction VETO
Bool_t Combinatorics::EarlyInteractionVeto(const CombEvent* event, const Double_t V0_cut_T[2][2], const Double_t AD_cut_T[2][2]) {

	const Double_t band = 500; // nanosec
	Bool_t filter_ok = kTRUE;

	// Interaction filtering (early events) based on timing veto
	if (  (IsInRange(event->fTreeData()->fADInfo.fTime[C], AD_cut_T[C][0] - band, AD_cut_T[C][0]) ) ||
		  (IsInRange(event->fTreeData()->fV0Info.fTime[C], V0_cut_T[C][0] - band, V0_cut_T[C][0]) ) ||
		  (IsInRange(event->fTreeData()->fV0Info.fTime[A], V0_cut_T[A][0] - band, V0_cut_T[A][0]) ) ||
		  (IsInRange(event->fTreeData()->fADInfo.fTime[A], AD_cut_T[A][0] - band, AD_cut_T[A][0]) ) ) {
	 	filter_ok = kFALSE;
	}
	return filter_ok;
}


// 1-dimensional Late interaction VETO
Bool_t Combinatorics::LateInteractionVeto(const CombEvent* event, const Double_t V0_cut_T[2][2], const Double_t AD_cut_T[2][2]) {

	const Double_t band = 500; // nanosec
	Bool_t filter_ok = kTRUE;

	// Interaction filtering (late events) based on timing veto
	if (  (IsInRange(event->fTreeData()->fADInfo.fTime[C], AD_cut_T[C][1] , AD_cut_T[C][1] + band) ) ||
		  (IsInRange(event->fTreeData()->fV0Info.fTime[C], V0_cut_T[C][1] , V0_cut_T[C][1] + band) ) ||
		  (IsInRange(event->fTreeData()->fV0Info.fTime[A], V0_cut_T[A][1] , V0_cut_T[A][1] + band) ) ||
		  (IsInRange(event->fTreeData()->fADInfo.fTime[A], AD_cut_T[A][1] , AD_cut_T[A][1] + band) ) ) {
	 	filter_ok = kFALSE;
	}
	return filter_ok;
}


// Calculate the number of FO bits in SPD inner and outer layers, on positive and negative eta
void Combinatorics::SPDBitsSplit(const CombEvent* event,
									UInt_t& nbits_inner_minus, UInt_t& nbits_inner_plus,
									UInt_t& nbits_outer_minus, UInt_t& nbits_outer_plus) {
	
	//printf("SPDBitsSplit:: \n");

	// SPD chips inner layer [0 ... 399]
	nbits_inner_minus = 0;  // eta ~ [-2.0, 0.0]
	nbits_inner_plus  = 0;  // eta ~ [0.0,  2.0]

	// SPD chips outer layer [400 ... 1199]
	nbits_outer_minus = 0;  // eta ~ [-1.4, 0.0]
	nbits_outer_plus  = 0;  // eta ~ [0.0,  1.4]
	
	UInt_t b = 0;
	UInt_t b_prev = 0;
	
	while (kTRUE) {

		b = event->fFiredChipMap()->FirstSetBit(b);

		if ((b == b_prev) && b > 0)
			break;

		// Special case, the start (0-th) bit
		if (b == 0 && !event->fFiredChipMap()->TestBitNumber(0))
			break;

		if (b >= 0 && b < 400) { // Inner layer
			if (event->fFiredChipMap()->TestBitNumber(b) == kTRUE && !IsNoisySPD(b)) {
				//		printf("X:b = %d \n", b);
				if (SPD_FO_z_pos_.at(b) < 0) {
					++nbits_inner_minus;
				} else {
					++nbits_inner_plus;
				}
			}
		}
		if (b >= 400 && b < 1200) { // Outer layer
			if (event->fFiredChipMap()->TestBitNumber(b) == kTRUE && !IsNoisySPD(b)) {
				//		printf("X:b = %d \n", b);
				if (SPD_FO_z_pos_.at(b) < 0) {
					++nbits_outer_minus;
				} else {
					++nbits_outer_plus;
				}
			}
		}
		b_prev = b;
		++b;
	}

	// THIS IS SIGNIFICANTLY MORE SLOW
	/*
	for (UInt_t b = 0; b < 400; ++b) {
		if (event->fFiredChipMap()->TestBitNumber(b) == kTRUE && !IsNoisySPD(b)) {
			printf("1:b = %d \n", b);
			if (SPD_FO_z_pos_.at(b) < 0) {
				++nbits_inner_minus;
			} else {
				++nbits_inner_plus;
			}
		}
	}
	
	
	for (UInt_t b = 400; b < 1200; ++b) {
		if (event->fFiredChipMap()->TestBitNumber(b) == kTRUE && !IsNoisySPD(b)) {
			printf("2:b = %d \n", b);
			if (SPD_FO_z_pos_.at(b) < 0) {
				++nbits_outer_minus;
			} else {
				++nbits_outer_plus;
			}
		}
	}
	*/
}


// Construct Boolean measurement vector at DETECTOR level
std::vector<Bool_t> Combinatorics::ConstructVectorDet(const CombEvent* event, 
	const Double_t V0_cut_T[2][2], const Double_t AD_cut_T[2][2], 
	const Double_t V0_cut_Q[2][2], const Double_t AD_cut_Q[2][2], 
	const UInt_t SPD_cut[2], const Double_t ZN_cut[2]) {

	std::vector<Bool_t> s(d_, kFALSE);

	// ZDN energy cut
	//s.at(ZDN_ind_[C]) = (event->fTreeData()->fZDCInfo.fZNEnergy[C] > ZN_cut[C]);

	// Own TIME CUTS as Decision + Charge cuts
	//s.at(AD_ind_[C])  = //TMath::Nint(event->fTreeData()->fADInfo.fDecisionOffline[C]) == AliTriggerAnalysis::kADBB && 
	s.at(AD_ind_[C])  = IsInRange(event->fTreeData()->fADInfo.fTime[C],   AD_cut_T[C][0], AD_cut_T[C][1]) &&
						IsInRange(event->fTreeData()->fADInfo.fCharge[C], AD_cut_Q[C][0], AD_cut_Q[C][1]);

	//s.at(V0_ind_[C])  = TMath::Nint(event->fTreeData()->fV0Info.fDecisionOffline[C]) == AliTriggerAnalysis::kV0BB && 
	s.at(V0_ind_[C])  = IsInRange(event->fTreeData()->fV0Info.fTime[C],   V0_cut_T[C][0], V0_cut_T[C][1]) &&
						IsInRange(event->fTreeData()->fV0Info.fCharge[C], V0_cut_Q[C][0], V0_cut_Q[C][1]);

	// SPD(C)- and SPD(A)+
	s.at(SPD_ind_[C]) = ((nbits_inner_[C] + nbits_outer_[C]) >= SPD_cut[C]); // At least 2 is usually crucial, due to noise
	s.at(SPD_ind_[A]) = ((nbits_inner_[A] + nbits_outer_[A]) >= SPD_cut[A]); // -||-

	// Own TIME CUTS as Decision + Charge cuts
	//s.at(V0_ind_[A])  = TMath::Nint(event->fTreeData()->fV0Info.fDecisionOffline[A]) == AliTriggerAnalysis::kV0BB && 
	s.at(V0_ind_[A])  = IsInRange(event->fTreeData()->fV0Info.fTime[A],   V0_cut_T[A][0], V0_cut_T[A][1]) &&
						IsInRange(event->fTreeData()->fV0Info.fCharge[A], V0_cut_Q[A][0], V0_cut_Q[A][1]);
						
	//s.at(AD_ind_[A])  = TMath::Nint(event->fTreeData()->fADInfo.fDecisionOffline[A]) == AliTriggerAnalysis::kADBB && 
	s.at(AD_ind_[A])  = IsInRange(event->fTreeData()->fADInfo.fTime[A],   AD_cut_T[A][0], AD_cut_T[A][1]) &&
						IsInRange(event->fTreeData()->fADInfo.fCharge[A], AD_cut_Q[A][0], AD_cut_Q[A][1]);
						

	// ZDN energy cut
	//s.at(ZDN_ind_[A]) = (event->fTreeData()->fZDCInfo.fZNEnergy[A] > ZN_cut[A]);

	return s;
}

// Construct Boolean measurement vector at GENERATOR / PARTICLE / FIDUCIAL ACCEPTANCE level
std::vector<Bool_t> Combinatorics::ConstructVectorGen(const CombEvent* event) {

	std::vector<Bool_t> s(d_, kFALSE);
	
	//s.at(ZDN_ind_[C]) = (event->fMCInfo()->ZDNC.NNeutral > 10000);

if  (FOLDING_MODE == 1) {  // Charged only

	s.at(AD_ind_[C])  = (event->fMCInfo()->ADC.NCharged > 0)  && (event->fMCInfo()->ADC.MaxPtCharged > PT_MIN);
	s.at(V0_ind_[C])  = ( (event->fMCInfo()->V0C.NCharged > 0) && (event->fMCInfo()->V0C.MaxPtCharged > PT_MIN));

	s.at(SPD_ind_[C]) = (event->fMCInfo()->SPDC.NCharged > 0) && (event->fMCInfo()->SPDC.MaxPtCharged > PT_MIN);
	s.at(SPD_ind_[A]) = (event->fMCInfo()->SPDA.NCharged > 0) && (event->fMCInfo()->SPDA.MaxPtCharged > PT_MIN);

	s.at(V0_ind_[A])  = ( (event->fMCInfo()->V0A.NCharged > 0) && (event->fMCInfo()->V0A.MaxPtCharged > PT_MIN));
	s.at(AD_ind_[A])  = (event->fMCInfo()->ADA.NCharged > 0)  && (event->fMCInfo()->ADA.MaxPtCharged > PT_MIN);

} else if
    (FOLDING_MODE == 2) {  // Charged + Neutral

	s.at(AD_ind_[C])  = ( (event->fMCInfo()->ADC.NCharged > 0) && (event->fMCInfo()->ADC.MaxPtCharged > PT_MIN) )   || ( (event->fMCInfo()->ADC.NNeutral > 0)  && (event->fMCInfo()->ADC.MaxPtNeutral > PT_MIN) );
	s.at(V0_ind_[C])  = ( (event->fMCInfo()->V0C.NCharged > 0) && (event->fMCInfo()->V0C.MaxPtCharged > PT_MIN) )   || ( (event->fMCInfo()->V0C.NNeutral > 0)  && (event->fMCInfo()->V0C.MaxPtNeutral > PT_MIN) );

	s.at(SPD_ind_[C]) = ( (event->fMCInfo()->SPDC.NCharged > 0) && (event->fMCInfo()->SPDC.MaxPtCharged > PT_MIN) ) || ( (event->fMCInfo()->SPDC.NNeutral > 0) && (event->fMCInfo()->SPDC.MaxPtNeutral > PT_MIN) );
	s.at(SPD_ind_[A]) = ( (event->fMCInfo()->SPDA.NCharged > 0) && (event->fMCInfo()->SPDA.MaxPtCharged > PT_MIN) ) || ( (event->fMCInfo()->SPDA.NNeutral > 0) && (event->fMCInfo()->SPDA.MaxPtNeutral > PT_MIN) );

	s.at(V0_ind_[A])  = ( (event->fMCInfo()->V0A.NCharged > 0) && (event->fMCInfo()->V0A.MaxPtCharged > PT_MIN) )   || ( (event->fMCInfo()->V0A.NNeutral > 0)  && (event->fMCInfo()->V0A.MaxPtNeutral > PT_MIN) );
	s.at(AD_ind_[A])  = ( (event->fMCInfo()->ADA.NCharged > 0) && (event->fMCInfo()->ADA.MaxPtCharged > PT_MIN) )   || ( (event->fMCInfo()->ADA.NNeutral > 0)  && (event->fMCInfo()->ADA.MaxPtNeutral > PT_MIN) );

} else if
	
    (FOLDING_MODE == 3) {  // Neutral only

	s.at(AD_ind_[C])  = (event->fMCInfo()->ADC.NNeutral > 0)  && (event->fMCInfo()->ADC.MaxPtNeutral > PT_MIN);
	s.at(V0_ind_[C])  = (event->fMCInfo()->V0C.NNeutral > 0)  && (event->fMCInfo()->V0C.MaxPtNeutral > PT_MIN);

	s.at(SPD_ind_[C]) = (event->fMCInfo()->SPDC.NNeutral > 0) && (event->fMCInfo()->SPDC.MaxPtNeutral > PT_MIN);
	s.at(SPD_ind_[A]) = (event->fMCInfo()->SPDA.NNeutral > 0) && (event->fMCInfo()->SPDA.MaxPtNeutral > PT_MIN);

	s.at(V0_ind_[A])  = (event->fMCInfo()->V0A.NNeutral > 0)  && (event->fMCInfo()->V0A.MaxPtNeutral > PT_MIN);
	s.at(AD_ind_[A])  = (event->fMCInfo()->ADA.NNeutral > 0)  && (event->fMCInfo()->ADA.MaxPtNeutral > PT_MIN);

} else {

	printf("Combinatorics::ConstructVectorGen:: Unknown folding matrix mode = %d! \n", FOLDING_MODE);
}

/*
	//printf("ADC minpt = %0.2f GeV \n",  event->fMCInfo()->ADC.MinPtCharged);
	
	s.at(AD_ind_[C])  = ( (event->fMCInfo()->ADC.NCharged > 0) & (event->fMCInfo()->ADC.MaxPtCharged > 0.05) );
	s.at(V0_ind_[C])  = ( (event->fMCInfo()->V0C.NCharged > 0) & (event->fMCInfo()->V0C.MaxPtCharged > 0.25) );
	
	s.at(SPD_ind_[C]) = ( (event->fMCInfo()->SPDC.NCharged > 2) & (event->fMCInfo()->SPDC.MaxPtCharged > 0.15) );
	s.at(SPD_ind_[A]) = ( (event->fMCInfo()->SPDA.NCharged > 2) & (event->fMCInfo()->SPDA.MaxPtCharged > 0.15) );
	
	s.at(V0_ind_[A])  = ( (event->fMCInfo()->V0A.NCharged > 0) & (event->fMCInfo()->V0A.MaxPtCharged > 0.25) );
	s.at(AD_ind_[A])  = ( (event->fMCInfo()->ADA.NCharged > 0) & (event->fMCInfo()->ADA.MaxPtCharged > 0.10) );
*/
	//s.at(ZDN_ind_[A]) = (event->fMCInfo()->ZDNA.NNeutral > 10000);
	return s;
}


void Combinatorics::CreateDirectories() {

	// Create output directory in a case (use \" \")
    gSystem->Exec(Form("mkdir \"./figures_xsec/%d/Acceptance\" -p",       fRunNumber_ ));
	//gSystem->Exec(Form("mkdir \"./figures_xsec/%d/Combinations\" -p",   fRunNumber_ ));
	gSystem->Exec(Form("mkdir \"./figures_xsec/%d/CrossSections\" -p",    fRunNumber_ ));
	gSystem->Exec(Form("mkdir \"./figures_xsec/%d/Matrix\" -p",  		  fRunNumber_ ));
	
	gSystem->Exec(Form("mkdir \"./figures_xsec/%d/Fit\" -p",		 	  fRunNumber_ ));
    gSystem->Exec(Form("mkdir \"./figures_xsec/%d/Unfolding\" -p", 	  	  fRunNumber_ ));
	gSystem->Exec(Form("mkdir \"./figures_xsec/%d/Detector\" -p", 	  	  fRunNumber_ ));
	gSystem->Exec(Form("mkdir \"./figures_xsec/%d/Cutflow\" -p", 		  fRunNumber_ ));
	gSystem->Exec(Form("mkdir \"./figures_xsec/%d/Ascii\" -p", 		  	  fRunNumber_ ));
	gSystem->Exec(Form("mkdir \"./figures_xsec/%d/Detector_1D\" -p",      fRunNumber_ ));

    gSystem->Exec(Form("mkdir \"./figures_xsec/%d/XC/%s\" -p", 		  	  fRunNumber_, fMCName_.Data() ));
    gSystem->Exec(Form("mkdir \"./figures_xsec/%d/SPD/%s\" -p", 		  fRunNumber_, fMCName_.Data() ));
    gSystem->Exec(Form("mkdir \"./figures_xsec/%d/Detector_2D/%s\" -p",   fRunNumber_, fMCName_.Data() ));

}

// Return true if the chip is noisy
Bool_t Combinatorics::IsNoisySPD(UInt_t chipnumber) {

	for (UInt_t i = 0; i < SPD_noisy.size(); ++i) {
		if (SPD_noisy.at(i) == chipnumber)
			return kTRUE; // is on the noisy list
	}

	return kFALSE;
}


// For finding out hot SPD fastOR chips
void Combinatorics::PrintHotSPD(UInt_t c) {

	// Print out the bin contents to find out "hot/noisy" SPD fastOR chips
	printf("Combinatorics::HotSPD:: combination = %d \n", c);
	if (c == 1) { // chosen combination
		for (Int_t i = 0; i < hSPDbit[c]->GetNbinsX(); ++i) {
			printf("bin = %d : %0.0f \n", i, hSPDbit[c]->GetBinContent(i + 1)); // +1 because ROOT starts bin indexing from 1
		}
	}
}


// Read data / MC from TTree
void Combinatorics::ReadTree() {


	printf("----------------------------------------------------------------\n");
	printf("Combinatorics::ReadTree:: %s \n", fMCName_.Data());

	if (fIsMC_ && GENERATOR_LEVEL) {
		printf("\n!!! GENERATOR level analysis mode on (GENERATOR_LEVEL == true) !!!\n\n");
	}


	// Init unfolding object (only for MC, though)
	NewUnfoldResponse();
	

	// Basic objects
	AliAnalysisTaskDiffCrossSectionsMM::TreeData* fTreeData   = new AliAnalysisTaskDiffCrossSectionsMM::TreeData();
	TBits* fFastOrMap 					                      = new TBits();
	TBits* fFiredChipMap 					  				  = new TBits();
	AliESDVertex* fVertexSPD 				  				  = new AliESDVertex();
	AliAnalysisTaskDiffCrossSectionsMM::MCInfo* fMCInfo 	  = new AliAnalysisTaskDiffCrossSectionsMM::MCInfo();

	// Read input
	TFile* fFile = new TFile(fFilename_);
	TDirectory* fDir = (TDirectory*)fFile->Get("results.root");
	TTree* fTE;

	// Construct the trigger name
	TString trigger_name;

	if (fRunNumber_ == 226062) {
		trigger_name = Form("TEC0SMB-%s-NOPF-ALLNOTRD_CINT10-%s-NOPF-ALLNOTRD", trdata_.BCMask.c_str(), trdata_.BCMask.c_str());
	} else if (fRunNumber_ == 234039) {
		trigger_name = Form("TECINT11-%s-NOPF-ALLNOTRD", trdata_.BCMask.c_str());
	} else if (fRunNumber_ == 274593 || fRunNumber_ == 274594 || fRunNumber_ == 274595) { 
		trigger_name = Form("TECINT11-%s-NOPF-CENTNOTRD", trdata_.BCMask.c_str());	
	} else {
		printf("Combinatorics::ReadTree:: Unknown input RUN %d !\n", fRunNumber_);
		return;
	}
	printf("Trigger: %s \n", trigger_name.Data());


	// Connect the TTree
	if (fIsMC_) { // MC Tree
   		fTE = (TTree*) fDir->Get("TE");
   	} else {      // DATA Tree
   		fTE = (TTree*) fDir->Get(trigger_name.Data());
   	}
   	if (!fTE) {
   		printf("Combinatorics::ReadTree:: Problem reading the TTree input!");
   		return;
   	}


	FILE* fp;
	// OUTPUT .csv for external analysis
	TString output_file( Form("./figures_xsec/%d/Ascii/%s_%s.csv", fRunNumber_, trigger_name.Data(), fMCName_.Data() ));
	fp = fopen(output_file, "w");
	if (!fp) {
		printf("Combinatorics::ReadTree:: Error, could not open %s output! \n", output_file.Data());
		return;
	}
	fprintf(fp, "Detector level ID (int, 0...2^N-1), Generator level ID, observable 1, ..., observable N (double), event class (MC) (int) \n");


	FILE* fpSPD;
	// OUTPUT .csv for external analysis
	TString output_file2( Form("./figures_xsec/%d/Ascii/%s_%s_SPDFiredChip.csv", fRunNumber_, trigger_name.Data(), fMCName_.Data() ));
	fpSPD = fopen(output_file2, "w");
	if (!fpSPD) {
		printf("Combinatorics::ReadTree:: Error, could not open %s output! \n", output_file.Data());
		return;
	}

   	// Print tree info (DEBUG)
   	//fTE->Print();

	// Set branchings
	fTE->SetBranchAddress("TreeData",     &fTreeData);
	fTE->SetBranchAddress("FastOrMap",    &fFastOrMap);
	fTE->SetBranchAddress("FiredChipMap", &fFiredChipMap);
	fTE->SetBranchAddress("VertexSPD",    &fVertexSPD);

	if (fIsMC_) {
		fTE->SetBranchAddress("MCInfo",   &fMCInfo);
	}

	printf("ROOT-file input: %s \n", fFilename_.Data());

	// Read the number of entries
	Long64_t nentries = fTE->GetEntries();
	printf("Number of entries in the TTree: %llu, ", nentries);

	// MAXIMUM NUMBER OF EVENTS (FOR SPEED/DEBUG)
	if (!fIsMC_) {
		nentries = std::round(nentries * MAXEVENTS_DATA) <= nentries ? std::round(nentries * MAXEVENTS_DATA) : nentries;
		printf("Maximum number of events limit set (%0.2f): %llu \n\n", MAXEVENTS_DATA, nentries);
	}  else {
		nentries = std::round(nentries * MAXEVENTS_MC) <= nentries ? std::round(nentries * MAXEVENTS_MC) : nentries;
		printf("Maximum number of events limit set (%0.2f): %llu \n\n", MAXEVENTS_MC, nentries);
	}
	// Create TTree code
	//fTE->MakeCode("fastCode.cxx");


	// Init CutFlow object (the results will be in an ascii output)
	CutFlow cuts_(nentries, 5, Form("%s", trigger_name.Data()) );
	cuts_.nameCut(0, "Online trigger", kTRUE);
	cuts_.nameCut(1, "!Beam-Gas VETO", kTRUE);
	cuts_.nameCut(2, "!Early Interactions (timing)", kTRUE);
	cuts_.nameCut(3, "!Late Interactions (timing)", kFALSE);
	cuts_.nameCut(4, "Offline minimum conditions", kTRUE);
	
	
	// EVENT LOOP speed test
  	TStopwatch watch;
  	watch.Start();

	// EVENT LOOP
	Double_t counter = 0.1; // For printing out the progression
	printf("0 "); fflush(stdout);

	for (UInt_t eventnr = 0; eventnr < nentries; ++eventnr) {

		// Read the event
		fTE->GetEntry(eventnr);

		// Event object
		CombEvent* event;

		// Create event object
		if (!fIsMC_) {
			event = new CombEvent(fTreeData, fFastOrMap, fFiredChipMap, fVertexSPD);
		} else {
			event = new CombEvent(fTreeData, fFastOrMap, fFiredChipMap, fVertexSPD, fMCInfo);
		}

		// Event "nano-object"
		evec ev;

		std::vector<Bool_t> s(d_, kFALSE);     // Create empty signal vector at Detector level
		std::vector<Bool_t> s_gen(d_, kFALSE); // Create empty signal vector at Generator level


		// Monte Carlo
		if (fIsMC_) {

			// 1. Get Scattering Process ID
			ev.proc = event->fMCInfo()->fEventType;

				// Exclude Central-Diffractive events completely
				if (SKIP_CD && ev.proc == 3) { 
					continue;
				}

			// 2. Get Diffractive system Masses at generator level
			ev.M2_L = pow2(event->fMCInfo()->fDiffSys.Mass[0]);
			ev.M2_R = pow2(event->fMCInfo()->fDiffSys.Mass[1]);

			// 3. Generate measurement vector at generator level
			s_gen = ConstructVectorGen(event);
			ev.c_gen = Vec2Ind(s_gen);
		}


		// *************************************************************************
		// This has been carefully constructed, do not play with this!
		// CUTFLOW START -->>

		// THIS SPD calculation is FIRST ALWAYS and done only once!
		SPDBitsSplit(event, nbits_inner_[C], nbits_inner_[A], nbits_outer_[C], nbits_outer_[A]);

		//printf("i[C] = %d, i[A] = %d, o[C] = %d, o[A] = %d \n", nbits_inner_[C], nbits_inner_[A], nbits_outer_[C], nbits_outer_[A]);

		Bool_t event_OK = kTRUE;

		// 0. ONLINE trigger
		if ( !cuts_.cut(eventnr, 0, OnlineTrigger(event)) ) {
			event_OK = kFALSE; // @@@ CUT FLOW !!
		}

		// 1. BEAM-GAS
		if ( !cuts_.cut(eventnr, 1, BGVeto(event)) )  {
			event_OK = kFALSE; // @@@ CUT FLOW !!
		}

		// 2. Interaction
		if ( !cuts_.cut(eventnr, 2, EarlyInteractionVeto(event, V0_cut_T, AD_cut_T)) )  {
			event_OK = kFALSE; // @@@ CUT FLOW !!
		}

		// 3. Interaction
		if ( !cuts_.cut(eventnr, 3, LateInteractionVeto(event, V0_cut_T, AD_cut_T)) )  {
			event_OK = kFALSE; // @@@ CUT FLOW !!
		}

		// 4. Minimum-bias OFFLINE conditions
		s  = ConstructVectorDet(event, V0_cut_T, AD_cut_T, V0_cut_Q, AD_cut_Q, SPD_cut, ZN_cut);
		ev.c  = VecOper::Vec2Ind(s); // Vector to index
		const Bool_t minimumbias = (ev.c == 0) ? kFALSE : kTRUE;
		if ( !cuts_.cut(eventnr, 4, minimumbias) )  {
			event_OK = kFALSE; // @@@ CUT FLOW !!
		}

		// Check does this event pass
		if (event_OK == kFALSE) {
			ev.c = 0; // We set it to zero
		}
		// <<-- CUTFLOW END
		// *************************************************************************

		// Counter for visual reasons
		if ( (eventnr /(double)nentries) > counter) {
			printf("%0.0f ", counter * 100); fflush(stdout);
			counter += 0.1;
		}


		// #### DEBUG ONLY #### Pure generator level studies
		if (fIsMC_ && GENERATOR_LEVEL) {
			ev.c = ev.c_gen;
			s = s_gen;
			event_OK = (ev.c_gen != 0) ? kTRUE : kFALSE;
		}

		//PrintData(event);

		// Fill nano-object output
		events.push_back(ev);

		// Fill graphical and csv output
		FillHistograms(event, ev);
		FillCSVOutput(event, ev, fp);
		FillSPDAscii(event, eventnr, fpSPD);
		FillGapFlow(event, ev);


		// Delete event object
		delete event;
	}
	printf("100 %% :: "); // at the end of event loop
	watch.Stop();
	printf("Reading finished in %0.1f sec! \n", watch.RealTime() );


	// Close ascii output
	fclose(fp);
	fclose(fpSPD);

	if (GAPFLOW_ON) {
		PrintGapFlowMatrix(gapflow1, "gapflow1");
		PrintGapFlowMatrix(gapflow2, "gapflow2");
		PrintGapFlowMatrix(gapflow3, "gapflow3");
	}


	// CUTFLOW
	std::string output_filename( Form("./figures_xsec/%d/Cutflow/CutFlow_%s.out", fRunNumber_, fMCName_.Data() ));
	cuts_.printFlow(output_filename, "w");
	cuts_.printCorrelation(output_filename, "a+");
	std::vector<Int_t> list = cuts_.getPass();


	// ** Finally generate default model with no modifications **
	GenerateModel(0.0, 1.0);


	// Delete and close the file input
	delete fTE;
	delete fDir;
	fFile->Close();
	delete fFile;

	// Delete temporary variables from heap
	delete fTreeData;
	delete fFastOrMap;
	delete fFiredChipMap;
	delete fMCInfo;


}


// Construct event statistics for MC
void Combinatorics::GenerateModel(Double_t POMERON_DELTA, Double_t XI_MAX) {

	// Call RESET ALWAYS!!
	EmptyData();

	// Generator truth modification parameters
	if (fIsMC_ && !fMCName_.Contains("EPOS-LHC")) {
		printf("\nCombinatorics::GenerateModel:: \n");
		printf(" + POMERON_DELTA = %0.3f with re-weight with dsigma/d(M^2)^(1+POMERON_DELTA) \n", POMERON_DELTA);
		printf(" + Maximum XI = %0.3E (Minimum <DeltaY> = %0.1f) \n", XI_MAX, VecOper::xi2deltaY(XI_MAX) );
		printf(" + Double Diffraction Cutoff Mode = %d (0 = product kinematics, 1 = separate kinematics)\n\n", DD_XIMAX_MODE);
	}


	// Loop over events
	for (UInt_t i = 0; i < events.size(); ++i) {

		// Take the event
		evec ev = events.at(i);

		// MC
		if (fIsMC_) {

			// MC mass cut-off and re-weight
			if (!MassCutOff(ev, XI_MAX))
				continue; // Event does not fullfill hard criteria

			// Calculate event weight for MC reweighting
			ev.weight = CalculateMCReWeight(ev, POMERON_DELTA);

			// Total process counts
			++fMCTotalProcessCount_.at(ev.proc);

			// GENERATOR level likelihood density for each scattering process
			F_gen_.at(ev.c_gen).at(ev.proc) += ev.weight;

			// DETECTOR level -||-
			F_.at(ev.c).at(ev.proc) += ev.weight;

			// Exception is EPOS, which is never weighted, so construct x here directly
			if (fMCName_.Contains("EPOS-LHC")) {
				++x_.at(ev.c);
				++x_gen_.at(ev.c_gen);
			}

		// DATA
		} else {
			if (ev.c != 0) { // Do not increase 0-bin, not an event passing low-level cuts
				++x_.at(ev.c);
				hxDet->Fill(ev.c);
			}
		}
	}

	// -------------------------------------------------------------------
	// Calculate process density functions (probability distributions)
	// Scale event rates to cross sections (only MC, data is processed later)

	if (fIsMC_) {
		CalculateF();           // 1. Construct Likelihood Densities
		x_cor_ = x_; 		    // 2. For MC, corrected == no corrections
		ConstructUnfolding();   // 3. Construct unfolding response functions
		VdmScale();             // 4. Scale to physical units (just to compare with data)
	}
}

// Empty data
void Combinatorics::EmptyData() {

	// Empty constructed event vectors
	for (UInt_t i = 0; i < x_.size(); ++i) {
		x_.at(i)     = 0;
		x_cor_.at(i) = 0;
		x_gen_.at(i) = 0;
		//x_unfolded_.at(i) = 0; // not this one (why not?)
	}

	// Unfolding
	hxDet->Reset();
	hxGen->Reset();
	hxEmpty->Reset();

	if (fIsMC_) {

		// Empty the process counts
		for (UInt_t j = 0; j < C_; ++j) {
			fMCTotalProcessCount_.at(j)   = 0;
		}

		// Empty the class density matrix
		for (Int_t i = 0; i < N_; ++i) {
			for (UInt_t j = 0; j < C_; ++j) {
				F_.at(i).at(j)     = 0;
				F_gen_.at(i).at(j) = 0;
			}
		}
	}
}


// Check diffractive mass hard cutoff routine
Bool_t Combinatorics::MassCutOff(const evec& ev, Double_t XI_MAX) {

	// --> Normal thing continue now
	if (!fMCName_.Contains("EPOS-LHC")) {

		// 2. Diffractive Mass cutoff procedure >>
		if ((ev.proc == 0 || ev.proc == 1 || ev.proc == 2) ) {

			// Single Diffraction (SD)
			if (ev.proc == 0 || ev.proc == 1) { 
				if ( ev.M2_L > XI_MAX * pow2(SQRTS) || ev.M2_R > XI_MAX * pow2(SQRTS) ) {
					return kFALSE; // Skip the event completely
				}
			}

			// Double Diffraction (DD)
			if (ev.proc == 2) {
				if (DD_XIMAX_MODE == 0) { // Combined xi limit (the traditional)
					const Double_t s0 = pow2(MP); // Normalization scale (GeV)^2 ~ proton mass ^2
					if ((ev.M2_L * ev.M2_R) > (XI_MAX * pow2(SQRTS) * s0) ) {
						return kFALSE; // Skip the event completely
					}
				}
				else if (DD_XIMAX_MODE == 1) { // Individual xi limit
					if (ev.M2_L > XI_MAX * pow2(SQRTS) || ev.M2_R > XI_MAX * pow2(SQRTS) ) {
						return kFALSE; // Skip the event completely
					}
				}
			}

			// Sanity check (caused by floating point precision problems in AliROOT with Pythia8)
			// Also at high mass, the calculation of the invariant mass may have sometimes
			// problems with Phojet (no mother-daughter tree available)
			const double M2_min = pow2(MP + 0.14); // ~proton + pi0 threshold
			if ( (ev.proc == 0 && ev.M2_L < M2_min) ||
				 (ev.proc == 1 && ev.M2_R < M2_min) || (ev.proc == 2 && (ev.M2_L < M2_min || ev.M2_R < M2_min) ) ) {
				//printf("Combinatorics::ReadEvent(): ev.proc = %d, M_L = %0.2f GeV, M_R = %0.3f GeV \n", ev.proc, std::sqrt(M2_L+EPS), std::sqrt(M2_R+EPS));
				return kFALSE; // Skip the event completely
			}
		}
	}

	return kTRUE;
}


// Calculate Diffractive events re-weighting
Double_t Combinatorics::CalculateMCReWeight(const evec& ev, Double_t POMERON_DELTA) {

	Double_t weight = 1.0; // Default

	// SDL weighting
	if (ev.proc == 0) {
		weight = std::pow(ev.M2_L + EPS, -POMERON_DELTA);

		//printf("Combinatorics::CalculateMCReWeight:: M2_L = %0.1f, M2_R = %0.1f, ev.weight = %0.5f \n", M2_L, M2_R, ev.weight);

	// SDR weighting
	} else if (ev.proc == 1) {
		weight = std::pow(ev.M2_R + EPS, -POMERON_DELTA);

	// DD
	} else if (ev.proc == 2) {
		// Left and Right system weighting
		Double_t weight_L = std::pow(ev.M2_L + EPS, -POMERON_DELTA);
		Double_t weight_R = std::pow(ev.M2_R + EPS, -POMERON_DELTA);

		// Total ev.weight
		weight = weight_L * weight_R;
	}

	return weight;
}



void Combinatorics::FillSPDAscii(const CombEvent* event, Int_t eventnr, FILE* fp) {
	
	// Save SPD information to outputfile
	if (eventnr < 30000) { // Limit maximum event count
		for (UInt_t b = 0; b < 1200; ++b) {
			
			if ( event->fFiredChipMap()->TestBitNumber(b) == kTRUE ) {
				fprintf(fp, "1,");
			} else {
				fprintf(fp, "0,");
			}
		}
		fprintf(fp, "\n");
	}
	
}


// Ascii output
void Combinatorics::FillCSVOutput(const CombEvent* event, const evec& ev, FILE* fp) {

	// Combination out
	fprintf(fp,"%d,%d,", ev.c, ev.c_gen);

	// Observables out
	for (Int_t i = 0; i < d_; ++i) {
		fprintf(fp, "%0.3f,", GetObservable(event, i));
	}
	if (fIsMC_) {
		const Double_t ETASAFE = 20; // Safety for spurious, because eta has range [-inf,inf]

		// Minimum particle within fiducial window
		Double_t min = 1e16; 
		if (-ETASAFE < event->fMCInfo()->MinEtaVisibleCharged && event->fMCInfo()->MinEtaVisibleCharged < ETASAFE) {
			min = event->fMCInfo()->MinEtaVisibleCharged < min ? event->fMCInfo()->MinEtaVisibleCharged : min;
		}
		/*
		if (-ETASAFE < event->fMCInfo()->MinEtaVisibleNeutral && event->fMCInfo()->MinEtaVisibleNeutral < ETASAFE) {
			min = event->fMCInfo()->MinEtaVisibleNeutral < min ? event->fMCInfo()->MinEtaVisibleNeutral : min;
		}
		*/
		// Maximum particle within fiducial window
		Double_t max = -1e16;
		if (-ETASAFE < event->fMCInfo()->MaxEtaVisibleCharged  && event->fMCInfo()->MaxEtaVisibleCharged < ETASAFE) {
			max = event->fMCInfo()->MaxEtaVisibleCharged > max ? event->fMCInfo()->MaxEtaVisibleCharged : max;
		}
		/*
		if (-ETASAFE < event->fMCInfo()->MaxEtaVisibleNeutral && event->fMCInfo()->MaxEtaVisibleNeutral < ETASAFE) {
			max = event->fMCInfo()->MaxEtaVisibleNeutral > max ? event->fMCInfo()->MaxEtaVisibleNeutral : max;
		}
		*/
		fprintf(fp, "%d,%0.3f,%0.3f,%0.3f\n", ev.proc, min,max, event->fMCInfo()->MaxGapVisibleCharged); // Event class (for Pythia,Phojet)
	} else {
		fprintf(fp, "0,0,0,0\n");
	}


}


// "GapFlow" .csv output
void Combinatorics::FillGapFlow(const CombEvent* event, const evec& ev) {

	// GAPFLOW Analysis
	if (GAPFLOW_ON && ev.c != 0) {

		// Calculate running mean only for events which have passed all the minimum cuts
		++N_ev;

		SPDC_sum += nbits_inner_[C] + nbits_outer_[C];
		SPDA_sum += nbits_inner_[A] + nbits_outer_[A];
		V0C_sum  += event->fTreeData()->fV0Info.fCharge[C];
		V0A_sum  += event->fTreeData()->fV0Info.fCharge[A];
		ADC_sum  += event->fTreeData()->fADInfo.fCharge[C];
		ADA_sum  += event->fTreeData()->fADInfo.fCharge[A];
		ZDC_sum  += event->fTreeData()->fZDCInfo.fZNEnergy[C];
		ZDA_sum  += event->fTreeData()->fZDCInfo.fZNEnergy[A];

		SPDC_sum2 += pow2(nbits_inner_[C] + nbits_outer_[C]);
		SPDA_sum2 += pow2(nbits_inner_[A]  + nbits_outer_[A]);
		V0C_sum2  += pow2(event->fTreeData()->fV0Info.fCharge[C]);
		V0A_sum2  += pow2(event->fTreeData()->fV0Info.fCharge[A]);
		ADC_sum2  += pow2(event->fTreeData()->fADInfo.fCharge[C]);
		ADA_sum2  += pow2(event->fTreeData()->fADInfo.fCharge[A]);
		ZDC_sum2  += pow2(event->fTreeData()->fZDCInfo.fZNEnergy[C]);
		ZDA_sum2  += pow2(event->fTreeData()->fZDCInfo.fZNEnergy[A]);


		if (N_ev > 10000) { // Minimum number of events so that the Running Mean Calculation is ok

			// Loop over thresholds
			for (UInt_t t = 0; t < GAPFLOW_N; ++t) {

				// Calculate new scale
				const Double_t cutscale = ( t / (Double_t) GAPFLOW_N) * GAPFLOW_MAXSCALE;

			    // New local cuts
				const Double_t AD_cut_Q__[2][2] = {{cutscale * ADC_sum/N_ev, 1e12},
				                         	      { cutscale * ADA_sum/N_ev, 1e12}};
				const Double_t V0_cut_Q__[2][2] = {{cutscale * V0C_sum/N_ev, 1e12},
				                         	      { cutscale * V0A_sum/N_ev, 1e12}};
         	    
     	        const UInt_t SPDC = cutscale * SPDC_sum / N_ev;
     	        const UInt_t SPDA = cutscale * SPDA_sum / N_ev;
			   	const UInt_t  SPD_cut__[2] 		= {SPDC, SPDA}; 									  // C-side and A-side >= number of fastOR
			   	const Double_t ZN_cut__[2] 		= {cutscale * ZDC_sum/N_ev, cutscale * ZDA_sum/N_ev}; // C-side and A-side Analog-to-digital units (A.U.)

				{ // Type 1 gapflow (all detectors (rapidity intervals) running)

					/*
				   	printf("Cutscale = %0.2f [ADC,V0C,SPDC,SPDA,V0A,ADA] = [%0.1f, %0.1f, %d, %d, %0.1f, %0.1f] \n", 
				   		cutscale, AD_cut_Q__[C][0], V0_cut_Q__[C][0], SPD_cut__[C], SPD_cut__[A], V0_cut_Q__[A][0], AD_cut_Q__[A][0]);
					*/

				   	// Construct Boolean vector
					std::vector<Bool_t> s_this = 
						ConstructVectorDet(event, V0_cut_T, AD_cut_T, V0_cut_Q__, AD_cut_Q__, SPD_cut__, ZN_cut__);

					// Increase measurement
					++gapflow1.at(t).at( Vec2Ind(s_this) );
				}
				{ // Type 2 gapflow (SPD fixed and the rest running)

				   	// Construct Boolean vector
					std::vector<Bool_t> s_this = 
						ConstructVectorDet(event, V0_cut_T, AD_cut_T, V0_cut_Q__, AD_cut_Q__, SPD_cut, ZN_cut__);

					// Increase measurement
					++gapflow2.at(t).at( Vec2Ind(s_this) );
				}
				{ // Type 3 gapflow (SPD running and the rest fixed)

				   	// Construct Boolean vector
					std::vector<Bool_t> s_this = 
						ConstructVectorDet(event, V0_cut_T, AD_cut_T, V0_cut_Q, AD_cut_Q, SPD_cut__, ZN_cut);

					// Increase measurement
					++gapflow3.at(t).at( Vec2Ind(s_this) );
				}
			}
		}
	}
}


// Basic histogram output
void Combinatorics::FillHistograms(const CombEvent* event, const evec& ev) {


	// Detector level Boolean vector
	std::vector<Bool_t> s = VecOper::Ind2Vec(ev.c, d_);

	if (HISTOGRAMS_ON) {

		if (fIsMC_) {

			// Generator level Boolean vector
			std::vector<Bool_t> s_gen = VecOper::Ind2Vec(ev.c_gen, d_);


			// For Folding Matrix Visualization

			// Gray code ordering
			/*
			UInt_t GC_c     = Binary2Gray(c);
			UInt_t GC_c_gen = Binary2Gray(c_gen);
			*/

			// Folding Matrix Visualizaiton
			M->Fill(ev.c, ev.c_gen);
			
			// Variables to do the shifting and flip in the rapidity edge space
			const Double_t Ddiff = (ev.proc == 0) ? std::log(SQRTS/MP) : -std::log(SQRTS/MP);
			const Double_t Dsign = (ev.proc == 0) ? 1.0 : -1.0;

			// Fill diffractive masses
			if (ev.proc == 0 || ev.proc == 1) { // SDL or SDR

				const Double_t mass   = event->fMCInfo()->fDiffSys.Mass[ev.proc] + EPS;
				const Double_t mass2  = pow2(mass);

				hDiffMassSDLowM_gene[ev.proc]->Fill(mass2, ev.weight);
				hDiffMassSD_gene[ev.proc]->Fill(mass2, ev.weight);
				hDiffDeltaYSD_gene[ev.proc]->Fill(Dsign*TMath::Log( mass2 / std::pow(SQRTS,2) + EPS ) + Ddiff, ev.weight);

				// Loop over C and A sides
				for (Int_t k = 0; k < 2; ++k) {

					// 1D
					UInt_t SPDbits = (k == 0) ? nbits_inner_[C] + nbits_outer_[C] : nbits_inner_[A]  + nbits_outer_[A];

					// Mass
					hAD_M_avgCharge[ev.proc][k]->Fill(std::log10(mass + 1e-5), event->fTreeData()->fADInfo.fCharge[k]);
					hV0_M_avgCharge[ev.proc][k]->Fill(std::log10(mass + 1e-5), event->fTreeData()->fV0Info.fCharge[k]);
					hSPD_M_avgFO[ev.proc][k]->Fill(std::log10(mass + 1e-5), SPDbits);
					hZDN_M_avgCharge[ev.proc][k]->Fill(std::log10(mass + 1e-5), event->fTreeData()->fZDCInfo.fZNEnergy[k]);
					hZDP_M_avgCharge[ev.proc][k]->Fill(std::log10(mass + 1e-5), event->fTreeData()->fZDCInfo.fZPEnergy[k]);

					h2AD_M_Charge[ev.proc][k]->Fill(mass, event->fTreeData()->fADInfo.fCharge[k]);
					h2V0_M_Charge[ev.proc][k]->Fill(mass, event->fTreeData()->fV0Info.fCharge[k]);
					h2SPD_M_FO[ev.proc][k]->Fill(mass, SPDbits);
					h2ZDN_M_Charge[ev.proc][k]->Fill(mass, event->fTreeData()->fZDCInfo.fZNEnergy[k]);
					h2ZDP_M_Charge[ev.proc][k]->Fill(mass, event->fTreeData()->fZDCInfo.fZPEnergy[k]);
				}


				// Particles within fiducial domain
				// C-side
				hAD_N_avgCharge[ev.proc][C]->Fill(event->fMCInfo()->ADC.NCharged + event->fMCInfo()->ADC.NNeutral, event->fTreeData()->fADInfo.fCharge[C]);
				
				// Charged only events
				if (event->fMCInfo()->ADC.NNeutral == 0 && event->fMCInfo()->ADC.NCharged > 0) {
					hAD_Nch_avgCharge[ev.proc][C]->Fill(event->fMCInfo()->ADC.NCharged, event->fTreeData()->fADInfo.fCharge[C]);
				}
				// Neutral only events
				if (event->fMCInfo()->ADC.NNeutral > 0 && event->fMCInfo()->ADC.NCharged == 0) {
					hAD_Nn_avgCharge[ev.proc][C]->Fill(event->fMCInfo()->ADC.NNeutral, event->fTreeData()->fADInfo.fCharge[C]);
				}
				
				hV0_N_avgCharge[ev.proc][C]->Fill(event->fMCInfo()->V0C.NCharged + event->fMCInfo()->V0C.NNeutral, event->fTreeData()->fV0Info.fCharge[C]);

				if (event->fMCInfo()->V0C.NNeutral == 0 && event->fMCInfo()->V0C.NCharged > 0) {
					hV0_Nch_avgCharge[ev.proc][C]->Fill(event->fMCInfo()->V0C.NCharged, event->fTreeData()->fV0Info.fCharge[C]);
				}

				if (event->fMCInfo()->V0C.NNeutral > 0 && event->fMCInfo()->V0C.NCharged == 0) {
					hV0_Nn_avgCharge[ev.proc][C]->Fill(event->fMCInfo()->V0C.NNeutral, event->fTreeData()->fV0Info.fCharge[C]);
				}

				hSPD_N_avgFO[ev.proc][C]->Fill(event->fMCInfo()->SPDC.NCharged + event->fMCInfo()->SPDC.NNeutral, nbits_inner_[C] + nbits_outer_[C]);
				
				if (event->fMCInfo()->SPDC.NNeutral == 0 && event->fMCInfo()->SPDC.NCharged > 0) {
					hSPD_Nch_avgFO[ev.proc][C]->Fill(event->fMCInfo()->SPDC.NCharged, nbits_inner_[C] + nbits_outer_[C]);
				}
				if (event->fMCInfo()->SPDC.NNeutral > 0 && event->fMCInfo()->SPDC.NCharged == 0) {
					hSPD_Nn_avgFO[ev.proc][C]->Fill(event->fMCInfo()->SPDC.NNeutral, nbits_inner_[C] + nbits_outer_[C]);
				}

				hZDN_Nn_avgCharge[ev.proc][C]->Fill(event->fMCInfo()->ZDNC.NNeutral, event->fTreeData()->fZDCInfo.fZNEnergy[C]);
				
				// A-side
				hAD_N_avgCharge[ev.proc][A]->Fill(event->fMCInfo()->ADA.NCharged + event->fMCInfo()->ADA.NNeutral, event->fTreeData()->fADInfo.fCharge[A]);
				
				if (event->fMCInfo()->ADA.NNeutral == 0 && event->fMCInfo()->ADA.NCharged > 0) {
					hAD_Nch_avgCharge[ev.proc][A]->Fill(event->fMCInfo()->ADA.NCharged, event->fTreeData()->fADInfo.fCharge[A]);
				}
				if (event->fMCInfo()->ADA.NNeutral > 0 && event->fMCInfo()->ADA.NCharged == 0) {
					hAD_Nn_avgCharge[ev.proc][A]->Fill(event->fMCInfo()->ADA.NNeutral, event->fTreeData()->fADInfo.fCharge[A]);
				}

				hV0_N_avgCharge[ev.proc][A]->Fill(event->fMCInfo()->V0A.NCharged + event->fMCInfo()->V0A.NNeutral, event->fTreeData()->fV0Info.fCharge[A]);
				
				if (event->fMCInfo()->V0A.NNeutral == 0 && event->fMCInfo()->V0A.NCharged > 0) {
					hV0_Nch_avgCharge[ev.proc][A]->Fill(event->fMCInfo()->V0A.NCharged, event->fTreeData()->fV0Info.fCharge[A]);
				}						
				if (event->fMCInfo()->V0A.NNeutral > 0 && event->fMCInfo()->V0A.NCharged == 0) {
					hV0_Nn_avgCharge[ev.proc][A]->Fill(event->fMCInfo()->V0A.NNeutral, event->fTreeData()->fV0Info.fCharge[A]);
				}

				hSPD_N_avgFO[ev.proc][A]->Fill(event->fMCInfo()->SPDA.NCharged + event->fMCInfo()->SPDA.NNeutral, nbits_inner_[A] + nbits_outer_[A]);
				
				if (event->fMCInfo()->SPDA.NNeutral == 0 && event->fMCInfo()->SPDA.NCharged > 0) {
					hSPD_Nch_avgFO[ev.proc][A]->Fill(event->fMCInfo()->SPDA.NCharged, nbits_inner_[A] + nbits_outer_[A]);
				}
				if (event->fMCInfo()->SPDA.NNeutral > 0 && event->fMCInfo()->SPDA.NCharged == 0) {
					hSPD_Nn_avgFO[ev.proc][A]->Fill(event->fMCInfo()->SPDA.NNeutral, nbits_inner_[A] + nbits_outer_[A]);
				}

				hZDN_Nn_avgCharge[ev.proc][A]->Fill(event->fMCInfo()->ZDNA.NNeutral, event->fTreeData()->fZDCInfo.fZNEnergy[A]);
			}

			if (ev.proc  == 2) { // DD
				h2DiffMassDD_gene->Fill(pow(event->fMCInfo()->fDiffSys.Mass[0],2), pow(event->fMCInfo()->fDiffSys.Mass[1],2), ev.weight);
				h1DiffDeltaYDD_gene->Fill(-TMath::Log( (ev.M2_L*ev.M2_R)/(pow2(MP)*pow2(SQRTS)) + EPS ), ev.weight);
			}

			if (ev.proc == 0 || ev.proc == 1) { // SDL or SDR
				// Loop over detectors
				for (Int_t k = 0; k < d_; ++k) {

					if (s_gen.at(k)) { // Generator level seen
						hDiffMassSD_seen_GEN[k][ev.proc]->Fill(pow(event->fMCInfo()->fDiffSys.Mass[ev.proc],2), ev.weight);
						hDiffDeltaYSD_seen_GEN[k][ev.proc]->Fill(Dsign*TMath::Log( std::pow(event->fMCInfo()->fDiffSys.Mass[ev.proc], 2) / pow2(SQRTS) + EPS ) + Ddiff, ev.weight);
					
					}
					if (s.at(k)) { 	   // Detector level seen
						hDiffMassSD_seen_DET[k][ev.proc]->Fill(pow(event->fMCInfo()->fDiffSys.Mass[ev.proc],2), ev.weight);
						hDiffDeltaYSD_seen_DET[k][ev.proc]->Fill(Dsign*TMath::Log( std::pow(event->fMCInfo()->fDiffSys.Mass[ev.proc], 2) / pow2(SQRTS) + EPS ) + Ddiff, ev.weight);
					
					}
				}
			}

			if (ev.proc  == 2) { // DD
				// Loop over detectors
				for (Int_t k = 0; k < d_; ++k) {
					if (s_gen.at(k)) { 	// Generator level seen
						h2DiffMassDD_seen_GEN[k]->Fill(pow2(event->fMCInfo()->fDiffSys.Mass[0]), pow2(event->fMCInfo()->fDiffSys.Mass[1]), ev.weight);
						h1DiffDeltaYDD_seen_GEN[k]->Fill(-TMath::Log( (ev.M2_L*ev.M2_R)/(pow2(MP)*pow2(SQRTS)) + EPS), ev.weight);	
					}
					if (s.at(k)) { 	    // Detector level seen
						h2DiffMassDD_seen_DET[k]->Fill(pow2(event->fMCInfo()->fDiffSys.Mass[0]), pow2(event->fMCInfo()->fDiffSys.Mass[1]), ev.weight);
						h1DiffDeltaYDD_seen_DET[k]->Fill(-TMath::Log( (ev.M2_L*ev.M2_R)/(pow2(MP)*pow2(SQRTS)) + EPS), ev.weight);
					}
				}
			}
		}

		// Loop over the SPD chips
		Double_t num_fired = 0;
		for (UInt_t b = 0; b < 1200; ++b) {

			if ( event->fFiredChipMap()->TestBitNumber(b) == kTRUE ) {
				 hSPDbit[ev.c]->Fill(b); // Fill if true
			 	 //h2->Fill(i, b);
			 	 ++num_fired;
			}
		}

		// Tracklets
		hSPDTR[ev.c]->Fill(event->fTreeData()->fEventInfo.fnTrklet);

		// ZDC, AD, V0 charge and time 1D histograms & 2D histograms
		for (UInt_t k = 0; k < 2; ++k) {

			if (k == 0) { // C-side
				hSPDFO[ev.c][k]->Fill(nbits_inner_[C] + nbits_outer_[C]);
				h2SPDFOTR[ev.c][k]->Fill(event->fTreeData()->fEventInfo.fnTrklet, nbits_inner_[C] + nbits_outer_[C]);
			} else {      // A-side
				hSPDFO[ev.c][k]->Fill(nbits_inner_[A] + nbits_outer_[A]);
				h2SPDFOTR[ev.c][k]->Fill(event->fTreeData()->fEventInfo.fnTrklet, nbits_inner_[A] + nbits_outer_[A]);		
			}

			// ZDN
			hZDN[ev.c][k]->Fill(event->fTreeData()->fZDCInfo.fZNEnergy[k]);
			hZDP[ev.c][k]->Fill(event->fTreeData()->fZDCInfo.fZPEnergy[k]);

			// AD
			hADCharge[ev.c][k]->Fill(event->fTreeData()->fADInfo.fCharge[k]);
			hADTime[ev.c][k]->Fill(event->fTreeData()->fADInfo.fTime[k]);

			// V0
			hV0Charge[ev.c][k]->Fill(event->fTreeData()->fV0Info.fCharge[k]);
			hV0Time[ev.c][k]->Fill(event->fTreeData()->fV0Info.fTime[k]);


			// 2D
			h2ADCT[ev.c][k]->Fill(event->fTreeData()->fADInfo.fTime[k], event->fTreeData()->fADInfo.fCharge[k]);
			h2V0CT[ev.c][k]->Fill(event->fTreeData()->fV0Info.fTime[k], event->fTreeData()->fV0Info.fCharge[k]);
		}


		// -----------------------------------------------------------
		// 2D-cross cross-correlations

		// First observable
		for (Int_t i = 0; i < d_; ++i) {

			Double_t OBS_i = GetObservable(event, i);

			// Second observable
			for (Int_t j = i; j < d_; ++j) {

				Double_t OBS_j = GetObservable(event, j);
				h2XC[ev.c][i][j]->Fill(OBS_i, OBS_j);
			}
		}

	} // Histograms end here

}


// Return observables
Double_t Combinatorics::GetObservable(const CombEvent* event, UInt_t index) {

	if 		(index == 0) {
		return event->fTreeData()->fADInfo.fCharge[C];
	}
	else if (index == 1) {
		return event->fTreeData()->fV0Info.fCharge[C];
	}
	else if (index == 2) {
		return nbits_inner_[C] + nbits_outer_[C];
	}
	else if (index == 3) {
		return nbits_inner_[A] + nbits_outer_[A];
	}
	else if (index == 4) {
		return event->fTreeData()->fV0Info.fCharge[A];
	}
	else if (index == 5){
		return event->fTreeData()->fADInfo.fCharge[A];
	}
	else {
		printf("Combinatorics::GetObservable:: Unknown index %d \n", index);
		return 0;
	}
}


// Print out GapFlow analysis to ascii
void Combinatorics::PrintGapFlowMatrix(const std::vector<std::vector<Double_t>>& matrix, std::string filename) {

	TString output_file( Form("./figures_xsec/%d/Ascii/%s_%s.out", fRunNumber_, filename.c_str(), fMCName_.Data() ));

	FILE* fp;
	fp = fopen(output_file, "w");
	for (UInt_t i = 0; i < matrix.size(); ++i) {
		for (UInt_t j = 0; j < matrix.at(0).size(); ++j) {
			fprintf(fp, "%0.0f ", matrix.at(i).at(j));
		}
		if (i < matrix.size()-1) {
			fprintf(fp, "\n");
		}
	}

	// Close Filepointer
	fclose(fp);	
}


// Return the data source name
TString Combinatorics::GetName() {
	return fMCName_;
}


// Return total class/process fractions in a vector
std::vector<Double_t> Combinatorics::GetTotalProcessCount(Bool_t normalize) {

	if (normalize) { // Normalize to sum to 1
		std::vector<Double_t> output(C_, 0);

		// Sum the number of events
		Double_t totsum = 0;
		for (UInt_t i = 0; i < fMCTotalProcessCount_.size(); ++i) {
			totsum += fMCTotalProcessCount_.at(i);
		}

		// Normalize
		for (UInt_t i = 0; i < fMCTotalProcessCount_.size(); ++i) {
			output.at(i) = fMCTotalProcessCount_.at(i) / (totsum + EPS);
		}
		return output;

	} else {
		return fMCTotalProcessCount_;
	}
}
