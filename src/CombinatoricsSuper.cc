// CombinatoricsSuper analysis class for Combinatorial Diffractive Cross Sections
//
// mikael.mieskolainen@cern.ch, 2018
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.


// ROOT Headers
#include "TF1.h"
#include "TArrayD.h"
#include "TGraph2D.h"
#include "TH2F.h"
#include "THStack.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TMatrixD.h"
#include "TVector3.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TGaxis.h"
#include "TMinuit.h"


// Own headers
#include "AliAnalysisTaskDiffCrossSectionsMM.h"
#include "CombinatoricsSuper.h"
#include "Combinatorics.h"
#include "CombEvent.h"
#include "VecOper.h"

// Libraries
#include "csv.h"


using namespace VecOper;


// ***********************************************************************
// Default analysis parameters

const Double_t MC_SIGMA_INEL = 80.0;  // N.B. only for reference (80.0 mb by default)
                                      // this affects only pure MC numbers
									  // -> no single effect on data analysis (easily tested by changing it)

// These set the analysis default without parameter scans
const Double_t DEFAULT_POMERON_DELTA = 0.085;
const Double_t DEFAULT_XI_MAX = 0.05;

// Fixed (default) number of unfold-iterations
Int_t UNFOLD_ITER = 5;

// Fixed number of EM-iterations in cross section fits (at least 50 is usually enough)
Int_t N_EM_ITER = 50;

// (DELTA,XI) scans
Bool_t SCAN_PARAMETERS = kTRUE;
Bool_t MINUIT_ON       = kTRUE;
Bool_t VERBOSE_ON      = kFALSE;
Int_t  SCAN_ND         = 30; // Discretization

// ***********************************************************************


// "Globals" for Minuit
static std::vector<Double_t> x_boot;
static std::vector<std::vector<Double_t> > F;


// USER sets these by hand here
Bool_t CombinatoricsSuper::Initialize(const Int_t RUN) {

	printf("CombinatoricsSuper::Initialize:: \n\n");

    std::vector<TString> run_path;
    std::vector<Int_t> run_num;
    std::vector<TString> name;

    std::vector<Double_t> Lscale;
    std::vector<Double_t> LscaleE;
    
    std::vector<Bool_t> is_MC;
    std::vector<Bool_t> weightMC;


    // See e.g. http://alice-analysis.web.cern.ch/sites/alice-analysis.web.cern.ch/files/documents/Analysis/FedericoJD.pdf
    //
    // LHC Beam 1 is the clock-wise beam (from ALICE to sectors 2-3), https://cds.cern.ch/journal/CERNBulletin/2015/12/News%20Articles/2001078
    //     Beam 2 is the anti-clock-wise beam
    //
    //
    // C-side in ALICE is negative pseudorapidities, beam 1 comes from A-side towards C-side -> Mask A-Empty uses beam 1
    // A-side in ALICE is positive pseudorapidities, beam 2 comes from C-side towards A-side -> Mask C-Empty uses beam 2
    // 
    // https://alice-servotech.web.cern.ch/ALICE-ServoTech/HELP_DCDB-SVT/Help_Files/ALICE-INT-2003-038.pdf

	/*
	
	// ALICE muon arm is on C-side (Gex side)

	ALICE local coordinate system definition (a right handed-Cartesian coordinate system)

	x-axis  
	  - perpendicular  to  the  mean  beam  direction,  aligned  with  the  local  horizontal  and  
	pointing to the accelerator centre. Positive x is from the point of origin toward the accelerator  
	centre  (visual  aid:  Saleve  mountain),  negative  x  is  from  the  point  of  origin  outward  (visual  
	aid: Jura mountain)

	y-axis
	  -  perpendicular  to  the  x  axis  and  the  mean  local  beam  direction,  pointing  upward.    
	Positive  y  is  from  the  point  of  origin  upward,  negative  y  is  from  the  point  of  origin  
	downward

	z-axis
	  -  parallel  to  the  mean  beam  direction.  POSITIVE  z  is  from  the  point  of  origin  toward  
	RB24  (visual  aid:  the  town  of  Bellegarde / Thoiry),  NEGATIVE  z  is  from  the  point  of  origin  toward  
	RB26 (visual aid: the town of Gex / Nyon). In ALICE the muon arm is at negative z
	*/

	// Beam mask "B","A","C","E" and triggerdata maps
	std::map<std::string, TriggerData> trdata;

    if (RUN == 226062) {

    	// LHC Filling Scheme (<SPACING>_<NBb>_<IP1/5>_<IP2>_<IP8>_<code>)
    	// Multi_39b_37_15_15_4bpi11inj
        
        run_path   = {"LHC15f_pass2_DiffCS", "LHC16a2a2_DiffCS", "LHC16a2b2_DiffCS", "LHC16a2c2_DiffCS", "LHC16a2d2_plus_DiffCS"};
        run_num    = {RUN, RUN, RUN, RUN, RUN};
        name       = {Form("Data-%d", RUN), "Pythia-6_(16a2a2)", "Pythia-8-MBR_(16a2b2)", "Phojet_(16a2c2)", "EPOS-LHC_(16a2d2_plus)"};

        Lscale     = {57.8, MC_SIGMA_INEL, MC_SIGMA_INEL, MC_SIGMA_INEL, MC_SIGMA_INEL}; //  V0-AND vdM (mb) result for data, _total inelastic_ for MC
        LscaleE    = {1.2, 0, 0, 0, 0}; // Error on vdM (mb)
        is_MC      = {kFALSE, kTRUE, kTRUE, kTRUE, kTRUE};
        
        // Beam mask scalers: [1/triggerdownscale] x [relative beam intensity compensation]
        // From ALICE e-logbook (set here manually!)

        scaleA = (1/0.1) * (1.157 / 1.776); // (1/triggerdownscale) * Intensity of interacting bunches in beam 1 / Intensity of Non-Interacting bunches beam in 1
        scaleC = (1/0.1) * (1.204 / 1.970); // (1/triggerdownscale) * Intensity of interacting bunches in beam 2 / Intensity of Non-Interacting bunches beam in 2
        scaleE = (1/0.1) * (15/2.0);        // (1/triggerdownscale) * (Number of interacting bunches / Number of Empty bunches)
    }

    else if (RUN == 234039) {

    	// LHC Filling Scheme (<SPACING>_<NBb>_<IP1/5>_<IP2>_<IP8>_<code>)
    	// Multi_51b_8_27_16_4bpi14inj_alt

        run_path   = {"LHC15h_pass1_DiffCS", "LHC16a2a2_DiffCS", "LHC16a2b2_DiffCS", "LHC16a2c2_DiffCS", "LHC16a2d2_DiffCS"};
        run_num    = {RUN, RUN, RUN, RUN, RUN};
        name       = {Form("Data-%d", RUN), "Pythia-6_(16a2a2)", "Pythia-8-MBR_(16a2b2)", "Phojet_(16a2c2)", "EPOS-LHC_(16a2d2)"};

        Lscale     = {57.8, MC_SIGMA_INEL, MC_SIGMA_INEL, MC_SIGMA_INEL, MC_SIGMA_INEL}; // Cross Section [mb], V0-AND vdM result for data, _total inelastic_ for MC
        LscaleE    = {1.2, 0, 0, 0, 0};
        is_MC      = {kFALSE, kTRUE, kTRUE, kTRUE, kTRUE};

        // From ALICE e-logbook (set here manually!)
        scaleA = (1/0.1) * (2.297 / 2.110); // 
        scaleC = (1/0.1) * (2.307 / 2.059); // 
        scaleE = (1/0.1) * (27/2.0);        // 
    }

    else if (RUN == 274593 || RUN == 274594 || RUN == 274595) {

    	// LHC Filling Scheme (<SPACING>_<NBb>_<IP1/5>_<IP2>_<IP8>_<code>)
    	// Multi_57b_56b_25_20_24_4bpi_15inj

        run_path   = {"LHC17j_pass1_DiffCS", "LHC17h7a_DiffCS", "LHC17h7b_DiffCS"};
        run_num    = {RUN, RUN, RUN};
        name       = {Form("Data-%d", RUN), "Pythia-6_(17h7a)", "Phojet_(17h7b)"};

        Lscale     = {57.8, MC_SIGMA_INEL, MC_SIGMA_INEL}; // Cross Section [mb], V0-AND vdM result for data, _total inelastic_ for MC
        LscaleE    = {1.2, 0, 0}; // uncertainty
        is_MC      = {kFALSE, kTRUE, kTRUE};


    } else {
    	printf("CombinatoricsSuper::Initialize:: Unknown input RUN %d !\n", RUN);
    	return kFALSE;
    }


    // ----------------------------------------------------------------------------------
    // Defined CSV-reader with X columns
	io::CSVReader<6> in(Form("./logbookdata/%d.csv", RUN));

	in.read_header(io::ignore_extra_column, "Name", "BCMask", "LMa", "LMb", "L0a", "L0b");
	std::string Name;
	std::string BCMask;
	double LMa      = 0.0;
	double LMb      = 0.0;
	double L0a      = 0.0;
	double L0b      = 0.0;

	while (in.read_row(Name, BCMask, LMa, LMb, L0a, L0b)){
		
		TriggerData data;
		data.Name   = Name;
		data.BCMask = BCMask;
		data.LMb    = LMb;
		data.LMa    = LMa;
		data.L0b    = L0b;
		data.L0a    = L0a;
		data.L0aL0b = L0a/L0b;
		data.L0bLMb = L0b/LMb;

		// Insert pair
		trdata.insert(std::make_pair(BCMask, data));

		data.Print();
	}

	// Read bunch intensity (Start of run, End of run, Average)
	io::CSVReader<4> inbeam(Form("./logbookdata/%d_beam.csv", RUN));
	double SOR = 0.0;
	double EOR = 0.0;
	double Avg = 0.0;
	inbeam.read_header(io::ignore_extra_column, "Name", "SOR", "EOR", "Avg");

	double IB1  = 0.0;
	double NIB1 = 0.0;
	double IB2  = 0.0;
	double NIB2 = 0.0;

	printf("\n");
	printf("Name \t\t\t\t SOR \t\t EOR \t\t Average intensity \n");
	while (inbeam.read_row(Name, SOR, EOR, Avg)) {

		if (Name.compare("Interacting Bunches Beam 1") == 0) {
			IB1  = Avg;
		}
		if (Name.compare("Non Interacting Bunches Beam 1") == 0) {
			NIB1 = Avg;
		}
		if (Name.compare("Interacting Bunches Beam 2") == 0) {
			IB2  = Avg;
		}
		if (Name.compare("Non Interacting Bunches Beam 2") == 0) {
			NIB2 = Avg;
		}
		printf("%s \t %0.2E \t %0.2E \t %0.2E \n", Name.c_str(), SOR, EOR, Avg);
	}
	printf("\n");

	// Beam based
	scaleA = (IB1/NIB1)*(1/0.5);
	scaleC = (IB2/NIB2)*(1/0.5);
	scaleE = (20.0/2.0);

	printf("Calculation 1: a = %0.3f, c = %0.3f, e = %0.3f \n", scaleA, scaleC, scaleE);

	// Calculate scale factors
	scaleA = (trdata["B"].LMb / trdata["A"].LMa) * (trdata["B"].L0aL0b / trdata["A"].L0aL0b);
	scaleC = (trdata["B"].LMb / trdata["C"].LMa) * (trdata["B"].L0aL0b / trdata["C"].L0aL0b);
	scaleE = (trdata["B"].LMb / trdata["E"].LMa) * (trdata["B"].L0aL0b / trdata["E"].L0aL0b);

	printf("Calculation 2: a = %0.3f, c = %0.3f, e = %0.3f \n", scaleA, scaleC, scaleE);

    // ----------------------------------------------------------------------------------

 	// Number of datasets
	const UInt_t N = run_path.size();

	// BEAM-BEAM
	Bool_t isMC = kFALSE;
	Combinatorics* cB = new Combinatorics(Form("%s/%s/%d/AnalysisResults.root", base_path_.Data(), run_path[0].Data(), run_num[0]), name[0], isMC, Lscale[0], LscaleE[0], run_num[0], trdata["B"]);

	cB->CreateDirectories();
	cB->ReadTree();
	cB->Printx();

	// Residual Beam-Gas substraction correction

	// EMPTY - A-side BEAM
	Combinatorics* cA = new Combinatorics(Form("%s/%s/%d/AnalysisResults.root", 
		base_path_.Data(), run_path[0].Data(), run_num[0]), "Data_A", isMC, 1.0, 0.0, run_num[0], trdata["A"]);
	cA->ReadTree();
	//cA->Printx();

	// C-side BEAM - EMPTY
	Combinatorics* cC = new Combinatorics(Form("%s/%s/%d/AnalysisResults.root", 
		base_path_.Data(), run_path[0].Data(), run_num[0]), "Data_C", isMC, 1.0, 0.0, run_num[0], trdata["C"]);
	cC->ReadTree();
	//cC->Printx();

	// EMPTY - EMPTY
	Combinatorics* cE = new Combinatorics(Form("%s/%s/%d/AnalysisResults.root", 
		base_path_.Data(), run_path[0].Data(), run_num[0]), "Data_E", isMC, 1.0, 0.0, run_num[0], trdata["E"]);
	cE->ReadTree();
	//cE->Printx();
	
	// Get vectors out
	std::vector<Double_t> xA = cA->GetX(); // A-side beam
	std::vector<Double_t> xC = cC->GetX(); // C-side beam
	std::vector<Double_t> xE = cE->GetX(); // Empty-Empty

	// Delete obsolete
	delete cA;
	delete cC;
	delete cE;
	
	if (!BG_substraction) {
		printf("Beam-Gas correction OFF \n");
		scaleA = 0;
		scaleC = 0;
		scaleE = 0;
	}

	// Set A/C/E trigger rates and correction factors
	cB->SetBeamGasX(xA, xC, xE, scaleA, scaleC, scaleE);

	// Now generate Bootstrap Sample
	cB->GenerateBootStrap();

	// Correct pileup
	cB->CorrectPileup();

	// Print out to file the corrected rates
	cB->Printx();

	// Plotting
	cB->Plot();
	cB->PlotCodingScheme();
	
	// Add object to the SuperClass
	AddSource(cB);

	// Here read MC, start from index 1
	for (UInt_t i = 1; i < N; ++i) {

		printf("-------------------------------------------------------------\n");
		printf("%s, Run %d \n", name[i].Data(), run_num[i]);
		printf("-------------------------------------------------------------\n");

		TriggerData trnull;
		Combinatorics* cObj = new Combinatorics(Form("%s/%s/%d/AnalysisResults.root", base_path_.Data(), run_path[i].Data(), run_num[i]), name[i], is_MC[i], Lscale[i], LscaleE[i], run_num[i], trnull);

		cObj->CreateDirectories();
		cObj->ReadTree();

		// Now generate Bootstrap Sample
		cObj->GenerateBootStrap();

		// Basic plots
		cObj->Plot();
		cObj->PlotCodingScheme();		

		// Add object to the SuperClass
		AddSource(cObj);
		printf("\n");
	}


    // Combined plots
    PlotAllMatrix();
    PlotAll1D();

	return kTRUE;
}


// Constructors
CombinatoricsSuper::CombinatoricsSuper(TString base_path) {

	printf("CombinatoricsSuper::CombinatoricsSuper:: \n");

	// Data base path
	base_path_ = base_path;

	// SETUP COMMON COLORS here
	const std::vector<Int_t> markers = {20, 26, 24, 25, 32, 27, 28};
	const std::vector<Int_t> colors  = {1, 46, 8, 9, kMagenta+3, 30, 41};
	markers_ = markers;
	colors_  = colors;
}


// Beam-Gas Correct 1D histograms
void CombinatoricsSuper::CorrectBGHist1(TH1F* hB, TH1F* hA, TH1F* hC, TH1F* hE, Double_t sA, Double_t sC, Double_t sE) {

	hB->Add(hA, -sA);
	hB->Add(hC, -sC);
	hB->Add(hE, -sE);
}

// Beam-Gas Correct 2D histograms
void CombinatoricsSuper::CorrectBGHist2(TH2F* hB, TH2F* hA, TH2F* hC, TH2F* hE, Double_t sA, Double_t sC, Double_t sE) {

	hB->Add(hA, -sA);
	hB->Add(hC, -sC);
	hB->Add(hE, -sE);
}

// Destructor
CombinatoricsSuper::~CombinatoricsSuper() {

	/*
	// Free memory allocations
	for (UInt_t i = 0; i < sources_.size(); ++i) {
		delete sources_.at(i);
	}
	*/
}

// Add Data or MC source
void CombinatoricsSuper::AddSource(Combinatorics* source) {
	
	sources_.push_back(source);

}


// Unfolding MC cross-test
Double_t CombinatoricsSuper::Unfold(Int_t Input_index, Int_t Model_index) {
  
  printf("CombinatoricsSuper::Unfold::\n");
  
  
  printf("\n\n=================================== UNFOLD ====================================\n");
  printf("Test input: %s, Folding model: %s \n", sources_.at(Input_index)->GetName().Data(),
  										         sources_.at(Model_index)->GetName().Data());
  printf("FOLDING_MODE = %d \n", FOLDING_MODE);
  if (FOLDING_MODE == 1) {
  	printf("-> Geometric eta windows; Charged particles with pT > %0.02f GeV \n", PT_MIN);
  }
  if (FOLDING_MODE == 2) {
  	printf("-> Geometric eta windows; Charged or Neutral particles with pT > %0.02f GeV \n", PT_MIN);
  }
  if (FOLDING_MODE == 3) {
  	printf("-> Geometric eta windows; Neutral particles with pT > %0.02f GeV \n", PT_MIN);
  }
  printf("FIXED Number of Unfolding iterations: %d \n\n", UNFOLD_ITER);
  printf("-------------------------------------------------------------------------------\n");
  	

  // Get the generator level truth and detector level measured from MC1
  TH1D* hMeas  = sources_.at(Input_index)->hxDet; // Measured
  TH1D* hTrue  = sources_.at(Input_index)->hxGen; // MC generator (particle level) fiducial truth
  TH1D* hEmpty = sources_.at(Input_index)->hxEmpty;
  
  std::vector<Double_t> x_count = sources_.at(Input_index)->GetX();
  
  
  // Get the folding matrix representation (training) from Model
  RooUnfoldResponse* response = sources_.at(Model_index)->response_;
  
  // Unfolded data
  TH1D* hReco;

  printf("Before unfolding: \n");
  printf("- Total number of visible events: %d \n", (int) hMeas->GetEntries());
  printf("- Visible inelastic:      %0.2f (mb) \n", sources_.at(Input_index)->GetVisSigmaInel() );


  Double_t sigma_inel_tot_unfolded = 0.0;
  Double_t sigma_inel_fid_unfolded = 0.0;


  // Last one is the fixed
  std::vector<Int_t> regularization = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25, UNFOLD_ITER};


  // Loop over the regularization parameter for systematic understanding
  for (UInt_t iter = 0; iter < regularization.size(); ++iter) {

  	printf("%d ", regularization.at(iter));

	// Only EM/Bayes will work, due to formal rank deficiency of the problem
	RooUnfoldBayes   unfold (response, hMeas, regularization.at(iter));    // Expectation-Maximization of Bayes-rule based
  	//RooUnfoldBinByBin unfold (response, hMeas);
	
  	//const Int_t tau = 5; // Regularization parameter (integer)
  	//RooUnfoldSvd     unfold (response, hMeas, tau);   	   // Singular Value Decomposition based 
  	//RooUnfoldTUnfold unfold (response, hMeas);  		 	   // TUnfold algorithm

	// == UNFOLD ==
  	hReco = (TH1D*) unfold.Hreco();

  	// Print
  	//unfold.PrintTable(std::cout, hTrue);

	// Export unfolded result to the Combinatorics object
	std::vector<Double_t> x(hReco->GetNbinsX(), 0);
	for (Int_t i = 0; i < hReco->GetNbinsX(); ++i) {
		x.at(i) = hReco->GetBinContent(i + 1); // +1 because ROOT starts bin indexing from 1
	}

	Double_t x0 = x.at(0); // Take it here, EXTRACTION_LEVEL 2 and 3 levels use this differently
	sources_.at(Input_index)->SetXUnfolded(x);
	
  	double unfold_ratio = hReco->GetSumOfWeights() / hMeas->GetEntries();

  	// Set unfolded + extrapolated total
  	sigma_inel_tot_unfolded = sources_.at(Input_index)->GetVisSigmaInel() * unfold_ratio;
  	sources_.at(Input_index)->SetTotSigmaInelUnfolded(sigma_inel_tot_unfolded);

  	// Set unfolded fiducial
  	sigma_inel_fid_unfolded = sigma_inel_tot_unfolded * vsum(x)/(vsum(x) + x0);
	sources_.at(Input_index)->SetFidSigmaInelUnfolded(sigma_inel_fid_unfolded);

	// Call output ascii print
	sources_.at(Input_index)->Printx_unfolded(regularization.at(iter), Model_index);		

		// Statistics
		if (iter == regularization.size()-1) {
		  	printf("\nAfter unfolding:  \n");
		  	printf("- Total number of events: %d \n", (int) hReco->GetSumOfWeights());
		  	printf("- Total number of fiducial events: %d \n", (int) (hReco->GetSumOfWeights() - x.at(0)) );
		  	printf("- Unfolded and extrapolated total inelastic: %0.2f (mb) \n", sigma_inel_tot_unfolded);
		  	printf("- Unfolded fiducial inelastic: %0.2f (mb) \n", 				 sigma_inel_fid_unfolded);
		  	printf("\n");


		  	// Obtain the bootstrapped data
			std::vector<std::vector<Double_t> > BOOTXCorrected = sources_.at(Input_index)->GetBOOTXCorrected();

			// This we calculate now
		  	std::vector<std::vector<Double_t> > BOOTXUnfolded(x.size(), std::vector<Double_t>(N_BOOTSTRAP, 0));

			// Loop over bootstrap samples
			for (UInt_t j = 0; j < BOOTXCorrected.at(0).size(); ++j) {

				std::vector<Double_t> x_cor = VecOper::getcolvec(BOOTXCorrected, j);

				// 1. CONSTRUCT INPUT 
				hEmpty->Reset();
				for (UInt_t c = 0; c < x_cor.size(); ++c) {
					for (UInt_t events = 0; events < std::round(x_cor.at(c)); ++events) {
						hEmpty->Fill(c);
					}
				}
				
				// 2. UNFOLD
			  	RooUnfoldBayes   unfold (response, hEmpty, regularization.at(iter));    // Expectation-Maximization of Bayes-rule based
			  	hReco = (TH1D*) unfold.Hreco();

				// 3. CONSTRUCT OUTPUT
				std::vector<Double_t> unfolded(hReco->GetNbinsX(), 0);
				for (Int_t i = 0; i < hReco->GetNbinsX(); ++i) {
					unfolded.at(i) = hReco->GetBinContent(i + 1); // +1 because ROOT starts bin indexing from 1
				}

				// Fill it as a column in the matrix
				for (UInt_t i = 0; i < unfolded.size(); ++i) {
					BOOTXUnfolded.at(i).at(j) = unfolded.at(i);
				}
			}

			// Save unfolded bootstrapped sample
			sources_.at(Input_index)->SetBOOTXUnfolded(BOOTXUnfolded);

  		} // final round

	} // loop
	
  	// *****************************************************************************
  	// Save unfolding plot
	TCanvas* c = new TCanvas("c", "canvas", 400, 300);
  	c->cd();
  	c->SetLogy();

  	// Colors
  	hTrue->SetLineColor(8);    // Green 
  	hTrue->SetLineWidth(2.0);

  	hMeas->SetLineColor(9);    // Blue

  	hReco->SetLineColor(46);   // Read
  	hReco->SetMarkerStyle(21); // Solid square
  	hReco->SetMarkerSize(0.3);
  	hReco->SetMarkerColor(46);

  	// Axis labels and title
	TString title_str = Form("Unfold_Input_%s_Model_%s.pdf",
  						sources_.at(Input_index)->GetName().Data(),
  						sources_.at(Model_index)->GetName().Data() );

	TString titlestr2 = Form("Unfold:: Input: %s, Model: %s", 
  						sources_.at(Input_index)->GetName().Data(),
  						sources_.at(Model_index)->GetName().Data() );

	// Axis labels
	hTrue->SetTitle(titlestr2.Data());
	hTrue->GetXaxis()->SetTitle("#Omega_{k} (combination)");
	hTrue->GetYaxis()->SetTitle("Counts");
    //
  	hMeas->SetTitle(titlestr2.Data());
	hMeas->GetXaxis()->SetTitle("#Omega_{k} (combination)");
	hMeas->GetYaxis()->SetTitle("Counts");

	// Set X-axis tick divisions
	const Int_t n1 = 14;  // Primary divisions
	const Int_t n2 = 5;   // Secondary divisions
	const Int_t n = n1 + 100*n2;

	hTrue->GetXaxis()->SetNdivisions(n);
	hTrue->SetMinimum(0.1);  // Along y-axis
    //
	hMeas->GetXaxis()->SetNdivisions(n);
	hMeas->SetMinimum(0.1);  // Along y-axis

	// Draw histograms
  	if (sources_.at(Input_index)->fIsMC_) {
		hTrue->Draw("SAME"); // True (only in MC)
	}
  	hMeas->Draw("SAME");     // Measured
  	hReco->Draw("SAME");	 // Reconstructed (unfolded) AS THE LAST (ON TOP)

	/*
	A2->SetTitle("exponential axis");
	A2->SetLabelSize(0.03);
	A2->SetTitleSize(0.03);
	A2->SetTitleOffset(1.2);
	*/

  	// Add legend
  	std::vector<std::string> names = {"Unfolded fiducial", "Measured visible", "True fiducial"};
  	std::vector<TH1*> histo; histo.push_back(hReco); histo.push_back(hMeas); histo.push_back(hTrue);
  	
	Double_t xl1 = 0.65; Double_t yl1 = 0.75; Double_t xl2 = xl1 + 0.20; Double_t yl2 = yl1 + 0.125;

	TLegend* leg = new TLegend(xl1,yl1,xl2,yl2);
	leg->AddEntry(hReco, names.at(0).c_str());
	leg->AddEntry(hMeas, names.at(1).c_str());
	if (sources_.at(Input_index)->fIsMC_) {
		leg->AddEntry(hTrue, names.at(2).c_str());
	}
	leg->Draw();

  	c->SaveAs(Form("./figures_xsec/%d/Unfolding/%s", sources_.at(Input_index)->GetRunNumber(), title_str.Data() ));

  	delete c;
  	delete leg;
  	// *****************************************************************************

  	// Call output print
  	sources_.at(Input_index)->PrintCombinatorics(
 									sources_.at(Input_index)->GetXUnfolded(), 
  									"fiducial (unfolded)",
  									sources_.at(Input_index)->GetFidSigmaInelUnfolded());

	return 1.0;
}


// Plot all detector distributions
void CombinatoricsSuper::PlotAll1D() {

	if (!HISTOGRAMS_ON) { return; }

	printf("CombinatoricsSuper::PlotAll1D:: \n");

	const Int_t d_ = sources_.at(0)->d_;
	UInt_t ncomb   = std::pow(2,d_);

   	std::vector<TCanvas*> can(ncomb, NULL);

   	for (UInt_t k = 0; k < ncomb; ++k) {
		can[k] = new TCanvas(Form("can%d",k), Form("hDet_%d",k), 7500, 2000);

		if (d_ == 6) {
			// +2 from ZDN + ZDP on sides
			can[k]->Divide(d_ + 2, 2, 0.01, 0.001); // (horizontal, vertical, x-margin, y-margin) boxes
		}
		if (d_ == 8) {
			// +2 from ZDN + ZDP on sides
			can[k]->Divide(d_, 2, 0.01, 0.001); // (horizontal, vertical, x-margin, y-margin) boxes
		}
	}
	
	// Loop over sources in inverse order (to have data on top)
	for (Int_t i = sources_.size() - 1; i > -1; --i) {
		sources_.at(i)->Plot1D(can, colors_.at(i), markers_.at(i));
	}

	// Create .pdf files
	for (UInt_t k = 0; k < ncomb; ++k) {
		can[k]->SaveAs(Form("./figures_xsec/%d/Detector_1D/Vector_ID%03u.pdf", sources_.at(0)->GetRunNumber(), k));		
	}

	// Delete canvases
	for (UInt_t k = 0; k < can.size(); ++k) {
		delete can[k];
	}

	// --------------------------------------------------------
	Int_t Ndet = 5; // FIXED Number of detector (ZDN,ZDP,AD,V0,SPD)
   	std::vector<TCanvas*> can2(Ndet, NULL);

   	for (Int_t k = 0; k < Ndet; ++k) {
		can2[k] = new TCanvas(Form("can2%d",k), Form("hDet_%d",k), 5000, 4000);
		can2[k]->Divide(2,2); // (horizontal, vertical, x-margin, y-margin) boxes
	}
	
	// Loop over sources in inverse order (to have data on top)
	for (Int_t i = sources_.size() - 1; i > -1; --i) {
		sources_.at(i)->Plot1DA(can2, colors_.at(i), markers_.at(i));
	}

	// Create .pdf files
	for (Int_t k = 0; k < Ndet; ++k) {
		can2[k]->SaveAs(Form("./figures_xsec/%d/Detector/Detector_M_%03u.pdf", sources_.at(0)->GetRunNumber(), k));		
	}

	// Delete canvases
	for (UInt_t k = 0; k < can2.size(); ++k) {
		delete can2[k];
	}	

	// --------------------------------------------------------
	
	Ndet = 4; // Number of detector (ZDP,AD,V0,SPD)
   	std::vector<TCanvas*> can3(Ndet, NULL);

   	for (Int_t k = 0; k < Ndet; ++k) {
		can3[k] = new TCanvas(Form("can2%d",k), Form("hDet_%d",k), 5000, 1600);
		can3[k]->Divide(6,2,0.001, 0.001); // (horizontal, vertical, x-margin, y-margin) boxes
	}
	
	// Loop over sources in inverse order (to have data on top)
	for (Int_t i = sources_.size() - 1; i > -1; --i) {
		sources_.at(i)->Plot1DB(can3, colors_.at(i), markers_.at(i));
	}

	// Create .pdf files
	for (Int_t k = 0; k < Ndet; ++k) {
		can3[k]->SaveAs(Form("./figures_xsec/%d/Detector/Detector_N_%03u.pdf", sources_.at(0)->GetRunNumber(), k));		
	}

	// Delete canvases
	for (UInt_t k = 0; k < can3.size(); ++k) {
		delete can3[k];
	}
}


// Make ratio coding matrices
void CombinatoricsSuper::PlotAllMatrix() {

if (!HISTOGRAMS_ON) { return; }

	printf("CombinatoricsSuper::PlotAllMatrix:: \n");

	const Int_t d = sources_.at(0)->d_;
	const Int_t N = std::pow(2,d);

	std::vector<UInt_t> ID(std::pow(2,d),0);
	for (UInt_t i = 0; i < ID.size()-1; ++i) {
		ID.at(i) = N-i-1;
	}

	TCanvas* c1 = new TCanvas("cMatrix","Coding matrices",700,5000);
	//c1->SetGrid();
	c1->SetLeftMargin(0.12);
	c1->SetTopMargin(0.04);
	c1->SetBottomMargin(0.04);

   	TH2F* h2data = sources_.at(0)->h2a;

	std::vector<Double_t> x_data = sources_.at(0)->GetX();
	x_data.at(0) = 0;
	x_data = nvec(x_data);

	// Loop over MC sources
	for (UInt_t i = 1; i < sources_.size(); ++i) {

		TH2F* h2mc = sources_.at(i)->h2a;
		
		// Ratio
		h2mc->Divide(h2data);

		c1->cd(1);
   		h2mc->SetMarkerSize(0.8);
		h2mc->LabelsDeflate("X");
		h2mc->LabelsDeflate("Y");
		h2mc->GetZaxis()->SetRangeUser(0.0, 2.0);

		h2mc->SetTitle(Form("Vector space B^{%d}, Ratio %s/Data;;", d, sources_.at(i)->GetName().Data()));
		h2mc->Draw("TEXT COLZ");


		// Change y-axis text labels to ratios between MC/DATA rates
		std::vector<Double_t> x_mc = sources_.at(i)->GetX();
		x_mc.at(0) = 0;
		x_mc = nvec(x_mc);

		for (Int_t k = 1; k < N+1; ++k) {
			Double_t ratio = TMath::Abs(x_data.at(ID[k-1]) - 0) < EPS ? 0 : x_mc.at(ID[k-1]) / x_data.at(ID[k-1]);
			h2mc->GetYaxis()->SetBinLabel(k, Form("R_{%02d} = %0.2f", ID[k-1], ratio));
		}

		// Add upper text box
		//TPaveText* pt = new TPaveText(-0.003558751,32.38788,4.994464,33.74872);
   		//pt->AddText("#LTcharge#GT    #LTcharge#GT    #LTtracklets#GT    #LTcharge#GT    #LTcharge#GT");
   		//pt->Draw("SAME");

		c1->SaveAs(Form("./figures_xsec/%d/Matrix/Ratio_%s_over_Data.pdf",  
			sources_.at(0)->GetRunNumber(), sources_.at(i)->GetName().Data()) );
	}
	delete c1;
}


// Method for creating canvas with "ratio plot" pad1 above, pad2 (ratio)
//
//  EXAMPLE:
//
//	TCanvas* c3 = new TCanvas("c3","c3",800,800);
//	Double_t ratio = 0.3;
//  Double_t epsilon = 0.0001;
//	TPad* pad1 = new TPad("pad1", "pad1", 0, ratio-epsilon, 1, 1);
//	TPad* pad2 = new TPad("pad2", "pad2", 0, 0, 1, ratio*(1-epsilon) );
//	
//	ratioplot(c3, pad1, pad2);

void CombinatoricsSuper::ratioplot(TCanvas*& c, TPad*& pad1, TPad*& pad2) {

   // Upper plot will be in pad1
   pad1->SetBottomMargin(0.015);    // Upper and lower plot are NOT joined
   //pad1->SetGridx();            	// Vertical grid
   pad1->SetFillColor(0);
   pad1->SetFrameFillColor(0);
   pad1->Draw();             		// Draw the upper pad: pad1

   // Do not draw the Y axis label on the upper plot and redraw a small
   // axis instead, in order to avoid the first label (0) to be clipped.

   //Double_t wmin = 0.0;
   //Double_t wmax = 0.15;
   c->cd();

   // Lower plot will be in pad2
   pad2->SetTopMargin(0);
   pad2->SetFillColor(0);
   pad2->SetFrameFillColor(0);
   pad2->SetBottomMargin(0.4);
   pad2->SetGridx();     			// vertical grid
   pad2->Draw();
}


// Cross Section estimation using
// Expectation-Maximization (Maximum Marginal Likelihood) iteration
void CombinatoricsSuper::EstimateEM(Int_t data_index, Int_t MC_index) {

	// IMPORTANT; CLEAR HERE FIRST!
	printf("CombinatoricsSuper::EstimateEM:: \n");

	if (data_index > (Int_t) (sources_.size() - 1) || data_index < 0) {
		printf("Invalid data_index value %d \n", data_index);
		return;
	}
	if (MC_index > (Int_t) (sources_.size() - 1) || MC_index < 0) {
		printf("Invalid MC_index value %d \n", MC_index);
		return;
	}

	printf("Clearing XSlevel1/2/3 objects \n");
	XSlevel1.clear();
	XSlevel2.clear();
	XSlevel3.clear();

	printf("\n\n======================= Max Likelihood-ESTIMATION =============================\n");
	printf("Test input: %s, MC model: %s \n", sources_.at(data_index)->GetName().Data(),
										      sources_.at(MC_index)->GetName().Data());
  	printf("-------------------------------------------------------------------------------\n");

	{ // Prepare Extraction level 1
		BSmatrix1_ = sources_.at(data_index)->GetBOOTXCorrected();
	}
	{ // Prepare Extraction level 2 (0-bin EXCLUDED)
		BSmatrix2_ = sources_.at(data_index)->GetBOOTXUnfolded();

		// Loop over bootstrap samples
		// Make sure 0-bin is 0
		for (UInt_t j = 0; j < BSmatrix2_.at(0).size(); ++j) {
			BSmatrix2_.at(0).at(j) = 0;
		}
	}
	{ // Prepare Extraction level 3 (0-bin INCLUDED)
		BSmatrix3_ = sources_.at(data_index)->GetBOOTXUnfolded();
	}
	
	// Differential mass distribution reweighting parameter range
	// linspace(a,b,n)
	std::vector<Double_t> delta_range  = linspace(0.0, 0.15,  SCAN_ND);
	std::vector<Double_t> deltaY_range = linspace(1e-3,  7,   SCAN_ND); // do not put exactly 0 lower bound -> unphysical
	
	// If parameter scanning off -> Use default values of (DELTA,XIMAX)
	if (!SCAN_PARAMETERS) {

		std::vector<Double_t> delta_range_   = {DEFAULT_POMERON_DELTA};              delta_range  = delta_range_;
		std::vector<Double_t> deltaY_range_  = {VecOper::xi2deltaY(DEFAULT_XI_MAX)}; deltaY_range = deltaY_range_;
	}
	
	// Create xi-loop vector
	std::vector<Double_t> ximax_range(deltaY_range.size(), 0);
	for (UInt_t i = 0; i < deltaY_range.size(); ++i) {
		ximax_range.at(i) = VecOper::deltaY2xi(deltaY_range.at(i));
	}

	// Output "metrics" of the analysis for 3 different EXTRACTION LEVELS
	std::vector<std::vector<std::vector<Double_t> > > LogL_values_mean(3, std::vector<std::vector<Double_t> >( delta_range.size(), std::vector<Double_t> (ximax_range.size(), 0.0)));
	std::vector<std::vector<std::vector<Double_t> > > KL_values_mean   = LogL_values_mean;
	std::vector<std::vector<std::vector<Double_t> > > KS_values_mean   = LogL_values_mean;
	std::vector<std::vector<std::vector<Double_t> > > CHI2_values_mean = LogL_values_mean;


	// DELTA parameter grid scan loop
	for (UInt_t i = 0; i < delta_range.size(); ++i) {

		std::vector<std::vector<CrossSection> > XI1;
		std::vector<std::vector<CrossSection> > XI2;
		std::vector<std::vector<CrossSection> > XI3;

		// XI_max parameter loop
		for (UInt_t j = 0; j < ximax_range.size(); ++j) {

			// ** GENERATE NEW RE-WEIGHTED MC for each (DELTA,XIMAX) pair **
			sources_.at(MC_index)->GenerateModel(delta_range.at(i), ximax_range.at(j));
			
			std::vector<std::vector<CrossSection> > xs_level;

			// Loop over the extraction levels
			for (UInt_t level = 1; level <= 3; ++level) {

				// Call the subroutine which does the Maximum Likelihood Fit
				std::vector<CrossSection> xs;
				EMsub(xs, data_index, MC_index, kFALSE, level);

				// Find the median values
				//std::vector<UInt_t> sortind = vsortind(KL_values);
				//std::reverse(sortind.begin(), sortind.end()); // sort to ascending order
				//ind = 0.5*(KL_values.size()-1);

				// We use just the 0-th process (SDL), they all have common/same metrics (combined fit)
				LogL_values_mean[level-1][i][j] = xs.at(0).LogL;
				KL_values_mean[level-1][i][j]   = xs.at(0).KL;
				KS_values_mean[level-1][i][j]   = xs.at(0).KS;
				CHI2_values_mean[level-1][i][j] = xs.at(0).CHI2;

				// Save (DELTA_P, XI_MAX) values for output prints etc.
				for (UInt_t process = 0; process < xs.size(); ++process) {
					xs.at(process).DELTA  = delta_range.at(i);
					xs.at(process).XI_MAX = ximax_range.at(j);
				}
				// Save this level
				xs_level.push_back( xs );
			}

			XI1.push_back(xs_level.at(0));
			XI2.push_back(xs_level.at(1));
			XI3.push_back(xs_level.at(2));
		}

		XSlevel1.push_back(XI1);
		XSlevel2.push_back(XI2);
		XSlevel3.push_back(XI3);
	}

	// -------------------------------------------------------------------
	// Call finally the EMsub once more with the best mean KL divergence to save it

	if (SCAN_PARAMETERS) {

		for (UInt_t level = 1; level <= 3; ++level) {

			Double_t min_KL = 1e32;
			UInt_t   min_ind[2] = {0, 0};

			// Find the best value
			for (UInt_t i = 0; i < delta_range.size(); ++i) {
				for (UInt_t j = 0; j < ximax_range.size(); ++j) {
					if (KL_values_mean[level-1][i][j] < min_KL) {
						min_KL = KL_values_mean[level-1][i][j];
						min_ind[0] = i;
						min_ind[1] = j;
					}
				}
			}

			printf("CombinatoricsSuper:: OPTIMAL VALUES of 2D-grid scan: delta = %0.4f, ximax = %0.4f \n", 
					delta_range.at(min_ind[0]), ximax_range.at(min_ind[1]));

			// ** GENERATE NEW RE-WEIGHTED MC for each (DELTA,XIMAX) pair **
			sources_.at(MC_index)->GenerateModel(delta_range.at(min_ind[0]), ximax_range.at(min_ind[1]));
			
			std::vector<CrossSection> xs;
			EMsub(xs, data_index, MC_index, kTRUE, level);
		}

	// -------------------------------------------------------------------
	// Write out the data

	// Loop over extraction levels
	for (UInt_t level = 1; level <= 3; ++level) {

	   	TCanvas* ccc = new TCanvas("c","MC-fit",0,0,700,600);
	   	TGraph2D* dt = new TGraph2D();
	   	
	   	Int_t N = 0;
		for (UInt_t i = 0; i < delta_range.size(); ++i) {
			for (UInt_t j = 0; j < ximax_range.size(); ++j) {

				dt->SetPoint(N, delta_range.at(i), deltaY_range.at(j),  KL_values_mean[level-1][i][j] );
				++N;
			}
		}

		dt->SetTitle("MC re-weighting fit (KL-divergence);Mass distribution re-weight parameter #Delta;Minimum <#DeltaY> cutoff");

	   	//gStyle->SetPalette(55);
	   	dt->Draw("COLZ");

		ccc->SaveAs(Form("./figures_xsec/%d/Fit/Fit_level_%d_Input_%s_Model_%s.pdf", sources_.at(0)->GetRunNumber(),
		       level, sources_.at(data_index)->GetName().Data(), sources_.at(MC_index)->GetName().Data()));

		delete dt;
		delete ccc;

		// ----------------------------------------------------------------------------
		// Write fit results to an ascii file
		{
		TString output_file( Form("./figures_xsec/%d/Fit/Fit_level_%d_Input_%s_Model_%s.csv", sources_.at(0)->GetRunNumber(),
		       level, sources_.at(data_index)->GetName().Data(), sources_.at(MC_index)->GetName().Data()) );

		FILE* fp;
		fp = fopen(output_file, "w");
		fprintf(fp, "#Fit point, Mass re-weight DELTA, Min <DeltaY> cutoff, <neg log(L)>, <KL-divergence>, <KS-error>, <CHI2> \n");
	   	N = 0;
		for (UInt_t i = 0; i < delta_range.size(); ++i) {
			for (UInt_t j = 0; j < ximax_range.size(); ++j) {

				fprintf(fp, "%d, %0.5f, %0.5f, %0.2f, %0.9f, %0.9f, %0.2f \n", 
					N, delta_range.at(i), deltaY_range.at(j),
					LogL_values_mean[level-1][i][j], KL_values_mean[level-1][i][j], KS_values_mean[level-1][i][j], CHI2_values_mean[level-1][i][j] );
				++N;
			}
		}

		// Close Filepointer
		fclose(fp);
		}


		{
		TString output_file( Form("./figures_xsec/%d/CrossSections/Extraction_level_%d_Input_%s_Model_%s.csv", sources_.at(0)->GetRunNumber(),
		       level, sources_.at(data_index)->GetName().Data(), sources_.at(MC_index)->GetName().Data()) );

		FILE* fp;
		fp = fopen(output_file, "w");

		if (level == 1)
			fprintf(fp, "#Process, Visible XS-value, stat. uncert, lumi uncert., MC efficiency*acceptance, negative LogL, KL-div, KS-err, CHI2, Mass re-weight DELTA, Min <DeltaY> cutoff \n");
		
		if (level == 2)
			fprintf(fp, "#Process, Fiducial XS-value, stat. uncert, lumi uncert., MC fiducial acceptance, negative LogL, KL-div, KS-err, CHI2, Mass re-weight DELTA, Min <DeltaY> cutoff \n");

		if (level == 3)
			fprintf(fp, "#Process, Total XS-value, stat. uncert, lumi uncert., MC fiducial acceptance, negative LogL, KL-div, KS-err, CHI2, Mass re-weight DELTA, Min <DeltaY> cutoff \n");

		for (UInt_t i = 0; i < delta_range.size(); ++i) {
			for (UInt_t j = 0; j < ximax_range.size(); ++j) {
				for (UInt_t k = 0; k < XSlevel1.at(i).at(j).size(); ++k) {

					CrossSection xs;
					if (level == 1)
						xs = XSlevel1.at(i).at(j).at(k);
					if (level == 2)
						xs = XSlevel2.at(i).at(j).at(k);
					if (level == 3)
						xs = XSlevel3.at(i).at(j).at(k);
					
					// Print out to the file
					fprintf(fp, "%d, %0.5f, %0.5f, %0.5f, %0.5f, %0.5f, %0.5f, %0.5f, %0.5f, %0.4f, %0.4f \n", 
								k, xs.value, xs.stat, xs.lumi, xs.eff, xs.LogL, xs.KL, xs.KS, xs.CHI2, xs.DELTA, VecOper::xi2deltaY(xs.XI_MAX));
				}
			}
		}

		// Close Filepointer
		fclose(fp);
		}

	}

	}

	// ----------------------------------------------------------------------------
	/*
	// PLOT optimization metric for different delta parameter values
	TCanvas* ccc = new TCanvas(Form("ccc%s_%s", sources_.at(data_index)->GetName().Data(), sources_.at(MC_index)->GetName().Data()),"DeltaSearch", 400, 300);
	
	// Median
	TGraph* deltagraph = new TGraph(delta_range.size(), &(delta_range[0]), &(KL_values_median[0]));
	deltagraph->SetLineStyle(1);
	deltagraph->GetYaxis()->SetTitleOffset(1.4); 
	deltagraph->SetTitle(Form("Input: %s, Model: %s;Parameter #Delta; KL(P|Q)", sources_.at(data_index)->GetName().Data(), sources_.at(MC_index)->GetName().Data()));
	deltagraph->GetXaxis()->SetLimits(delta_range.at(0), delta_range.at(delta_range.size()-1));

	// Lower bound
	TGraph* deltagraph_low = new TGraph(delta_range.size(), &(delta_range[0]), &(KL_values_low_1sigma[0]));
	deltagraph_low->SetLineStyle(7);
	deltagraph_low->SetLineColor(2);

	// Upper bound
	TGraph* deltagraph_upp = new TGraph(delta_range.size(), &(delta_range[0]), &(KL_values_upp_1sigma[0]));
	deltagraph_upp->SetLineStyle(7);
	deltagraph_upp->SetLineColor(2);

	ccc->cd();
	deltagraph->Draw("AC");
	deltagraph_low->Draw("C");
	deltagraph_upp->Draw("C");


	ccc->SaveAs(Form("./figures_xsec/%d/CrossSections/M2DistReweight_Input_%s_Model_%s.pdf", sources_.at(0)->GetRunNumber(),
	       sources_.at(data_index)->GetName().Data(), sources_.at(MC_index)->GetName().Data()));

	delete deltagraph;
	delete ccc;
	*/

}

// Minuit negative log-likelihood
void CombinatoricsSuper::negLogLfunc(int& npar, double* gin, double& f, double* par, int iflag) {

	// Collect process probabilities
	std::vector<Double_t> p(F.at(0).size(), 0);

	for (UInt_t i = 0; i < F.at(0).size(); ++i) {
		p.at(i) = par[i];
	}

	// Evaluate negative log-likelihood
	Double_t negLogL = logLmultimix(x_boot, p, F);

	//printf("-log(L) = %0.1f \n\n", negLogL);
	f = negLogL;
}

// EM subroutine
std::vector<Double_t> CombinatoricsSuper::EMsub(std::vector<CrossSection>& xs, Int_t data_index, Int_t MC_index, Bool_t final_round, UInt_t EXTRACTION_LEVEL) {

	printf("========================================================================\n");

	const Int_t C_ = sources_.at(0)->C_; // The number of process classes

	// Vector space dimension
	const Int_t d_ = sources_.at(0)->d_;
	const UInt_t NDF = (std::pow(2,d_)-1) - C_ - 2;


	// Integrated efficiency x acceptance (or only acceptance)
	TMatrixD effVec(C_, 1); // 1,1,1,... by default

	// Visible 
	if (EXTRACTION_LEVEL == 1) {
		// Read the MC density matrix
		F = sources_.at(MC_index)->GetF();
	}
	// Unfolded fiducial
	if (EXTRACTION_LEVEL == 2) {
		// Read the MC density matrix
		F = sources_.at(MC_index)->GetFGen();
	}
	// Extrapolated and unfolded
	if (EXTRACTION_LEVEL == 3) {
		// Read the MC density matrix
		F = sources_.at(MC_index)->GetFGen();
	}
	
	// Construct the MC efficiencies
	for (Int_t j = 0; j < C_; ++j) {    	// Loop over processes

		effVec(j,0) = 1.0 - F.at(0).at(j); 	// Save the efficiency

		// IMPORTANT, make sure that the 0-bin is nullified out of the density matrix,
		// at the level of visible (1), or fiducial (2)
		if (EXTRACTION_LEVEL == 1 || EXTRACTION_LEVEL == 2) {
			F.at(0).at(j) = 0; 		   		// Then nullify
		}

		// Normalize to probability (binned likelihood) vector
		Double_t sum = EPS;
		for (UInt_t i = 0; i < pow(2,d_); ++i) {
			sum += F.at(i).at(j);
		}
		for (UInt_t i = 0; i < pow(2,d_); ++i) {
			F.at(i).at(j) /= sum;
		}

		if (!(sum > EPS)) {
			effVec(j,0) = 0.0; // This process is null (MC does not include it)
		}
	}

	if (EXTRACTION_LEVEL == 1) {
		printf("EXTRACTION_LEVEL == 1 :: MC Model Process ACCEPTANCE x EFFICIENCY (from MC Generator + ALICE simulation): \n");
		effVec.Print();
	}
	if (EXTRACTION_LEVEL == 2) {
		printf("EXTRACTION_LEVEL == 2 :: MC Model Process ACCEPTANCES (from MC Generator): \n");
		effVec.Print();
	}
	if (EXTRACTION_LEVEL == 3) {
		printf("EXTRACTION_LEVEL == 3 :: \n");
	}

	TRandom3* r = new TRandom3();


	// Construct Diagonal Priors Matrix with uniform initialization
/*		std::vector<std::vector<Double_t> > P(pow(2,d_), std::vector<Double_t>(C_, 0.0));

	for (Int_t i = 0; i < C_; ++i) {
		for (Int_t j = 0; j < C_; ++j) {
			if (i == j)
				P.at(i).at(j) = 1.0 / (Double_t) C_; // Uniform init
		}
	}
*/

	// Create [2^d x C_] class density matrix
	TMatrixD Fmat(std::pow(2,d_), C_);

	for (Int_t i = 0; i < std::pow(2,d_); ++i) {
		for (Int_t j = 0; j < C_; ++j) {
			Fmat[i][j] = F.at(i).at(j);
		}
	}

	TMatrixD Kmat(std::pow(2,d_), C_);
	TMatrixD KmatT(C_, std::pow(2,d_));
	
	// Create [C_ x C_] diagonal priors matrix with random initialization
	TMatrixD Pmat(C_, C_);
	TMatrixD Pvec(C_,1);

	Double_t rsum = 0;
	for (Int_t i = 0; i < C_; ++i) {
		for (Int_t j = 0; j < C_; ++j) {
			if (i == j) {

				Pmat[i][j] = r->Rndm();
				rsum += Pmat[i][j];

			} else { Pmat[i][j] = 0.0; }
		}
	}	
	// Normalize
	for (Int_t i = 0; i < C_; ++i) {
		for (Int_t j = 0; j < C_; ++j) {
			if (i == j) {

				Pmat[i][j] /= rsum;

			} else { Pmat[i][j] = 0.0; }
		}
	}

	// -------------------------------------------------------------------
	// Expectation-Maximization (EM) iteration
	if (VERBOSE_ON) {
		printf("MC Process Likelihood Density Matrix (each column with sum 1) (2^N x C) \n");
		Fmat.Print();
	}

	printf("EM-iteration:: N_EM = %d \n", N_EM_ITER);

	TMatrixD Xvec(std::pow(2,d_),1);

	// -----------------------------------------------------------------
	// Bootstrap loop
	std::vector<Double_t> KL(N_BOOTSTRAP, 0);
	std::vector<Double_t> KS(N_BOOTSTRAP, 0);
	std::vector<Double_t> chi2(N_BOOTSTRAP, 0);
	std::vector<Double_t> logL(N_BOOTSTRAP, 0);


	// Extracted process probabilities
	std::vector<std::vector<Double_t>> C_PROB(N_BOOTSTRAP, std::vector<Double_t> (C_, 0.0));


	// WE NEGLECT statistical uncertainty of Monte Carlo (Pythia etc.) productions here
	// -> make sure you have enough statistics for all the combinations of interest.

	// This procedure works at the level of event counts.

	// Bootstrapping loop for statistical uncertanties
	for (UInt_t BOOT_iter = 0; BOOT_iter < N_BOOTSTRAP; ++BOOT_iter) {

		if ((BOOT_iter % 10) == 0) {
			printf("."); fflush(stdout);
		}

		Double_t KL_this = 0;
		Double_t KS_this = 0;
		Double_t chi2_this = 0;
		Double_t logL_this = 0;


		// Draw the bootstrapped event rate vector (counts)
		if (EXTRACTION_LEVEL == 1) {
			x_boot = getcolvec(BSmatrix1_, BOOT_iter);
		}
		if (EXTRACTION_LEVEL == 2) {
			x_boot = getcolvec(BSmatrix2_, BOOT_iter);
		}
		if (EXTRACTION_LEVEL == 3) {
			x_boot = getcolvec(BSmatrix3_, BOOT_iter);
		}

		// Collect [2^N x 1] event rate vector
		for (UInt_t i = 0; i < x_boot.size(); ++i) {
			Xvec[i][0] = x_boot.at(i);
		}

		// Randomize process fractions -> Random initial quess
		rsum = 0;
		for (Int_t i = 0; i < C_; ++i) {
			for (Int_t j = 0; j < C_; ++j) {
				if (i == j) {

					Pmat[i][j] = r->Rndm();
					rsum += Pmat[i][j];

				} else { Pmat[i][j] = 0.0; }
			}
		}	
		// Normalize to event counts (sum of each prosess count = sum of all events)
		for (Int_t i = 0; i < C_; ++i) {
			for (Int_t j = 0; j < C_; ++j) {
				if (i == j) {

					Pmat[i][j] /= (rsum * vsum(x_boot));

				} else { Pmat[i][j] = 0.0; }
			}
		}

		// The actual iterative EM-loop
		for (Int_t EM_iter = 0; EM_iter < N_EM_ITER; ++EM_iter) {

			Kmat = Fmat * Pmat; // Likelihood matrix X Prior process count matrix

			// Matrix Transpose
			for (Int_t i = 0; i < std::pow(2,d_); ++i) {
				for (Int_t j = 0; j < C_; ++j) {
					KmatT[j][i] = Kmat[i][j];
				}
			}

			// Normalize posterior probabilities
			for (Int_t j = 0; j < std::pow(2,d_); ++j) {
				Double_t sum = 1e-12;
				for (Int_t i = 0; i < C_; ++i) {
					sum += KmatT[i][j];
				}
				for (Int_t i = 0; i < C_; ++i) {
					KmatT[i][j] /= sum;
				}
			}

			// New (prior) process event counts
			Pvec = KmatT * Xvec;

			// Update diagonal matrix
			for (Int_t i = 0; i < C_; ++i) {
				Pmat[i][i] = Pvec[i][0];
			}


			// Get process probabilities with sum = 1
			std::vector<Double_t> p(C_, 0);
			for (Int_t s = 0; s < C_; ++s) {
				p.at(s) = Pvec[s][0];
			}
			p = nvec(p);

			// Construct MC synthesized measurement
			Double_t scale = vsum(x_boot);
			std::vector<Double_t> x_hat = GenerateXHat(F, p, scale);

			// Data combinatorial prob. distribution (normalize sum to 1)
			std::vector<Double_t> x_n = nvec(x_boot);

			// Reweighted MC combinatorial prob. distribution (normalize sum to 1)
			std::vector<Double_t> x_hat_n = nvec(x_hat);


			// Measures
			KL_this   = calcKL(x_n, x_hat_n);       // Kullback-Leibler divergence [data, model]
			KS_this   = calcKS(x_n, x_hat_n);       // Kolmogorov-Smirnov 		   [data, model]
			chi2_this = calchi2(x_boot, x_hat);     // Chi^2 count of events 	   [data, model]

			/*if ((EM_iter % 5) == 0 && (BOOT_iter == 0 || BOOT_iter == (N_BOOTSTRAP - 1)))
				printf("Bootstrap %d/%d, EM %2d/%2d : logL: %0.1f, KL-div: %0.4E, KS: %0.4E, chi^2/NDF: %0.1f / %d \n", 
					BOOT_iter+1, N_BOOTSTRAP, EM_iter+1, N_EM_ITER, logL_this, KL_this, KS_this, chi2_this, NDF );
					*/

		} // EM-iteration loop


		// Get process probabilities with sum = 1
		std::vector<Double_t> p(C_, 0);
		for (Int_t s = 0; s < C_; ++s) {
			p.at(s) = Pvec[s][0];
		}
		p = nvec(p);


		// ******************************************************************************
		if (MINUIT_ON) {

			// Call MINUIT to finalize the optimization -> Make sure we are the minimum of negative log-likelihood
			const int NPARAM = C_;

			// Init TMinuit
			TMinuit* gMinuit = new TMinuit(NPARAM);  //initialize TMinuit with a maximum of 5 params
			gMinuit->SetFCN(negLogLfunc);

	        // Set Print Level
	        // -1 no output
	        // 1 standard output
	        gMinuit->SetPrintLevel(-1);


			double arglist[10];
			int ierflg = 0;

			arglist[0] = 1;
			gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);


			// Set starting values and step sizes for parameters
			const double stepsize = 1e-3;
			
			for (UInt_t i = 0; i < p.size(); ++i) {
				gMinuit->mnparm(i, Form("P%d", i),  p[i], stepsize, 0.0, 1.0, ierflg);
			}

			// Fix the central diffraction
			if (SKIP_CD) {
				gMinuit->FixParameter(3);
			}

			// Now ready for minimization step

			// First simplex 
			//gMinuit->mnexcm("SIMPLEX", arglist, 2, ierflg);

			arglist[0] = 500;
			arglist[1] = 1.0;

			// Then migrad
			gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

			// Print results
			/*
			double amin,edm,errdef;
			int nvpar,nparx,icstat;
			gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
			*/

			// Collect parameters out to p[i]
			Double_t empty = 0.0;
			for (UInt_t i = 0; i < p.size(); ++i) {
				gMinuit->GetParameter(i, p[i], empty);
			}

			delete gMinuit;

		}

		// Negative Log-likelihood
		logL_this = logLmultimix(x_boot, p, F); 

		// Get values out of this bootstrapped sample
		for (Int_t i = 0; i < C_; ++i) {
			C_PROB.at(BOOT_iter).at(i) = p.at(i);
		}

		// "Metrics"
		KL.at(BOOT_iter)   = KL_this;
		KS.at(BOOT_iter)   = KS_this;
		chi2.at(BOOT_iter) = chi2_this;
		logL.at(BOOT_iter) = logL_this;

	} // BOOTSTRAP loop

	if (MINUIT_ON) {
		printf("Fit executed with EM + MINUIT finalization \n");
	} else {
		printf("Fit executed with EM \n");
	}


	if (VERBOSE_ON) {

		// Get process probabilities with sum = 1
		std::vector<Double_t> p(C_, 0);
		for (Int_t s = 0; s < C_; ++s) {
			p.at(s) = Pvec[s][0];
		}
		p = nvec(p);

		// Construct MC synthesized measurement
		Double_t scale = vsum(x_boot);
		std::vector<Double_t> x_hat = GenerateXHat(F, p, scale);

		printf("\na = Data, b = MC :: Component statistics after Max Likelihood-fit \n");
		PrintXDiff(x_boot, x_hat);
	}

	// Average quantities
	Double_t LogL_ = vsum(logL)   / double(logL.size());
	Double_t KL_   = vsum(KL)     / double(KL.size());
	Double_t KS_   = vsum(KS)     / double(KS.size());
	Double_t CHI2_ = vsum(chi2)   / double(chi2.size());


	printf("\n");
	printf("<negative Log(L)> = %0.1f, <KL-divergence> = %0.5f, <KS-error> = %0.5f, <CHI^2>/NDF = %0.1f / %d \n", 
			LogL_, KL_, KS_, CHI2_, NDF);


	// Normalization to physical units (mb) with +-1 sigma uncertainties
 	// one could also get exact confidence intervals here from the bootstrap sample.
	std::vector<Double_t> process_xsec(C_, 0.0);
	std::vector<Double_t> process_xsec_stat_uncertainty(C_, 0.0);
	std::vector<Double_t> process_xsec_lumi_uncertainty(C_, 0.0);


	Double_t	normalization       = sources_.at(data_index)->GetVisSigmaInel(); 
	Double_t	normalization_error = sources_.at(data_index)->GetVisSigmaInelError();
	Double_t    relative_error      = normalization_error / normalization;

	if (EXTRACTION_LEVEL == 1) {
		// Already put there ^
	}
	if (EXTRACTION_LEVEL == 2) {
		normalization       = sources_.at(data_index)->GetFidSigmaInelUnfolded();
	}
	if (EXTRACTION_LEVEL == 3) {
		normalization       = sources_.at(data_index)->GetTotSigmaInelUnfolded();
	}


	// *** Calculate mean value ***
	for (Int_t process = 0; process < C_; ++process) {
		for (UInt_t k = 0; k < N_BOOTSTRAP; ++k) {
			process_xsec.at(process) += C_PROB.at(k).at(process);
		}
		process_xsec.at(process) /= (double) N_BOOTSTRAP; // 1/N

		// Normalize to mb
		process_xsec.at(process)             *= normalization;
	}
	// *** Calculate statistical uncertainty standard deviation ***
	for (Int_t process = 0; process < C_; ++process) {
		for (UInt_t k = 0; k < N_BOOTSTRAP; ++k) {
			process_xsec_stat_uncertainty.at(process) += std::pow(process_xsec.at(process) - C_PROB.at(k).at(process)*normalization, 2);
		}
		process_xsec_stat_uncertainty.at(process) /= (double) N_BOOTSTRAP;

		// Standard deviation
		process_xsec_stat_uncertainty.at(process) = std::sqrt( process_xsec_stat_uncertainty.at(process) );
	}
	// *** Calculate luminosity uncertainty standard deviation ***
	for (Int_t process = 0; process < C_; ++process) {

		// Standard deviation
		process_xsec_lumi_uncertainty.at(process) = process_xsec.at(process) * relative_error;
	}



	if (EXTRACTION_LEVEL == 1) {
		printf("\nEXTRACTION_LEVEL = 1 :: Extraction fit in VISIBLE (detector) space: \n");
		printf("\nExtracted VISIBLE process cross sections (mb): \n");

		for (Int_t i = 0; i < C_; ++i) {
			printf("%2d | %7.3f +- %0.3f (stat) +- %0.3f (lumi) \n",
				i, process_xsec.at(i), process_xsec_stat_uncertainty.at(i), process_xsec_lumi_uncertainty.at(i));
		}
		printf("\nSigma_inel^{visible}: %0.3f +- %0.3f (stat) +- %0.3f (lumi) mb\n\n", 
			sources_.at(data_index)->GetVisSigmaInel(), vsum(process_xsec_stat_uncertainty), sources_.at(data_index)->GetVisSigmaInel() * relative_error );

		// --------------------------------------------------------------------------
		// MC based acceptance x efficiency correction (-> extrapolation to low diffractive masses)
		printf("\nMC based integrated ACCEPTANCE x EFFICIENCY extrapolation to TOTAL cross sections (mb): \n");
		std::vector<Double_t> process_xsec_ext(C_, 0.0);
		std::vector<Double_t> process_xsec_ext_stat_uncertainty(C_, 0.0);
		std::vector<Double_t> process_xsec_ext_lumi_uncertainty(C_, 0.0);

		for (Int_t i = 0; i < C_; ++i) {
			process_xsec_ext.at(i) 			   		= process_xsec.at(i) / (effVec[i][0] + 1e-12);
			process_xsec_ext_stat_uncertainty.at(i) = process_xsec_stat_uncertainty.at(i) / (effVec[i][0] + 1e-12);
			process_xsec_ext_lumi_uncertainty.at(i) = process_xsec_lumi_uncertainty.at(i) / (effVec[i][0] + 1e-12);
		}
		Double_t DATA_inel_ext = vsum(process_xsec_ext);

		for (Int_t i = 0; i < C_; ++i) {
			printf("%2d | %7.3f +- %0.3f (stat) +- %0.3f (lumi) \n",
				i, process_xsec_ext.at(i), process_xsec_ext_lumi_uncertainty.at(i), process_xsec_ext_lumi_uncertainty.at(i));
		}
		printf("\nSigma_inel^{total}: %0.3f +- %0.3f (stat) +- %0.3f (lumi) mb \n\n",
			DATA_inel_ext, vnorm2(process_xsec_ext_stat_uncertainty), DATA_inel_ext * relative_error);
	}


	if (EXTRACTION_LEVEL == 2) {
		printf("\nEXTRACTION_LEVEL = 2 :: Extraction fit in FIDUCIAL (unfolded) space : \n");
		printf("\nExtracted FIDUCIAL process cross sections (mb): \n");
		for (Int_t i = 0; i < C_; ++i) {
			printf("%2d | %7.3f +- %0.3f (stat) +- %0.3f (lumi) \n",
				i, process_xsec.at(i), process_xsec_stat_uncertainty.at(i), process_xsec_lumi_uncertainty.at(i));
		}
		printf("\nSigma_inel^{fiducial unfolded}: %0.3f +- %0.3f (stat) +- %0.3f (lumi) mb \n\n",
			sources_.at(data_index)->GetFidSigmaInelUnfolded(), vsum(process_xsec_stat_uncertainty), sources_.at(data_index)->GetFidSigmaInelUnfolded() * relative_error );

		// --------------------------------------------------------------------------
		// MC based acceptance correction (-> extrapolation to low diffractive masses)
		printf("\nMC based integrated ACCEPTANCE extrapolation to TOTAL cross sections (mb): \n");
		std::vector<Double_t> process_xsec_ext(C_, 0.0);
		std::vector<Double_t> process_xsec_ext_stat_uncertainty(C_, 0.0);
		std::vector<Double_t> process_xsec_ext_lumi_uncertainty(C_, 0.0);

		for (Int_t i = 0; i < C_; ++i) {
			process_xsec_ext.at(i) 			        = process_xsec.at(i) / (effVec[i][0] + 1e-12);
			process_xsec_ext_stat_uncertainty.at(i) = process_xsec_stat_uncertainty.at(i) / (effVec[i][0] + 1e-12);
			process_xsec_ext_lumi_uncertainty.at(i) = process_xsec_lumi_uncertainty.at(i) / (effVec[i][0] + 1e-12);
		}
		Double_t DATA_inel_ext = vsum(process_xsec_ext);

		for (Int_t i = 0; i < C_; ++i) {
			printf("%2d | %7.3f +- %0.3f (stat) +- %0.3f (lumi) \n", 
				i, process_xsec_ext.at(i), process_xsec_ext_stat_uncertainty.at(i), process_xsec_ext_lumi_uncertainty.at(i));
		}
		printf("\nSigma_inel^{total estimate}: %0.3f +- %0.3f (stat) +- %0.3f (lumi) mb \n\n", 
			DATA_inel_ext, vnorm2(process_xsec_ext_stat_uncertainty), DATA_inel_ext * relative_error);
	}


	if (EXTRACTION_LEVEL == 3) {
		printf("\nEXTRACTION_LEVEL = 3 :: Extraction fit in TOTAL (unfolded+extrapolated) space : \n");
		printf("\nExtracted TOTAL process cross sections (mb): \n");
		for (Int_t i = 0; i < C_; ++i) {
			printf("%2d | %7.3f +- %0.3f (stat) +- %0.3f (lumi) \n", 
				i, process_xsec.at(i), process_xsec_stat_uncertainty.at(i), process_xsec_lumi_uncertainty.at(i));
		}
		printf("\nSigma_inel^{total unfolded+extrapoled}: %0.3f +- %0.3f (stat) +- %0.3f (lumi) mb \n\n", 
			sources_.at(data_index)->GetTotSigmaInelUnfolded(), vsum(process_xsec_stat_uncertainty), sources_.at(data_index)->GetTotSigmaInelUnfolded() * relative_error);
	}


	// Add Cross Sections to the object
	xs.clear();
	for (Int_t i = 0; i < C_; ++i) {

		CrossSection xsec;
		xsec.value = process_xsec.at(i);              	   // Value
		xsec.stat  = process_xsec_stat_uncertainty.at(i);  // Statistical uncertainty
		xsec.lumi  = process_xsec_lumi_uncertainty.at(i);  // Luminosity uncertainty
		xsec.eff   = effVec[i][0];                    	   // MC efficiency

		// Fit measures
		xsec.LogL  = LogL_;
		xsec.KL    = KL_;
		xsec.KS    = KS_;
		xsec.CHI2  = CHI2_;

		xs.push_back(xsec);
	}


	// --------------------------------------------------------------------------
	// Total inelastic extrapolation visualization

	if (final_round) {
	
	// Vectors
	std::vector<Double_t> inel_vals; // Total inelastic values
	std::vector<Double_t> L_vals;    // Lagrangian cost for the total inelastic

	std::vector<Double_t> tot_SDL;
	std::vector<Double_t> tot_SDR;
	std::vector<Double_t> tot_DD;
	std::vector<Double_t> tot_CD;
	std::vector<Double_t> tot_ND;


	// Maximum value of the loop
	const Double_t sigma_inel_max = 95;


	// Total inelastic loop
	for (Double_t param_tot = sources_.at(data_index)->GetVisSigmaInel(); 
		 	param_tot <= sigma_inel_max; param_tot += 0.01) {

		// ---------------------------------------------------------------
		// "Variational loop"

		// Now extract the efficiency and subsequantly the total inelastic,
		// under the solid constraint: P(SDL_tot) = P(SDR_tot). This routine
		// is plain logic.

		// Technical parameter ("Floating point control", 0 gives floating point problems)
		const Double_t lambda = 0.05; 

		// Boundary condition (minimum efficiency for single diffraction (left or right))
		const Double_t SD_eff_min = 0.5;

		Double_t min_L = 1e99;

		Double_t best_eff[5] = {0}; // Save here the efficiencies

		// These are fixed
		best_eff[3] = effVec[3][0] + 1e-12; // CD process efficiency from MC
		best_eff[4] = effVec[4][0]; // SCAN_ND process efficiency from MC (almost ~100%)

		// The efficiency range from SD_eff_min to 0.99
		for (Double_t eff_SDL = SD_eff_min; eff_SDL < 0.99; eff_SDL += 0.005) {
			for (Double_t eff_SDR = SD_eff_min; eff_SDR < 0.99; eff_SDR += 0.005) {

				//    P(A ∪ B) = P(A) + P(B) - P(A ∩ B)), then assuming A and B are independent gives by definition
				// => P(A ∪ B) = P(A) + P(B) - P(A)P(B)
				Double_t eff_DD = (eff_SDL + eff_SDR) - (eff_SDL*eff_SDR);

				// Main term
				Double_t term1 = pow(Pvec[0][0]/eff_SDL + Pvec[1][0]/eff_SDR + Pvec[2][0]/eff_DD + Pvec[3][0]/best_eff[3] + Pvec[4][0]/best_eff[4] - param_tot, 2);

				// Lambda-term
				Double_t term2 = pow( (Pvec[0][0]/eff_SDL - Pvec[1][0]/eff_SDR), 2);

				//cout << term1 << " " << term2 << endl;

				Double_t lagrangian = term1 + lambda*term2;
				lagrangian /= pow(param_tot,2);

	            if (lagrangian < min_L) {
	            	min_L = lagrangian;
	            	best_eff[0] = eff_SDL;
	            	best_eff[1] = eff_SDR;
	            	best_eff[2] = eff_DD;
	            }
			}
		}

		//printf("Lagrangian x 10^4 = %0.06f (eff_SDL = %0.3f, eff_SDR = %0.3f) [sigma = %0.2f %0.2f %0.2f %0.2f mb] Total Inelastic: %0.2f \n", 
		//		min_L*1e4, best_eff[0], best_eff[1], Pvec[0][0]/best_eff[0], Pvec[1][0]/best_eff[1], Pvec[2][0]/best_eff[2], Pvec[4][0]/best_eff[4], param_tot);

		inel_vals.push_back(param_tot);
		L_vals.push_back(min_L);

		//printf("Sigma_inel = %0.2f mb, Best efficiencies: %0.2f %0.2f %0.2f Eff_SDL/Eff_SDR ratio: %0.2f\n", param_tot, best_eff[0], 
		//	best_eff[1], best_eff[2], best_eff[0]/(best_eff[1]+1e-12));

		tot_SDL.push_back(Pvec[0][0]/best_eff[0]);
		tot_SDR.push_back(Pvec[1][0]/best_eff[1]);
		tot_DD.push_back( Pvec[2][0]/best_eff[2]);
		tot_CD.push_back( Pvec[3][0]/best_eff[3]);
		tot_ND.push_back( Pvec[4][0]/best_eff[4]);
	}

} // final_round

	return KL;
}
