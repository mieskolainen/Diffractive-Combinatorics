// CombinatoricsSuper analysis class for Combinatorial Diffractive Cross Sections
//
// mikael.mieskolainen@cern.ch, 2018
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.


#ifndef COMBINATORICSSUPER_H
#define COMBINATORICSSUPER_H


// C++ headers
#include <iostream>
#include <vector>

// ROOT headers
#include <TObject.h>
#include <TString.h>
#include <TBits.h>
#include <TClonesArray.h>
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TArrayD.h"
#include "THnSparse.h"

// Own
#include "AliAnalysisTaskDiffCrossSectionsMM.h"
#include "Combinatorics.h"
#include "CombEvent.h"


extern Int_t N_EM_ITER;
extern Int_t UNFOLD_ITER;
extern Bool_t SCAN_PARAMETERS;
extern Bool_t VERBOSE_ON;
extern Bool_t MINUIT_ON;
extern Int_t  SCAN_ND;


// Integrated Cross Section object
class CrossSection {

public:
	CrossSection() {}
	~CrossSection() {}

	// Name of the process
	std::string name;

	// Central value, statistical uncertainty, luminosity uncertainty
	Double_t value   = 0.0;
	Double_t stat    = 0.0;
	Double_t lumi    = 0.0;

	// Systematic uncertainties and their names
	std::vector<Double_t> syst;
	std::vector<std::string> syst_name;

	// MC efficiency
	Double_t eff 	 = 0.0;

	// MC parameters
	Double_t DELTA   = 0.0; // Pomeron delta
	Double_t XI_MAX  = 0.0; // Maximum rapidity

	// Fit measures
	Double_t LogL    = 0.0;
	Double_t KL      = 0.0;
	Double_t KS      = 0.0;
	Double_t CHI2    = 0.0;

	void Print() {
		printf("%s: %0.2f +- %0.2f (stat) +- %0.2f (lumi) \n", name.c_str(), value, stat, lumi);
	}
};



class CombinatoricsSuper {

public:
	
	// Constructor and destructor
	CombinatoricsSuper(TString base_path);
	~CombinatoricsSuper();
	
	// Initialization
	Bool_t Initialize(const Int_t RUN);
	
	// Add Data or MC sources
	void AddSource(Combinatorics* source);
	
	// Estimate cross sections
	static void negLogLfunc(int& npar, double* gin, double& f, double* par, int iflag);
	void EstimateEM(Int_t data_index, Int_t MC_index);
	std::vector<Double_t> EMsub(std::vector<CrossSection>& xs, Int_t data_index, Int_t MC_index, Bool_t final_round, UInt_t EXTRACTION_LEVEL);
	
	// Unfolding of combinatorial rates
	Double_t Unfold(Int_t Input_index, Int_t Model_index);

	// Combined plots
	void PrintCrossSections();
	void PlotAll1D();
	void PlotAllMatrix();
	void ratioplot(TCanvas*& c, TPad*& pad1, TPad*& pad2);
	
	// Beam-Gas correction for the histograms
	void CorrectBGHist1(TH1F* hB, TH1F* hA, TH1F* hC, TH1F* hE, Double_t sA, Double_t sC, Double_t sE);
	void CorrectBGHist2(TH2F* hB, TH2F* hA, TH2F* hC, TH2F* hE, Double_t sA, Double_t sC, Double_t sE);

	// Beam mask scalers: [1/triggerdownscale] x [relative beam intensity compensation]
	Bool_t BG_substraction = kTRUE;

	Double_t scaleA = 0.0;   // Empty x A-side beam
	Double_t scaleC = 0.0;   // C-side beam x Empty
	Double_t scaleE = 0.0;   // Empty x Empty


private:

	// Data base path
	TString base_path_;

	// BOOTstrapped sample for 3 different extraction levels
	std::vector<std::vector<Double_t> > BSmatrix1_;
	std::vector<std::vector<Double_t> > BSmatrix2_;
	std::vector<std::vector<Double_t> > BSmatrix3_;


	// Cross Sections for 3 different extraction levels
	std::vector<std::vector<std::vector<CrossSection> > > XSlevel1;
	std::vector<std::vector<std::vector<CrossSection> > > XSlevel2;
	std::vector<std::vector<std::vector<CrossSection> > > XSlevel3;


	// Plot standards
	std::vector<Int_t> markers_;
	std::vector<Int_t> colors_;
	
	// Data and MC sources on "equal" footing
	std::vector<Combinatorics*> sources_;
	
		
	//ClassDef(CombinatoricsSuper,1);        	   // ROOT system integration
};

#endif
