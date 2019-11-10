// Main program for Combinatorial Diffractive Cross Section analysis
// -----------------------------------------------------------------------
//
// The classes are compiled. User might need to slightly modify classes
// to incorporate spesific runs, i.e., include run numbers, trigger masks etc.
// under CombinatoricsSuper::Initialize
// 
// -----------------------------------------------------------------------
// Reading TTree directly (for a quick debug):
// For example:
// TE->Draw("","(fEventInfo.fClassMask & (1<<20)) != 0 && fADInfo.fDecisionOffline[0] == 1
// && fV0Info.fDecisionOffline[0] == 0 && fV0Info.fDecisionOffline[1] == 0 && fADInfo.fDecisionOffline[1] == 0", "GOFF")
//
//
// - Fit based on also on the multiplicity information (1.average), (2.histogram shape)
//   how to take care of gain invariance?
//
// valgrind --tool=callgrind ./(Your binary)
// generates file "callgrind.out.x",
// then use "kcachegrind" to read it
//
// use -pg switch with g++
// -----------------------------------------------------------------------
// 
// mikael.mieskolainen@cern.ch, 2019
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.


// C++
#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>

// ROOT
#include "TROOT.h"
#include "TStyle.h"
#include "TColor.h"
#include "TProfile.h"
#include "TObject.h"
#include "TSystem.h"
#include "TStopwatch.h"

// Own
#include "CombinatoricsSuper.h"
#include "Combinatorics.h"


// Prototypes
void SetROOTStyle();
void SetPlotStyle();


// Main
int main(void) {

  printf("\n\n\nCombinatorial Diffractive Cross Section Analysis \n");
  printf("mikael.mieskolainen@cern.ch, 2019 \n\n\n");
  
  // Supress enum MsgLevel { DEBUG=0, INFO=1, PROGRESS=2, WARNING=3, ERROR=4, FATAL=5 };
  gErrorIgnoreLevel = kWarning;

  SetROOTStyle();
  SetPlotStyle();

  TStopwatch watch;
  watch.Start();

  // ********************************************************************

  // Database working directory
  const TString base_path   = "/home/user/cernbox/ALICE/diffxsdata";

  // Set runs
  std::vector<UInt_t> runs = {274593};//, 274594, 274595};

  // Histogramming on
  HISTOGRAMS_ON = kFALSE;

  // VERBOSE PRINTING
  VERBOSE_ON = kFALSE;

  // Pure generator level study (SKIPS GEANT SIMUATION, ONLY FOR DEBUG!)
  // GENERATOR_LEVEL = kFALSE;

  // Set Gapflow variables
  GAPFLOW_ON = kFALSE;
  GAPFLOW_N  = 100;         // Discretization
  GAPFLOW_MAXSCALE = 10.0;  // Maximum scale

  // Set maximum fraction [0 ... 1] of data or mc (for speed)
  MAXEVENTS_DATA = 0.01;
  MAXEVENTS_MC   = 0.01;

  // Skip central diffraction
  SKIP_CD = kTRUE;

  // Folding matrix construction mode (1 = Charged only, 2 == Charged+Neutral, 3 = Neutral only)
  FOLDING_MODE = 1;
  PT_MIN = 0.050; // At least one particle within fiducial window above this (GeV)

  // Double Diffraction kinematics cutoff mode (0 standard, 1 separate)
  DD_XIMAX_MODE = 0;

  // Number of bootstrap samples
  N_BOOTSTRAP = 100;
  FASTSTAT = kTRUE;

  // Fixed (default) number of unfold-iterations
  UNFOLD_ITER = 5;

  // Fixed number of EM-iterations in cross section fits (at least 50 is usually enough)
  N_EM_ITER = 50;
  
  // (POMERON DELTA, XIMAX) scans
  SCAN_PARAMETERS = kTRUE;
  MINUIT_ON       = kFALSE;
  VERBOSE_ON      = kFALSE;
  SCAN_ND         = 6;      // Discretization

  // ********************************************************************

  // Loop over runs
  for (UInt_t i = 0; i < runs.size(); ++i) {

    // Superclass object
    CombinatoricsSuper* cSuperObj = new CombinatoricsSuper(base_path);

    // Initialize and Run
    cSuperObj->Initialize(runs[i]);

    // HERE BE CAREFUL, BECAUSE AFTER REFITTING HAS BEEN DONE, THE MODEL
    // DEFINITIONS ARE DIFFERENT. THIS NEEDS TO BE FIXED TO MORE ROBUST!

    // ===================================================================
    // EM cross section extraction
    //
    // For 274593, 274594, 274595
    // 0 = data, 1 = Pythia, 2 = Phojet (See CombinatoricsSuper::Initialize)

    for (UInt_t input = 0; input <= 0; ++input) {
      for (UInt_t model = 1; model <= 2; ++model) {

          cSuperObj->Unfold(input, model);
          cSuperObj->EstimateEM(input, model);
      }
    }
    // ===================================================================

    delete cSuperObj;
  }

  watch.Stop();
  printf("\nAnalysis:: Finished in %0.1f sec! \n", watch.RealTime() );

  // printf("Combining pdfs... \n");
  // gSystem->Exec("/home/user/cernbox/ALICE/offlineclass/combinepdfs.sh");

  return EXIT_SUCCESS;
}


// Global Style Setup
void SetROOTStyle() {

  gStyle->SetOptStat(0); // Statistics BOX OFF with argument 0
  gStyle->SetTitleSize(0.0475,"t"); // Title with "t" (or anything else than xyz)
  gStyle->SetStatY(1.0);
  gStyle->SetStatX(1.0);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.09);
}


// Set "nice" 2D-plot style
// Read here more about problems with the Rainbow
void SetPlotStyle() {

  // Set Smooth color gradients
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

  // Black-Red palette
  //gStyle->SetPalette(56); // 53 for inverted

  // Number of decimals in text in TH2 plots
  gStyle->SetPaintTextFormat("4.2f");

  gStyle->SetTitleOffset(1.6,"x");  //X-axis title offset from axis
  gStyle->SetTitleOffset(1.6,"y");  //X-axis title offset from axis
  gStyle->SetTitleSize(0.03,"x");   //X-axis title size
  gStyle->SetTitleSize(0.03,"y");
  gStyle->SetTitleSize(0.03,"z");
  gStyle->SetLabelOffset(0.025);

}
