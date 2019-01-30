// GENERAL COMMENT; Member variables needs to be (are) initialized per event in functions below,
// the class constructor does it only once!!
//
// mikael.mieskolainen@cern.ch, 2018
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.
//
// This class is based on ALICE collaboration skeleton code

#include <memory>
#include <iostream>

// ROOT
#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <TParticle.h>

// AliROOT
#include "AliLog.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliESDVertex.h"
#include "AliStack.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliRawEventHeaderBase.h"
#include "AliESDVZERO.h"
#include "AliESDAD.h"
#include "AliESDZDC.h"

#include "AliAnalysisTaskDiffCrossSectionsMM.h"

ClassImp(AliAnalysisTaskDiffCrossSectionsMM)
ClassImp(AliAnalysisTaskDiffCrossSectionsMM::TreeData)
ClassImp(AliAnalysisTaskDiffCrossSectionsMM::MCInfo)


// These are constant, so being global is just fine
const Double_t proton_mass = 0.938272; // Proton mass
const Int_t pdgDiff_p = 9902210;       // Diffractive proton (pseudosystem) Pythia6/8 PDG code
const Int_t pdgProton = 2212;          // Proton PDG code

// *********************************************************************
// Fiducial measurement definitions

// Generator level (=RIVET level) fiducial eta (From ALICE performance paper & AD note)
// For SPD, we take the inner (0) layer acceptance
const Double_t ZDNC_eta[2] = {-1e9, -8.8};
const Double_t ADC_eta[2]  = {-7.0, -4.9};
const Double_t V0C_eta[2]  = {-3.7, -1.7};
const Double_t SPDC_eta[2] = {-2.0, 0.0};
const Double_t SPDA_eta[2] = {0.0, 2.0};
const Double_t V0A_eta[2]  = {2.8, 5.1};
const Double_t ADA_eta[2]  = {4.8, 6.3};
const Double_t ZDNA_eta[2] = {8.8, 1e9};

// Generator level (=RIVET level) minimum pT cutoff in GeV
// These are based on a separate simulation.
const Double_t ZDNC_minPt  = 0.0;
const Double_t ADC_minPt   = 0.0;
const Double_t V0C_minPt   = 0.0;
const Double_t SPDC_minPt  = 0.0;
const Double_t SPDA_minPt  = 0.0;
const Double_t V0A_minPt   = 0.0;
const Double_t ADA_minPt   = 0.0;
const Double_t ZDNA_minPt  = 0.0;
// *********************************************************************

// -----------------------------------------------------------------------
// MC Generator level information collection
void AliAnalysisTaskDiffCrossSectionsMM::MCInfo::Fill(const AliMCEvent* mcEvent, TString mcType) {
  fEventType = kInvalid;
  if (!mcEvent)
    AliFatal("NULL == mcEvent");

  const AliGenEventHeader* h(mcEvent->GenEventHeader());
  if (!h)
    AliFatal("NULL == h");


  // Init diffractive system(s) object (IMPORTANT)
  DiffSystem empty;
  fDiffSys = empty;


  if (mcType.Contains("Pythia")) {
    const AliGenPythiaEventHeader* ph = dynamic_cast<const AliGenPythiaEventHeader* >(h);
    if (!ph)
      AliFatal("NULL == ph");
    switch (ph->ProcessType()) {
    case 92:
      fEventType = kSDR;
      break;
    case 93:
      fEventType = kSDL;
      break;
    case 94:
      fEventType = kDD;
      break;
    case 91:
      fEventType = kElastic;
      break;
    default:
      fEventType = kND;
   }
  }
  if (mcType.Contains("Pythia8")) {
    const AliGenPythiaEventHeader* ph = dynamic_cast<const AliGenPythiaEventHeader* >(h);
    if (!ph)
      AliFatal("NULL == ph");
    switch (ph->ProcessType()) {
    case 103:
      fEventType = kSDR;
      break;
    case 104:
      fEventType = kSDL;
      break;
    case 105:
      fEventType = kDD;
      break;
    case 106:
      fEventType = kCD;
      break;
    case 102:
      fEventType = kElastic;
      break;
    default:
      fEventType = kND;
    }
  }
  if (mcType.Contains("PHOJET")) {
    const AliGenDPMjetEventHeader* ph = dynamic_cast<const AliGenDPMjetEventHeader* >(h);
    if (!ph)
      AliFatal("NULL == ph");
    switch (ph->ProcessType()) {
    case  5: 
      fEventType = kSDR;
      break;
    case  6: 
      fEventType = kSDL;
      break;
    case  7: 
      fEventType = kDD;
      break;
    case  4: 
      fEventType = kCD;
      break;
    case  2: 
      fEventType = kElastic;
      break;
    default:
      fEventType = kND;
    }    
  }

  //printf("Process type = %0.0f \n", fEventType);

  // Get AliStack of particles
  AliStack* stack  = dynamic_cast<AliStack*> (const_cast<AliMCEvent*>(mcEvent)->Stack());
  if (!stack)
    AliFatal("NULL == stack");

  // Check Energy-Momentun conservation
  EMCheck(stack);

  // Find out the initial state protons
  TLorentzVector p_plus(0,0,0,0);  // Proton on positive direction
  TLorentzVector p_minus(0,0,0,0); // negative direction
  CalcInitialState(p_plus, p_minus, stack);


  if (mcType.Contains("Pythia") || mcType.Contains("Pythia8")) {
    if (fEventType == kSDL || fEventType == kSDR || fEventType == kDD){
      CalcPseudoSystem(stack); // Pseudoparticle method 
      CalcSDslashDD(p_plus, p_minus, stack); // Mother particle backtrack method
    }
  }

  // Phojet in AliROOT does not have full mother-daughter tree information,
  // at least turned on and saved by default thus different methods..

  if (mcType.Contains("PHOJET")) {
    if (fEventType == kSDL || fEventType == kSDR) {
      CalcSD(p_plus, p_minus, stack); // Forward proton method
    }
    if (fEventType == kDD) {
      CalcAssignmentSDslashDD(p_plus, p_minus, stack); // Rapidity ordering method
    }
  }
  
  // Find pseudorapidity gaps
  FindEtaGaps(stack);
  
  //printf("\n");
}


void AliAnalysisTaskDiffCrossSectionsMM::MCInfo::CalcInitialState(TLorentzVector& p_plus, TLorentzVector& p_minus, AliStack* stack) {

  Double_t y_max_initial = -1e9;   // -inf!
  Double_t y_min_initial =  1e9;   // +inf!
  TLorentzVector v(0,0,0,0);       // Temp 4-vector

  // Loop over particle stack
  for (Int_t i = 0, n = stack->GetNprimary(); i < n; ++i) {

    const TParticle* p = const_cast<AliStack*>(stack)->Particle(i);
    if (!p) // Not a valid particle
      continue;

    // Choose initial state protons by choosing the protons with the largest +-z-momentum
    if (p->GetPdgCode() == pdgProton && p->GetStatusCode() != 1) {

      // *******************************************************************
      ReCalc4Vector(v, p); // Re-calculate 4-momentum (due to floating point precision losses in AliROOT)
      // *******************************************************************

      if (v.Pz() > y_max_initial) {
        y_max_initial = v.Pz();
        p_plus = v;
      }
      if (v.Pz() < y_min_initial) {
        y_min_initial = v.Pz();
        p_minus = v;
      }
    }
  }

}


void AliAnalysisTaskDiffCrossSectionsMM::MCInfo::CalcPseudoSystem(AliStack* stack) {

  // ---------------------------------------------------------------------
  // Pseudoparticle method of Pythia
  // THIS IS FOR PYTHIA 6 and 8 
  // (N.B.! Pythia 8 with ALIROOT shows corrupted results. Problem with floating point numbers
  // in AliGenPythiaPlus:: class)

  TLorentzVector v(0,0,0,0);

  // Single Diffraction
  if (fEventType == kSDL || fEventType == kSDR) {

    // Loop over particle stack
    for (Int_t i = 0, n = stack->GetNprimary(); i < n; ++i) {

      const TParticle* p = const_cast<AliStack*>(stack)->Particle(i);
      if (!p) // Not a valid particle
        continue;

      // If we found the pseudoparticle
      if (p->GetPdgCode() == pdgDiff_p) {

        p->Momentum(v); // Here we do not re-calculate the 4-vector (impossible)

        if (v.Pz() < 0) { // Diffractive pseudosystem with negative z-momentum => SDL
          fDiffSys.Mass_PDG[0] = v.M();
          fDiffSys.Mass_PDG[1] = proton_mass;
        } else {          // => SDR
          fDiffSys.Mass_PDG[0] = proton_mass;
          fDiffSys.Mass_PDG[1] = v.M();
        }
        break; // No need to loop anymore, we need to find only one pseudoparticle
      }
    }
  }

  // Double Diffraction
  if (fEventType == kDD) {

    // Loop over particle stack
    for (Int_t i = 0, n = stack->GetNprimary(); i < n; ++i) {

      const TParticle* p = const_cast<AliStack*>(stack)->Particle(i);
      if (!p) // Not a valid particle
        continue;

      // If we found the pseudoparticle
      if (p->GetPdgCode() == pdgDiff_p) {

        p->Momentum(v); // Here we do not re-calculate the 4-vector (impossible)

        if (v.Pz() < 0) {    // Diffractive pseudosyst. with neg. z-momentum => left DD system
          fDiffSys.Mass_PDG[0] = v.M();
        } // Note that no else here, the difference to SD. There are two (2) pseudoparticles.
        if (v.Pz() > 0) {    // => Right DD system
          fDiffSys.Mass_PDG[1] = v.M();
        }
      }
    }
  }
  
  //printf("PDG:: Mass[0] = %0.3f, Mass[1] = %0.3f\n", fDiffSys.Mass_PDG[0], fDiffSys.Mass_PDG[1]);
}


void AliAnalysisTaskDiffCrossSectionsMM::MCInfo::EMCheck(AliStack* stack) {

  // ---------------------------------------------------------------------
  // Here we compare GetStatusCode() == 1 versus 
  //                     IsPhysicalPrimary() | IsSecondaryFromWeakDecay()
  //
  // Check energy-momentum conservation
  // Loop over particle stack

  TLorentzVector vSC_sum(0,0,0,0);
  TLorentzVector vSC(0,0,0,0);
  for (Int_t i = 0, n = stack->GetNprimary(); i < n; ++i) {

    const TParticle* p = const_cast<AliStack*>(stack)->Particle(i);
    if (!p) // Not a valid particle
      continue;

      // Accept only generator tagged final states
      if ( p->GetStatusCode() != 1 )
        continue;

      // *******************************************************************
      ReCalc4Vector(vSC, p); // Re-calculate 4-momentum (due to floating point precision losses in AliROOT)
      // *******************************************************************

      vSC_sum += vSC;
  }
  EM_check_SC = vSC_sum.M(); // should be \sqrt(s) for each event, by strict 4-momentum conservation


  TLorentzVector vPP_sum(0,0,0,0);
  TLorentzVector vPP(0,0,0,0);
  for (Int_t i = 0, n = stack->GetNprimary(); i < n; ++i) {

    const TParticle* p = const_cast<AliStack*>(stack)->Particle(i);
    if (!p) // Not a valid particle
      continue;

      // Accept "Physical Primaries" or "Secondaries from Weak Decays",
      // these should carry all energy-momentum of the final states
      if ( !(stack->IsPhysicalPrimary(i) || stack->IsSecondaryFromWeakDecay(i)) )
        continue;

      // *******************************************************************
      ReCalc4Vector(vPP, p); // Re-calculate 4-momentum (due to floating point precision losses in AliROOT)
      // *******************************************************************

      vPP_sum += vPP;
  }
  EM_check_PP = vPP_sum.M(); // should be \sqrt(s) for each event, by strict 4-momentum conservation

}

// Find pseudorapiditygap variables
void AliAnalysisTaskDiffCrossSectionsMM::MCInfo::FindEtaGaps(AliStack* stack) {


  std::vector<UInt_t> ind;           // Final state indices

  // Loop over particle stack for final states and collect the indices to a vector
  for (Int_t i = 0, n = stack->GetNprimary(); i < n; ++i) {

    //printf("Particle %d/%d \n", i, n);

    const TParticle* p = const_cast<AliStack*>(stack)->Particle(i);
    if (!p) // Not a valid particle
      continue;

    // Generator tagged final state
    if (p->GetStatusCode() != 1)
      continue;

    /*
    // AliROOT assignment
    if (!stack->IsPhysicalPrimary(i))
      continue;
    */

    ind.push_back(i);
  }

  // Calculate pseudorapidities of the final states
  std::vector<Double_t> eta_charged;
  std::vector<Double_t> eta_neutral;
  std::vector<Double_t> eta_all;  

  TLorentzVector v(0,0,0,0);

  for (UInt_t i = 0; i < ind.size(); ++i) {

    const TParticle* p = const_cast<AliStack*>(stack)->Particle( ind.at(i) ); // <- ind.at(i)

    // Get PDG id information
    TParticlePDG* pdgPart = p->GetPDG();
    Int_t pdgCode_ABS = TMath::Abs(p->GetPdgCode());
    
    // First skip neutrinos (pdg.lbl.gov/2007/reviews/montecarlorpp.pdf)
    if (pdgCode_ABS == 12 || pdgCode_ABS == 14 || pdgCode_ABS == 16)
      continue;

    // *******************************************************************
    ReCalc4Vector(v, p); // Re-calculate 4-momentum (due to floating point precision losses in AliROOT)
    // *******************************************************************
    Double_t pseudorapidity = v.Eta();
    

    if (ADC_eta[0] < pseudorapidity && pseudorapidity < ADA_eta[1]) { // AD gives our fiducial min/max edges  

      // Charged
      if (pdgPart->Charge() != 0) {
          eta_charged.push_back( pseudorapidity );
      // Neutral    
      } else {
          eta_neutral.push_back( pseudorapidity );
      }
      // All
      eta_all.push_back( pseudorapidity );
    }
  }

  // Then get pseudorapidity ordering eta_1 < eta_2 < ... < eta_n

  // ----------------------------------------------------------------
  // We take all particles with no pt cutoffs, i.e. pt_min > 0

  // CHARGED
  if (eta_charged.size() > 0) {
    GapCalculus(eta_charged, MinEtaVisibleCharged, MaxEtaVisibleCharged, MaxGapVisibleCharged);
  } else {
    MinEtaVisibleCharged =  1024;
    MaxEtaVisibleCharged = -1024;
    MaxGapVisibleCharged = -1024;
  }

  // NEUTRAL
  if (eta_neutral.size() > 0) {
    GapCalculus(eta_neutral, MinEtaVisibleNeutral, MaxEtaVisibleNeutral, MaxGapVisibleNeutral);
  } else {
    MinEtaVisibleNeutral =  1024;
    MaxEtaVisibleNeutral = -1024;
    MaxGapVisibleNeutral = -1024;
  }

  // ALL
  if (eta_all.size() > 0) {
    Float_t nullA = 0;
    Float_t nullB = 0;
    GapCalculus(eta_all, nullA, nullB, MaxGapVisibleAll);
  } else {
    MaxGapVisibleAll = -1024;
  }


  //printf("MinEtaVisibleCharged = %0.2f, MaxEtaVisibleCharged = %0.2f, MaxGapVisibleCharged = %0.2f \n", MinEtaVisibleCharged, MaxEtaVisibleCharged, MaxGapVisibleCharged);

}


// Do the gap counting. Input as the vector of eta values
void AliAnalysisTaskDiffCrossSectionsMM::MCInfo::GapCalculus(const std::vector<Double_t>& eta, Float_t& MinEta, Float_t& MaxEta, Float_t& MaxGap) {

    std::vector<UInt_t> eta_order = vsortind(eta); // ind vector indices now in the pseudorapidity ordering

    // TYPE I: Gap-Boundary conditions fixed edge
    MinEta = eta.at(eta_order.at(0));
    MaxEta = eta.at(eta_order.at(eta.size()-1));

    // TYPE II: Gap-Boundary conditions floating (maximum absolute difference between two particle over pseudorapidity)
    Double_t max = -1024;
    for (UInt_t i = 0; i < eta.size() - 1; ++i) { // - 1 guarantees that we have at least two particles
      Double_t diff = std::abs( eta.at(eta_order.at(i)) - eta.at(eta_order.at(i+1)) );
      max = (diff > max) ? diff : max;
    }
    MaxGap = max;
}


// CORRECT the energy component due to floating point precision dilution (for example in AliGenPythiaPlus)
void AliAnalysisTaskDiffCrossSectionsMM::MCInfo::ReCalc4Vector(TLorentzVector& v, const TParticle* p) {

  Double_t mass = 0;
  TParticlePDG* pdgPart = p->GetPDG();
  p->Momentum(v);
  mass = pdgPart->Mass();
  
  if (mass < 0 || mass > 200) { // Something funny going on
    AliWarning(Form("TParticlePDG with PDG-ID = %d gives mass = %0.3f GeV! \n", pdgPart->PdgCode(), mass));
    mass = 0; // Set to zero
  }

  v.SetXYZM(v.Px(), v.Py(), v.Pz(), mass); // Re-calculate
}


void AliAnalysisTaskDiffCrossSectionsMM::MCInfo::CalcSD(const TLorentzVector& p_plus, const TLorentzVector& p_minus, AliStack* stack) {

  // ---------------------------------------------------------------------
  // THIS IS FOR PYTHIA 6, 8, PHOJET
  // N.B. The selection of initial state protons is based on momentum z-component
  // 
  // This routine is verified with Pythia 6 to give identical results to pseudoparticle method
  //
  // SD system mass by finding the forward proton and utilizing 4-momentum conservation

  const Int_t sign = (fEventType == kSDL ? 1 : -1); // The surviving proton momentum sign

  TLorentzVector v(0,0,0,0);       // temp
  Double_t y_abs_max = 0;          // zero!
  Double_t y = 0;

  // Final state variables
  TLorentzVector f_proton;         // Forward proton

  // Loop over particle stack
  for (Int_t i = 0, n = stack->GetNprimary(); i < n; ++i) {

    const TParticle* p = const_cast<AliStack*>(stack)->Particle(i);
    if (!p) // Not a valid particle
      continue;

    // From here, accept only generator tagged final states, i.e., reject initial and intermediate states
    if (p->GetStatusCode() != 1)  
      continue;
    
    // *******************************************************************
    ReCalc4Vector(v, p); // Re-calculate 4-momentum (due to floating point precision losses in AliROOT)
    // *******************************************************************

    y = v.Pz();

    // For SDL or SDR, the surviving forward proton has the largest absolute p_z of the final states
    // and also matches the side we are searching for
    if ( (p->GetPdgCode() == pdgProton) && (TMath::Abs(y) > y_abs_max) && (TMath::Sign(1,(Int_t)y) == sign) ) {
      y_abs_max = TMath::Abs(y);
      f_proton  = v;
    }
  }

  // Calculate the diffractive system invariant mass by 4-momentum conservation: 
  // p_1 + p_2 = p_A + p_B <=> p_B = p_1 + p_2 - p_A => M_B^2 = (p_1 + p_2 - p_A)^2
  TLorentzVector sum_v(0,0,0,0);
  sum_v = p_plus + p_minus - f_proton;

  if (f_proton.Pz() > 0) {       // Surviving forward proton with positive p_z => SDL
    TLorentzVector diff = p_plus - f_proton;
    fDiffSys.t = diff.M2(); // Mandelstam t
    
    fDiffSys.Mass[0] = sum_v.M();
    fDiffSys.Mass[1] = proton_mass;
  } else {                       // => SDR
    TLorentzVector diff = p_minus - f_proton;
    fDiffSys.t = diff.M2(); // Mandelstam t

    fDiffSys.Mass[0] = proton_mass;
    fDiffSys.Mass[1] = sum_v.M();
  }

  //printf("SDD:: Mass[0] = %0.3f, Mass[1] = %0.3f \n", fDiffSys.Mass[0], fDiffSys.Mass[1]);
}

void AliAnalysisTaskDiffCrossSectionsMM::MCInfo::CalcAssignmentSDslashDD(const TLorentzVector& p_plus, const TLorentzVector& p_minus, AliStack* stack) {

  // THIS IS MAINLY FOR PHOJET
  // Single and Double diffraction by final state rapidity gap (optimality)


  TLorentzVector v(0,0,0,0);         // temp vector

  // Final state systems
  TLorentzVector f_minus(0,0,0,0);   // Dissociated system on negative z-direction
  TLorentzVector f_plus(0,0,0,0);    // Dissociated system on positive z-direction
  TLorentzVector zerovec(0,0,0,0);   // zero 4-vector

  std::vector<UInt_t> ind;           // Final state indices

  // Loop over particle stack for final states and collect the indices to a vector
  for (Int_t i = 0, n = stack->GetNprimary(); i < n; ++i) {

    //printf("Particle %d/%d \n", i, n);

    const TParticle* p = const_cast<AliStack*>(stack)->Particle(i);
    if (!p) // Not a valid particle
      continue;

    // Generator tagged final state
    if (p->GetStatusCode() != 1)
      continue;

    /*
    // AliROOT assignment
    if (!stack->IsPhysicalPrimary(i))
      continue;
    */

    ind.push_back(i);
  }

  // First calculate rapidities of the final states
  std::vector<Double_t> y(ind.size(), 0);

  for (UInt_t i = 0; i < ind.size(); ++i) {

    const TParticle* p = const_cast<AliStack*>(stack)->Particle( ind.at(i) ); // <- ind.at(i)

    // *******************************************************************
    ReCalc4Vector(v, p); // Re-calculate 4-momentum (due to floating point precision losses in AliROOT)
    // *******************************************************************

    y.at(i) = v.Rapidity();
  }

  // Then get rapidity ordering rapidity y_1 < y_2 < ... < y_n
  std::vector<UInt_t> y_order(ind.size(), 0);
  y_order = vsortind(y); // ind vector indices now in the rapidity ordering


  // Plot ordered rapidities
  for (UInt_t i = 0; i < y.size(); ++i) {
    const TParticle* p = const_cast<AliStack*>(stack)->Particle( ind.at( y_order.at(i) ) ); // <- ind.at( y_order.at(i) )

    // *******************************************************************
    ReCalc4Vector(v, p); // Re-calculate 4-momentum (due to floating point precision losses in AliROOT)
    // *******************************************************************

    //printf("Rapidity [%d] = %0.6f \n", i, v.Rapidity() );
    //cout << "Rapidity [" << i << "] = " << v.Rapidity() << endl;
  }


  // Assignment boolean vectors
  std::vector<Bool_t> x_minus(ind.size(), kFALSE);
  std::vector<Bool_t> x_minus_optimal(ind.size(), kFALSE);

  Double_t opt_sol = -1e30;

  // Do the assigment over rapidity
  for (UInt_t i = 0; i < ind.size() - 1; ++i) {

    x_minus.at(i) = kTRUE;

    // Calculate sum of 4-vectors
    f_minus.SetXYZM(0,0,0,0);
    f_plus.SetXYZM(0,0,0,0);
    for (UInt_t k = 0; k < ind.size(); ++k) {

      const TParticle* p = const_cast<AliStack*>(stack)->Particle( ind.at( y_order.at(k) ) );  // <- ind.at( y_order.at(i) )

      // *******************************************************************
      ReCalc4Vector(v, p); // Re-calculate 4-momentum (due to floating point precision losses in AliROOT)
      // *******************************************************************

      if (x_minus.at(k)) {
        f_minus += v;
      } else {
        f_plus  += v;
      }
    }

    // Calculate rapidity gap
    Double_t DY = -TMath::Log( (f_minus.M2() * f_plus.M2())  / ( pow(proton_mass,2) * pow(13000,2) ) + 1e-12 );

    // Physical constraints put here (gap must be between 0 ... -log(m_p^2 / s)))
    if (DY < 0 || DY > (-TMath::Log(pow(proton_mass,2) / pow(13000,2))) || f_minus.M() < (proton_mass-0.01) || f_plus.M() < (proton_mass-0.01) ) {
          //printf("DeltaY = %0.16f i = %d (NOT VALID) \n", DY, i);
      continue; // not valid
    }

    // Check optimality
    //printf("DeltaY = %0.16f i = %d \n", DY, i);
    if (DY > opt_sol) {
      opt_sol = DY;
      x_minus_optimal = x_minus;
    }
  }

/*
  // Plot the optimal vector
  for (UInt_t i = 0; i < ind.size(); ++i) {
    printf("X sequence [%d] = %d \n", i, (UInt_t)x_minus_optimal.at(i));
  }
*/

  // ---------------------------------------------------------------------
  // Final assignment based on the optimal solution
  f_minus.SetXYZM(0,0,0,0);
  f_plus.SetXYZM(0,0,0,0);
  for (UInt_t k = 0; k < ind.size(); ++k) {

    const TParticle* p = const_cast<AliStack*>(stack)->Particle( ind.at( y_order.at(k) ) );

    // *******************************************************************
    ReCalc4Vector(v, p); // Re-calculate 4-momentum (due to floating point precision losses in AliROOT)
    // *******************************************************************

    if (x_minus_optimal.at(k)) {
      f_minus += v;
    } else {
      f_plus  += v;
    }
  }

  // Calculate the invariant masses
  fDiffSys.Mass[0] = f_minus.M(); // Left system of DD
  fDiffSys.Mass[1] = f_plus.M();  // Right system of DD

  //printf("ASS:: Mass[0] = %0.3f, Mass[1] = %0.3f \n", fDiffSys.Mass[0], fDiffSys.Mass[1]);

  // Mandelstam t
  TLorentzVector diff = p_minus - f_minus;
  fDiffSys.t = diff.M2(); // Mandelstam t

  // Check it
  TLorentzVector diff2 = p_plus - f_plus;
  Double_t t = diff2.M2();
  if (TMath::Abs(fDiffSys.t - diff2.M2()) > 1e-6) {
  //  AliWarning(Form("Event type = %0.0f : Kinematics accuracy O(1e-6) error (t- = %0.12f vs t+ = %0.12f)", fEventType, fDiffSys.t, t));
  }

}


// Return reordered vector according to given index order
void AliAnalysisTaskDiffCrossSectionsMM::MCInfo::vsort(std::vector<Double_t>& v, const std::vector<UInt_t>& sort_ind) {

  std::vector<Double_t> reordered;

  for (UInt_t i = 0; i < v.size(); ++i) {
    reordered.push_back( v.at(sort_ind.at(i)) );
  }
  for (UInt_t i = 0; i < v.size(); ++i) {
    v.at(i) = reordered.at(i);
  }
}


// Return reordered indices in ascending order
// std::vector <Double_t> v = {666, 23, 884, 483}; // input data example
std::vector<UInt_t> AliAnalysisTaskDiffCrossSectionsMM::MCInfo::vsortind(const std::vector<Double_t>& v) {

  std::multimap <Double_t, UInt_t> m; // mapping from value to its index
  for (auto it = v.begin(); it != v.end(); ++it) {
      m.insert(std::make_pair(*it, it - v.begin()));
  }

  // reordered indices in ascending order
  std::vector<UInt_t> s;

  for (auto it = m.begin(); it != m.end(); ++it) {
    s.push_back( it->second );
  }

  return s;
}


void AliAnalysisTaskDiffCrossSectionsMM::MCInfo::CalcSDslashDD(const TLorentzVector& p_plus, const TLorentzVector& p_minus, AliStack* stack) {

  // Generic method for both Single and Double diffraction by mother particle backtracking
  //

  TLorentzVector v(0,0,0,0);         // temp vector

  // Final state systems
  TLorentzVector f_minus(0,0,0,0);   // Dissociated system on negative z-direction
  TLorentzVector f_plus(0,0,0,0);    // Dissociated system on positive z-direction

  std::vector<Double_t> PtCharged[2];
  std::vector<Double_t> PtNeutral[2];
  std::vector<Double_t> EtaCharged[2];
  std::vector<Double_t> EtaNeutral[2];

  Int_t counter = 0;                 // Backtrack counter for safety
  const Int_t btr_max = 1000;        // Safety max (if the last particle never returns -1)

  // Loop over particle stack for final states
  for (Int_t i = 0, n = stack->GetNprimary(); i < n; ++i) {

    //printf("Particle %d/%d \n", i, n);

    const TParticle* p = const_cast<AliStack*>(stack)->Particle(i);
    if (!p) // Not a valid particle
      continue;

    // From here, accept only generated tagged final states, i.e., reject initial and intermediate states
    if (p->GetStatusCode() != 1)
      continue;

    // *******************************************************************
    ReCalc4Vector(v, p); // Re-calculate 4-momentum (due to floating point precision losses in AliROOT)
    // *******************************************************************

    // Get mother particle index
    Int_t mother_ind = p->GetFirstMother();
    if (mother_ind == -1)
      continue; // This must initial state, or intermediate state

    // Now loop back to the grandmother, i.e., to the initial state proton 1 or 2
    Double_t pZ_fm = 0;
    do {
      TParticle* p_fm = (TParticle*)stack->Particle(mother_ind);
      pZ_fm = p_fm->Pz();
      mother_ind = p_fm->GetFirstMother(); // gives index -1 when we hit the initial state proton

      //printf("Mother particle stack index: [%d] with PDG code: %d and Pz = %0.3f \n", mother_ind, p_fm->GetPdgCode(), p_fm->Pz());
      ++counter;
    } while (mother_ind != -1 && (counter != btr_max));

    if (counter == btr_max) {
      AliWarning("Faulty mother-daughter index tree in AliStack!");
      continue;
    }

    // Finally do the classification
    Int_t side = 0;
    if (pZ_fm < 0) {   // The mother proton is with negative z-momentum
      f_minus += v;
      side = 0;
    } else {           // The mother proton is with positive z-momentum
      f_plus  += v;
      side = 1;
    }

    TParticlePDG* pdgPart = p->GetPDG();
    
    // Charged
    if (pdgPart->Charge() != 0) {
      PtCharged[side].push_back(v.Pt());
      EtaCharged[side].push_back(v.Eta());

    // Neutral
    } else { 
      PtNeutral[side].push_back(v.Pt());
      EtaNeutral[side].push_back(v.Eta());
    }
    counter = 0;
  }
  
  // Calculate multiplicity, min/max eta & mean pT values
  for (UInt_t i = 0; i < 2; ++i) {

    fDiffSys.NCharged[i]      = PtCharged[i].size();
    fDiffSys.NNeutral[i]      = PtNeutral[i].size();

    fDiffSys.MeanPtCharged[i] = VectorMean(PtCharged[i]);
    fDiffSys.MeanPtNeutral[i] = VectorMean(PtNeutral[i]);

    fDiffSys.MinEtaCharged[i] = VectorMin(EtaCharged[i]);
    fDiffSys.MinEtaNeutral[i] = VectorMin(EtaNeutral[i]);

    fDiffSys.MaxEtaCharged[i] = VectorMax(EtaCharged[i]);
    fDiffSys.MaxEtaNeutral[i] = VectorMax(EtaNeutral[i]);
  }

  // Calculate the invariant masses
  fDiffSys.Mass[0] = f_minus.M(); // Left system of DD
  fDiffSys.Mass[1] = f_plus.M();  // Right system of DD

  //printf("MOT:: Mass[0] = %0.3f, Mass[1] = %0.3f \n", fDiffSys.Mass[0], fDiffSys.Mass[1]);


  // Mandelstam t
  TLorentzVector diff = p_minus - f_minus;
  fDiffSys.t = diff.M2(); // Mandelstam t

  // Now use the Mandelstam t as a numerical accuracy test of kinematics

  // Cross check it (THIS COULD BE AN ALGORITHMIC BASIS FOR DD SELECTION IN MC)
  TLorentzVector diff2 = p_plus - f_plus;

  Double_t order = TMath::Abs( (TMath::Abs(diff.M2()) - TMath::Abs(diff2.M2()) )) / TMath::Abs(diff.M2() + 1e-12);
  if (order > 1e-2) { // Threshold
    AliWarning(Form("Event type = %0.0f, Kinematics accuracy O(%0.0e) warning (t- = %0.7f vs t+ = %0.7f)", fEventType, order, diff.M2(), diff2.M2()));
  }

}

// Calculate mean of a std::vector
Double_t AliAnalysisTaskDiffCrossSectionsMM::MCInfo::VectorMean(const std::vector<Double_t>& vec) {

  Double_t sum = 0;
  if (vec.size() == 0) 
    return -1e9;

  for (UInt_t i = 0; i < vec.size(); ++i) {
    sum += vec.at(i);
  }
  return sum / ( static_cast<Double_t>(vec.size()) );
}

// Calculate the minimum value of a std::vector
Double_t AliAnalysisTaskDiffCrossSectionsMM::MCInfo::VectorMin(const std::vector<Double_t>& vec) {

  Double_t min = 1e9; // inf

  for (UInt_t i = 0; i < vec.size(); ++i) {
    if (vec.at(i) < min) {
      min = vec.at(i);
    }
  }
  return min;
}

// Calculate the maximimum value of a std::vector
Double_t AliAnalysisTaskDiffCrossSectionsMM::MCInfo::VectorMax(const std::vector<Double_t>& vec) {
  
  Double_t max = -1e9; // -inf

  for (UInt_t i = 0; i < vec.size(); ++i) {
    if (vec.at(i) > max) {
      max = vec.at(i);
    }
  }
  return max;
}

// -----------------------------------------------------------------------
// MC Fiducial Combinatorial truth for unfolding and control

void AliAnalysisTaskDiffCrossSectionsMM::MCInfo::FillComb(const AliMCEvent* mcEvent) {

  if (!mcEvent)
    AliFatal("NULL == mcEvent");

  // Get AliStack of particles
  AliStack* stack  = dynamic_cast<AliStack*> (const_cast<AliMCEvent*>(mcEvent)->Stack());
  if (!stack)
    AliFatal("NULL == stack");

  // These are initialized here by first creating a new object, 
  // which is initialized properly by the constructor of DetGenLevel::
  DetGenLevel zerostruct;

  ZDNC = zerostruct;
  ADC  = zerostruct;
  V0C  = zerostruct;
  SPDC = zerostruct;
  SPDA = zerostruct;
  V0A  = zerostruct;
  ADA  = zerostruct;
  ZDNA = zerostruct;

  // Loop over particle stack for final states
  for (Int_t i = 0, n = stack->GetNprimary(); i < n; ++i) {

    const TParticle* p = const_cast<AliStack*>(stack)->Particle(i);
    if (!p) // Not a valid particle
      continue;

    // Generator level final state
    if (p->GetStatusCode() != 1)
      continue;

/*
    // AliROOT assignment
    if ( !(stack->IsPhysicalPrimary(i) || stack->IsSecondaryFromWeakDecay(i)) )
      continue;
*/

    // Get PDG id information
    Int_t pdgCode_ABS = TMath::Abs(p->GetPdgCode());
    
    // First skip neutrinos (pdg.lbl.gov/2007/reviews/montecarlorpp.pdf)
    if (pdgCode_ABS == 12 || pdgCode_ABS == 14 || pdgCode_ABS == 16)
      continue;
    
    // Now check if the particle is inside the fiducial acceptance region ->
    CheckFiducial(p, ZDNC, ZDNC_eta, ZDNC_minPt);
    CheckFiducial(p, ADC, ADC_eta, ADC_minPt);
    CheckFiducial(p, V0C, V0C_eta, V0C_minPt);
    CheckFiducial(p, SPDC, SPDC_eta, SPDC_minPt);
    CheckFiducial(p, SPDA, SPDA_eta, SPDA_minPt);
    CheckFiducial(p, V0A, V0A_eta, V0A_minPt);
    CheckFiducial(p, ADA, ADA_eta, ADA_minPt);
    CheckFiducial(p, ZDNA, ZDNA_eta, ZDNA_minPt);
  }
}

// -----------------------------------------------------------------------
// Check if the particle belongs to the fiducial region
void AliAnalysisTaskDiffCrossSectionsMM::MCInfo::CheckFiducial(const TParticle* p, DetGenLevel& det, const Double_t eta_range[], const Double_t minPt) {


  TParticlePDG* pdgPart = p->GetPDG();

  // Check fiducial eta window
  if ( !((eta_range[0] < p->Eta()) && (p->Eta() < eta_range[1])) )
    return;

  // Check pt threshold
  if (p->Pt() < minPt)
    return;

  // Charged
  if (pdgPart->Charge() != 0) {

      // The previous numbers
      Double_t N_prev      = (Double_t) det.NCharged;

      // New mean, min and max
      ++det.NCharged;
      det.MeanPtCharged  = (det.MeanPtCharged * N_prev + p->Pt()) / (N_prev + 1);
      det.MinPtCharged   = p->Pt() < det.MinPtCharged ? p->Pt() : det.MinPtCharged;
      det.MaxPtCharged   = p->Pt() > det.MaxPtCharged ? p->Pt() : det.MaxPtCharged;

  // Neutral
  } else {

      // The previous numbers
      Double_t N_prev      = (Double_t) det.NNeutral;

      // New mean, min and max
      ++det.NNeutral;
      det.MeanPtNeutral  = (det.MeanPtNeutral * N_prev + p->Pt()) / (N_prev + 1);
      det.MinPtNeutral   = p->Pt() < det.MinPtNeutral ? p->Pt() : det.MinPtNeutral;
      det.MaxPtNeutral   = p->Pt() > det.MaxPtNeutral ? p->Pt() : det.MaxPtNeutral;
  }
}


// -----------------------------------------------------------------------
// Event information fill
void AliAnalysisTaskDiffCrossSectionsMM::EventInfo::Fill(const AliESDEvent* esdEvent) {
  const AliESDHeader* esdHeader = esdEvent->GetHeader();
  if (NULL == esdHeader) // this is already dealt with in UserExec
    return;

  fClassMask       = esdHeader->GetTriggerMask();
  fClassMaskNext50 = esdHeader->GetTriggerMaskNext50();
  
  fRunNumber       = esdEvent->GetRunNumber();
  
  fL0Inputs        = esdHeader->GetL0TriggerInputs();
  fL1Inputs        = esdHeader->GetL1TriggerInputs();
  fL2Inputs        = esdHeader->GetL2TriggerInputs();
  
  fBCID            = esdHeader->GetBunchCrossNumber();
  fOrbitID         = esdHeader->GetOrbitNumber();
  fPeriod          = esdHeader->GetPeriodNumber();
  fTimeStamp       = esdHeader->GetTimeStamp();
}

// -----------------------------------------------------------------------
// Vertex information fill
void AliAnalysisTaskDiffCrossSectionsMM::VtxInfo::Fill(const AliESDVertex* vtx) {
  if (vtx) {
    fZ      = vtx->GetZ();
    fNcontr = vtx->GetNContributors();
    
  } else {
    fZ      =  0.0f;
    fNcontr = -4;
    AliErrorClass("NULL == vtx");
  }
}

// -----------------------------------------------------------------------
// AD/V0 Fill AD
void AliAnalysisTaskDiffCrossSectionsMM::ADV0::FillAD(const AliESDEvent* esdEvent, AliTriggerAnalysis& trigAna) {

  // Trigger analysis using AliTriggerAnalysis:: class
  fDecisionOnline[0]  = trigAna.ADTrigger(esdEvent, AliTriggerAnalysis::kCSide, kFALSE);
  fDecisionOnline[1]  = trigAna.ADTrigger(esdEvent, AliTriggerAnalysis::kASide, kFALSE);

  fDecisionOffline[0] = trigAna.ADTrigger(esdEvent, AliTriggerAnalysis::kCSide, kTRUE);
  fDecisionOffline[1] = trigAna.ADTrigger(esdEvent, AliTriggerAnalysis::kASide, kTRUE);

  AliESDAD* esdAD = esdEvent->GetADData();
  if (NULL == esdAD) {
    FillInvalid();
    return;
  }

  // Time
  fTime[0] = esdAD->GetADCTime();
  fTime[1] = esdAD->GetADATime();

  // Beam-Beam/Gas flags
  fBB[0] = fBB[1] = fBG[0] = fBG[1] = 0;
  for (Int_t ch = 0; ch < 4; ++ch) {
    fBB[0] += (esdAD->GetBBFlag(ch  ) && esdAD->GetBBFlag(ch+ 4));
    fBB[1] += (esdAD->GetBBFlag(ch+8) && esdAD->GetBBFlag(ch+12));
    fBG[0] += (esdAD->GetBGFlag(ch  ) && esdAD->GetBGFlag(ch+ 4));
    fBG[1] += (esdAD->GetBGFlag(ch+8) && esdAD->GetBGFlag(ch+12));
  }
  
  // Sum of charges
  fCharge[0] = fCharge[1] = 0.0;
  for (Int_t ch = 0; ch < 16; ++ch)
    fCharge[ch/8] += esdAD->GetMultiplicity(ch);
}

// -----------------------------------------------------------------------
// AD/V0 Fill V0
void AliAnalysisTaskDiffCrossSectionsMM::ADV0::FillV0(const AliESDEvent* esdEvent, AliTriggerAnalysis& trigAna) {

  // Trigger analysis using AliTriggerAnalysis:: class
  fDecisionOnline[0]  = trigAna.V0Trigger(esdEvent, AliTriggerAnalysis::kCSide, kFALSE);
  fDecisionOnline[1]  = trigAna.V0Trigger(esdEvent, AliTriggerAnalysis::kASide, kFALSE);

  fDecisionOffline[0] = trigAna.V0Trigger(esdEvent, AliTriggerAnalysis::kCSide, kTRUE);
  fDecisionOffline[1] = trigAna.V0Trigger(esdEvent, AliTriggerAnalysis::kASide, kTRUE);

  AliESDVZERO* esdV0 = esdEvent->GetVZEROData();
  if (NULL == esdV0) {
    FillInvalid();
    return;
  }

  // Time
  fTime[0] = esdV0->GetV0CTime();
  fTime[1] = esdV0->GetV0ATime();

  // Beam-Beam(Gas) flags
  fBB[0] = fBB[1] = fBG[0] = fBG[1] = 0;
  for (Int_t ch = 0; ch < 64; ++ch) {
    fBB[ch/32] += esdV0->GetBBFlag(ch);
    fBG[ch/32] += esdV0->GetBGFlag(ch);
  }

  // Sum of charges
  fCharge[0] = fCharge[1] = 0.0;
  for (Int_t ch=0; ch<64; ++ch)
    fCharge[ch/32] += esdV0->GetMultiplicity(ch);
}

// AD/V0 Fill invalid
void AliAnalysisTaskDiffCrossSectionsMM::ADV0::FillInvalid() {
  fTime[0] = fTime[1] = -10240.0f;
  fCharge[0] = fCharge[1] = fBB[0] = fBG[0] = fBB[1] = fBG[1] = -1.0f;
}

// -----------------------------------------------------------------------
// ZDC information
void AliAnalysisTaskDiffCrossSectionsMM::ZDC::FillZDC(const AliESDEvent* esdEvent, AliTriggerAnalysis& trigAna) {
  
  AliESDZDC* esdZDC = esdEvent->GetESDZDC();
  if (NULL == esdZDC) {
    FillInvalid();
    return;
  }

  // Neutron hit
  fZNHit[0] = esdZDC->IsZNChit();
  fZNHit[1] = esdZDC->IsZNAhit();

  // Proton hit
  fZPHit[0] = esdZDC->IsZPChit();
  fZPHit[1] = esdZDC->IsZPAhit();

  // Neutron energy
  fZNEnergy[0] = esdZDC->GetZNCEnergy();
  fZNEnergy[1] = esdZDC->GetZNAEnergy();
  
  // Proton energy
  fZPEnergy[0] = esdZDC->GetZPCEnergy();
  fZPEnergy[1] = esdZDC->GetZPAEnergy();
  
}

// ZDC Fill invalid
void AliAnalysisTaskDiffCrossSectionsMM::ZDC::FillInvalid() {
  fZNHit[0] = fZNHit[1] = fZPHit[0] = fZPHit[1] = kFALSE;
  fZNEnergy[0] = fZNEnergy[1] = fZPEnergy[0] = fZPEnergy[1] = -1.0f;
}


// -----------------------------------------------------------------------
// Constructor
AliAnalysisTaskDiffCrossSectionsMM::AliAnalysisTaskDiffCrossSectionsMM(const char* name)
  : AliAnalysisTaskSE(name)
  , fIsMC(kFALSE)
  , fMCType("")
  , fTriggerSelection("")
  , fTriggerAnalysis()
  , fAnalysisUtils()
  , fTE(NULL)
  , fFastOrMap()
  , fFiredChipMap()
  , fVertexSPD()
  , fTreeData()
  , fMCInfo()
{  
  fTriggerAnalysis.SetAnalyzeMC(fIsMC);

  DefineOutput(1, TTree::Class());
}

// -----------------------------------------------------------------------
// Destructor
AliAnalysisTaskDiffCrossSectionsMM::~AliAnalysisTaskDiffCrossSectionsMM() {
  const AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
  if (NULL != man && man->GetAnalysisType() == AliAnalysisManager::kProofAnalysis)
    return;

  if (NULL != fTE)
    delete fTE;
  fTE = NULL;
}

// -----------------------------------------------------------------------
// Set branches of ROOT tree
void AliAnalysisTaskDiffCrossSectionsMM::SetBranches(TTree* t) {
  t->Branch("TreeData",     &fTreeData);
  t->Branch("FastOrMap",    &fFastOrMap,    32000, 0);
  t->Branch("FiredChipMap", &fFiredChipMap, 32000, 0);
  t->Branch("VertexSPD",    &fVertexSPD,    32000, 0);
  if (fIsMC)
    t->Branch("MCInfo", &fMCInfo);
}

// -----------------------------------------------------------------------
// Create the objects of the analysis
void AliAnalysisTaskDiffCrossSectionsMM::UserCreateOutputObjects() {
  TDirectory* owd = gDirectory;
  OpenFile(1);
  fTE = new TTree;
  fTE->SetName(GetTreeName());
  SetBranches(fTE);
  PostData(1, fTE);
  owd->cd();
}

// -----------------------------------------------------------------------
void AliAnalysisTaskDiffCrossSectionsMM::NotifyRun() {
}

// -----------------------------------------------------------------------
// The event loop
void AliAnalysisTaskDiffCrossSectionsMM::UserExec(Option_t*) {

  AliVEvent* event = InputEvent();
  if (NULL == event) {
    AliFatal("NULL == event");
    return;
  }

  AliESDEvent* esdEvent = dynamic_cast<AliESDEvent*>(InputEvent());
  if (NULL == esdEvent) {
    AliFatal("NULL == esdEvent");
    return;
  }
  
  const AliAnalysisManager* man(AliAnalysisManager::GetAnalysisManager());
  if (NULL == man) {
    AliFatal("NULL == man");
    return;
  }

  AliESDInputHandler* inputHandler(dynamic_cast<AliESDInputHandler*>(man->GetInputEventHandler()));  
  if (NULL == inputHandler) {
    AliFatal("NULL == inputHandler");
    return;
  }

  const AliESDHeader* esdHeader = esdEvent->GetHeader();
  if (NULL == esdHeader) {
    AliFatal("NULL == esdHeader");
    return;
  }

  if (kFALSE == fIsMC && esdHeader->GetEventType() != AliRawEventHeaderBase::kPhysicsEvent)
    return;

  // Get event multiplicity information
  const AliMultiplicity* mult = esdEvent->GetMultiplicity();
  if (NULL == mult) {
    AliFatal("NULL == mult");
    return;
  }

  if (!fIsMC) {
    std::unique_ptr<const TObjArray> split(fTriggerSelection.Tokenize("|"));
    Bool_t selected = kFALSE;
    for (Int_t i = 0, n = split->GetEntries(); i < n && !selected; ++i)
      selected = esdEvent->GetFiredTriggerClasses().Contains(split->At(i)->GetName());
    
    AliInfo(Form("selected: %d %s", selected, esdEvent->GetFiredTriggerClasses().Data()));
    if (!selected)
      return;
  }
  
  // Fill event data
  fTreeData.fEventInfo.Fill(esdEvent);

  // Fill aux information
  fTreeData.fPhysSelBits              = inputHandler->IsEventSelected();
  fTreeData.fIsIncompleteDAQ          = esdEvent->IsIncompleteDAQ();
  fTreeData.fIsSPDClusterVsTrackletBG = fAnalysisUtils.IsSPDClusterVsTrackletBG(esdEvent);
  fTreeData.fVtxInfo.Fill(esdEvent->GetPrimaryVertexSPD());

  // Fill V0, AD and ZDC information
  fTreeData.fV0Info.FillV0(esdEvent, fTriggerAnalysis);
  fTreeData.fADInfo.FillAD(esdEvent, fTriggerAnalysis);
  fTreeData.fZDCInfo.FillZDC(esdEvent, fTriggerAnalysis);
  
  // Fill SPD FastOR and FiredChip maps
  fFastOrMap    = mult->GetFastOrFiredChips();
  fFiredChipMap = mult->GetFiredChipMap();

  // Fill SPD tracklets and clusters per layer
  fTreeData.fEventInfo.fnTrklet         = mult->GetNumberOfTracklets();
  fTreeData.fEventInfo.fnSPDClusters[0] = mult->GetNumberOfITSClusters(0);
  fTreeData.fEventInfo.fnSPDClusters[1] = mult->GetNumberOfITSClusters(1);

  // MC ONLY
  if (fIsMC) {
    // Fill Regge style info
    fMCInfo.Fill(MCEvent(), fMCType);

    // Fill combinatorics (Rivet style)
    fMCInfo.FillComb(MCEvent());
  }

  // Final step
  fTE->Fill();
  PostData(1, fTE);
}

// -----------------------------------------------------------------------
// Called at the end
void AliAnalysisTaskDiffCrossSectionsMM::Terminate(Option_t*) {
  fTE  = dynamic_cast<TTree*>(GetOutputData(1));
  if (NULL == fTE)
    Error("Terminate","fTE is not available");
}

