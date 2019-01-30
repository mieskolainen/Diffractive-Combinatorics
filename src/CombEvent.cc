// Event object class for Combinatorial Diffractive Cross Sections
//
// mikael.mieskolainen@cern.ch, 2018
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.


// AliROOT headers
//#include "AliPID.h"
//#include "AliPIDResponse.h"
//#include "AliESDVertex.h"

// ROOT headers
#include "TArrayD.h"
#include "TH2F.h"
#include "THStack.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TTree.h"

// Own headers
#include "AliAnalysisTaskDiffCrossSectionsMM.h"
#include "CombEvent.h"


//ClassImp(CombEvent) // ROOT

// Constructor for DATA
CombEvent::CombEvent(AliAnalysisTaskDiffCrossSectionsMM::TreeData* fTreeData, TBits* fFastOrMap, TBits* fFiredChipMap, AliESDVertex* fVertexSPD) {
	//Initialize();

	// Create copies of objects pointer by pointers
	fTreeData_     = fTreeData;
	fFastOrMap_    = fFastOrMap;
	fFiredChipMap_ = fFiredChipMap;
	fVertexSPD_    = fVertexSPD;

}

// Constructor for MC
CombEvent::CombEvent(AliAnalysisTaskDiffCrossSectionsMM::TreeData* fTreeData, TBits* fFastOrMap, 
	TBits* fFiredChipMap, AliESDVertex* fVertexSPD, AliAnalysisTaskDiffCrossSectionsMM::MCInfo* fMCInfo) {
	//Initialize();
	
	// Take the pointers
	fTreeData_     = fTreeData;
	fFastOrMap_    = fFastOrMap;
	fFiredChipMap_ = fFiredChipMap;
	fVertexSPD_    = fVertexSPD;

	fMCInfo_       = fMCInfo;

}

CombEvent::CombEvent() {
	//Initialize();
}

void CombEvent::Initialize() {

}

// Destructor
CombEvent::~CombEvent() {

}

// Return references for efficient access
const AliAnalysisTaskDiffCrossSectionsMM::TreeData* CombEvent::fTreeData() const { 
	return fTreeData_;
}
const AliAnalysisTaskDiffCrossSectionsMM::MCInfo* CombEvent::fMCInfo() const {
 	return fMCInfo_; 
}

const TBits* CombEvent::fFastOrMap() const { 
	return fFastOrMap_; 
}
const TBits* CombEvent::fFiredChipMap() const { 
	return fFiredChipMap_; 
}
const AliESDVertex* CombEvent::fVertexSPD() const { 
	return fVertexSPD_; 
}

// Run information
Int_t CombEvent::runNumber() const {
	return fTreeData_->fEventInfo.fRunNumber;
}

