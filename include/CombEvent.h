// Event object class for Combinatorial Diffractive Cross Sections
//
// mikael.mieskolainen@cern.ch, 2018
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.


#ifndef COMBEVENT_CLASS
#define COMBEVENT_CLASS

// C++ headers
#include <iostream>
#include <vector>

// ROOT headers
#include "TFile.h"
#include "TTree.h"

// Own headers
#include "AliAnalysisTaskDiffCrossSectionsMM.h"


// Triggerdata object
struct TriggerData {

	std::string Name;
	std::string BCMask;
	double LMa      = 0.0; // Counts
	double LMb      = 0.0; // Counts
	double L0a      = 0.0;
	double L0b      = 0.0;
	double L0aL0b   = 0.0; // L0a / L0b ratio
	double L0bLMb   = 0.0; // L0b / LMb ratio

	void Print() {
		printf("TriggerData::Print::\n .Name = %s, .BCMask = %s, .LMb = %0.0f, .LMa = %0.0f, .L0b = %0.0f, .L0a = %0.0f, .L0aL0b = %0.2f, .L0bLMb = %0.3E\n", 
			Name.c_str(), BCMask.c_str(), LMb, LMa, L0b, L0a, L0aL0b, L0bLMb);
	}
};


class CombEvent {

public:
	
	// Constructors and destructor
    CombEvent(AliAnalysisTaskDiffCrossSectionsMM::TreeData* fTreeData, TBits* fFastOrMap, TBits* fFiredChipMap, AliESDVertex* fVertexSPD);
    CombEvent(AliAnalysisTaskDiffCrossSectionsMM::TreeData* fTreeData, TBits* fFastOrMap, TBits* fFiredChipMap, AliESDVertex* fVertexSPD, AliAnalysisTaskDiffCrossSectionsMM::MCInfo* fMCInfo);
	CombEvent();
	~CombEvent();

	// References for efficient object access (they are const, not able to modify data)
	const AliAnalysisTaskDiffCrossSectionsMM::TreeData* fTreeData() const;
	const AliAnalysisTaskDiffCrossSectionsMM::MCInfo* fMCInfo() const;

	const TBits* fFastOrMap() const;
	const TBits* fFiredChipMap() const;
	const AliESDVertex* fVertexSPD() const;

	void Initialize();

	// Get event information
	Int_t runNumber() const;

private:

	AliAnalysisTaskDiffCrossSectionsMM::TreeData* fTreeData_;
	AliAnalysisTaskDiffCrossSectionsMM::MCInfo* fMCInfo_;

	TBits* fFastOrMap_;
	TBits* fFiredChipMap_;
	AliESDVertex* fVertexSPD_;

	// ---------------------------------------------------------------
	// No Copy Constructor (right now, no need to copy events)
	// - by making privates with no implementation below
	//CombEvent(const CombEvent& foo);
	//CombEvent& operator=(const CombEvent& foo);

	//ClassDef(CombEvent, 1);       // ROOT system integration
};

#endif
