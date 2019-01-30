// Cutflow object class
//
// mikael.mieskolainen@cern.ch, 2018
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef CUTFLOW_H
#define CUTFLOW_H

// C++
#include <iostream>
#include <vector>
#include <cmath>


class CutFlow {

public:

	// Constructor and destructor
	CutFlow(int N_events, int N_cuts, std::string input_name);
    ~CutFlow();

	// Call cut, return true if cut is pass, false if not
	bool cut(int event, int cut_index, bool cut_result);
	// Print out the linear cut flow
	void printFlow(const std::string& output_filename, const char* write_mode);

	// Print out the cut correlation matrix
	void printCorrelation(const std::string& output_filename, const char* write_mode);
	// Add name for a cut
	void nameCut(int cut_index, const std::string& str, bool status);
	
	// Get name of the cut
	std::string getName(int cut_index);
	// Get vector of event IDs which pass the cuts
	std::vector<int> getPass();

private:

	// Cut matrix [N_events x N_cuts]
	// - Cut passed     = 1
	// - Cut not passed = 0
	// - Not evaluated  = -1 (init value)
	std::vector<std::vector<int> > CMat_;
	
	// Total data sample name;
	std::string input_name_;

	// Name of the cuts
	std::vector<std::string> N_names_;

	// Status of the cut (on or off)
	std::vector<bool> status_;

};

#endif // CUTFLOW_H
