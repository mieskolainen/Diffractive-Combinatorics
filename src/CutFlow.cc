// CutFlow object class. Class evaluates and keeps track of Cut & Count, 
// and evaluates correlations between cuts. Both "serial" and "parallel" flow.
// -----------------------------------------------------------------------
//
// Example of usage:
//
//	CutFlow cuts_(nentries, 4, "TriggerXYZ-123");
//	cuts_.nameCut(0, "Online trigger", true);
//	cuts_.nameCut(1, "!Beam-Gas", true);
//	cuts_.nameCut(2, "!Satellites", true);
//	cuts_.nameCut(3, "Offline minimum conditions", true);
//
//
//  EVENT LOOP ...
//
//		bool event_OK = true;
//
//		// 1. ONLINE trigger
//		if ( !cuts_.cut(eventnr, 0, OnlineTrigger(event)) ) {
//			event_OK = false; // @@@ CUT FLOW !!
//		}
//
//		// 2. BEAM-GAS
//		if ( !cuts_.cut(eventnr, 1, BGVeto(event)) )  {
//			event_OK = false; // @@@ CUT FLOW !!
//		}
//
//      ...
//  ...
//
//	cuts_.printFlow(output_filename, "w");
//	cuts_.printCorrelation(output_filename, "a+");
//
//	std::vector<int> list = cuts_.getPass();
//
//
// TODO: make a memory efficient version of this for very large event statistics.
//
//
// mikael.mieskolainen@cern.ch, 2018
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.


// C++
#include <iostream>
#include <vector>
#include <cmath>


// Own headers
#include "CutFlow.h"


// Constructors
CutFlow::CutFlow(int N_events, int N_cuts, std::string input_name) {

	// Init cut vector and status of cuts
	const int init_value = -1;
	std::vector<std::vector<int> > temp(N_cuts, std::vector<int>(N_events, init_value) );
	std::vector<bool> temp2(N_cuts, false);

    // Init cut names
	std::vector<std::string> temp3(N_cuts);
	
	CMat_ = temp;
	status_ = temp2;
	N_names_ = temp3;

	// Input event set name
	input_name_ = input_name;
}


// Destructor
CutFlow::~CutFlow() {
}


//Call cut, return true if cut is pass, false if not
bool CutFlow::cut(int event, int cut_index, bool cut_result) {

	// Save the cut result
	CMat_.at(cut_index).at(event) = cut_result;

	// Now check if the cut was on, and save the result
	bool cut_return = false;

	if (status_.at(cut_index) == true) { // CUT is on
		cut_return = cut_result;
	} else {
		cut_return = true; 			  // <CUT off> = always pass
	}

	return cut_return;
}


// Name and set status of the cut
void CutFlow::nameCut(int cut_index, const std::string& str, bool status) {

	N_names_.at(cut_index) = str;
	status_.at(cut_index) = status;
}

// Get name of the cut
std::string CutFlow::getName(int cut_index) {

	return N_names_.at(cut_index);
}


// Print out cut correlation matrix
// corr(X,Y) = cov(X,Y) / sigma_x * sigma_y

void CutFlow::printCorrelation(const std::string& output_filename, const char* write_mode) {

	FILE* fp;
	fp = fopen(output_filename.c_str(), write_mode );

	fprintf(fp, "--------------------------------------------------------------------------\n");
	fprintf(fp, "CutFlow:: printCorrelation() <Parallel flow> results: \n");
	fprintf(fp, "--------------------------------------------------------------------------\n");

	int N_cuts = CMat_.size();
	int N_events = CMat_.at(1).size();

	fprintf(fp, "Input  : %d\t %s\n", N_events, input_name_.c_str());

	// -------------------------------------------------------------------
	// Calculate mean values

	std::vector<double> mean(N_cuts, 0.0);

	// Loop over all cuts
	for (int x = 0; x < N_cuts; ++x) {

		// Loop over all events
		for (int i = 0; i < N_events; ++i) {

			mean.at(x) += (double) CMat_.at(x).at(i);
		}
		mean.at(x) /= std::max(1.0, (double) N_events); // 1 / n  (max to avoid division by zero)
	}
	// -------------------------------------------------------------------
	fprintf(fp, "Mean:  [" );
	for (int x = 0; x < N_cuts; ++x) {
		fprintf(fp, "%0.2f ", mean.at(x));
	}
	fprintf(fp, "]\n");


	// -------------------------------------------------------------------
	// Calculate the variable standard deviations sigma_x

	std::vector<double> sigma(N_cuts, 0.0);

	// Loop over all cuts
	for (int x = 0; x < N_cuts; ++x) {

		// Loop over all events
		for (int i = 0; i < N_events; ++i) {

			sigma.at(x) += pow( CMat_.at(x).at(i) - mean.at(x), 2);
		}
		sigma.at(x) /= (double) N_events;     // 1 / n
		sigma.at(x) = std::sqrt(sigma.at(x)); // sqrt()
	}
	fprintf(fp, "Sigma: [" );
	for (int x = 0; x < N_cuts; ++x) {
		fprintf(fp, "%0.2f ", sigma.at(x));
	}
	fprintf(fp, "]\n\n");


	// -------------------------------------------------------------------
	// Calculate the covariance matrix cov(X,Y) and correlation matrix corr(X,Y)

	// Init [N_cuts x N_cuts] matrix elements with 0
	std::vector<std::vector<double> > cov_xy(N_cuts, std::vector<double>(N_cuts, 0.0) );
	std::vector<std::vector<double> >   r_xy(N_cuts, std::vector<double>(N_cuts, 0.0) );


	// Loop over all cuts x = [rows], y = [columns]
	for (int x = 0; x < N_cuts; ++x) {
		for (int y = 0; y < N_cuts; ++y) {

			// Loop over all events
			for (int i = 0; i < N_events; ++i) {

				cov_xy[x][y] += ( CMat_.at(x).at(i) - mean.at(x)) * ( CMat_.at(y).at(i) - mean.at(y));
			}

			cov_xy[x][y] /= (double) N_events; 							    	// 1 / n
			  r_xy[x][y] = cov_xy[x][y] / ( sigma.at(x) * sigma.at(y) + 1e-9 ); // 1 / s_xs_y
		}
	}

	// -------------------------------------------------------------------
	// Print out the correlation matrix corr(X,Y)

	fprintf(fp, "(X,Y) ");
	for (int x = 0; x < N_cuts; ++x) {
		fprintf(fp, "[%d]     ", x);
	}
	fprintf(fp,"\n");
	for (int x = 0; x < N_cuts; ++x) {
		fprintf(fp,"[%d]   ", x);
		for (int y = 0;  y < N_cuts; ++y) {

			double corr = r_xy[x][y];

			// if else, for visual printing reasons
			if (copysign(1.0, corr) == 1) {
				fprintf(fp, "%0.2f    ", corr);
			} else {
				fprintf(fp, "%0.2f   ",  corr);
			}
		}
		fprintf(fp, "\n");
	}

	// Close filepointer
	fclose(fp);
}


// Return the passed events list
//
std::vector<int> CutFlow::getPass() {

	std::vector<int> list;

	int N_cuts = CMat_.size();
	int N_events = CMat_.at(1).size();

	// Loop over all events
	for (int j = 0; j < N_events; ++j) {

		bool all_passed = true;
		// Loop over all cuts
		for (int i = 0; i < N_cuts; ++i) {		

			// It cut was turned ON
			if (status_.at(i) == true) {
				// If cut evaluated to one, i.e., passed
				if (CMat_.at(i).at(j) == 1) {
					// OK
				} else {
					all_passed = false;
					break; // Break inner loop
				}
			}
		}
		if (all_passed)
			list.push_back(j); // Push event ID
	}

	printf("CutFlow::getPass() %d / %d events passed \n", (int) list.size(), (int) N_events);

	return list;
}

// Print out cuts
// This function prints out a serial cut flow, i.e., first cut to cut breaks the flow

void CutFlow::printFlow(const std::string& output_filename, const char* write_mode) {

	FILE* fp;
	fp = fopen(output_filename.c_str(), write_mode);
	
	fprintf(fp, "--------------------------------------------------------------------------\n");
	fprintf(fp, "CutFlow:: printFlow() <Serial flow> results: \n");
	fprintf(fp, "--------------------------------------------------------------------------\n");
	
	int N_cuts   = CMat_.size();
	int N_events = CMat_.at(1).size();

	fprintf(fp, "Input  : %7d (1.00)          %s\n\n", N_events, input_name_.c_str()); 

	// Loop over all events
	std::vector<int> cut_flow(N_cuts,0);

	for (int i = 0; i < N_events; ++i) {

		// Loop over all cuts
		for (int c = 0; c < N_cuts; ++c) {

			// If cut was turned ON
			if (status_.at(c) == true) {

				// If cut evaluated to zero
				if (CMat_.at(c).at(i) == 0) {
					++cut_flow.at(c);
					break; // Linear cutflow, so first cut == 0 breaks flow
				}
			}
		}
	}

	// Loop over all cuts
	int sum = N_events;
	for (int c = 0; c < N_cuts; ++c) {

		// If cut was turned ON
		if (status_.at(c) == true) {
			sum -= cut_flow.at(c);
			fprintf(fp, "Cut [%d]: %7d (%0.2f) %7d  (%s) \n", c, sum, sum / (double) N_events, -cut_flow.at(c), N_names_.at(c).c_str() );
		} else {
			fprintf(fp, "Cut [%d]:                  <OFF>  (%s) \n", c, N_names_.at(c).c_str() );
		}
	}

	int totsumm = 0;
	for (int c = 0; c < N_cuts; ++c) {
		totsumm += cut_flow.at(c);
	}

	fprintf(fp, "--------------------------------------------------------------------------\n");
	fprintf(fp, "Output : %7d (%0.2f) \n\n", N_events - totsumm, (N_events - totsumm) / (double) N_events);

	// close filepointer
	fclose(fp);
}




