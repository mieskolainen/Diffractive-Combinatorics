// Pileup inverse class
//
// mikael.mieskolainen@cern.ch, 2018
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef PILEUP_H
#define PILEUP_H

// C++
#include <iostream>
#include <vector>
#include <cmath>
#include <map>


class PileUp {

public:

	// Constructor and destructor
	PileUp();
	~PileUp();

	// Operator map
	std::vector<double> MapCounts(const std::string& type, const std::vector<double>& y, double R);

	// Return Poisson mu
	double GetMu(double R);

	// Verbose printing
	void SetVerbose(bool verb) {
		VERBOSE = verb;
	}

	// Verbose printing
	void SetPositivity(bool verb) {
		ENFORCE_POSITIVITY = verb;
	}


private:

	bool VERBOSE = true; 			// Verbose print
	bool ENFORCE_POSITIVITY = true; // Require positive definite result
	
	
	// Inversion matrices
	std::map<int, std::vector<std::vector<double> > > A;
	std::map<int, std::vector<std::vector<double> > > invA;

	// Matrix-Vector product
	std::vector<double> MatVecProd(const std::vector<std::vector<double> >& M, const std::vector<double>& x);

	std::vector<double> Map(const std::vector<double>& );

	void PrintComparison(const std::vector<double>& a, const std::vector<double>& b);

};

#endif // PILEUP_H
