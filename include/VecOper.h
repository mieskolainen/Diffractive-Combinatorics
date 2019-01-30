// Vector, statistics and histogram operations in a namespace
//
// mikael.mieskolainen@cern.ch, 2018
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef VECOPER_H
#define VECOPER_H

// C++ headers
#include <iostream>
#include <vector>
#include <cmath>
#include <utility>

// ROOT
#include "TH1D.h"
#include "TLegend.h"


namespace VecOper {


	// Trigger-Mask statistics
	void TriggerMaskStatistics(const std::vector<double>& xB,
							   const std::vector<double>& xA,
							   const std::vector<double>& xC,
							   const std::vector<double>& xE,
							   double a,
							   double c,
							   double e);

	// Beam-Gas substraction
	std::vector<double> BGSubstract(const std::vector<double>& xB,
								  	  const std::vector<double>& xA,
								  	  const std::vector<double>& xC,
								  	  const std::vector<double>& xE,
								  	  double a,
								  	  double c,
								  	  double e);

	// Construct MC synthesized combinatorial measurement vector
	std::vector<double> GenerateXHat(const std::vector<std::vector<double> >& F, const std::vector<double>& P, double scale);

	// Print difference statistics between two combinatorial vectors
	void PrintXDiff(const std::vector<double>& a, const std::vector<double>& b);
	
	double pow2(double x);

	// Diffractive xi to average rapidity gap <DeltaY>
	double xi2deltaY(double xi);

	// Average rapidity gap <DeltaY> to xi = M^2/s
	double deltaY2xi(double deltaY);
	
	// Histogram operations
	TLegend* AddLegend(std::vector<TH1*> histo, std::vector<std::string> names);
	void NormHist(TH1F* h);

	// Statistical errors
	double GetRatioError(double x, double x_e, double y, double y_e, double cov);
	std::vector<double> GetBinomError(const std::vector<double>& x_count);

	// Vector operations
	std::vector<double> nvec(const std::vector<double>& a);
	std::vector<double> scalevec(const std::vector<double>& a, double scale);
	double vsum(const std::vector<double>& v);
	double vnorm2(const std::vector<double>& v);

	// Bootstrapped samples
	std::vector<std::vector<double> > CreateBootStrapSample(const std::vector<double>& prob, uint N_events, uint N_sample_size, bool FASTMODE);
	std::vector<double> BootStrap(const std::vector<double>& prob, uint N_samples);
	std::vector<double> BootStrapFast(const std::vector<double>& prob, uint N_events);

	// Get and create vectors
	std::vector<double> getcolvec(const std::vector<std::vector<double> >& matrix, uint column_number);
	std::vector<double> linspace(double start, double end, uint n);
	void setcolvec(std::vector<std::vector<double> >& matrix, const std::vector<double>& vec, uint column_number);
		
	// Sort vector elements
	std::vector<uint> vsortind(const std::vector<double>& v);
	void vsort(std::vector<double>& v, const std::vector<uint>& sort_ind);

	// Log-likelihoods
	double logLmultimix(const std::vector<double>& n, const std::vector<double>& p, const std::vector<std::vector<double> >& F);

	// Comparison "metrics"
	double calchi2(const std::vector<double>& obs, const std::vector<double>& exp);
	double calcKS(const std::vector<double>& a, const std::vector<double>& b);
	double calcKL(const std::vector<double>& p, const std::vector<double>& q);
	double calcH(const std::vector<double>& p);

	// Test functions
	bool IsInRange(double value, double min, double max);

	// Bit operations
	uint Gray2Binary(uint number);
	uint Binary2Gray(uint number);

	std::vector<std::vector<int> > ConstructB(uint d);
	std::vector<bool> Ind2Vec(uint ind, uint d);
	std::vector<uint> LRsequence(uint dim);
	int Vec2Ind(std::vector<bool> vec);
	void PrintBits(size_t const size, void const * const ptr);

	// Subspace selection
	void Subspace(const std::vector<uint>& x, const std::vector<int>& subind, double scale);

	// Math functions
	double Binomial(int k, int n, double p);
	int    Factorial(int n);
}

#endif
