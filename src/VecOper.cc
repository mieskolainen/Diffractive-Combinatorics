// Vector, statistics and histogram operations in a namespace
//
// mikael.mieskolainen@cern.ch, 2018
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.


// C++ headers
#include <iostream>
#include <vector>
#include <cmath>
#include <utility>

// ROOT headers
#include "TMath.h"
#include "TRandom3.h"
#include "TH1.h"

// ROOT headers
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH2.h"
#include "TArrayD.h"
#include <TObject.h>
#include <TBits.h>
#include <TClonesArray.h>


// Own headers
#include "VecOper.h"


namespace VecOper {

// We need this for Bootstrap routines, for example
// Do not delete this, it is within namespace -> scope reasonably safe!
TRandom3* r = new TRandom3(123); // seed with 123

// Print out Beam-Empty statistics
void TriggerMaskStatistics(const std::vector<double>& xB,
						   const std::vector<double>& xA,
				 	  	   const std::vector<double>& xC,
				  	  	   const std::vector<double>& xE,
								  	  			double a,
								  	  			double c,
								  	  			double e) {

	printf("\n\n");
	printf("VecOper::TriggerMaskStatistics:: BEAM-GAS combinations: \n");
	printf("Scale factors: a = %0.3f, c = %0.3f, e = %0.3f \n", a, c, e);

	// Do the substraction
	std::vector<double> x_corrected = BGSubstract(xB, xA, xC, xE, a, c, e);

	// Only printing from here
	for (uint i = 0; i < x_corrected.size(); ++i) {

		// Get binary expansion (vector) for this k = 0...2^N-1
		std::vector<Bool_t> binvec = Ind2Vec(i, std::log2(xB.size()));
		printf("%3d [", i);
		for (uint bit = 0; bit < binvec.size(); ++bit) {
			printf("%d", (int) binvec.at(bit));
		}

		printf("] %-7.0f = [%-7.0f] - a[%-7.0f] - c[%-7.0f] + 2e[%-7.0f] ratio %0.3f \n",
			x_corrected.at(i),
			xB.at(i),
			xA.at(i),
			xC.at(i),
			xE.at(i),
			x_corrected.at(i) / (double) xB.at(i));

		// THIS SHOULD BE DONE IN THE FUTURE, BACKGROUND SUBTRACTION OF HISTOGRAMS

		// Histograms for each combination
   		//CorrectBGHist1(cB->hSPDbit[i],  cA->hSPDbit[i],  cC->hSPDbit[i],  cE->hSPDbit[i],      scaleA, scaleC, scaleE);
   		//CorrectBGHist1(cB->hSPDFO[i],   cA->hSPDFO[i],   cC->hSPDFO[i],   cE->hSPDFO[i],       scaleA, scaleC, scaleE);
   		//CorrectBGHist1(cB->hSPDTR[i],   cA->hSPDTR[i],   cC->hSPDTR[i],   cE->hSPDTR[i],       scaleA, scaleC, scaleE);

   		// 2D
   		//CorrectBGHist2(cB->h2SPDFOTR[i], cA->h2SPDFOTR[i], cC->h2SPDFOTR[i], cE->h2SPDFOTR[i], scaleA, scaleC, scaleE);

   		// C and A-sides
		//for (uint j = 0; j < 2; ++j) {

		//	CorrectBGHist1(cB->hADCharge[i][j], cA->hADCharge[i][j], cC->hADCharge[i][j], cE->hADCharge[i][j], scaleA, scaleC, scaleE);
		//	CorrectBGHist1(cB->hADTime[i][j],   cA->hADTime[i][j],   cC->hADTime[i][j],   cE->hADTime[i][j],   scaleA, scaleC, scaleE);

		//	CorrectBGHist1(cB->hV0Charge[i][j], cA->hV0Charge[i][j], cC->hV0Charge[i][j], cE->hV0Charge[i][j], scaleA, scaleC, scaleE);
		//	CorrectBGHist1(cB->hV0Time[i][j],   cA->hV0Time[i][j],   cC->hV0Time[i][j],   cE->hV0Time[i][j],   scaleA, scaleC, scaleE);

			// 2D
		//	CorrectBGHist2(cB->h2ADCT[i][j],    cA->h2ADCT[i][j],    cC->h2ADCT[i][j],    cE->h2ADCT[i][j],    scaleA, scaleC, scaleE);
		//	CorrectBGHist2(cB->h2V0CT[i][j],    cA->h2V0CT[i][j],    cC->h2V0CT[i][j],    cE->h2V0CT[i][j],    scaleA, scaleC, scaleE);
		//}
	}
	printf("\n");
}




// Make a statistcal Beam-Gas substraction
std::vector<double> BGSubstract(const std::vector<double>& xB,
								  const std::vector<double>& xA,
								  const std::vector<double>& xC,
								  const std::vector<double>& xE,
								  double a,
								  double c,
								  double e) {

	std::vector<double> y(xB.size(), 0.0);

	// Loop over vector components
	for (uint i = 0; i < y.size(); ++i) {
		
		// Beam-Gas correction: Y = B - aA - cC + 2eE
		// ********************************************************************************************
		// See the simple logic of this here: https://twiki.cern.ch/twiki/pub/ALICE/PWG1EvSelDocumentation/BeamGas.pdf
		// The positive sign on E-mask is due to double counting sign flip.

		y.at(i) = std::round(xB.at(i) - a*xA.at(i) - c*xC.at(i) + 2*e*xE.at(i));

		// Do not let event count go below zero (physical boundary condition)
		y.at(i) = (y.at(i) < 0) ? y.at(i) = 0 : y.at(i);
	}

	return y;
}

// Construct MC synthesized combinatorial measurement vector
// F is 2^N x K process likelihood matrix
// P is K x 1 process fraction matrix
// scale is a scaling factor (events, cross sections (barn) etc.)
std::vector<double> GenerateXHat(const std::vector<std::vector<double> >& F, const std::vector<double>& P, double scale) {

	// Construct Monte Carlo synthesized measurement vector
	std::vector<double> x_hat(F.size(), 0.0);

	// Loop over components, 2^N
	for (uint i = 0; i < x_hat.size(); ++i) {

		// Loop over processes, SDL, SDR, ...
		for (uint j = 0; j < P.size(); ++j) {
			x_hat.at(i) += F.at(i).at(j) * P.at(j) * scale;
		}
	}

	return x_hat;
}


// Print out difference statistics component by components between two combinatorial vectors
void PrintXDiff(const std::vector<double>& a, const std::vector<double>& b) {

	const double SAFE_EPS = 1; // 1 event
	
	printf("\n");
	printf("                  a/b          a-b          (a-b)^2/b    (a-b)^2/a \n");

	std::vector<double> val1(a.size(), 0.0);
	std::vector<double> val2(a.size(), 0.0);
	std::vector<double> val3(a.size(), 0.0);
	std::vector<double> val4(a.size(), 0.0);

	for (uint i = 0; i < a.size(); ++i) {

		// Get binary expansion (vector) for this i = 0...2^N-1
		std::vector<Bool_t> binvec = Ind2Vec(i, std::log2(a.size()));
		printf("%3d [", i);
		for (uint bit = 0; bit < binvec.size(); ++bit) {
			printf("%d", (int) binvec.at(bit));
		}

		val1.at(i) = a.at(i)/(b.at(i) + SAFE_EPS);
		val2.at(i) = a.at(i) - b.at(i);
		val3.at(i) = pow2(a.at(i) - b.at(i))/(b.at(i) + SAFE_EPS);
		val4.at(i) = pow2(a.at(i) - b.at(i))/(a.at(i) + SAFE_EPS);

		printf("]:     %-12.2f %-12.1f %-12.1f %-12.1f \n", val1.at(i), val2.at(i), val3.at(i), val4.at(i) );
	}
	double N = a.size();
	printf("\n");
	printf("Mean:             %-12.2f %-12.1f %-12.1f %-12.1f \n", vsum(val1)/N, vsum(val2)/N, vsum(val3)/N, vsum(val4)/N);

}



double pow2(double x) {
	return x*x;
}

// Diffractive xi to average rapidity gap <DeltaY>, xi here by definition M^2/s
double xi2deltaY(double xi)     { return -std::log(xi + 1e-12); }

// Average rapidity gap <DeltaY> to diffractive xi=M^2/s
double deltaY2xi(double deltaY) { return  std::exp(-deltaY); }

// Add legend
TLegend* AddLegend(std::vector<TH1*> histo, std::vector<std::string> names) {

	double xl1 = 0.05, yl1 = 0.75, xl2 = xl1 + 0.3, yl2 = yl1 + 0.125;

	TLegend* leg = new TLegend(xl1,yl1,xl2,yl2);

	for (uint i = 0; i < names.size(); ++i) {
		leg->AddEntry(histo.at(i), names.at(i).c_str());
	}
	leg->Draw();
	return leg;
}


// Normalize ROOT histogram
void NormHist(TH1F* h) {

	double norm = h->GetEntries();
	if (norm) h->Scale(1.0/norm);
}


// Count ratio X/Y error by Taylor expansion
double GetRatioError(double x, double x_e, double y, double y_e, double cov) {

	const double SAFE_EPS = 1e-12;

	// Var(X/Y) ~= (E[x]/E[y])^2 ( Var(x)/E[x]^2 + Var(y)/E[y]^2 - 2Cov(x,y)/(E[x]E[y]))
	double var = pow2(x / (y + SAFE_EPS)) * ( pow2(x_e) / (pow2(x) + SAFE_EPS) + 
				   pow2(y_e)/(pow2(x) + SAFE_EPS) - 2*cov/(x*y + SAFE_EPS));

	return std::sqrt(var);
}



// Get binomial / multinomial errors (The simple formula)
std::vector<double> GetBinomError(const std::vector<double>& x_count) {

	// First get normalized version
	std::vector<double> x_prob = VecOper::nvec(x_count);

	// Count the total number of events (trials)
	double N = VecOper::vsum(x_count);

	// Create error vector
	std::vector<double> x_err(x_count.size(), 0);

	for (uint i = 0; i < x_err.size(); ++i) {
		x_err.at(i) = std::sqrt( x_prob.at(i)*(1 - x_prob.at(i)) * N );
	}

	return x_err;
}



// Create bootstrapped sample of vectors
std::vector<std::vector<double> > CreateBootStrapSample(const std::vector<double>& prob, uint N_events, uint N_sample_size, Bool_t FASTMODE) {

	printf("CreateBootStrapSample:: N = %d, events = %d \n", N_sample_size, N_events);

	// Init matrix of size [prob.size() x N_sample_size] with zeros
	std::vector<std::vector<double> > xrix(prob.size(), std::vector<double>(N_sample_size, 0));

	double counter = 0.1; // For printing out the progression
	printf("0 "); fflush(stdout);

	// Fast Poisson sampling (approximation) (kTRUE), Slow multinomial sampling (accurate) (kFALSE)
	// Makes a statistically significant difference with bins with low statistics, but speed difference is huge.

	// Now create data set
	for (uint j = 0; j < N_sample_size; ++j) {

		if ( (j /(double)N_sample_size) > counter) {
			printf("%0.0f ", counter * 100); fflush(stdout);
			counter += 0.1;
		}

		// Create one bootstrapped vector of rates, seed with different j
		std::vector<double> vec = (FASTMODE == kTRUE) ? 
			VecOper::BootStrapFast(prob, N_events) : VecOper::BootStrap(prob, N_events);
		
		// Fill it as a column in the matrix
		for (uint i = 0; i < vec.size(); ++i) {
			xrix.at(i).at(j) = vec.at(i);
		}

		//if (j % 10 == 0)
		//	printf("Sample %d/%d \n", j+1, N_sample_size);
	}
	printf("100 %%\n"); // at the end of event loop
	
	return xrix;
}

// Create one bootstrapped vector
std::vector<double> BootStrap(const std::vector<double>& prob, uint N_events) {

	std::vector<double> vec(prob.size(), 0);
	uint events = 0;

	// Draw bootstrap samples
	while (events < N_events) {

		// Draw the bin, for example bins 0...7, use Integer(8)
		uint bin = r->Integer(prob.size());

		// Draw second random number betwen 0...1, if less than probability, then fill it
		if (r->Rndm() < prob.at(bin)) {
			++vec.at(bin);
			++events;
		}
	}
	return vec;
}


// Create one bootstrapped vector fast with Poisson approximation
std::vector<double> BootStrapFast(const std::vector<double>& prob, uint N_events) {

	std::vector<double> vec(prob.size(), 0);

	for (uint i = 0; i < vec.size(); ++i) {

		// Poisson mean
		double lambda = prob.at(i) * N_events;

		// Draw Poisson distributed random number
		vec.at(i) = r->Poisson(lambda);
	}
	return vec;
}


// Return column vector of a matrix
std::vector<double> getcolvec(const std::vector<std::vector<double> >& matrix, uint column_number) {

	const uint number_of_rows = matrix.size();
	std::vector<double> vec(number_of_rows, 0);

	// Loop over rows of the particular column
	for (uint i = 0; i < number_of_rows; ++i) {

		vec.at(i) = matrix.at(i).at(column_number);
	}

	return vec;
}

// Create MATLAB style linearly spaced vector
std::vector<double> linspace(double start, double end, uint n) {

	std::vector<double> vec(n, 0);
	double step = (end - start) / (n-1);

	for (uint i = 0; i < n; ++i) {
		vec.at(i) = start + step*i;
	}

	return vec;
}

// Set column vector
void setcolvec(std::vector<std::vector<double> >& matrix, const std::vector<double>& vec, uint column_number) {

	// Fill it as a column in the matrix
	for (UInt_t i = 0; i < vec.size(); ++i) {
		matrix.at(i).at(column_number) = vec.at(i);
	}
}


// Return reordered vector according to given index order
void vsort(std::vector<double>& v, const std::vector<uint>& sort_ind) {

	std::vector<double> reordered;

	for (uint i = 0; i < v.size(); ++i) {
		reordered.push_back( v.at(sort_ind.at(i)) );
	}
	for (uint i = 0; i < v.size(); ++i) {
		v.at(i) = reordered.at(i);
	}

}


// Return reordered indices in descending order
// std::vector <double> v = {666, 23, 884, 483}; // input data example
std::vector<uint> vsortind(const std::vector<double>& v) {
	
	std::multimap <double, uint> m; // mapping from value to its index
	for (auto it = v.begin(); it != v.end(); ++it) {
    	m.insert(std::make_pair(*it, it - v.begin()));
	}

		// reordered indices in ascending order
	std::vector<uint> s;

	for (auto it = m.begin(); it != m.end(); ++it) {
		s.push_back( it->second );
	}

	// Turn into descending order
	std::reverse(s.begin(), s.end());

	return s;
}


// Return normalized vector
std::vector<double> nvec(const std::vector<double>& a) {

	double sum = VecOper::vsum(a);
	std::vector<double> vec;

	for (uint i = 0; i < a.size(); ++i) {
		vec.push_back(a.at(i) / sum);
	}

	return vec;
}

// Return scaled vector
std::vector<double> scalevec(const std::vector<double>& a, double scale) {

	std::vector<double> vec;
	for (uint i = 0; i < a.size(); ++i) {
		vec.push_back(a.at(i) * scale);
	}

	return vec;
}


// Return sum of vector elements
double vsum(const std::vector<double>& v) {

	double sum = 1e-12;
	for (uint i = 0; i < v.size(); ++i) {
		sum += v.at(i);
	}

	return sum;
}

// Return vector L2-norm (Euclidian)
double vnorm2(const std::vector<double>& v) {

	double sum = 1e-12;
	for (uint i = 0; i < v.size(); ++i) {
		sum += v.at(i)*v.at(i);
	}

	return sqrt(sum);
}


// Return NEGATIVE log-likelihood of a multinomial distribution with mixture process probabilities p_s
// ln(L) = \sum_{k=1}^K n_k ln( \sum_s^S p_s q_k^{(s)} ), with constraint \sum_{s=1}^S p_s = 1
// 
// where 
//
// q_k^{s} is the probability of process s falling into the bin q_k (model), s = 1 ... S
// n_k is the observed number of events in the bin k of mixture process (measured) 
// 
// K is the number of bins (2^N for combinatorics)
// S is the number of processes (in mixture superposition, e.g. SDL,SDR,DD,ND)
//
//
// below q_k^(s) is the K x S process likelihood matrix F
//
// See discussion in e.g. https://arxiv.org/pdf/0803.2711.pdf
//
// This is non-extended likelihood, i.e., total number of events N = \sum_k n_k is not considered here.
//
// Remember that chi^2 fits are easily biased in normalization type fits!
// Extended Maximum Likelihood needs to be used for fitting the normalization.

double logLmultimix(const std::vector<double>& n, const std::vector<double>& p, const std::vector<std::vector<double> >& F) {

	double logL = 0.0;
	const double SAFE_EPS = 1e-12;

	// Make sure the constraint sum = 1 holds
	std::vector<double> Pn = nvec(p);

	// Loop over measured bins
	for (uint k = 0; k < n.size(); ++k) {

		double summ = 0.0;
		// Loop over scattering processes
		for (uint s = 0; s < Pn.size(); ++s) {
			summ += Pn.at(s) * F.at(k).at(s);
		}
		if (summ > SAFE_EPS)
			logL += n.at(k) * std::log(summ);
	}

	return (-1.0)*logL;
}


// Return bin by bin Chi^2 between two count vectors (or histogram vectors)
//
// Example:
//
// obs is the observed number of events per vector bin/element (measurement)
// exp is the expected number of events per vector bin/element (H_0 null hypothesis ~ model)
//
// The test assumes n->infty limit, then the limiting distribution is chi^2 distributed
// with k - 1 degrees of freedom
double calchi2(const std::vector<double>& obs, const std::vector<double>& exp) {

	double chi2 = 0;
	const double SAFE_EPS = 1e-12;

	for (uint k = 0; k < obs.size(); ++k) {
		if (std::abs(exp.at(k)) > SAFE_EPS) {
			chi2 += pow2(obs.at(k) - exp.at(k)) / exp.at(k);
		}
	}
	return chi2;
}


// Return Kullback-Leibler divergence (relation to likelihoods)
// D_KL(P,Q) ) = -sum_i P(i)log(Q(i)) + sum_i (P(i)Log(P(i)) = H(P,Q) - H(P),
// where H(P,Q) is the cross entropy of P and Q, H(P) is the entropy of P.
//
// P is the data and Q is the (theory) model probability distributions
double calcKL(const std::vector<double>& p, const std::vector<double>& q) {

	const double SAFE_EPS = 1e-12;
	double KL = 0;

	for (uint i = 0; i < p.size(); ++i) {

		if (std::abs(q.at(i) - 0) > SAFE_EPS && std::abs(p.at(i) - 0) > SAFE_EPS) { // 0-bins are skipped
			KL += p.at(i) * std::log(p.at(i) / q.at(i));
		}
	}

	return KL;
}


// Return Kolmogorov-Smirnov error between two discrete probability distributions
double calcKS(const std::vector<double>& x, const std::vector<double>& y) {

	double KS = -1e9;

	// First calculate cumulative distributions
	std::vector<double> xF(x.size(), 0);
	std::vector<double> yF(y.size(), 0);

	for (uint i = 0; i < x.size()-1; ++i) {

		double cumsum_x = 0;
		double cumsum_y = 0;
		for (uint j = 0; j < i+1; ++j) {
			cumsum_x += x.at(j);
			cumsum_y += y.at(j);
		}

		xF.at(i) = cumsum_x;
		yF.at(i) = cumsum_y;
	}
	xF.at(xF.size()-1) = 1; // Integral (sum) to one by definition
	yF.at(yF.size()-1) = 1;

	// Now evaluate Max(|F_X(i) - F_Y(i)|)
	for (uint i = 0; i < xF.size(); ++i) {

		double eval = std::abs(xF.at(i) - yF.at(i));
		if (eval > KS) {
			KS = eval;
		}
	}

	return KS;
}

// Calculate Shannon entropy (in bits as is the most common -> log2() )
double calcH(const std::vector<double>& p) {

	const double SAFE_EPS = 1e-12;
	double H = 0;

	for (uint i = 0; i < p.size(); ++i) {
		if (p.at(i) > SAFE_EPS) { // Not defined for 0 bins
			H -= std::log2(p.at(i))*p.at(i);
		}
	}

	return H;
}

// Construct binary reflect Gray code (BRGC) 
// which is the Hamiltonian Path on a d-dim unit hypercube (2^d Boolean vector space)
// The operator ^ is Exclusive OR (XOR) and the operator >> is bit shift right
uint Binary2Gray(uint number) {

   return number ^ (number >> 1);
}


// Convert BRGC (binary reflect Gray code) to a binary number
// Gray code is obtained as XOR (^) with all more significant bits
uint Gray2Binary(uint number) {

  for (uint mask = number >> 1; mask != 0; mask = mask >> 1) {
    number = number ^ mask;
  }
  return number;
}


// Construct binary matrix in normal binary order
std::vector<std::vector<int> > ConstructB(uint d) {

	//printf("ConstructB() \n");

	std::vector<std::vector<int> > B;

	uint N = std::pow(2,d);

	// First initialize
	std::vector<int> rvec(d, 0);
	for (uint i = 0; i < N; ++i) {
        B.push_back(rvec);
    }

	// Now fill
	for (uint i = 0; i < N; ++i) {

		// Get binary expansion (vector) for this i = 0...2^d-1
		std::vector<Bool_t> binvec = VecOper::Ind2Vec(i, d);

		for (uint j = 0; j < d; ++j) {
			B.at(i).at(j) = binvec.at(j); 
		}
	}

	return B;
}


// OEIS.org A030109 sequence of left-right bit reversed sequence
// Input: dim = Boolean vector space dimension
std::vector<uint> LRsequence(uint dim) {

	// The sequence
	std::vector<uint> seq(std::pow(2,dim), 0);

	for (uint i = 0; i < seq.size(); ++i) {

		// Write i in binary with dimension d
		std::vector<Bool_t> vec = VecOper::Ind2Vec(i, dim);

		// Reverse the bits
		reverse(vec.begin(), vec.end());

		// Turn back to index representation
		seq.at(i) = VecOper::Vec2Ind(vec);
	}

	/* // DEBUG
	for (uint i = 0; i < seq.size(); ++i) {
		printf("%u: %u \n", i, seq.at(i) );
	}
	printf("\n");*/

	return seq;
}


// Boolean vector to index representation
// [the normal order: 0 ~ 00, 1 ~ 01, 2 ~ 10, 3 ~ 11 ...]
int Vec2Ind(std::vector<Bool_t> vec) {

	// Swap bit direction (comment this for the reverse ordering)
	reverse(vec.begin(), vec.end());

    int retval = 0;
    int i = 0;

    for (std::vector<Bool_t>::iterator it = vec.begin(); it != vec.end(); ++it, ++i) {
        if (*it) {
            retval |= 1 << i;
        }
    }
    return retval;
}


// Index to Boolean vector representation
// [the normal order: 0 ~ 00, 1 ~ 01, 2 ~ 10, 3 ~ 11 ...]
// Input: x = index representation
//        d = Boolean vector space dimension
std::vector<Bool_t> Ind2Vec(uint ind, uint d) {

	std::vector<Bool_t> ret;

	while (ind) {
		if (ind&1)
			ret.push_back(1);
		else
			ret.push_back(0);

		ind >>= 1;
	}
	reverse(ret.begin(),ret.end());

    // Now will the d-dimensional bit vector starting from rightmost bit
    std::vector<Bool_t> final(d, 0);
    for (uint i = 0; i < ret.size(); ++i) {
    	final.at(d-1 - i) = ret.at(ret.size()-1 - i);
    }

	return final;
}

// Assumes big endian
// Example of how to use: int i = 6; PrintBits(sizeof(i), &i);
// >> 00000000000000000000000000000110
void PrintBits(size_t const size, void const * const ptr) {

    unsigned char* b = (unsigned char*) ptr;
    unsigned char byte;

    // Loop over bytes, e.g. for int32 the size is 4
    for (int i = size-1; i >= 0; i--) {
    	// Loop over bits of a byte
        for (int j = 7; j >= 0; j--) {
            byte = (b[i] >> j) & 1;
            printf("%u", byte);
        }
    }
    puts("");
}

// Selector function
Bool_t IsInRange(double value, double min, double max) {

	if (value >= min && value <= max) {
		return kTRUE;
	} else {
		return kFALSE;
	}

}

// Print subspace combinations
void Subspace(const std::vector<uint>& x, const std::vector<int>& subind, double scale) {

  // Sum over the full space
  uint N_events = 0;
  for (uint i = 0; i < x.size(); ++i) {
    N_events += x.at(i);
  }

  // Number of events in this particular subspace
  uint N_visible = 0;

  // Loop over subspace combinations
  for (uint i = 0; i < std::pow(2, subind.size()); ++i) {

    std::vector<bool> subvec = Ind2Vec(i, subind.size());

    uint N_this = 0;

    // Find these subspace combinations from the fullspace
    for (uint i = 0; i < x.size(); ++i) {

      // Create combination vector, check condition
      std::vector<bool> fullvec = Ind2Vec(i, std::log(x.size()) / std::log(2));
      std::vector<bool> booleans(subind.size(), false);
        
      // Go through dimensions
      uint ok = 0;
      for (uint z = 0; z < subind.size(); ++z) {
        if (fullvec[subind[z]] == subvec[z])
          ++ok;        
      }
      if (ok == subind.size()) // All ok
        N_this += x.at(i);
    }

    // Print ID and [binary vector]
    printf("%2d ", i);
    std::vector<bool> bvec = Ind2Vec(i, subind.size());
    printf("[");
    for (uint k = 0; k < bvec.size(); ++k) {
      printf("%d", (int)bvec.at(k));
    }
    printf("] ");

    // Print rate, Simple Poisson error -> change to binomial
    printf("%0.3f +- %0.3f", N_this / (double) N_events * scale, std::sqrt(N_this) / (double) N_events * scale); 
    printf("\n");

    if (i != 0)
      N_visible += N_this;
  }
  // Simple Poisson error -> change to binomial
  int ncomb = std::pow(2, subind.size()) - 1;
  printf("Sum [1:%d] = %0.3f +- %0.3f \n", ncomb, N_visible / (double) N_events * scale, std::sqrt(N_visible) / (double) N_events * scale);
  printf("\n");
}

// Binomial probability, k success in n trials, with p being probability of success
// n!/(k!*(n-k)!) p^k (1-p)^(n-k)

double Binomial(int k, int n, double p) {

	double fact = Factorial(n)/(Factorial(k)*Factorial(n-k) + 1e-12);
	double binprob = fact * std::pow(p,k) * std::pow(1-p, n-k);

	return binprob;
}

// n!
int Factorial(int n) {

    int ret = 1;
    for (int i = 1; i <= n; ++i){
        ret *= i;
    }
    return ret;
}


} // Namespace ends
