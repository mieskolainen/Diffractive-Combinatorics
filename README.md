# Diffractive-Combinatorics
https://arxiv.org/abs/1910.06279

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## High Energy Diffraction at the LHC via Combinatorial (incidence) Algebra

<img width="600px" src="figs/comb1.jpg">

### Framework with novel constructions and algorithms:

- Unfolded multidimensional 'fiducial vector cross sections'
- Multivariate Maximum Likelihood fit extraction of the effective Pomeron intercept 1 + &Delta;
- Multivariate Maximum Likelihood fit extraction of the diffractive cross sections (SD, DD, ND) and the maximum diffractive mass dependence (minimum average rapidity gap size <&Delta;Y>)
- Möbius inversion of pileup
- Multidimensional algebraic hierarchy of detector level distributions (raw multiplicity etc.) - highly constraining comparison with MC models, probing AGK rules / multipomeron cuts / soft multiparton interactions

For reference, you may cite the paper above.

mikael.mieskolainen@cern.ch, 2019


#### External dependencies:

- <a href="https://alice-doc.github.io/alice-analysis-tutorial/building/">ALICE experiment software</a> (FYI: algorithms here are fully experiment generic)
- <a href="http://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html">ROOUnfold</a> (in /libs)
- <a href="https://github.com/ben-strasser/fast-cpp-csv-parser">csv parser</a> (in /libs, BSD-3 license)

In general, one needs to modify Makefile for your local installation.
