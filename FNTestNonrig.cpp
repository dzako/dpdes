/*
 * FNTest.cpp
 *
 * Fitz Hugh Nagumo example rigorous integration
 *
 *
 *  Created on: Mar 31, 2016
 *      Author: Jacek Cyranka
 *
 *
 */

#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <cstring>
#include "config.h"
#include "capd/intervals/Interval.hpp"

#define PI 3.1415926535897932384626433832795
typedef double Double;
typedef Double Interval;


// MOCK FUNCTIONS NEEDED FOR NONRIGOROUS INTEGRATOR TO WORK
Double diam( Double i ){
  return 0.;
}
Double power( Double i, Double e ){
  return pow(i, e);
}
void setLeftBound( Double& i, Double b){
  i = b;
}
void setRightBound( Double& i, Double b){
  i = b;
}



#include "ComplexScalar.h"
#include "Odd.h"
#include "Even.h"
#include "Index.h"
#include "Coefficients.h"

#include "FFT.h"
#include "Equations.h"
#include "FFTDynSys.h"

#include "Pair.h"
#include "FirstOrderJet.h"


#include "capd/alglib/ap.h"
#include "capd/alglib/hessenberg.h"
#include "capd/alglib/hsschur.h"
#include "capd/alglib/lib.h"

#define _D 0


typedef capd::vectalg::Matrix<Interval,_D,_D> IntervalMatrix;
typedef capd::vectalg::Vector<Interval, _D> IntervalVector;
typedef capd::vectalg::MaxNorm<IntervalVector, IntervalMatrix> MaxNorm;
typedef capd::vectalg::EuclNorm<IntervalVector, IntervalMatrix> EuclNorm;
typedef capd::jaco::ComplexScalar<Interval> ComplexScalar;
typedef capd::jaco::ComplexDerivativePair<ComplexScalar> ComplexDerivativePair;

//Integrators 1D
//Basic Fad Taylor 1D
typedef capd::jaco::Index1D Index1D;
typedef capd::jaco::MaximumNorm<Index1D> MaximumNorm;
//typedef capd::jaco::Odd<Interval, Index1D, MaximumNorm> OddSubspace;
typedef capd::jaco::Even<Interval, Index1D, MaximumNorm, IntervalMatrix> EvenSubspace;
typedef capd::jaco::RealPolynomialBound<ComplexScalar, Index1D, 0, EvenSubspace> RealPolynomialBound;
//typedef capd::jaco::RealPolynomialBound<ComplexScalar, Index1D, 0> RealPolynomialBound;
typedef capd::jaco::Burgers<RealPolynomialBound> Burgers;
typedef capd::jaco::KS<RealPolynomialBound> KS;
typedef capd::jaco::GL<RealPolynomialBound> GL;
typedef capd::jaco::SH<RealPolynomialBound> SH;
typedef capd::jaco::FN<RealPolynomialBound> FN;

//FFT1D
typedef capd::jaco::FFT1D<ComplexScalar, ComplexScalar, 0, 0, RealPolynomialBound> FFT1D;

typedef capd::jaco::DPDE3<FN, FFT1D, 0> FNDPDE; ///set equation here

typedef capd::jaco::FFTBasicDynSys<FNDPDE, 0> FFTDynSys;

typedef capd::jaco::FirstOrderJet<ComplexDerivativePair, 0> FOJ1D;

typedef capd::jaco::ComplexPolyBdJetOptimized<Interval, FOJ1D, Index1D, 0, EvenSubspace> ModesContainer;

typedef capd::jaco::FFT1D<FOJ1D, ComplexScalar, 0, 0, ModesContainer> JetFFT1D;

typedef capd::jaco::DPDE3<FN, JetFFT1D, 0> JetFNDPDE; ///set equation here


void basicTest(int testNumber, int approach){
  int n = 15, //change FOJ1D stack dimension
      m = 45,
      order = 5;
  Interval nu,
           step(0.0001);

  ///1.initialize FOJ, the dimension of the first order jet and the initial condition type defined by a DPDEContainer instance
  ///should be provided, see FOJ1D::initialize description
  
  clock_t start, end;

  
  int STEPS;
  capd::auxil::OutputStream log(std::cout, false, true);
  std::stringstream ss;

  if(testNumber == 0){
    ///3.choose the algorithm type, see enum AlgorithmType in FFTDynSys.h file
    nu = _NU_;
    STEPS = 100000;
    n = 41;
    step = 0.001;

    ss << "nonrig_FNproj_";
  }

  if(testNumber > 0){
    std::cerr << "Test with nr > 0 not implemented!";
    exit(1);
  }

  ///begin FOJ initialization for FFT integrator
  capd::jaco::DPDEContainer container;
  ///2.set here the subspace of the initial condition e.g. setToRealValuedOdd means that the initial condition is real valued odd
  container.setToRealValuedEven();
  FOJ1D::initialize(ModesContainer::modes2arraySizeStatic(n), container);
  ///end FOJ initialization

  ///3.choose the algorithm type, see enum AlgorithmType in FFTDynSys.h file
  FFTDynSys dynsys( n, m, step, order, PI, nu);
  
  ss  << "nu_" << rightBound(nu);

  ss  << "_eps_" << rightBound(nu);

  log.logfile(ss.str().c_str(), true); log.log = true;

  Index1D idx;
  RealPolynomialBound u_0(n), enclosure(n);

  idx[0] = 0;
  u_0[idx].re = 1 ;
  idx[0] = 1;
  u_0[idx].re =  1 ;

  idx[0] = n/2 + 3;
  u_0[idx].re = 1 ;
    
  Interval max;

  IntervalVector u = u_0;

  start = clock();
  int i;

  for(i=0; i < STEPS; ++i){

    dynsys.move(u, u);
    if(i % 10 == 0){
      log << u << "\n";
    }
  }

  end = clock();
  log << "time: "<<(end-start)/CLOCKS_PER_SEC<<"\n";
}




int main(int argc, char * argv[]){
  setLoggers();
  if(argc != 4){
      std::cerr << "Number of arguments wrong.\nUsage: ./SHTest [projection|inclusion] test_number approach_number\n" << 
                   "if \"projection\" is specified then the Galerkin projection is integrated, when \"inclusion\" is specified the " <<
                   "whole system is integrated (using the Lohner algorithm for inclusions)" <<
                   "\ntest_number is the test index, which corresponds to the index provided in Figure 1.21," <<
                   "\napproach_number is the approach type, 0 - the direct approach, 1 - the FFT approach, 2 - the FFT approach, but the first normalized derivative is calculated directly, in order to avoid blow-ups.\n";
    }else{
      if(strcmp(argv[1], "projection") == 0){
        basicTest(atoi(argv[2]), atoi(argv[3]));
      }else{
        if(strcmp(argv[1], "inclusion") == 0){
          //diffInclTest(atoi(argv[2]), atoi(argv[3]));
        }else{
          std::cerr << "Number of arguments too small.\nUsage: ./SHTest [projection|inclusion] test_number approach_number\n" << 
                   "if \"projection\" is specified then the Galerkin projection is integrated, when \"inclusion\" is specified the " <<
                   "whole system is integrated (using the Lohner algorithm for inclusions)" <<
                   "\ntest_number is the test index, which corresponds to the index provided in Figure 1.21," <<
                   "\napproach_number is the approach type, 0 - the direct approach, 1 - the FFT approach, 2 - the FFT approach, but the first normalized derivative is calculated directly, in order to avoid blow-ups.\n";
        }
      }    
    }
  FOJ1D::destroy();
}

