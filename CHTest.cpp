/*
 * CahnHillardTest.cpp
 *
 *  Created on: Sep 23, 2013
 *      Author: cyranka
 */


#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <cstring>
//#include "capd/filib/Interval.h"
#include "capd/intervals/Interval.hpp"
#include "config.h"

#define DOUBLE double
#if __FILIB__
  typedef capd::filib::Interval<DOUBLE> Interval;
  #define PI Interval::pi()
#else
  typedef capd::rounding::IntRounding IntRounding;
  typedef capd::intervals::Interval<DOUBLE> Interval;
  #define PI Interval::pi()
//  typedef DOUBLE Interval;
//  #define PI 3.1415926535897932384626433832795
#endif

#include "capd/dynsys/BasicFadTaylor.h"
#include "capd/dynsys/FadTaylor.h"

#include "ComplexScalar.h"
#include "FirstOrderJet.h"
#include "Odd.h"
#include "Index.h"
#include "Coefficients.h"
#include "Pair.h"

#include "Odd.h"
#include "capd/dynset/C0Rect2Set.hpp"
#include "capd/dynset/C0Rect2RSet.hpp"

#include "FFT.h"
#include "Equations.h"
#include "FFTDynSys.h"
#include "PolyBd.h"

#include "DPDEInclusionCW.h"
#include "InclRect2Set.hpp"

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
typedef capd::jaco::MaximumNorm<Index1D> Norm1D;
typedef capd::jaco::Odd<Interval, Index1D, Norm1D, IntervalMatrix> OddSubspace;
typedef capd::jaco::RealPolynomialBound<ComplexScalar, Index1D, 0, OddSubspace> RealPolynomialBound;
//typedef capd::jaco::RealPolynomialBound<ComplexScalar, Index1D, 0> RealPolynomialBound;
typedef capd::jaco::Burgers<RealPolynomialBound> Burgers;
typedef capd::jaco::KS<RealPolynomialBound> KS;

typedef capd::dynset::C0Rect2RSet<IntervalMatrix> C0Set;
//FFT1D

typedef capd::jaco::FFT1D<ComplexScalar, ComplexScalar, 0, 0, RealPolynomialBound> FFT1D;

typedef capd::jaco::DPDE2<KS, FFT1D, 0> DPDE; ///set equation here

typedef capd::jaco::FirstOrderJet<ComplexDerivativePair, 0> FOJ1D;

typedef capd::jaco::ComplexPolyBdJetOptimized<Interval, FOJ1D, Index1D, 0, OddSubspace> ModesContainer;

typedef capd::jaco::FFT1D<FOJ1D, ComplexScalar, 0, 0, ModesContainer> JetFFT1D;

typedef capd::jaco::DPDE2<KS, JetFFT1D, 0> JetDPDE; ///set equation here

typedef capd::jaco::FFTTaylorDynSys<DPDE, JetDPDE, 0> FFTDynSys;

typedef FFTDynSys::ModesContainerType ModesContainer1D;

typedef capd::jaco::DPDEInclusionCW<DPDE, FFTDynSys> DPDEInclusionCW3;
typedef capd::diffIncl2::InclRect2Set<IntervalMatrix, RealPolynomialBound> InclRect2Set;

void calculateDiams(const IntervalVector& v, Interval& maxD){
  int i;
  Interval d;
  maxD = 0;
  generalDebug << "diams: ";
  for(i=0; i < v.size(); i++){
    d = diam(v[i]);
    if(d > maxD) maxD = d;
    generalDebug << d << ", ";
  }
  generalDebug << "\n";
}

void inclPOTest(int testNumber, int approach){
  int m = 23, //change FOJ1D stack dimension
      M = 69,
      dftPts = 48,
      dftPts2 = 150,
      order = 23;
  Interval nu(0.032),
           step(0.0011);

  ///1.initialize FOJ, the dimension of the first order jet and the initial condition type defined by a DPDEContainer instance
  ///should be provided, see FOJ1D::initialize description

  ///begin FOJ initialization for FFT integrator
  capd::jaco::DPDEContainer container;
  ///2.set here the subspace of the initial condition e.g. setToRealValuedOdd means that the initial condition is real valued odd
  container.setToRealValuedOdd();
  FOJ1D::initialize(ModesContainer1D::modes2arraySizeStatic(m), container);
  ///end FOJ initialization

  clock_t start, end;

  DPDEInclusionCW3 diffIncl(m, dftPts, M, dftPts2, PI, nu, order, step, MaxNorm());

  double diam;
  capd::auxil::OutputStream log(std::cout, false, true);
  std::stringstream ss;
  if(testNumber == 0){
    ///3.choose the algorithm type, see enum AlgorithmType in FFTDynSys.h file
    diam = 1e-08;
    ss << "test1_KSincl_";
  }
  if(testNumber == 1){
    ///3.choose the algorithm type, see enum AlgorithmType in FFTDynSys.h file
    diam = 1e-07;
    ss << "test2_KSincl_";
  }
  if(testNumber == 2){
    ///3.choose the algorithm type, see enum AlgorithmType in FFTDynSys.h file
    diam = 1e-06;
    ss << "test3_KSincl_";
  }
  if(testNumber == 3){
    ///3.choose the algorithm type, see enum AlgorithmType in FFTDynSys.h file
    diam = 1e-05;
    ss << "test4_KSincl_";
  }
  if(testNumber == 4){
    diam = 1e-04;
    ss << "test5_KSincl_";
  }

  ss  << "diam_" << diam;
  if(approach == 0)
    ss << "_direct.txt";
  if(approach == 1)
    ss << "_fft.txt";
  if(approach == 2)
    ss << "_fftbutforddirect.txt";
  log.logfile(ss.str().c_str(), true); log.log = true;
  time_t rawtime;
  time ( &rawtime );
  log << "The current local time is: " << ctime (&rawtime) << "\n";
  log << "the whole infinite dimensional system is being integrated (the Lohner algorithm for differential inclusions is used).\n";
  if(approach == 0){
    diffIncl.getDynamicalSystem().setJetDynSysAlgorithmType(capd::jaco::direct);
    log << "Using the direct approach\n";
  }else{
    if(approach == 1){
      diffIncl.getDynamicalSystem().setJetDynSysAlgorithmType(capd::jaco::FFT);
      log << "Using the FFT approach\n";
    }else{
      diffIncl.getDynamicalSystem().setJetDynSysAlgorithmType(capd::jaco::FFTButFirstOrderDirect);
      log << "Using the FFT approach, but the first order normalized derivative is calculated directly to avoid a blowup\n";
    }
  }

  Interval unit(-1, 1);
  Interval r(-diam/2, diam/2);
  Index1D idx;
  RealPolynomialBound u_0(m, M, container), enclosure(m);


  log << "initial condition (infinite dimensional):\n" << u_0 << "\nitegration started, the output below is the Galerkin projection of the set at each timestep\n";
  double retTime = 0.1;
  int i;
  InclRect2Set set(u_0);
  IntervalMatrix mx(m, m);
  Interval max;
  int STEPS = 100;
           //= 4092;
           //=1;

  start = clock();
  for(i=0; i < STEPS; ++i){
    std::cout << "step #" << i << " ";
    set.move(diffIncl);
    log << (IntervalVector)set << "\n";
    if(i * step > retTime){
      std::cout << "Return time passed, full revolution of the orbit\n";
      break;
    }

  }
  log << "\nset at the end (infinite dimensional): " << set.getPerturbationParams() << "\n";
  calculateDiams((IntervalVector)set, max);
  log << "max diameter of the set at the end: " << max << "\n";
  end = clock();
  log << "time: "<<(end-start)/CLOCKS_PER_SEC<<"\n";
}

void basicPOTest(int testNumber, int approach){
  int m = 23, //change FOJ1D stack dimension
      M = 48,
      order = 5;
  Interval nu(0.032),
           step(0.00015);

  ///1.initialize FOJ, the dimension of the first order jet and the initial condition type defined by a DPDEContainer instance
  ///should be provided, see FOJ1D::initialize description

  ///begin FOJ initialization for FFT integrator
  capd::jaco::DPDEContainer container;
  ///2.set here the subspace of the initial condition e.g. setToRealValuedOdd means that the initial condition is real valued odd
  container.setToRealValuedOdd();
  FOJ1D::initialize(ModesContainer1D::modes2arraySizeStatic(m), container);
  ///end FOJ initialization

  clock_t start, end;

  FFTDynSys dynsys( m, M, step, order, PI, nu);

  double diam;
  capd::auxil::OutputStream log(std::cout, false, true);
  std::stringstream ss;
  if(testNumber == 0){
    ///3.choose the algorithm type, see enum AlgorithmType in FFTDynSys.h file
    diam = 1e-08;
    ss << "test1_KSproj_";
  }
  if(testNumber == 1){
    ///3.choose the algorithm type, see enum AlgorithmType in FFTDynSys.h file
    diam = 1e-07;
    ss << "test2_KSproj_";
  }
  if(testNumber == 2){
    ///3.choose the algorithm type, see enum AlgorithmType in FFTDynSys.h file
    diam = 1e-06;
    ss << "test3_KSproj_";
  }
  if(testNumber == 3){
    ///3.choose the algorithm type, see enum AlgorithmType in FFTDynSys.h file
    diam = 1e-05;
    ss << "test4_KSproj_";
  }
  if(testNumber == 4){
    diam = 1e-04;
    ss << "test5_KSproj_";
  }

  ss  << "diam_" << diam;
  if(approach == 0)
    ss << "_direct.txt";
  if(approach == 1)
    ss << "_fft.txt";
  if(approach == 2)
    ss << "_fftbutforddirect.txt";
  log.logfile(ss.str().c_str(), true); log.log = true;
  time_t rawtime;
  time ( &rawtime );
  log << "The current local time is: " << ctime (&rawtime) << "\n";
  log << "only a Galerkin projection is being integrated.\n";
  if(approach == 0){
    dynsys.setJetDynSysAlgorithmType(capd::jaco::direct);
    log << "Using the direct approach\n";
  }else{
    if(approach == 1){
      dynsys.setJetDynSysAlgorithmType(capd::jaco::FFT);
      log << "Using the FFT approach\n";
    }else{
      dynsys.setJetDynSysAlgorithmType(capd::jaco::FFTButFirstOrderDirect);
      log << "Using the FFT approach, but first order is calculated directly to avoid a blowup\n";
    }
  }


  Interval r(-diam/2, diam/2);
  Index1D idx;
  RealPolynomialBound u_0(m), enclosure(m);


  double retTime = 0.1;
  int i;
  C0Set set(u_0, 1);
  IntervalMatrix mx(m, m);
  Interval max;
//  int STEPS //= 100;
//         = 4092;
         // = 4000;
           //=1;

  start = clock();
  log << "initial condition (finite dimensional):\n" << u_0 << "\nitegration started, the output below is the Galerkin projection of the set at each timestep\n";
  for(i=0; true; ++i){
//    std::cout << "step #" << i << " ";
    set.move(dynsys);
    log << (IntervalVector)set << "\n";
    if(i * step > retTime){
      std::cout << "Return time passed, full revolution of the orbit\n";
      break;
    }
  }
  log << "\nset at the end: " << (IntervalVector)set << "\n";
  calculateDiams((IntervalVector)set, max);
  log << "max diameter of the set at the end: " << max << "\n";
  end = clock();
  log << "time: "<<(end-start)/CLOCKS_PER_SEC<<"\n";
}

int main(int argc, char * argv[]){
//  MpFloat::setDefaultPrecision(100);
  setLoggers();
  if(argc != 4){
    std::cerr << "Number of arguments wrong.\nUsage: ./CHTest [projection|inclusion] test_number approach_number\n" <<
                 "if \"projection\" is specified then the Galerkin projection is integrated, when \"inclusion\" is specified the " <<
                 "whole system is integrated (using the Lohner algorithm for inclusions)" <<
                 "\n test_number is the test index, which corresponds to the index provided in Figure 1.18 and 1.19," <<
                 "\napproach_number is the approach type, 0 - the direct approach, 1 - the FFT approach, 2 - the FFT approach, but the first normalized derivative is calculated directly, in order to avoid blow-ups.\n";
  }else{
    if(strcmp(argv[1], "projection")==0){
      basicPOTest(atoi(argv[2]), atoi(argv[3]));
    }else{
      if(strcmp(argv[1], "inclusion")==0){
        inclPOTest(atoi(argv[2]), atoi(argv[3]));
      }else{
        std::cerr << "Number of arguments too small.\nUsage: ./CHTest [projection|inclusion] test_number approach_number\n" <<
                 "if \"projection\" is specified then the Galerkin projection is integrated, when \"inclusion\" is specified the " <<
                 "whole system is integrated (using the Lohner algorithm for inclusions)" <<
                 "\n test_number is the test index, which corresponds to the index provided in Figure 1.18 and 1.19," <<
                 "\napproach_number is the approach type, 0 - the direct approach, 1 - the FFT approach, 2 - the FFT approach, but the first normalized derivative is calculated directly, in order to avoid blow-ups.\n";
      }
    }
  }
  FOJ1D::destroy();
}

