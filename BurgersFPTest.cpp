/*
 * BurgersFPTest.cpp
 *
 *  Created on: Sep 14, 2011
 *      Author: cyranka
 */
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <cstring>
//#include "capd/filib/Interval.h"
#include "capd/intervals/Interval.hpp"
#include "config.h"

//#include "capd/intervals/MpInterval.h"
//#include "capd/multiPrec/MpReal.h"
//typedef capd::multiPrec::MpReal MpFloat;


#define DOUBLE long double
#if __FILIB__
  typedef capd::filib::Interval<DOUBLE> Interval;
#else
  typedef capd::intervals::Interval<DOUBLE> Interval;
  #define PI Interval::pi()
//  typedef DOUBLE Interval;
//  #define PI 3.1415926535897932384626433832795
#endif

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

#define _D 0
#define _m 3
#define _ni 5
#define _order 3


typedef capd::vectalg::Matrix<Interval,_D,_D> IntervalMatrix;
typedef capd::vectalg::Vector<Interval, _D> IntervalVector;
typedef capd::vectalg::EuclNorm<IntervalVector, IntervalMatrix> EuclNorm;
typedef capd::jaco::ComplexScalar<Interval> ComplexScalar;
typedef capd::jaco::ComplexDerivativePair<ComplexScalar> ComplexDerivativePair;

//Integrators 1D
//Basic Fad Taylor 1D
typedef capd::jaco::Index1D Index1D;
typedef capd::jaco::MaximumNorm<Index1D> MaximumNorm;
//typedef capd::jaco::Odd<Interval, Index1D, MaximumNorm> OddSubspace;
typedef capd::jaco::RealPolynomialBound<ComplexScalar, Index1D, 0> RealPolynomialBound;
//typedef capd::jaco::RealPolynomialBound<ComplexScalar, Index1D, 0> RealPolynomialBound;
typedef capd::jaco::Burgers<RealPolynomialBound> Burgers;
typedef capd::jaco::KS<RealPolynomialBound> KS;

typedef capd::jaco::MaximumNorm<Index1D> Norm1D;
//typedef capd::jaco::Real<Interval, Index1D, Norm1D> Real1D;
typedef capd::jaco::Odd<Interval, Index1D, Norm1D> Real1D;
typedef capd::dynset::C0Rect2RSet<IntervalMatrix> C0Set;
//FFT1D

typedef capd::jaco::FFT1D<ComplexScalar, ComplexScalar, 0, 0, RealPolynomialBound> FFT1D;

typedef capd::jaco::DPDE2<Burgers, FFT1D, 0> BurgersDPDE; ///set equation here

typedef capd::jaco::FirstOrderJet<ComplexDerivativePair, 0> FOJ1D;

typedef capd::jaco::ComplexPolyBdJetOptimized<Interval, FOJ1D, Index1D, 0> ModesContainer;

typedef capd::jaco::FFT1D<FOJ1D, ComplexScalar, 0, 0, ModesContainer> JetFFT1D;

typedef capd::jaco::DPDE2<Burgers, JetFFT1D, 0> JetBurgersDPDE; ///set equation here

typedef capd::jaco::FFTTaylorDynSys<BurgersDPDE, JetBurgersDPDE, 0> FFTDynSys;

typedef FFTDynSys::ModesContainerType ModesContainer1D;

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

void basicFPtest(int testNumber, int approach){
  int n = 19, //change FOJ1D stack dimension
      m = 40,
      order = 6;
  Interval nu(0.1), //this was chosen empirically such that set will decrease by 10% at the time 1
           step(0.0001);

  ///begin FOJ initialization for FFT integrator
  capd::jaco::DPDEContainer container;
  ///2.set here the subspace of the initial condition e.g. setToRealValuedOdd means that the initial condition is real valued odd
  container.setToRealValued();
  FOJ1D::initialize(ModesContainer1D::modes2arraySizeStatic(m), container);
  ///end FOJ initialization

  clock_t start, end;

  FFTDynSys dynsys( n, m, step, order, PI, nu);

  Index1D idx;
  RealPolynomialBound u_0(n), enclosure(n);
  double diam;  
  capd::auxil::OutputStream log(std::cout, false, true);
  std::stringstream ss;
  if(testNumber == 0){
    ///3.choose the algorithm type, see enum AlgorithmType in FFTDynSys.h file
    diam = 1e-08;
    ss << "test1_Bproj_";
  }
  if(testNumber == 1){
    ///3.choose the algorithm type, see enum AlgorithmType in FFTDynSys.h file
    diam = 1e-07;
    ss << "test2_Bproj_";
  }
  if(testNumber == 2){
    ///3.choose the algorithm type, see enum AlgorithmType in FFTDynSys.h file
    diam = 1e-06;
    ss << "test3_Bproj_";
  }
  if(testNumber == 3){
    ///3.choose the algorithm type, see enum AlgorithmType in FFTDynSys.h file
    diam = 1e-05;
    ss << "test4_Bproj_";
  }
  if(testNumber == 4){
    diam = 1e-04;
    ss << "test5_Bproj_";
  }
  if(testNumber == 5){
    diam = 1e-03;
    ss << "test6_Bproj_";
  }
  if(testNumber == 6){
    diam = 1e-02;
    ss << "test7_Bproj_";
  }
  Interval r(-diam/2, diam/2);
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
  log << "Only a Galerkin projection is being integrated.\n";
  if(approach == 0){
    dynsys.setJetDynSysAlgorithmType(capd::jaco::direct);
    log << "Using the direct approach\n";
  }else{
    if(approach == 1){
      dynsys.setJetDynSysAlgorithmType(capd::jaco::FFT);
      log << "Using the FFT approach\n";
    }else{
      dynsys.setJetDynSysAlgorithmType(capd::jaco::FFTButFirstOrderDirect);
      log << "Using the FFT approach, but the first order normalized derivative is calculated directly to avoid a blowup\n";
    }
  }
  
  int i;
  for(i=1; i <= n; ++i){
    u_0.set(Index1D(i), ComplexScalar(r, r));
  } 
  
  //  //end init mo
  C0Set set(u_0, 1);
  IntervalMatrix mx(m, m);
  Interval max;
  set = C0Set(u_0, 1);
  start = clock();
  log << "initial condition (finite dimensional):\n" << u_0 << "\nitegration started, the output below is the Galerkin projection of the set at each timestep\n";
  for(i=0; i < 10000; ++i){    
    set.move(dynsys);
    log << (IntervalVector)set << "\n";
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
    std::cerr << "Number of arguments wrong.\nUsage: ./BurgersFPTest projection test_number approach_number\n" << 
                 "\"projection\" is specified then the Galerkin projection is integrated, " <<
                 "\ntest_number is the test index, which corresponds to the index provided in Figure 1.16," <<
                 "\napproach_number is the approach type, 0 - the direct approach, 1 - the FFT approach, 2 - the FFT approach, but the first normalized derivative is calculated directly, in order to avoid blow-ups.\n";
  }else{
    if(strcmp(argv[1], "projection")==0){
      basicFPtest(atoi(argv[2]), atoi(argv[3]));
    }else{
      if(strcmp(argv[1], "inclusion")==0){
//        inclPOTest(atoi(argv[2]), atoi(argv[3]));
//      }else{
        std::cerr << "This test is not implemented.\nUsage: ./BurgersFPTest projection test_number approach_number\n" << 
                 "\"projection\" is specified then the Galerkin projection is integrated, " <<
                 "\ntest_number is the test index, which corresponds to the index provided in Figure 1.16," <<
                 "\napproach_number is the approach type, 0 - the direct approach, 1 - the FFT approach, 2 - the FFT approach, but the first normalized derivative is calculated directly, in order to avoid blow-ups.\n";
      }
    }    
  }
  FOJ1D::destroy();
}

