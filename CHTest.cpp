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

#include "config.h"
//#include "capd/intervals/Interval.hpp"

bool COUNT_OPERATIONS = false;
#include "intervals/Interval.hpp" // for operations count

#define DOUBLE long double
#if __FILIB__
  typedef capd::filib::Interval<DOUBLE> Interval;
#else
  typedef capd::intervals::Interval<DOUBLE> Interval;
  #define PI Interval::pi()
//  typedef DOUBLE Interval;
//  #define PI 3.1415926535897932384626433832795
#endif

// MOCK FUNCTIONS NEEDED FOR NONRIGOROUS INTEGRATOR TO WORK
void setLeftBound( Interval& i, DOUBLE b){
  i.setLeftBound(b);
}
void setRightBound( Interval& i, DOUBLE b){
  i.setRightBound(b);
}

long long unsigned TOTAL_ITER;

#include "capd/dynsys/BasicFadTaylor.h"
#include "capd/dynsys/FadTaylor.h"

#include "ComplexScalar.h"
#include "FirstOrderJet.h"
#include "Odd.h"
#include "Even.h"
#include "Index.h"
#include "Coefficients.h"
#include "Pair.h"

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
typedef capd::jaco::MaximumNorm<Index1D> MaximumNorm;
//typedef capd::jaco::Odd<Interval, Index1D, MaximumNorm> OddSubspace;
typedef capd::jaco::Even<Interval, Index1D, MaximumNorm, IntervalMatrix> EvenSubspace;
typedef capd::jaco::RealPolynomialBound<ComplexScalar, Index1D, 0, EvenSubspace> RealPolynomialBound;
//typedef capd::jaco::RealPolynomialBound<ComplexScalar, Index1D, 0> RealPolynomialBound;
typedef capd::jaco::Burgers<RealPolynomialBound> Burgers;
typedef capd::jaco::KS<RealPolynomialBound> KS;
typedef capd::jaco::GL<RealPolynomialBound> GL;
typedef capd::jaco::SH<RealPolynomialBound> SH;
typedef capd::jaco::CH<RealPolynomialBound> CH;
typedef capd::jaco::DBCP<RealPolynomialBound> DBCP;

typedef capd::dynset::C0Rect2RSet<IntervalMatrix> C0Set;
//FFT1D

typedef capd::jaco::FFT1D<ComplexScalar, ComplexScalar, 0, 0, RealPolynomialBound> FFT1D;

typedef capd::jaco::DPDE3<DBCP, FFT1D, 0> SHDPDE; ///set equation here

typedef capd::jaco::FirstOrderJet<ComplexDerivativePair, 0> FOJ1D;

typedef capd::jaco::ComplexPolyBdJetOptimized<Interval, FOJ1D, Index1D, 0, EvenSubspace> ModesContainer;

typedef capd::jaco::FFT1D<FOJ1D, ComplexScalar, 0, 0, ModesContainer> JetFFT1D;

typedef capd::jaco::DPDE3<DBCP, JetFFT1D, 0> JetSHDPDE; ///set equation here

typedef capd::jaco::FFTTaylorDynSys<SHDPDE, JetSHDPDE, 0> FFTDynSys;

typedef FFTDynSys::ModesContainerType ModesContainer1D;

typedef capd::jaco::DPDEInclusionCW<SHDPDE, FFTDynSys> DPDEInclusionCW3;
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


void basicTest(int testNumber, int approach){
  int n = 22, //change FOJ1D stack dimension
      m = 45,
      order = 7;
  Interval nu,
           step(0.0001);

  ///1.initialize FOJ, the dimension of the first order jet and the initial condition type defined by a DPDEContainer instance
  ///should be provided, see FOJ1D::initialize description

  clock_t start, end;

  int STEPS;
  capd::auxil::OutputStream log(std::cout, false, true);
  std::stringstream ss;

  ///begin FOJ initialization for FFT integrator
  capd::jaco::DPDEContainer container;
  ///2.set here the subspace of the initial condition e.g. setToRealValuedOdd means that the initial condition is real valued odd
  container.setToRealValuedEven();
  FOJ1D::initialize(ModesContainer1D::modes2arraySizeStatic(n), container);
  ///end FOJ initialization


  Interval r(-1e-11, 1e-11);
  Index1D idx;
  RealPolynomialBound u_0(n), enclosure(n);

  if(testNumber == 0){
    //nu = 3.546241427;
    nu = 3.546241427;

    STEPS = 10000;
    step = 0.00007;

    srand(time(0));
    for(int i = 1; i < 5; i++){
      if( rand() % 2 == 0 )
        u_0[ Index1D(i) ].re = 1e-10 + r;
      else
        u_0[ Index1D(i) ].re = -1e-10 + r;
    }
    /*u_0[Index1D(3)] = 0.250095;
    u_0[Index1D(9)] = -0.0030652;
    u_0[Index1D(15)] = 3.72671e-05;*/

    /*u_0 = r;
    u_0[Index1D(0)] = 0.;
    u_0[Index1D(3)] = 0.250095 + r;
    u_0[Index1D(9)] = -0.0030652 + r;
    u_0[Index1D(15)] = 3.72671e-05 + r;*/

    /*u_0[Index1D(0)] = 0.;
    u_0[Index1D(2)] = 0.330288 + r;
    u_0[Index1D(6)] = -0.0161974 + r;
    u_0[Index1D(10)] = 0.000751484 + r;
    u_0[Index1D(14)] = -3.57331e-05 + r;*/

    ss << "test1_CHproj_";
  }
  if(testNumber == 1){

    nu = 8;
    STEPS = 200000;
    step = 0.00001;

    for(int i = 1; i < 5; i++){
      if( rand() % 2 == 0 )
        u_0[ Index1D(i) ].re = 1e-10 + r;
      else
        u_0[ Index1D(i) ].re = -1e-10 + r;
    }

    /*u_0 = r;
    u_0[Index1D(0)] = 0.;
    u_0[Index1D(3)] = 0.250095 + r;
    u_0[Index1D(9)] = -0.0030652 + r;
    u_0[Index1D(15)] = 3.72671e-05 + r;*/

    /*u_0[Index1D(0)] = 0.;
    u_0[Index1D(2)] = 0.330288 + r;
    u_0[Index1D(6)] = -0.0161974 + r;
    u_0[Index1D(10)] = 0.000751484 + r;
    u_0[Index1D(14)] = -3.57331e-05 + r;*/

    ss << "test2_CHproj_";
  }

  ss  << "nu_" << rightBound(nu);
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

  ///3.choose the algorithm type, see enum AlgorithmType in FFTDynSys.h file
  FFTDynSys dynsys( n, m, step, order, PI, nu);
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


  //  //end init mo
  C0Set set(u_0, 1);
  Interval max;

  start = clock();
  int i;
  log << "initial condition (finite dimensional):\n" << u_0 << "\nitegration started, the output below is the Galerkin projection of the set at each timestep\n";
  for(i=0; i < STEPS; ++i){
    set.move(dynsys);
    if(i % 100 == 0){
      //    std::cout << "step #" << i << " ";
      IntervalVector s = (IntervalVector)set;
      log << s << "\n";

      for(int j = 0; j < s.dimension(); j++){
        log << rightBound( diam(s[j]) ) << ", ";
      }
      log << "\n";
    }
  }
  log << "\nset at the end: " << (IntervalVector)set << "\n";
  calculateDiams((IntervalVector)set, max);
  log << "max diameter of the set at the end: " << max << "\n";
  end = clock();
  log << "time: "<<(end-start)/CLOCKS_PER_SEC<<"\n";
  log << "RRmultiplicationsSum=" << RRmultiplicationsSum << "\n";
  log << "RRadditionsSum=" << RRadditionsSum << "\n";
}

void diffInclTest(int testNumber, int approach){
  int m = 22, //change FOJ1D stack dimension
      M = 100,
      dftPts = 50,
      dftPts2 = 500,
      order = 7;
  Interval nu,
           step(0.0001);

  ///1.initialize FOJ, the dimension of the first order jet and the initial condition type defined by a DPDEContainer instance
  ///should be provided, see FOJ1D::initialize description

  clock_t start, end;

  int STEPS;
  capd::auxil::OutputStream log(std::cout, false, true);
  std::stringstream ss;

  ///begin FOJ initialization for FFT integrator
  capd::jaco::DPDEContainer container;
  ///2.set here the subspace of the initial condition e.g. setToRealValuedOdd means that the initial condition is real valued odd
  container.setToRealValuedEven();
  FOJ1D::initialize(ModesContainer1D::modes2arraySizeStatic(m), container);
  ///end FOJ initialization


  Interval r(-1e-11, 1e-11);
  Index1D idx;
  RealPolynomialBound u_0(m, M, container), enclosure(m);

  setC(u_0, 0);
  setS(u_0, 6);

  if(testNumber == 0){
    //nu = 3.546241427;
    nu = 3.546241427;

    STEPS = 10000;
    step = 0.00007;

    srand(time(0));
    for(int i = 1; i < 5; i++){
      if( rand() % 2 == 0 )
        u_0[ Index1D(i) ].re = 1e-10 + r;
      else
        u_0[ Index1D(i) ].re = -1e-10 + r;
    }
    /*u_0[Index1D(3)] = 0.250095;
    u_0[Index1D(9)] = -0.0030652;
    u_0[Index1D(15)] = 3.72671e-05;*/

    /*u_0 = r;
    u_0[Index1D(0)] = 0.;
    u_0[Index1D(3)] = 0.250095 + r;
    u_0[Index1D(9)] = -0.0030652 + r;
    u_0[Index1D(15)] = 3.72671e-05 + r;*/

    /*u_0[Index1D(0)] = 0.;
    u_0[Index1D(2)] = 0.330288 + r;
    u_0[Index1D(6)] = -0.0161974 + r;
    u_0[Index1D(10)] = 0.000751484 + r;
    u_0[Index1D(14)] = -3.57331e-05 + r;*/

    ss << "test1_CHincl_";
  }
  if(testNumber == 1){

    nu = 8;
    STEPS = 20000;
    step = 0.0001;

    for(int i = 1; i < 5; i++){
      if( rand() % 2 == 0 )
        u_0[ Index1D(i) ].re = 1e-10 + r;
      else
        u_0[ Index1D(i) ].re = -1e-10 + r;
    }

    /*u_0 = r;
    u_0[Index1D(0)] = 0.;
    u_0[Index1D(3)] = 0.250095 + r;
    u_0[Index1D(9)] = -0.0030652 + r;
    u_0[Index1D(15)] = 3.72671e-05 + r;*/

    /*u_0[Index1D(0)] = 0.;
    u_0[Index1D(2)] = 0.330288 + r;
    u_0[Index1D(6)] = -0.0161974 + r;
    u_0[Index1D(10)] = 0.000751484 + r;
    u_0[Index1D(14)] = -3.57331e-05 + r;*/

    ss << "test2_CHincl_";
  }


  DPDEInclusionCW3 diffIncl(m, dftPts, M, dftPts2, PI, nu, order, step, MaxNorm());

  ss  << "nu_" << rightBound(nu);
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
  log << "Taylor method order=" << order << ", constant time step=" << step << "\n";
  log << "Galerkin projection dimension m=" << m << ", M_{FFT}=" << dftPts << "\n";
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

  log << "initial condition (infinite dimensional):\n" << u_0 << "\nitegration started, the output below is the Galerkin projection of the set at each timestep\n";
  int i;
  //  //end init mo
  InclRect2Set set(u_0);
  Interval max;

  start = clock();
  for(i=0; i < STEPS; ++i){
//    std::cout << "step #" << i << " ";
    set.move(diffIncl);
    log << (IntervalVector)set << "\n";
  }
  log << "\nset at the end (infinite dimensional): " << set.getPerturbationParams() << "\n";
  calculateDiams((IntervalVector)set, max);
  log << "max diameter of the set at the end: " << max << "\n";
  end = clock();
  log << "time: "<<(end-start)/CLOCKS_PER_SEC<<"\n";
  log << "RRmultiplicationsSum=" << RRmultiplicationsSum << "\n";
  log << "RRadditionsSum=" << RRadditionsSum << "\n";
}

int main(int argc, char * argv[]){
  setLoggers();
  if(argc != 4){
      std::cerr << "Number of arguments wrong.\nUsage: ./CHTest [projection|inclusion] test_number approach_number\n" <<
                   "if \"projection\" is specified then the Galerkin projection is integrated, when \"inclusion\" is specified the " <<
                   "whole system is integrated (using the Lohner algorithm for inclusions)" <<
                   "\ntest_number is the test index, which corresponds to the index provided in Figure 1.21," <<
                   "\napproach_number is the approach type, 0 - the direct approach, 1 - the FFT approach, 2 - the FFT approach, but the first normalized derivative is calculated directly, in order to avoid blow-ups.\n";
    }else{
      if(strcmp(argv[1], "projection") == 0){
        basicTest(atoi(argv[2]), atoi(argv[3]));
      }else{
        if(strcmp(argv[1], "inclusion") == 0){
          diffInclTest(atoi(argv[2]), atoi(argv[3]));
        }else{
          std::cerr << "Number of arguments too small.\nUsage: ./CHTest [projection|inclusion] test_number approach_number\n" <<
                   "if \"projection\" is specified then the Galerkin projection is integrated, when \"inclusion\" is specified the " <<
                   "whole system is integrated (using the Lohner algorithm for inclusions)" <<
                   "\ntest_number is the test index, which corresponds to the index provided in Figure 1.21," <<
                   "\napproach_number is the approach type, 0 - the direct approach, 1 - the FFT approach, 2 - the FFT approach, but the first normalized derivative is calculated directly, in order to avoid blow-ups.\n";
        }
      }
    }
  FOJ1D::destroy();
}
