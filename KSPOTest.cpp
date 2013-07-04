/*
 * 
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
  
  //an another periodic orbit
  //nu=0.127 stb
//  idx[0] = 1;
//  u_0[dpde.mode2array(idx, 0)] =  Interval(0.20121)+r;
//  idx[0] = 2;
//  u_0[dpde.mode2array(idx, 0)] = Interval(1.28998)+r;
//  idx[0] = 3;
//  u_0[dpde.mode2array(idx, 0)] = Interval(0.20121)+r;
//  idx[0] = 4;
//  u_0[dpde.mode2array(idx, 0)] = Interval(-0.377866)+r;
//  idx[0] = 5;
//  u_0[dpde.mode2array(idx, 0)] = Interval(-0.0423095)+r;
//  idx[0] = 6;
//  u_0[dpde.mode2array(idx, 0)] = Interval(0.0431619)+r;
//  idx[0] = 7;
//  u_0[dpde.mode2array(idx, 0)] = Interval(0.00694169)+r;
//  idx[0] = 8;
//  u_0[dpde.mode2array(idx, 0)] = Interval(-0.00415743)+r;
//  idx[0] = 9;
//  u_0[dpde.mode2array(idx, 0)] = Interval(-0.000796969)+r;
//  idx[0] = 10;
//  u_0[dpde.mode2array(idx, 0)] = Interval(0.000331914)+r;
  idx[0] = 1;
  u_0[idx].im =  Interval(0.350668)+r;
  idx[0] = 2;
  u_0[idx].im = Interval(0.0252289)+r;
  idx[0] = 3;
  u_0[idx].im = Interval(0.350668)+r;
  idx[0] = 4;
  u_0[idx].im = Interval(-2.27674)+r;
  idx[0] = 5;
  u_0[idx].im = Interval(-1.11533)+r;
  idx[0] = 6;
  u_0[idx].im = Interval(-0.369307)+r;
  idx[0] = 7;
  u_0[idx].im = Interval(0.460388)+r;
  idx[0] = 8;
  u_0[idx].im = Interval(-0.460455)+r;
  idx[0] = 9;
  u_0[idx].im = Interval(-0.311502)+r;
  idx[0] = 10;
  u_0[idx].im = Interval(-0.14496)+r;
  idx[0] = 11;
  u_0[idx].im = Interval(0.0510489)+r;
  idx[0] = 12;
  u_0[idx].im = Interval(-0.021659)+r;
  idx[0] = 13;
  u_0[idx].im = Interval(-0.0341329)+r;
  idx[0] = 14;
  u_0[idx].im = Interval(-0.0261352)+r;
  idx[0] = 15;
  u_0[idx].im = Interval(0.00130753)+r;
  idx[0] = 16;
  u_0[idx].im = Interval(8.74462e-005)+r;
  idx[0] = 17;
  u_0[idx].im = Interval(-0.00211553)+r;
  idx[0] = 18;
  u_0[idx].im = Interval(-0.00289166)+r;
  idx[0] = 19;
  u_0[idx].im = Interval(-0.000500848)+r;
  idx[0] = 20;
  u_0[idx].im = Interval(3.35395e-005)+r;
  idx[0] = 21;
  u_0[idx].im = Interval(-4.42524e-005)+r;
  idx[0] = 22;
  u_0[idx].im = Interval(-0.000228317)+r; 
  idx[0] = 23;
  u_0[idx].im = Interval(-9.03847e-005)+r;
  u_0[24] = -8.072115e-05 + r;// 1.308558e-04*unit;   
  u_0[25] = 4.061472e-08 + r;// 9.198090e-05*unit;   
  u_0[26] = -1.771179e-05 + r;//3.554454e-05*unit; 
  u_0[27] = 2.160790e-08 + r;//2.349743e-05*unit;
  u_0[28] = -7.683240e-06 + r;//1.260077e-05*unit;   
  u_0[29] = 2.596207e-08 + r;//9.292410e-06*unit;
  u_0[30] = -2.010571e-06 + r;//3.919528e-06*unit;
  u_0[31] = 6.457350e-10 + r;//2.345435e-06*unit;
  u_0[32] = -6.027532e-07 + r;//1.234009e-06*unit;
  u_0[33] = -1.827155e-10 + r;//9.444884e-07*unit;
  u_0[34] = -1.912217e-07 + r;//4.313323e-07*unit;
  u_0[35] = 3.132762e-10 + r;//2.868533e-07*unit;
  u_0[36] = -4.913630e-08 + r;//1.266784e-07*unit;
  u_0[37] = -9.639931e-11 + r;//8.683696e-08*unit;
  u_0[38] = -1.765817e-08 + r;//4.377705e-08*unit;
  u_0[39] = 5.956034e-11 + r;//2.952584e-08*unit;
  u_0[40] = -4.319419e-09 + r;//1.356505e-08*unit;
  u_0[41] = -5.973398e-13 + r;//8.370721e-09*unit;
  u_0[42] = -1.504751e-09 + r;//4.220192e-09*unit;
  u_0[43] = 1.344236e-12 + r;//2.861523e-09*unit;
  u_0[44] = -3.968956e-10 + r;//1.380749e-09*unit;  
  u_0[45] = 1.644704e-12 + r;//8.409932e-10*unit;
  u_0[46] = -1.194526e-10 + r;//4.089145e-10*unit;
  u_0[47] = 5.600085e-14 + r;//2.657488e-10*unit;
  u_0[48] = -3.735948e-11 + r;//1.302601e-10*unit;
  u_0[49] = 6.923386e-14 + r;//8.303362e-11*unit;
  u_0[50] = -1.115027e-11 + r;//3.944118e-11*unit;
  u_0[51] = 7.625038e-15 + r;//2.587359e-11*unit;
  u_0[52] = -3.235955e-12 + r;//1.367747e-11*unit;
  u_0[53] = 5.355447e-15 + r;//9.556918e-12*unit;
  u_0[54] = -8.851600e-13 + r;//6.710554e-12*unit;
  u_0[55] = 8.157066e-15 + r;//7.432209e-12*unit;
  u_0[56] = -2.575738e-13 + r;//1.124376e-11*unit;
  u_0[57] = -1.227735e-14 + r;//1.289762e-11*unit;
  u_0[58] = -7.847187e-14 + r;//2.020907e-11*unit;
  u_0[59] = -9.994174e-16 + r;//2.716484e-11*unit;
  u_0[60] = -2.302107e-14 + r;//5.846583e-11*unit;
  u_0[61] = -1.752133e-15 + r;//7.916104e-11*unit;
  u_0[62] = -6.511992e-15 + r;//1.314991e-10*unit;
  u_0[63] = -3.185013e-16 + r;//1.656574e-10*unit;
  u_0[64] = -1.860478e-15 + r;//2.432627e-10*unit;
  u_0[65] = 1.817169e-16 + r;//3.177145e-10*unit;
  u_0[66] = -5.958659e-16 + r;//4.635467e-10*unit;
  u_0[67] = -3.186116e-17 + r;//4.507650e-10*unit;
  u_0[68] = -1.644389e-16 + r;//3.785445e-10*unit;
  u_0[69] = -8.100358e-18 + r;//3.550455e-10*unit;
  u_0.setToRealValuedOdd();
  setC(u_0, 5.039534e+15);
  setS(u_0, 12);
  log << "initial condition (infinite dimensional):\n" << u_0 << "\nitegration started, the output below is the Galerkin projection of the set at each timestep\n";
  double retTime = 4.093009e-001;
  int i;
  InclRect2Set set(u_0);
  IntervalMatrix mx(m, m);
  Interval max;
  int STEPS// = 100;
  		   = 4092;
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
//    if(i*step > retTime - 0.003){
//      generalDebug.log = true;
//      tailDebug.log = true;
//    }
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
  
  //an another periodic orbit
  //nu=0.127 stb
//  idx[0] = 1;
//  u_0[dpde.mode2array(idx, 0)] =  Interval(0.20121)+r;
//  idx[0] = 2;
//  u_0[dpde.mode2array(idx, 0)] = Interval(1.28998)+r;
//  idx[0] = 3;
//  u_0[dpde.mode2array(idx, 0)] = Interval(0.20121)+r;
//  idx[0] = 4;
//  u_0[dpde.mode2array(idx, 0)] = Interval(-0.377866)+r;
//  idx[0] = 5;
//  u_0[dpde.mode2array(idx, 0)] = Interval(-0.0423095)+r;
//  idx[0] = 6;
//  u_0[dpde.mode2array(idx, 0)] = Interval(0.0431619)+r;
//  idx[0] = 7;
//  u_0[dpde.mode2array(idx, 0)] = Interval(0.00694169)+r;
//  idx[0] = 8;
//  u_0[dpde.mode2array(idx, 0)] = Interval(-0.00415743)+r;
//  idx[0] = 9;
//  u_0[dpde.mode2array(idx, 0)] = Interval(-0.000796969)+r;
//  idx[0] = 10;
//  u_0[dpde.mode2array(idx, 0)] = Interval(0.000331914)+r;
  idx[0] = 1;
  u_0[idx].im =  Interval(0.350668)+r;
  idx[0] = 2;
  u_0[idx].im = Interval(0.0252289)+r;
  idx[0] = 3;
  u_0[idx].im = Interval(0.350668)+r;
  idx[0] = 4;
  u_0[idx].im = Interval(-2.27674)+r;
  idx[0] = 5;
  u_0[idx].im = Interval(-1.11533)+r;
  idx[0] = 6;
  u_0[idx].im = Interval(-0.369307)+r;
  idx[0] = 7;
  u_0[idx].im = Interval(0.460388)+r;
  idx[0] = 8;
  u_0[idx].im = Interval(-0.460455)+r;
  idx[0] = 9;
  u_0[idx].im = Interval(-0.311502)+r;
  idx[0] = 10;
  u_0[idx].im = Interval(-0.14496)+r;
  idx[0] = 11;
  u_0[idx].im = Interval(0.0510489)+r;
  idx[0] = 12;
  u_0[idx].im = Interval(-0.021659)+r;
  idx[0] = 13;
  u_0[idx].im = Interval(-0.0341329)+r;
  idx[0] = 14;
  u_0[idx].im = Interval(-0.0261352)+r;
  idx[0] = 15;
  u_0[idx].im = Interval(0.00130753)+r;
  idx[0] = 16;
  u_0[idx].im = Interval(8.74462e-005)+r;
  idx[0] = 17;
  u_0[idx].im = Interval(-0.00211553)+r;
  idx[0] = 18;
  u_0[idx].im = Interval(-0.00289166)+r;
  idx[0] = 19;
  u_0[idx].im = Interval(-0.000500848)+r;
  idx[0] = 20;
  u_0[idx].im = Interval(3.35395e-005)+r;
  idx[0] = 21;
  u_0[idx].im = Interval(-4.42524e-005)+r;
  idx[0] = 22;
  u_0[idx].im = Interval(-0.000228317)+r;
  idx[0] = 23;
  u_0[idx].im = Interval(-9.03847e-005)+r;
  
  double retTime = 4.093009e-001;
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
    std::cerr << "Number of arguments wrong.\nUsage: ./KSPOTest [projection|inclusion] test_number approach_number\n" << 
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
        std::cerr << "Number of arguments too small.\nUsage: ./KSPOTest [projection|inclusion] test_number approach_number\n" << 
                 "if \"projection\" is specified then the Galerkin projection is integrated, when \"inclusion\" is specified the " <<
                 "whole system is integrated (using the Lohner algorithm for inclusions)" <<
                 "\n test_number is the test index, which corresponds to the index provided in Figure 1.18 and 1.19," <<
                 "\napproach_number is the approach type, 0 - the direct approach, 1 - the FFT approach, 2 - the FFT approach, but the first normalized derivative is calculated directly, in order to avoid blow-ups.\n";
      }
    }    
  }
  FOJ1D::destroy();
}

