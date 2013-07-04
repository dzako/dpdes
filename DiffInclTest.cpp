/*
 * 2DIntegrationTests.cpp
 *
 *  Created on: Sep 14, 2011
 *      Author: cyranka
 */
#include <time.h>
#include <math.h>
#include <stdlib.h>

//#include "capd/filib/Interval.h"
#include "capd/intervals/Interval.hpp"
#include "config.h"


#define DOUBLE long double
#if __FILIB__
  typedef capd::filib::Interval<DOUBLE> Interval;
#else
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
#include "Even.h"
#include "capd/dynset/C0Rect2Set.hpp"
#include "capd/dynset/C0Rect2RSet.hpp"

#include "FFT.h"
#include "Equations.h"
#include "FFTDynSys.h"
#include "PolyBd.h"

#include "DPDEInclusionCW.h"
#include "InclRect2Set.hpp"

#define _D 0
#define _m 3
#define _ni 5
#define _order 3


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
typedef capd::jaco::Odd<Interval, Index1D, MaximumNorm> OddSubspace;
typedef capd::jaco::Even<Interval, Index1D, MaximumNorm> EvenSubspace;
typedef capd::jaco::RealPolynomialBound<ComplexScalar, Index1D, 0, EvenSubspace> RealPolynomialBound;
//typedef capd::jaco::RealPolynomialBound<ComplexScalar, Index1D, 0> RealPolynomialBound;
typedef capd::jaco::Burgers<RealPolynomialBound> Burgers;
typedef capd::jaco::KS<RealPolynomialBound> KS;
typedef capd::jaco::SH<RealPolynomialBound> SH;

typedef capd::jaco::MaximumNorm<Index1D> Norm1D;
//typedef capd::jaco::Real<Interval, Index1D, Norm1D> Real1D;
typedef capd::jaco::Odd<Interval, Index1D, Norm1D> Real1D;
typedef capd::dynset::C0Rect2RSet<IntervalMatrix> C0Set;
//FFT1D

typedef capd::jaco::FFT1D<ComplexScalar, ComplexScalar, 0, 0, RealPolynomialBound> FFT1D;
typedef capd::jaco::DPDE3<SH, FFT1D, 0> BurgersDPDE; ///set equation here
typedef capd::jaco::FirstOrderJet<ComplexDerivativePair, 0> FOJ1D;
typedef capd::jaco::ComplexPolyBdJetOptimized<Interval, FOJ1D, Index1D, 0, EvenSubspace> ModesContainer;
typedef capd::jaco::FFT1D<FOJ1D, ComplexScalar, 0, 0, ModesContainer> JetFFT1D;
typedef capd::jaco::DPDE3<SH, JetFFT1D, 0> JetBurgersDPDE; ///set equation here

typedef capd::jaco::FFTTaylorDynSys<BurgersDPDE, JetBurgersDPDE, 0> FFTDynSys;
typedef FFTDynSys::ModesContainerType ModesContainer1D;
typedef capd::jaco::DPDEInclusionCW<BurgersDPDE, FFTDynSys> DPDEInclusionCW3;
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

void basicTest(){
  int m = 23, //change FOJ1D stack dimension
      M = 47,
      dftPts = 80,
      dftPts2 = 100,
      order = 5;
  Interval nu(0.032),
           step(0.0001);

  ///begin FOJ initialization for FFT integrator
  capd::jaco::DPDEContainer container;
  container.setToRealValuedOdd();
  FOJ1D::initialize(ModesContainer1D::modes2arraySizeStatic(m), container);
  ///end FOJ initialization

  clock_t start, end;

//  BurgersDPDE burgersDpde(m, M, dftPts, dftPts2, nu, PI);
//  burgersDpde.useFFT = false;

  DPDEInclusionCW3 diffIncl(m, dftPts, M, dftPts2, PI, nu, order, step, MaxNorm());
  Interval r(-2e-05,2e-05);

  Index1D idx;

  RealPolynomialBound u_0(m, M), enclosure(m);

  //begin init modes
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

  u_0.setToRealValuedOdd();

  setC(u_0, 1e-05);
  setS(u_0, 6);

  diffIncl.getDynamicalSystem().setInitialConditionSubspace(u_0);

  diffIncl.getDynamicalSystem().setJetDynSysAlgorithmType(capd::jaco::FFT);
  int i;
  //  //end init mo
  InclRect2Set set(u_0);
  IntervalMatrix mx(m, m);
  Interval max;
  int STEPS = 100;
  std::cout << "INTEGRATION START\n";
  start = clock();
  for(i=0; i < STEPS; ++i){
    std::cout << "step #" << i << " ";
    set.move(diffIncl);
//    burgersDpde.printModes( (IntervalVector)set, generalDebug );
  }
  end = clock();
  std::cout << "clock: " << (end - start) << "\n";
  generalDebug << "dynsys: " << (IntervalVector)set << "\n";
  calculateDiams((IntervalVector)set, max);
  generalDebug << "maxDiam: " << max << "\n";

}

void basicTest2(){
  int m = 10, //change FOJ1D stack dimension
      M = 40,
      dftPts = 25,
      dftPts2 = 90,
      order = 5;
  Interval nu(25),
           step(0.0001);

  ///begin FOJ initialization for FFT integrator
  capd::jaco::DPDEContainer container;
  container.setToRealValuedEven();
  FOJ1D::initialize(ModesContainer1D::modes2arraySizeStatic(m), container);
  ///end FOJ initialization


//  BurgersDPDE burgersDpde(m, M, dftPts, dftPts2, nu, PI);
//  burgersDpde.useFFT = false;

  DPDEInclusionCW3 diffIncl(m, dftPts, M, dftPts2, PI, nu, order, step, MaxNorm());
  Interval r(-2e-05,2e-05);

  Index1D idx;

  RealPolynomialBound u_0(m, M), enclosure(m);

  //begin init modes
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


  u_0[Index1D(1)].re =  Interval(1.11401,1.32634);
  u_0[Index1D(3)].re = Interval(-0.0217897,-0.0109837);
  u_0[Index1D(5)].re = Interval(2.51306e-05,0.000192297);
  u_0[Index1D(7)].re = Interval(-1.50552e-06,3.09592e-07);
  u_0[Index1D(9)].re = Interval(-6.10386e-09,1.21895e-08);
  u_0[Index1D(11)].re = Interval(-1.55895e-10,1.26636e-10);
//  idx[0] = 2;
//  u_0[idx].im = Interval(0.0252289)+r;
//  idx[0] = 3;
//  u_0[idx].im = Interval(0.350668)+r;
//  idx[0] = 4;
//  u_0[idx].im = Interval(-2.27674)+r;
//  idx[0] = 5;
//  u_0[idx].im = Interval(-1.11533)+r;
//  idx[0] = 6;
//  u_0[idx].im = Interval(-0.369307)+r;
//  idx[0] = 7;
//  u_0[idx].im = Interval(0.460388)+r;
//  idx[0] = 8;
//  u_0[idx].im = Interval(-0.460455)+r;
//  idx[0] = 9;
//  u_0[idx].im = Interval(-0.311502)+r;
//  idx[0] = 10;
//  u_0[idx].im = Interval(-0.14496)+r;
//  idx[0] = 11;
//  u_0[idx].im = Interval(0.0510489)+r;
//  idx[0] = 12;
//  u_0[idx].im = Interval(-0.021659)+r;
//  idx[0] = 13;
//  u_0[idx].im = Interval(-0.0341329)+r;
//  idx[0] = 14;
//  u_0[idx].im = Interval(-0.0261352)+r;
//  idx[0] = 15;
//  u_0[idx].im = Interval(0.00130753)+r;
//  idx[0] = 16;
//  u_0[idx].im = Interval(8.74462e-005)+r;
//  idx[0] = 17;
//  u_0[idx].im = Interval(-0.00211553)+r;
//  idx[0] = 18;
//  u_0[idx].im = Interval(-0.00289166)+r;
//  idx[0] = 19;
//  u_0[idx].im = Interval(-0.000500848)+r;
//  idx[0] = 20;
//  u_0[idx].im = Interval(3.35395e-005)+r;
//  idx[0] = 21;
//  u_0[idx].im = Interval(-4.42524e-005)+r;
//  idx[0] = 22;
//  u_0[idx].im = Interval(-0.000228317)+r;
//  idx[0] = 23;
//  u_0[idx].im = Interval(-9.03847e-005)+r;


  setC(u_0, 0.);
  setS(u_0, 6);

  diffIncl.getDynamicalSystem().setJetDynSysAlgorithmType(capd::jaco::FFT);
  
  int i;
  //  //end init mo
  InclRect2Set set(u_0);
  IntervalMatrix mxFFT(ModesContainer1D::modes2arraySizeStatic(m), ModesContainer1D::modes2arraySizeStatic(m)),
                 mxDirect(ModesContainer1D::modes2arraySizeStatic(m), ModesContainer1D::modes2arraySizeStatic(m)),
                 mxR(ModesContainer1D::modes2arraySizeStatic(m), ModesContainer1D::modes2arraySizeStatic(m));
  Interval max;
  int STEPS = 1;
  std::cout << "INTEGRATION START\n";
  for(i=0; i < STEPS; ++i){
    std::cout << "step #" << i << " ";
    set.move(diffIncl);
    //generalDebug << (IntervalVector)set << "\n\n";
    diffIncl.getDynamicalSystem().currentD(mxFFT);
    generalDebug << "matrix FFT:\n" << mxFFT << "\n\n";
  }
  
  set = u_0;
  std::cout << "INTEGRATION START\n";
  for(i=0; i < STEPS; ++i){
    std::cout << "step #" << i << " ";
    set.move(diffIncl);
    //generalDebug << (IntervalVector)set << "\n\n";
    diffIncl.getDynamicalSystem().currentD(mxDirect);
    generalDebug << "matrix Direct:\n" << mxDirect << "\n\n";
  }
    
  generalDebug << "intersection(" << intersection(mxFFT, mxDirect, mxR) << "):\n" << mxR << "\n";
  
//  generalDebug << "dynsys: " << (IntervalVector)set << "\n";
//  calculateDiams((IntervalVector)set, max);
//  generalDebug << "maxDiam: " << max << "\n";

}

int main(){
//  MpFloat::setDefaultPrecision(100);
  setLoggers();
  basicTest2();
}
