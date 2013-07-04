/*  dPDEs - this program is an open research software performing rigorous integration in time of partial differential equations
    Copyright (C) 2010-2013  Jacek Cyranka

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    Please consult the webpage www.cyranka.net,
    or contact me on jcyranka@gmail.com for further details.
*/

/*
 * GLTest.cpp
 *
 *  Created on: Mar 27, 2012
 *      Author: cyranka
 */

#include <time.h>
#include <math.h>
#include <stdlib.h>

#include "capd/filib/Interval.h"
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
//typedef capd::jaco::Odd<Interval, Index1D, MaximumNorm> OddSubspace;
typedef capd::jaco::Even<Interval, Index1D, MaximumNorm> EvenSubspace;
typedef capd::jaco::RealPolynomialBound<ComplexScalar, Index1D, 0, EvenSubspace> RealPolynomialBound;
//typedef capd::jaco::RealPolynomialBound<ComplexScalar, Index1D, 0> RealPolynomialBound;
typedef capd::jaco::Burgers<RealPolynomialBound> Burgers;
typedef capd::jaco::KS<RealPolynomialBound> KS;
typedef capd::jaco::GL<RealPolynomialBound> GL;

typedef capd::jaco::MaximumNorm<Index1D> Norm1D;
//typedef capd::jaco::Real<Interval, Index1D, Norm1D> Real1D;
//typedef capd::jaco::Odd<Interval, Index1D, Norm1D> Real1D;
typedef capd::jaco::Even<Interval, Index1D, MaximumNorm> Real1D;
typedef capd::dynset::C0Rect2RSet<IntervalMatrix> C0Set;
//FFT1D

typedef capd::jaco::FFT1D<ComplexScalar, ComplexScalar, 0, 0, RealPolynomialBound> FFT1D;

typedef capd::jaco::GLDPDE<GL, FFT1D, 0> BurgersDPDE; ///set equation here

typedef capd::jaco::FirstOrderJet<ComplexDerivativePair, 0> FOJ1D;

typedef capd::jaco::ComplexPolyBdJetOptimized<Interval, FOJ1D, Index1D, 0, EvenSubspace> ModesContainer;

typedef capd::jaco::FFT1D<FOJ1D, ComplexScalar, 0, 0, ModesContainer> JetFFT1D;

typedef capd::jaco::GLDPDE<GL, JetFFT1D, 0> JetBurgersDPDE; ///set equation here


typedef capd::jaco::FFTTaylorDynSys<BurgersDPDE, JetBurgersDPDE, 0> DPDEDynSys1D;

typedef DPDEDynSys1D::ModesContainerType ModesContainer1D;


typedef capd::jaco::DPDEInclusionCW<BurgersDPDE, DPDEDynSys1D> DPDEInclusionCW3;
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

void basicPOTest(){
  int n = 23, //change FOJ1D stack dimension
      m = 48,
      order = 5;
  Interval nu(1),
           step(0.001);

  ///begin FOJ initialization for FFT integrator
  FOJ1D::dim = ModesContainer1D::modes2arraySizeStatic(n);
  capd::jaco::DPDEContainer container;
  container.setToRealValuedEven();
//  container.initialConditionSubspaceType = capd::jaco::odd;
  FOJ1D::initialCondition = &container;
  FOJ1D buffer;
  FOJ1D::initialized = 1;
  FOJ1D::buffer = &buffer;
  ///end FOJ initialization

  clock_t start, end;

  BurgersDPDE burgersDpde(n, 2*n+1, m, 2 * m, nu, PI);

  DPDEDynSys1D dynsys( n, m, step, order, PI, nu);

  Interval r(-2e-05,2e-05);

  Index1D idx;

  RealPolynomialBound u_0(n), enclosure(n);
  idx[0] = 1;
  u_0[idx].re =  Interval(1.)+r;
  idx[0] = 2;
  u_0[idx].re =  Interval(1.)+r;

  u_0.setToRealValuedEven();
  dynsys.setInitialConditionSubspace(u_0);

  int i;
  //  //end init mo
  C0Set set(u_0, 1);
  IntervalMatrix mx(m, m);
  Interval max;
  int STEPS = 1000;
         //= 4092;
           //=1;
//  start = clock();
//  for(i=0; i < STEPS; ++i){
//    std::cout << "step #" << i << " ";
//    set.move(basicTaylor);
////    dynsys.move(u_0, u_0);
////    generalDebug << u_0<<"\n";
//  }
//  end = clock();
//  generalDebug << "basicTaylor: " << (IntervalVector)set << "\n";
//  calculateDiams((IntervalVector)set, max);
//  generalDebug << "maxDiam: " << max << "\n";
//  std::cout << "time: " << (end-start)<<"\n\n";

//  generalDebug2.log = false;
  set = C0Set(u_0, 1);
  start = clock();
  for(i=0; i < STEPS; ++i){
    std::cout << "step #" << i << " ";
    set.move(dynsys);
    generalDebug << (IntervalVector)set << "\n";
  }
  generalDebug << "dynsys: " << (IntervalVector)set << "\n";
  calculateDiams((IntervalVector)set, max);
  generalDebug << "maxDiam: " << max << "\n";
  end = clock();
  std::cout << "time: "<<(end-start)<<"\n";
}

void diffInclTest(){
  int m = 23, //change FOJ1D stack dimension
      M = 100,
      dftPts = 48,
      dftPts2 = 128,
      order = 5;
  Interval nu(5),
           step(0.001);

  ///begin FOJ initialization for FFT integrator
  FOJ1D::dim = ModesContainer1D::modes2arraySizeStatic(m);
  capd::jaco::DPDEContainer container;
  container.setToRealValuedEven();
//  container.initialConditionSubspaceType = capd::jaco::odd;
  FOJ1D::initialCondition = &container;
  FOJ1D buffer;
  FOJ1D::initialized = 1;
  FOJ1D::buffer = &buffer;
  ///end FOJ initialization

  clock_t start, end;

  BurgersDPDE burgersDpde(m, M, dftPts, dftPts2, nu, PI);
  burgersDpde.useFFT = true;

  DPDEDynSys1D dynsys( m, dftPts, step, order, PI, nu);

  DPDEInclusionCW3 diffIncl(burgersDpde, order, step, MaxNorm());

  Interval r(-2e-05,2e-05);

  Index1D idx;

  RealPolynomialBound u_0(m, M, true), enclosure(m);
  idx[0] = 1;
  u_0[idx].re =  Interval(1.)+r;
  idx[0] = 2;
  u_0[idx].re =  Interval(1.)+r;

  u_0.setToRealValuedEven();
  setC(u_0, 1e-08);
  setS(u_0, 5);

  dynsys.setInitialConditionSubspace(u_0);

  int i;
  //  //end init mo
  InclRect2Set set(u_0, 1);
  IntervalMatrix mx(m, m);
  Interval max;
  int STEPS = 100;
         //= 4092;
           //=1;
//  start = clock();
//  for(i=0; i < STEPS; ++i){
//    std::cout << "step #" << i << " ";
//    set.move(basicTaylor);
////    dynsys.move(u_0, u_0);
////    generalDebug << u_0<<"\n";
//  }
//  end = clock();
//  generalDebug << "basicTaylor: " << (IntervalVector)set << "\n";
//  calculateDiams((IntervalVector)set, max);
//  generalDebug << "maxDiam: " << max << "\n";
//  std::cout << "time: " << (end-start)<<"\n\n";

//  generalDebug2.log = false;
  set = InclRect2Set(u_0);
  start = clock();
  for(i=0; i < STEPS; ++i){
    std::cout << "step #" << i << " ";
    set.move(diffIncl);
    generalDebug << (IntervalVector)set << "\n";
  }
  generalDebug << "dynsys: " << (IntervalVector)set << "\n";
  calculateDiams((IntervalVector)set, max);
  generalDebug << "maxDiam: " << max << "\n";
  end = clock();
  std::cout << "time: "<<(end-start)<<"\n";

}

int main(){
//  MpFloat::setDefaultPrecision(100);
  setLoggers();
  diffInclTest();
}

