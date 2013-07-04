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
 * 
 *
 *  Created on: May 11, 2013
 *      Author: cyranka
 */
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <cstring>
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

/**ATTENTION: THE PADDING METHOD MUST BE USED IN FFT 
 */
void inclPOTest(int testNumber, int approach){
  int m = 25, //change FOJ1D stack dimension
      M = 75,
      dftPts = 64,
      dftPts2 = 250,
      order = 5;
  
  
  //this are the parameter values considered in [CCP]. The orbit is on the chaotic attractor
  Interval nu(0.02991), 
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
    ss << "test1_KS_unst_incl_";
  }
  if(testNumber == 1){
    ///3.choose the algorithm type, see enum AlgorithmType in FFTDynSys.h file
    diam = 1e-07;
    ss << "test2_KS_unst_incl_";
  }
  if(testNumber == 2){
    ///3.choose the algorithm type, see enum AlgorithmType in FFTDynSys.h file
    diam = 1e-06;
    ss << "test3_KS_unst_incl_";
  }
  if(testNumber == 3){
    ///3.choose the algorithm type, see enum AlgorithmType in FFTDynSys.h file
    diam = 1e-05;
    ss << "test4_KS_unst_incl_";
  }
  if(testNumber == 4){
    diam = 1e-04;
    ss << "test5_KS_unst_incl_";
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
  
  Interval unit(-1, 1);  
  Interval r(-diam/2, diam/2);
  Index1D idx;
  RealPolynomialBound u_0(m, M, container), enclosure(m);
  
  idx[0] = 1;
  u_0[idx].im =  Interval(0.268929)+r;
  idx[0] = 2;
  u_0[idx].im = Interval(0.217826)+r;
  idx[0] = 3;
  u_0[idx].im = Interval(0.268929)+r;
  idx[0] = 4;
  u_0[idx].im = Interval(-1.73394)+r;
  idx[0] = 5;
  u_0[idx].im = Interval(-1.15828)+r;
  idx[0] = 6;
  u_0[idx].im = Interval(-0.745116)+r;
  idx[0] = 7;
  u_0[idx].im = Interval(0.451721)+r;
  idx[0] = 8;
  u_0[idx].im = Interval(-0.217102)+r;
  idx[0] = 9;
  u_0[idx].im = Interval(-0.263587)+r;
  idx[0] = 10;
  u_0[idx].im = Interval(-0.203271)+r;
  idx[0] = 11;
  u_0[idx].im = Interval(0.00102834)+r;
  idx[0] = 12;
  u_0[idx].im = Interval(-0.00070666)+r;
  idx[0] = 13;
  u_0[idx].im = Interval(-0.0109867)+r;
  idx[0] = 14;
  u_0[idx].im = Interval(-0.0267955)+r;
  idx[0] = 15;
  u_0[idx].im = Interval(-0.00765627)+r;
  idx[0] = 16;
  u_0[idx].im = Interval(-0.000913737)+r;
  idx[0] = 17;
  u_0[idx].im = Interval(0.000396578)+r;
  idx[0] = 18;
  u_0[idx].im = Interval(-0.00173248)+r;
  idx[0] = 19;
  u_0[idx].im = Interval(-0.00112683)+r;
  idx[0] = 20;
  u_0[idx].im = Interval(-0.000417885)+r;
  idx[0] = 21;
  u_0[idx].im = Interval(3.42432e-05)+r;
  idx[0] = 22;
  u_0[idx].im = Interval(-5.5054e-05)+r; 
  idx[0] = 23;
  u_0[idx].im = Interval(-8.2428e-05)+r;
  idx[0] = 24;
  u_0[idx].im = Interval(-5.72438e-05)+r;
  idx[0] = 25;
  u_0[idx].im = Interval(-9.69028e-06)+r;
  
  
  u_0[26] = -3.576929e-05 + r;
  u_0[27] = 1.398043e-07 + r;
  u_0[28] = -1.004207e-05 + r;
  u_0[29] = 5.541160e-08 + r;
  u_0[30] = -3.411275e-06 + r;
  u_0[31] = 9.474237e-09 + r;
  u_0[32] = -1.116238e-06 + r;
  u_0[33] = 4.014312e-09 + r;
  u_0[34] = -2.909364e-07 + r;
  u_0[35] = 1.640397e-09 + r;
  u_0[36] = -1.026859e-07 + r;
  u_0[37] = 4.643364e-10 + r;
  u_0[38] = -3.167789e-08 + r;
  u_0[39] = 1.921545e-10 + r;
  u_0[40] = -9.813174e-09 + r;
  u_0[41] = 5.510014e-11 + r;
  u_0[42] = -3.214223e-09 + r;
  u_0[43] = 1.729007e-11 + r;
  u_0[44] = -9.170883e-10 + r;
  u_0[45] = 5.777833e-12 + r;
  u_0[46] = -2.778709e-10 + r;
  u_0[47] = 1.481838e-12 + r;
  u_0[48] = -8.266239e-11 + r;
  u_0[49] = 4.982945e-13 + r;
  u_0[50] = -2.489663e-11 + r;
  u_0[51] = 1.586148e-13 + r;
  u_0[52] = -8.329328e-12 + r;
  u_0[53] = 5.179451e-14 + r;
  u_0[54] = -2.642210e-12 + r;
  u_0[55] = 1.759783e-14 + r;
  u_0[56] = -7.666569e-13 + r;
  u_0[57] = 5.332019e-15 + r;
  u_0[58] = -2.443670e-13 + r;
  u_0[59] = 2.052596e-15 + r;
  u_0[60] = -7.122830e-14 + r;
  u_0[61] = 3.418136e-15 + r;
  u_0[62] = -2.286871e-14 + r;
  u_0[63] = 2.683071e-15 + r;
  u_0[64] = -6.570857e-15 + r;
  u_0[65] = 7.687055e-16 + r;
  u_0[66] = -2.067558e-15 + r;
  u_0[67] = 2.268319e-16 + r;
  u_0[68] = -5.922949e-16 + r;
  u_0[69] = 2.318853e-17 + r;
  u_0[70] = -1.915800e-16 + r;
  u_0[71] = 1.963924e-17 + r;
  u_0[72] = -5.517140e-17 + r;
  u_0[73] = 2.401302e-18 + r;
  u_0[74] = -1.741311e-17 + r;
  u_0[75] = 5.526748e-19 + r;
  
  u_0.setToRealValuedOdd();
  setC(u_0, 9.227553e+15);
  setS(u_0, 12);
  log << "initial condition (infinite dimensional):\n" << u_0 << "\nitegration started, the output below is the Galerkin projection of the set at each timestep\n";
  double retTime = 4.490232e-01;
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
  int m = 25, //change FOJ1D stack dimension
      dftPts = 80,
      order = 5;
  Interval nu(0.02991),
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
  
  FFTDynSys dynsys( m, dftPts, step, order, PI, nu);
  
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
  
  idx[0] = 1;
  u_0[idx].im =  Interval(0.268929)+r;
  idx[0] = 2;
  u_0[idx].im = Interval(0.217826)+r;
  idx[0] = 3;
  u_0[idx].im = Interval(0.268929)+r;
  idx[0] = 4;
  u_0[idx].im = Interval(-1.73394)+r;
  idx[0] = 5;
  u_0[idx].im = Interval(-1.15828)+r;
  idx[0] = 6;
  u_0[idx].im = Interval(-0.745116)+r;
  idx[0] = 7;
  u_0[idx].im = Interval(0.451721)+r;
  idx[0] = 8;
  u_0[idx].im = Interval(-0.217102)+r;
  idx[0] = 9;
  u_0[idx].im = Interval(-0.263587)+r;
  idx[0] = 10;
  u_0[idx].im = Interval(-0.203271)+r;
  idx[0] = 11;
  u_0[idx].im = Interval(0.00102834)+r;
  idx[0] = 12;
  u_0[idx].im = Interval(-0.00070666)+r;
  idx[0] = 13;
  u_0[idx].im = Interval(-0.0109867)+r;
  idx[0] = 14;
  u_0[idx].im = Interval(-0.0267955)+r;
  idx[0] = 15;
  u_0[idx].im = Interval(-0.00765627)+r;
  idx[0] = 16;
  u_0[idx].im = Interval(-0.000913737)+r;
  idx[0] = 17;
  u_0[idx].im = Interval(0.000396578)+r;
  idx[0] = 18;
  u_0[idx].im = Interval(-0.00173248)+r;
  idx[0] = 19;
  u_0[idx].im = Interval(-0.00112683)+r;
  idx[0] = 20;
  u_0[idx].im = Interval(-0.000417885)+r;
  idx[0] = 21;
  u_0[idx].im = Interval(3.42432e-05)+r;
  idx[0] = 22;
  u_0[idx].im = Interval(-5.5054e-05)+r; 
  idx[0] = 23;
  u_0[idx].im = Interval(-8.2428e-05)+r;
  idx[0] = 24;
  u_0[idx].im = Interval(-5.72438e-05)+r;
  idx[0] = 25;
  u_0[idx].im = Interval(-9.69028e-06)+r;
  
  double retTime = 4.490232e-01;
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
    std::cerr << "Number of arguments wrong.\nUsage: ./KSUnstPOTest [projection|inclusion] test_number approach_number\n" << 
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
        std::cerr << "Number of arguments too small.\nUsage: ./KSUnstPOTest [projection|inclusion] test_number approach_number\n" << 
                 "if \"projection\" is specified then the Galerkin projection is integrated, when \"inclusion\" is specified the " <<
                 "whole system is integrated (using the Lohner algorithm for inclusions)" <<
                 "\n test_number is the test index, which corresponds to the index provided in Figure 1.18 and 1.19," <<
                 "\napproach_number is the approach type, 0 - the direct approach, 1 - the FFT approach, 2 - the FFT approach, but the first normalized derivative is calculated directly, in order to avoid blow-ups.\n";
      }
    }    
  }
  FOJ1D::destroy();
}