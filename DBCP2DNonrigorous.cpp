/*
 * DBCP2DNonrigorous.cpp
 *
 *  Created on: Oct 12, 2013
 *      Author: cyranka
 */

///this test probably does not make sense in sense of 2D dynamics of DBCP model
///- its purpose is to check the software genericity

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

#define _D 0

typedef capd::vectalg::Matrix<Interval,_D,_D> IntervalMatrix;
typedef capd::vectalg::Vector<Interval, _D> IntervalVector;
typedef capd::vectalg::MaxNorm<IntervalVector, IntervalMatrix> MaxNorm;
typedef capd::vectalg::EuclNorm<IntervalVector, IntervalMatrix> EuclNorm;
typedef capd::jaco::ComplexScalar<Interval> ComplexScalar;

//Integrators 2D

typedef capd::jaco::Index2D Index2D;
typedef capd::jaco::MaximumNorm<Index2D> MaximumNorm;
//typedef capd::jaco::Odd<Interval, Index1D, MaximumNorm> OddSubspace;
typedef capd::jaco::Even<Interval, Index2D, MaximumNorm, IntervalMatrix> EvenSubspace;
typedef capd::jaco::RealPolynomialBound<ComplexScalar, Index2D, 0, EvenSubspace> RealPolynomialBound;
//typedef capd::jaco::RealPolynomialBound<ComplexScalar, Index1D, 0> RealPolynomialBound;
typedef capd::jaco::Burgers<RealPolynomialBound> Burgers;
typedef capd::jaco::KS<RealPolynomialBound> KS;
typedef capd::jaco::GL<RealPolynomialBound> GL;
typedef capd::jaco::SH<RealPolynomialBound> SH;
typedef capd::jaco::CH<RealPolynomialBound> CH;
typedef capd::jaco::DBCP<RealPolynomialBound> DBCP;

//FFT1D
typedef capd::jaco::FFT2DOneComponent<ComplexScalar, ComplexScalar, MaximumNorm, 0, 0, RealPolynomialBound> FFT2D;
typedef capd::jaco::DPDE3<DBCP, FFT2D, 0> SHDPDE; ///set equation here

typedef capd::jaco::FFTBasicDynSys<SHDPDE, 0> FFTDynSys;


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
  int n = 5, //change FOJ1D stack dimension
      m = 30,
      order = 7;
  Interval nu,
           step(0.0001);

  ///1.initialize FOJ, the dimension of the first order jet and the initial condition type defined by a DPDEContainer instance
  ///should be provided, see FOJ1D::initialize description

  clock_t start, end;

  Index2D idx;
  RealPolynomialBound u_0(n), enclosure(n);

  for(int i=0; i <= 5; i++){
    for(int j=0; j <= 5; j++){
      idx[0] = i;
      idx[1] = j;
      if( rand() % 2 == 0 )
        u_0[ idx ].re = 1e-15 ;
      else
        u_0[ idx ].re = -1e-15 ;
    }
  }

  int STEPS;
  capd::auxil::OutputStream log(std::cout, false, true);
  std::stringstream ss;
  if(testNumber == 0){
    //a test case
    ///3.choose the algorithm type, see enum AlgorithmType in FFTDynSys.h file
    //nu = 6.085;
    nu = 3.546241427; //scaling 2 * pi
    //nu = 14.18496571; //scaling pi
    //nu = 56.73986284; //scaling pi/2
    STEPS = 10000000;
    n = 5;
    step = 0.00005;
    ss << "test1_CHproj_";
    srand(time(0));
    //rand() % 2;
    //a case for which the solution goes really close to one attr fixed point, but ends at another fixed point

  }
  if(testNumber == 1){
    //the left panel from DBCP paper p. 3691
    nu = 10.14;
    STEPS = 100000;
    n = 5;
    step = 0.00001;

    ss << "test2_CHproj_";
  }
  if(testNumber == 2){
    //the right panel from DBCP paper p. 3691
    //nu = 26.37; //(without scaling the eigenvalues
    nu = 6.59; //(with scaling)
    STEPS = 400000;
    n = 5;
    step = 0.00001;
    ss << "test3_CHproj_";

  }

  ///3.choose the algorithm type, see enum AlgorithmType in FFTDynSys.h file
  FFTDynSys dynsys( n, m, step, order, PI, nu);
  dynsys.changeStepAdaptively = true;

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

  Interval max;

  start = clock();
  int i;
  log << "initial condition (finite dimensional):\n" << u_0 << "\nitegration started, the output below is the Galerkin projection of the set at each timestep\n";
  IntervalVector u = u_0;

  for(i=0; i < STEPS; ++i){
//    std::cout << "step #" << i << " ";
    dynsys.move(u, u);
    if(i % 100 == 0)
      log << u << "\n";
    if(i % 3000 == 0){
      std::cout << dynsys.step << "\n";
      std::cout << u << "\n";
    }
  }

  log << "\nset at the end: " << u << "\n";
  calculateDiams(u, max);
  log << "max diameter of the set at the end: " << max << "\n";
  end = clock();
  log << "time: "<<(end-start)/CLOCKS_PER_SEC<<"\n";
}


int main(int argc, char * argv[]){
  setLoggers();
  if(argc != 3){
    std::cerr << "Number of arguments wrong.\nUsage: ./CHTest test_number approach_number\n" <<
                 "\ntest_number is the test index, which corresponds to the index provided in Figure 1.21," <<
                 "\napproach_number is the approach type, 0 - the direct approach, 1 - the FFT approach, 2 - the FFT approach, but the first normalized derivative is calculated directly, in order to avoid blow-ups.\n";
  }else{
    basicTest(atoi(argv[1]), atoi(argv[2]));
  }
}
