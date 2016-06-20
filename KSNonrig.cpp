/*
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
typedef capd::jaco::MaximumNorm<Index1D> Norm1D;
typedef capd::jaco::Odd<Interval, Index1D, Norm1D, IntervalMatrix> OddSubspace;
typedef capd::jaco::RealPolynomialBound<ComplexScalar, Index1D, 0, OddSubspace> RealPolynomialBound;
//typedef capd::jaco::RealPolynomialBound<ComplexScalar, Index1D, 0> RealPolynomialBound;
typedef capd::jaco::Burgers<RealPolynomialBound> Burgers;
typedef capd::jaco::KS<RealPolynomialBound> KS;

//FFT1D

typedef capd::jaco::FFT1D<ComplexScalar, ComplexScalar, 0, 0, RealPolynomialBound> FFT1D;

typedef capd::jaco::DPDE2<KS, FFT1D, 0> DPDE; ///set equation here

typedef capd::jaco::FirstOrderJet<ComplexDerivativePair, 0> FOJ1D;

typedef capd::jaco::ComplexPolyBdJetOptimized<Interval, FOJ1D, Index1D, 0, OddSubspace> ModesContainer;

typedef capd::jaco::FFT1D<FOJ1D, ComplexScalar, 0, 0, ModesContainer> JetFFT1D;

typedef capd::jaco::FFTBasicDynSys<DPDE, 0> FFTDynSys;

typedef FFTDynSys::ModesContainerType ModesContainer1D;

//define nonrigorous variational solver
typedef capd::jaco::FFT1D<FOJ1D, ComplexScalar, 0, 0, ModesContainer> JetFFT1DV;

typedef capd::jaco::DPDE2<KS, JetFFT1DV, 0> DPDEV;

typedef capd::jaco::FFTBasicDynSys<DPDEV, 0> FFTDynSysV;



void basicPOTest(int approach){
  int m = 24, //change FOJ1D stack dimension
      M = 50,
      order = 10;
  Interval nu(0.029),
           step(0.0005);

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
  ss << "nonrig";

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

  capd::auxil::OutputStream log2(std::cout, false, true);
  std::stringstream ss2;
  ss2 << "energy.txt";
  log2.logfile(ss2.str().c_str(), true);

  Index1D idx;
  RealPolynomialBound u_0(m), modes(m), enclosure(m);

  idx[0] = 1;
  u_0[idx].im =  Interval(-1.);


  u_0[Index1D(1)].im = -1.05447758550998e+00;

  u_0[Index1D(2)].im = 4.11248251743269e-02;
  u_0[Index1D(3)].im = 2.54658807658815e+00;
  u_0[Index1D(4)].im = -6.51634713012213e+00;
  u_0[Index1D(5)].im = -5.82913036831442e+00;
  u_0[Index1D(6)].im = -6.37282571173927e+00;
  u_0[Index1D(7)].im = 2.26228048880386e+00;
  u_0[Index1D(8)].im = -3.86319449452104e-01;
  u_0[Index1D(9)].im = -9.71389872638077e-01;
  u_0[Index1D(10)].im = -1.94372664384366e+00;
  u_0[Index1D(11)].im = -4.71133048041808e-01;
  u_0[Index1D(12)].im = -1.14834079057883e-01;
  u_0[Index1D(13)].im = 1.20946843961438e-01;
  u_0[Index1D(14)].im = -1.89938648397725e-01;
  u_0[Index1D(15)].im = -1.33270651763804e-01;
  u_0[Index1D(16)].im = -7.79934571868897e-02;
  u_0[Index1D(17)].im = 7.01513111997249e-03;
  u_0[Index1D(18)].im = -4.00665877855380e-03;
  u_0[Index1D(19)].im = -8.80237319681345e-03;
  u_0[Index1D(20)].im = -1.37632255607465e-02;
  //-4.11253107515612e-03 -8.86790222541490e-04 2.04892693474443e-04 -9.68097455425838e-04 -7.34556760395398e-04 -3.91107163515728e-04 -2.02993518318263e-05 -2.65398448306088e-05 -4.55921907375380e-05 -5.48925255262046e-05


  /*idx[0] = 1;
  u_0[idx].im =  Interval(0.350668);
  idx[0] = 2;
  u_0[idx].im = Interval(0.0252289);
  idx[0] = 3;
  u_0[idx].im = Interval(0.350668);
  idx[0] = 4;
  u_0[idx].im = Interval(-2.27674);
  idx[0] = 5;
  u_0[idx].im = Interval(-1.11533);
  idx[0] = 6;
  u_0[idx].im = Interval(-0.369307);
  idx[0] = 7;
  u_0[idx].im = Interval(0.460388);
  idx[0] = 8;
  u_0[idx].im = Interval(-0.460455);
  idx[0] = 9;
  u_0[idx].im = Interval(-0.311502);
  idx[0] = 10;
  u_0[idx].im = Interval(-0.14496);
  idx[0] = 11;
  u_0[idx].im = Interval(0.0510489);
  idx[0] = 12;
  u_0[idx].im = Interval(-0.021659);
  idx[0] = 13;
  u_0[idx].im = Interval(-0.0341329);
  idx[0] = 14;
  u_0[idx].im = Interval(-0.0261352);
  idx[0] = 15;
  u_0[idx].im = Interval(0.00130753);
  idx[0] = 16;
  u_0[idx].im = Interval(8.74462e-005);
  idx[0] = 17;
  u_0[idx].im = Interval(-0.00211553);
  idx[0] = 18;
  u_0[idx].im = Interval(-0.00289166);
  idx[0] = 19;
  u_0[idx].im = Interval(-0.000500848);
  idx[0] = 20;
  u_0[idx].im = Interval(3.35395e-005);*/
  //idx[0] = 21;
  //u_0[idx].im = Interval(-4.42524e-005);
  //idx[0] = 22;
  //u_0[idx].im = Interval(-0.000228317);
  //idx[0] = 23;
  //u_0[idx].im = Interval(-9.03847e-005);

  IntervalVector u = u_0;

  const int STEPS = 1000000;

  const int DFTP = 400;
  FFT1D fft(m, DFTP, PI);
  FFT1D::DFTGridType grid(DFTP);
  capd::auxil::OutputStream log3(std::cout, false, true);

  Interval time;
  for(int i=0; i < STEPS; ++i){
    dynsys.move(u, u);
    if(i % 10 == 0){
      Interval norm = 0;
      for(int j=1; j < m; j++){
        log << u[j] << " ";
        norm += u[j] * u[j];
      }
      time = i * step;

      modes = u;
      fft.fastTransform(modes, grid);
//      std::stringstream ss;
//      ss << "data/time" << time << "solution.txt";
//      log3.logfile(ss.str().c_str(), true); log3.log = true;
//      log3<< grid ;

      log2 << norm << "\n";
      log << "\n";
    }
  }

  end = clock();
  log << "time: "<<(end-start)/CLOCKS_PER_SEC<<"\n";
}

int main(int argc, char * argv[]){
//  MpFloat::setDefaultPrecision(100);
  setLoggers();

  basicPOTest( atoi(argv[1]) );

  FOJ1D::destroy();
}

