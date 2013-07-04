/*
 * PolyBd2DTests.cpp
 *
 *  Created on: Jun 27, 2013
 *      Author: cyranka
 */

#include <time.h>

//#include "capd/filib/Interval.h"
#include "capd/intervals/Interval.hpp"
#include "config.h"

#define DOUBLE long double
#if __FILIB__
  typedef capd::filib::Interval<DOUBLE> Interval;
#else
  typedef capd::intervals::Interval<DOUBLE> Interval;
#endif


#include "ComplexScalar.h"
#include "FirstOrderJet.h"

#include "FFT.h"
#include "Equations.h"
#include "Even.h"
#include "PolyBd.h"
#include "norms.h"
#include "capd/vectalg/Vector.h"

typedef capd::vectalg::Vector<Interval, 0> IntervalVector;
typedef capd::jaco::ComplexScalar<Interval> ComplexScalar;
typedef capd::jaco::Index2D Index;
typedef capd::jaco::MaximumNorm<Index> MaximumNorm;
typedef capd::jaco::Even<Interval, Index, MaximumNorm > EvenSubspace;
typedef capd::jaco::RealPolynomialBound<ComplexScalar, Index, 0> RealPolynomialBound;
typedef capd::jaco::RealPolynomialBound<ComplexScalar, Index, 0, EvenSubspace> RealPolynomialBoundEven;

typedef capd::jaco::Burgers<RealPolynomialBound> Burgers;

typedef capd::jaco::FFT2DOneComponent<ComplexScalar, ComplexScalar, MaximumNorm, 0, 0, RealPolynomialBound> FFT2D;
typedef capd::jaco::FFT2DOneComponent<ComplexScalar, ComplexScalar, MaximumNorm, 0, 0, RealPolynomialBoundEven> FFT2DEven;

typedef capd::jaco::DPDE2<Burgers, FFT2D, 0> DPDE2;

typedef FFT2D::DFTGridType DFT2DGridType;
typedef FFT2DEven::ModesContainerType Modes2DContainer;

void basicTest(){
  int n = 5,
      m = 10,
      N = 12,
      M = 25;

  RealPolynomialBound polybd(n, N),
                      polybd2(n, N),
                      r(n, N),
                      rPS(n, N),
                      c(n, N);

  FFT2D fft(N, M , Interval::pi());
  DFT2DGridType grid(M), s(M), sPS(M), grid2(M), gridPS(M), grid2PS(M);
  DPDE2 dpde2(n, N, m, M, 1, Interval::pi(), true);

  srand(time(0));
  int i, j;
  for(i=-N; i <= N; ++i){
    for(j=-N; j <= N; ++j){
      Index ind = Index(i, j);
      if(ind.upperHalfspace()){
        polybd[ind] = ComplexScalar(double(rand())/RAND_MAX, double(rand())/RAND_MAX);
        polybd2[ind] = ComplexScalar(double(rand())/RAND_MAX, double(rand())/RAND_MAX);
      }
    }
  }

  generalDebug << polybd << "\n\n";

  fft.fastTransform(polybd, gridPS, capd::jaco::phaseShiftRegular);
  fft.fastTransform(polybd, grid, capd::jaco::none);
  fft.fastTransform(polybd2, grid2PS, capd::jaco::phaseShiftRegular);
  fft.fastTransform(polybd2, grid2, capd::jaco::none);

  s.multiply(grid, grid2);
  sPS.multiply(gridPS, grid2PS);

  fft.fastInverseTransform(s, r, capd::jaco::none);
  fft.fastInverseTransform(sPS, rPS, capd::jaco::phaseShiftRegular);
  r *= 0.5;
  rPS *= 0.5;
  r += rPS;


  generalDebug << "r:\n" << r << "\n\n";

}

/*void polyBdTest(){
  int n = 42,
      m = 128,
      N = 85,
      M = 270;

  IntervalVector d1(10000), d2(10000), d3(10000);
  int i, j;
  srand(time(0));

  for(i=0; i < 10000; ++i){
    d1[i] = rand();
    d2[i] = rand();
  }
  clock_t end, start = clock();
  for(i=0; i < 10000; ++i)
    for(j=0; j < 10000; ++j){
      d2[j] * d1[i];
      d1[i] * d2[j];
    }
  end = clock();
  std::cout << end - start << "\n";


  RealPolynomialBound polybd(n, N),
                      polybd2(n, N),
                      r(n, N),
                      rPS(n, N),
                      c(n, 2 * N);

  FFT1D fft(N, M , Interval::pi());
  DFT1DGridType grid(M), s(M), sPS(M), grid2(M), gridPS(M), grid2PS(M);
  DPDE2 dpde(n, N, m, M, 1, Interval::pi(), true);


  for(i=1; i <= N; ++i){
    polybd[Index(i)] = ComplexScalar(double(rand())/RAND_MAX, double(rand())/RAND_MAX);
    polybd2[Index(i)] = ComplexScalar(double(rand())/RAND_MAX, double(rand())/RAND_MAX);
  }
  polybd.farTail.m_c = 1;
  polybd.farTail.m_s = 4;
  polybd2.farTail.m_c = 1;
  polybd2.farTail.m_s = 4;

  generalDebug << "polybd:\n" << polybd << "\n";
  generalDebug << "polybd2:\n" << polybd2 << "\n";

  dpde.CalculateNonlinearTermUsingFFT(polybd, r);

  generalDebug << "r:\n" << r << "\n";

  c = fft.calculateConvolution(polybd, polybd2);

  generalDebug << "c:\n" << c << "\n";

}*/



int main(){
  setLoggers();
  basicTest();

}
