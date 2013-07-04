/*
 * PolyBdTests.cpp
 *
 *  Created on: Dec 6, 2011
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
typedef capd::jaco::Index1D Index;
typedef capd::jaco::MaximumNorm<Index> MaximumNorm;
typedef capd::jaco::Even<Interval, Index, MaximumNorm > EvenSubspace;
typedef capd::jaco::RealPolynomialBound<ComplexScalar, Index, 0> RealPolynomialBound;
typedef capd::jaco::RealPolynomialBound<ComplexScalar, Index, 0, EvenSubspace> RealPolynomialBoundEven;

typedef capd::jaco::Burgers<RealPolynomialBound> Burgers;
typedef capd::jaco::GL<RealPolynomialBoundEven> GL;

typedef capd::jaco::FFT1D<ComplexScalar, ComplexScalar, 0, 0, RealPolynomialBound> FFT1D;
typedef capd::jaco::FFT1D<ComplexScalar, ComplexScalar, 0, 0, RealPolynomialBoundEven> FFT1DEven;

typedef capd::jaco::DPDE2<Burgers, FFT1D, 0> DPDE2;
typedef capd::jaco::DPDE3<GL, FFT1DEven, 0> GLPDE;

typedef FFT1D::DFTGridType DFT1DGridType;
typedef FFT1DEven::ModesContainerType Modes1DContainer;

void basicTest(){
  int n = 42,
      m = 90,
      N = 85,
      M = 200;

  RealPolynomialBound polybd(n, N),
                      polybd2(n, N),
                      r(n, N),
                      rPS(n, N),
                      c(n, N);

  FFT1D fft(N, M , Interval::pi());
  DFT1DGridType grid(M), s(M), sPS(M), grid2(M), gridPS(M), grid2PS(M);
  DPDE2 dpde2(n, N, m, M, 1, Interval::pi(), true);

  srand(time(0));
  int i;
  for(i=1; i <= N; ++i){
    polybd[Index(i)] = ComplexScalar(double(rand())/RAND_MAX, double(rand())/RAND_MAX);
    polybd2[Index(i)] = ComplexScalar(double(rand())/RAND_MAX, double(rand())/RAND_MAX);
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

  c = fft.calculateConvolution(polybd, polybd2);

  generalDebug << "r:\n" << r << "\n\n";
  generalDebug << "c:\n" << c << "\n\n";
}

void polyBdTest(){
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

}

Modes1DContainer calculateThirdDegConvolution(int n, const Modes1DContainer& u){
  int k, k_1;
  ComplexScalar first, second, r;
  Modes1DContainer s(n), t(2*n);
  for(k=0; k <= 2*n; k++){
      t[Index(k)] = ComplexScalar(0);
      for(k_1 = -n; k_1 <= n; k_1++){
        if(abs(k-k_1) <= n && abs(k_1) <= n){

            first = u[Index(k_1)];
            second = u[Index(k-k_1)];
            r = first * second;
            t[Index(k)] += r;
        }
    }
  }
  for(k=0; k <= n; k++){
      std::cout << "k=" << k << "\n";
      s[Index(k)] = ComplexScalar(0);
      for(k_1 = -2*n; k_1 <= 2*n; k_1++){
        if(abs(k-k_1) <= n && abs(k_1) <= 2*n){
            first = ((const Modes1DContainer&)t)[Index(k_1)];
            second = u[Index(k-k_1)];
            r = first * second;
            s[Index(k)] += r;
        }
    }
  }
  return s;
}

void polyBdTest2(){
  int n = 42,
      m = 128,
      N = 86,
      M = 270;


  RealPolynomialBoundEven polybd(n, N),
                      polybd2(n, N),
                      r(n, N),
                      rPS(n, N),
                      c(n, 2 * N);

  FFT1DEven fft(N, M , Interval::pi());
  DFT1DGridType grid(M), s(M), sPS(M), grid2(M), gridPS(M), grid2PS(M);
  GLPDE dpde(n, N, m, M, 1, Interval::pi(), true);

  int i;
  for(i=1; i <= n; ++i){
    polybd[Index(i)] = ComplexScalar(double(rand())/RAND_MAX, 0);
    polybd2[Index(i)] = ComplexScalar(double(rand())/RAND_MAX, 0);
  }
  polybd.farTail.m_c = 1e-5;
  polybd.farTail.m_s = 4;
  polybd2.farTail.m_c = 1e-5;
  polybd2.farTail.m_s = 4;

  polybd.setToRealValuedEven();
  dpde.useFFT = true;

  generalDebug << "polybd:\n" << polybd << "\n";
  generalDebug << "polybd2:\n" << polybd2 << "\n";

  dpde.useFFT = true;
  dpde.N(polybd, r);

  generalDebug << "N (fft):\n" << r << "\n";

  r.clean();
  dpde.useFFT = false;
  dpde.N(polybd, r);

  generalDebug << "N (direct):\n" << r << "\n";


  c = calculateThirdDegConvolution(N, polybd);

  generalDebug << "c:\n" << c << "\n";

}

int main(){
  setLoggers();
  polyBdTest2();

}
