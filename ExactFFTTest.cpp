/*
 * ExactFFTTests.cpp
 *
 *  Created on: Jun 25, 2013
 *      Author: cyranka
 */

#include <time.h>

//#include "capd/filib/Interval.h"
#include "capd/intervals/Interval.hpp"
#include "config.h"
#include "Pair.h"
#include "FirstOrderJet.h"

#define DOUBLE long double
#if __FILIB__
  typedef capd::filib::Interval<DOUBLE> Interval;
#else
  typedef capd::intervals::Interval<DOUBLE> Interval;
#endif

#include "Odd.h"
#include "ComplexScalar.h"
#include "FirstOrderJet.h"

#include "FFT.h"

#define M 0
#define N 0
#define _m 3
#define _ni 5
#define _order 4

typedef capd::vectalg::Matrix<Interval,N,N> IntervalMatrix;
typedef capd::vectalg::Vector<Interval, M> IntervalVec;
typedef capd::jaco::ComplexScalar<Interval> ComplexScalar;
typedef capd::jaco::FFT1D<ComplexScalar, ComplexScalar, N, M> FFT1D;
typedef FFT1D::ModesContainerType Modes1DContainer;
typedef FFT1D::VectorType Complex1DVector;
typedef FFT1D::IndexType Index1D;

typedef capd::jaco::ComplexDerivativePair<ComplexScalar> ComplexDerivativePair;
typedef capd::jaco::FirstOrderJet<ComplexDerivativePair, 0> FOJ1D;
typedef capd::jaco::MaximumNorm<Index1D> MaximumNorm;
typedef capd::jaco::Odd<Interval, Index1D, MaximumNorm, IntervalMatrix> OddSubspace;
typedef capd::jaco::ComplexPolyBdJetOptimized<Interval, FOJ1D, Index1D, 0, OddSubspace> OddModesContainer;
//typedef FFT1D JetFFT1D;
typedef capd::jaco::FFT1D<FOJ1D, ComplexScalar, N, M, OddModesContainer> JetFFT1D;
//typedef capd::jaco::FFT1D<FOJ1D, ComplexScalar, N, M> JetFFT1D;

typedef JetFFT1D::ModesContainerType Modes1DJetContainer;

typedef JetFFT1D::VectorType IntervalVector;

typedef capd::jaco::Index2D Index2D;
typedef capd::jaco::FFT2D<ComplexScalar, ComplexScalar, capd::jaco::MaximumNorm<Index2D>, N, M> FFT2D;
typedef capd::jaco::FFT2DOneComponent<ComplexScalar, ComplexScalar, capd::jaco::MaximumNorm<Index2D>, N, M> FFT2DOneComponent;
typedef FFT2D::ModesContainerType Modes2DContainer;
typedef Complex1DVector ComplexVector;
typedef Modes2DContainer::VectorType Complex2DVector;


typedef FFT2D::IndexRangeType IndexRange;

void printVector(const Modes1DContainer::VectorType& v){
  int i;
  for(i=0; i < v.dimension(); ++i)
    generalDebug << i << ". " << v[i] << "\n";
}


/**
 * for test purpose only. Calculates naively a convolution (using O(N^2) operations).
 */
Modes2DContainer calculateConvolution(int n, const Modes2DContainer& u, const Modes2DContainer& v) {
  Index2D k, k_1;
  ComplexScalar first, second, r;
  Modes2DContainer s(n, Modes2DContainer::modes2arraySizeStatic(n));
  IndexRange kr(Index2D::zero()),
             k1r(Index2D::zero());
  kr.setRange(0, capd::jaco::strong, n, capd::jaco::weak);
  k1r.setRange(0, capd::jaco::strong, n, capd::jaco::weak);
  k[0] = -n;
  k[1] = 0;
  for(; !k.limitReached(kr); k.inc(kr)){
      s.set(k, ComplexScalar(0));
      k1r.k = k;
      //for(k_1 = k/2, k_1.inc(k1r); !k_1.limitReached(k1r); k_1.inc(k1r)){//this is wrong loop, because all modes should be taken into account
      for(k_1[0] = -n, k_1[1] = -n; !k_1.limitReached(k1r); k_1.inc(k1r)){
        if(kr.squareNorm(k-k_1) <= n*n && kr.squareNorm(k_1) <= n*n){
          first = u[k_1];
          second = v[k-k_1];
          r =  first * second;
          s[k] +=r;
        }
      }
  }
  return s;
}


void fft2dOneComponentTest(){
  int m = 16,
      n = 5,
      d = Modes2DContainer::modes2arraySizeStatic(n);
  
  FFT2DOneComponent fft2d(n, m, Interval::pi());          

  Modes2DContainer u(n, d), v(n, d), r(n, d), r1(n, d), r2(n, d), r3(n, d), r4(n, d), c(n, d), mu(n, d), delta_u(n, d), abs_mu(n, d), mv(n, d), delta_v(n, d), abs_mv(n, d);
  
  FFT2DOneComponent::DFT2DGridType s1(m), s2(m), s3(m), s4(m), ru(m), rv(m), r_mu(m), r_delta_u(m), r_abs_mu(m), r_mv(m), r_delta_v(m), r_abs_mv(m);
  ComplexVector temp(m), tempR(m), pad(m);
  
  Index2D index1, index2, ind;
  index1[0] = 1;
  index1[1] = 0;
    
  u.set(index1, ComplexScalar(0, -1));
  v.set(index1, u[index1]);
  
  index2[0] = -1;
  index2[1] = 0;
  u.set(index2, ComplexScalar(0, 1));
  v.set(index2, u[index2]);

  int i,j;
  Interval diam = Interval(-1e-5, 1e-5);
  IndexRange range;
  range.setRange(0, capd::jaco::strong, n, capd::jaco::weak);

  for(ind = fft2d.firstModeIndex(range); !ind.limitReached(range); ind.inc(range, true)){
   // u.set(ind, ComplexScalar(double(rand())/RAND_MAX + diam, double(rand())/RAND_MAX + diam));
   // v.set(ind, ComplexScalar(double(rand())/RAND_MAX + diam, double(rand())/RAND_MAX + diam));
    u.set(ind, ComplexScalar(double(rand())/RAND_MAX + diam, 0));
    v.set(ind, ComplexScalar(double(rand())/RAND_MAX + diam, 0));
  }
  clock_t start, end;   
  
  start = clock();
    
  //std::cout << "conv t=" << end-start << "\n";
  /*for(i=0; i <= n; ++i)
    for(j=-n; j <= n; ++j){
      ind[0] = j; ind[1] = i;
      generalDebug << "["<<i<<"]["<<j<<"] "<<u[ind]<<"\n";
    }*/
  
  //executing 4 transforms instead of one in order to reduce the overestimations
  u.split(mu, delta_u);
  delta_u.abs_supremum(delta_u);
  mu.abs_supremum(abs_mu);
  v.split(mv, delta_v);
  delta_v.abs_supremum(delta_v);
  mv.abs_supremum(abs_mv);
  generalDebug << "abs_mu:\n" << abs_mu << "\n";
  generalDebug << "delta_u:\n" << delta_u << "\n";
  
  fft2d.transform(mu, r_mu);
  fft2d.transform(delta_u, r_delta_u);
  fft2d.transform(abs_mu, r_abs_mu);
  
  fft2d.transform(mv, r_mv);
  fft2d.transform(delta_v, r_delta_v);
  fft2d.transform(abs_mv, r_abs_mv);
  
  //generalDebug << "ru:\n" << ru << "\n";
  
  //ru *= Interval(1.) /Interval(m*m);
  
  //fft2d.inverseTransform(ru, r);
  
  s1.multiply(r_mu, r_mv);
  s2.multiply(r_abs_mu, r_delta_v);
  s3.multiply(r_delta_u, r_abs_mv);
  s4.multiply(r_delta_u, r_delta_v);   
  
  fft2d.inverseTransform(s1, r1);
  fft2d.inverseTransform(s2, r2);
  fft2d.inverseTransform(s3, r3);
  fft2d.inverseTransform(s4, r4);
  
  
  r = r2;
  r += r3;
  r += r4;
  r *= Interval(-1, 1);
  
  r += r1;
  
  end = clock();

  std::cout << "transform t=" << end-start << "\n";
  start = clock();
  c = calculateConvolution(n, u, v);
  end = clock();

  generalDebug << "c:\n";
  for(i=0; i <= n; ++i)
    for(j=-n; j <= n; ++j){
      ind[0] = j; ind[1] = i;
      generalDebug << "["<<i<<"]["<<j<<"] "<<c[ind]<<"\n";
    }
  generalDebug << "r:\n";
  //std::cout << "conv t=" << end-start << "\n";
  for(i=0; i <= n; ++i)
    for(j=-n; j <= n; ++j){
      ind[0] = j; ind[1] = i;
      generalDebug << "["<<i<<"]["<<j<<"] "<<r[ind]<<"\n";
    }
  
  generalDebug << "r subset c: " << r.subset(c) << "\n";

}


int main(int argc, char * argv[]){

  setLoggers();
  fftDebug.log = true;
  generalDebug2.log = true;
  
  fft2dOneComponentTest();
  
}
