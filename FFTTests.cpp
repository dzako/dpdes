/*
 * FFTTests.cpp
 *
 *  Created on: Oct 13, 2011
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
#include <cassert>

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
typedef capd::jaco::FFT2DOneComponent<ComplexScalar, ComplexScalar, N, M> FFT2DOneComponent;
typedef typename FFT2DOneComponent::ModesContainerType Modes2DContainer;
typedef Complex1DVector ComplexVector;
typedef Modes2DContainer::VectorType Complex2DVector;
typedef FFT2DOneComponent::IndexRangeType IndexRange;

typedef capd::jaco::Index2DTwoComponents Index2DTwoComponents;
typedef capd::jaco::FFT2D<ComplexScalar, ComplexScalar, N, M> FFT2DTwoComponents;
typedef typename FFT2DTwoComponents::ModesContainerType Modes2DContainerTwoComponents;
typedef Modes2DContainerTwoComponents::VectorType Complex2DVectorTwoComponents;
typedef FFT2DTwoComponents::IndexRangeType IndexRangeTwoComponents;

typedef capd::jaco::FFT2D<FOJ1D, ComplexScalar, N, M> JetFFT2DTwoComponents;
typedef typename JetFFT2DTwoComponents::ModesContainerType JetModes2DContainerTwoComponents;


void printVector(const Modes1DContainer::VectorType& v){
  int i;
  for(i=0; i < v.dimension(); ++i)
    generalDebug << i << ". " << v[i] << "\n";
}


/**
 * for test purpose only. Calculates naively a convolution (using O(N^2) operations).
 */
template< typename ModesContainerT>
ModesContainerT calculateConvolution(int n, const ModesContainerT& u, const ModesContainerT& v) {
  Index2D k, k_1;
  ComplexScalar first, second, r;
  ModesContainerT s(n, ModesContainerT::modes2arraySizeStatic(n));
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

template<typename ModesContainerT>
ModesContainerT calculateTwoComponentProduct(int n, const ModesContainerT& u, const ModesContainerT& grad1, const ModesContainerT& grad2){
  ModesContainerT r(n);

  IndexRangeTwoComponents r1(Index2DTwoComponents::zero()),
                          k1r(Index2DTwoComponents::zero());

  Index2DTwoComponents k, t, k_1, k_1_2, kmk1_2;

  r1.setRange(0, capd::jaco::strong, n, capd::jaco::weak);
  k1r.setRange(0, capd::jaco::strong, n, capd::jaco::weak);

  k[0] = -n;
  k[1] = 0;
  k.l = 0;
  for(; !k.limitReached(r1); k.inc(r1, true)){
    k1r.k = k;
    for(k_1[0] = -n, k_1[1] = -n, k_1.l = 0; !k_1.limitReached(k1r); k_1.inc(k1r, true)){
      if(abs((k-k_1)[0]) <= n && abs((k-k_1)[1]) <=n && abs(k_1[0]) <= n && abs(k_1[1]) <= n){
        k_1_2 = k_1;
        k_1_2.l = 1;
        kmk1_2 = k-k_1;
        kmk1_2.l = 1;
        r[k] += u[k_1] * grad1[k-k_1];
        r[k] += u[k_1_2] * grad1[kmk1_2];
      }

    }
  }

  k[0] = -n;
  k[1] = 0;
  k.l = 1;
  for(; !k.limitReached(r1); k.inc(r1, true)){
    k1r.k = k;
    for(k_1[0] = -n, k_1[1] = -n, k_1.l = 1; !k_1.limitReached(k1r); k_1.inc(k1r, true)){
      if(abs((k-k_1)[0]) <= n && abs((k-k_1)[1]) <=n && abs(k_1[0]) <= n && abs(k_1[1]) <= n){
        k_1_2 = k_1;
        k_1_2.l = 0;
        kmk1_2 = k-k_1;
        kmk1_2.l = 0;
        r[k] += u[k_1_2] * grad2[kmk1_2];
        r[k] += u[k_1] * grad2[k-k_1];
      }
    }
  }

  return r;

}


void fft2dOneComponentTest(){
  int m = 16,
      n = 5,
      d = Modes2DContainer::modes2arraySizeStatic(n);
  
  FFT2DOneComponent fft2d(n, m, Interval::pi());  
    
  
  FFT2DOneComponent::DFT2DGridType s(m), ru(m), rv(m);
  ComplexVector temp(m), tempR(m), pad(m);

  Modes2DContainer u(n, d), v(n, d), r(n, d), c(n, d);
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
    u.set(ind, ComplexScalar(double(rand())/RAND_MAX + diam, double(rand())/RAND_MAX + diam));
    v.set(ind, ComplexScalar(double(rand())/RAND_MAX + diam, double(rand())/RAND_MAX + diam));
  }
  clock_t start, end;   
  
  start = clock();
    
  //std::cout << "conv t=" << end-start << "\n";
  for(i=0; i <= n; ++i)
    for(j=-n; j <= n; ++j){
      ind[0] = j; ind[1] = i;
      generalDebug << "["<<i<<"]["<<j<<"] "<<u[ind]<<"\n";
    }
  
  
  fft2d.transform(u, ru);
  
  //generalDebug << "ru:\n" << ru << "\n";
  
  //ru *= Interval(1.) /Interval(m*m);
  
  //fft2d.inverseTransform(ru, r);
  
  
  fft2d.transform(v, rv);

  s.multiply(ru, rv);
  fft2d.inverseTransform(s, r);
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

  //assign lacking elements
  for(i=-n; i<=0; i++)
    c[Index2D(i, 0)] = r[Index2D(i, 0)];

  //direct result should be subset of r
  assert(c.subset(r));
  std::cout << "assertion c.subset(r) passed for ONE component FFT\n";
}


void fft2dTwoComponentTest(){
  int m = 16,
      n = 5,
      d = Modes2DContainerTwoComponents::modes2arraySizeStatic(n);

  FFT2DTwoComponents fft2d(n, m, Interval::pi());

  FFT2DTwoComponents::DFTGridType s(m), ru(m), rv(m), rGrad1(m), rGrad2(m), scp(m);

  ComplexVector temp(m), tempR(m), pad(m);

  Modes2DContainerTwoComponents u(n, d), v(n, d), r(n, d), c(n, d), grad1(n, d), grad2(n, d), projectedu(n, d);
  Index2DTwoComponents index1, index2, ind, ind_2;
  index1[0] = 1;
  index1[1] = 0;
  index1.l = 1;
  u.set(index1, ComplexScalar(0, -1));

  index2[0] = 0;
  index2[1] = 1;
  index2.l = 0;
  u.set(index2, ComplexScalar(0, -1));
  //these are easier test values

  int i,j;
  Interval diam = Interval(-1e-5, 1e-5);
  IndexRangeTwoComponents range;
  range.setRange(0, capd::jaco::strong, n, capd::jaco::weak);

  for(ind = fft2d.firstModeIndex(range); !ind.limitReached(range); ind.inc(range)){
    //u.set(ind, ComplexScalar(double(rand())/RAND_MAX + diam, double(rand())/RAND_MAX + diam));
  }

  fft2d.CalculateGradients(u, grad1, grad2);

  clock_t start, end;

  FOJ1D::switchToComplexValued();

  fft2d.fastTransform(u, ru); //calculates transforms of Pu, and its gradients
  fft2d.scalarProduct(ru, scp);
  fft2d.inverseExtendedTransform(scp, r, 0);
  fft2d.inverseExtendedTransform(scp, r, 1);
  fft2d.project(u, projectedu);

  c = calculateTwoComponentProduct(n, projectedu, grad1, grad2);

  generalDebug << "r:\n" << r ;
  generalDebug << "c:\n" << c ;

  assert(c.subset(r));
  std::cout << "assertion c.subset(r) passed for TWO component FFT\n";

}

void jet2DFFTTest(){
  int m = 15,
      n = 4,
      d = JetModes2DContainerTwoComponents::modes2arraySizeStatic(n);

  ///begin FOJ initialization for FFT integrator
  capd::jaco::DPDEContainer container;
  ///2.set here the subspace of the initial condition e.g. setToRealValuedOdd means that the initial condition is real valued odd
  container.setToRealValued();
  FOJ1D::initialize(JetModes2DContainerTwoComponents::modes2arraySizeStatic(n), container);
  ///end FOJ initialization

  JetFFT2DTwoComponents fft2d(n, m, Interval::pi());

  JetFFT2DTwoComponents::DFTGridType s(m), ru(m), rv(m), rGrad1(m), rGrad2(m), scp(m);

  ComplexVector temp(m), tempR(m), pad(m);

  Modes2DContainerTwoComponents polybd(n, d);

  JetModes2DContainerTwoComponents u(n, d), v(n, d), r(n, d), c(n, d), grad1(n, d), grad2(n, d), projectedu(n, d);
  Index2DTwoComponents index1, index2, ind, ind_2;
  index1[0] = 1;
  index1[1] = 0;
  index1.l = 1;
  polybd.set(index1, ComplexScalar(0, -1));

  index2[0] = 0;
  index2[1] = 1;
  index2.l = 0;
  polybd.set(index2, ComplexScalar(0, -1));
  //these are easier test values

  int i,j;
  Interval diam = Interval(0.);
  IndexRangeTwoComponents range;
  range.setRange(0, capd::jaco::strong, n, capd::jaco::weak);

  for(ind = fft2d.firstModeIndex(range); !ind.limitReached(range); ind.inc(range)){
    polybd.set(ind, ComplexScalar(double(rand())/RAND_MAX + diam, double(rand())/RAND_MAX + diam));
  }
  u = polybd;

  fft2d.CalculateGradients(u, grad1, grad2);

  clock_t start, end;

  fft2d.fastTransform(u, ru); //calculates transforms of Pu, and its gradients
  fft2d.scalarProduct(ru, scp);
  fft2d.inverseExtendedTransform(scp, r, 0);
  fft2d.inverseExtendedTransform(scp, r, 1);

  fft2d.project(u, projectedu);

  c = calculateTwoComponentProduct(n, projectedu, grad1, grad2);

  generalDebug << "r:\n" << r ;
  generalDebug << "c:\n" << c ;

//  assert(c.subset(r));
  std::cout << "XXX\n";
}


int main(int argc, char * argv[]){

  setLoggers();
  fftDebug.log = true;
  generalDebug.log = true;
  generalDebug2.log = true;
  

  //fft2dOneComponentTest();
  //fft2dTwoComponentTest();
  
  jet2DFFTTest();

}
