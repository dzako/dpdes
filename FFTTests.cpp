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
typedef FFT2D::ModesContainerType Modes2DContainer;
typedef Complex1DVector ComplexVector;
typedef Modes2DContainer::VectorType Complex2DVector;


typedef FFT2D::IndexRangeType IndexRange;

void printVector(const Modes1DContainer::VectorType& v){
  int i;
  for(i=0; i < v.dimension(); ++i)
    generalDebug << i << ". " << v[i] << "\n";
}

void pickRandomSet(int n, DOUBLE diam, DOUBLE rateOfDecay, Modes1DContainer& even, Modes1DContainer& odd){
  int i;
  int sign;
  DOUBLE v;
  srand(time(0));
  for(i = 1; i <= n; ++i){
    sign = (rand() % 2 == 0 ? 1 : -1);
    v = double(rand()) / RAND_MAX * (1 / power(i, rateOfDecay));
    even.set(Index1D(i), ComplexScalar(sign*v + Interval(-diam / 2, diam / 2), 0));
    odd.set(Index1D(i), ComplexScalar(0, sign*v + Interval(-diam / 2, diam / 2)));
  }
  even.setToRealValuedEven();
  odd.setToRealValuedOdd();
}

void pickRandomSet(int n, DOUBLE diam, DOUBLE rateOfDecay, Modes1DJetContainer& even, Modes1DJetContainer& odd){
  int i;
  int sign;
  DOUBLE v;
  srand(time(0));
  for(i = 1; i <= n; ++i){
    sign = (rand() % 2 == 0 ? 1 : -1);
    v = double(rand()) / RAND_MAX * (1 / power(i, rateOfDecay));
    even.set(Index1D(i), ComplexScalar(sign*v + Interval(-diam / 2, diam / 2), 0));
    even.setVariationalPartToId(Index1D(i), Index1D(i));
    odd.set(Index1D(i), ComplexScalar(0, sign*v + Interval(-diam / 2, diam / 2)));
    odd.setVariationalPartToId(Index1D(i), Index1D(i));
  }
  even.setToRealValuedEven();
  odd.setToRealValuedOdd();
}

Interval calculateMaxDiffBetweenDiams(const Modes1DContainer& mc1, const Modes1DContainer& mc2, int& coord){
  int i;
  coord = 0;
  DOUBLE maxDiam = 0, t;
  for(i=0; i < mc1.dimension(); i++){
    if((t = (rightBound(diam(mc1.m_upperHalfspace[i].re)) - rightBound(diam(mc2.m_upperHalfspace[i].re))) > 0 ?
        (rightBound(diam(mc1.m_upperHalfspace[i].re)) - rightBound(diam(mc2.m_upperHalfspace[i].re))) :
        - (rightBound(diam(mc1.m_upperHalfspace[i].re)) - rightBound(diam(mc2.m_upperHalfspace[i].re)))) >= maxDiam){
      coord = i;
      maxDiam = t;
    }
  }
  return maxDiam;
}

const int PROJECTIONS = 9;
int n[PROJECTIONS];
int m1[PROJECTIONS];
int m2[PROJECTIONS];
int m3[PROJECTIONS];
int m4[PROJECTIONS];
int m5[PROJECTIONS];

void initN(){
  n[0] = 5;   n[1] = 7;   n[2] = 11;  n[3] = 21;  n[4] = 31;  n[5] = 39;  n[6] = 49;   n[7] = 59;   n[8] = 63;
  m1[0] = 12; m1[1] = 15; m1[2] = 24; m1[3] = 45; m1[4] = 64; m1[5] = 81; m1[6] = 100; m1[7] = 120; m1[8] = 128;
  m2[0] = 18; m2[1] = 24; m2[2] = 36; m2[3] = 72; m2[4] = 96; m2[5] = 120; m2[6] = 150; m2[7] = 180; m2[8] = 192;
  m3[0] = 24; m3[1] = 32; m3[2] = 48; m3[3] = 90; m3[4] = 128; m3[5] = 160; m3[6] = 200; m3[7] = 240; m3[8] = 256;
  m4[0] = 50; m4[1] = 48; m4[2] = 120; m4[3] = 144; m4[4] = 320; m4[5] = 400; m4[6] = 500; m4[7] = 360; m4[8] = 640;
  m5[0] = 50; m5[1] = 64; m5[2] = 120; m5[3] = 192; m5[4] = 320; m5[5] = 400; m5[6] = 500; m5[7] = 480; m5[8] = 640;
}
const int DIAMS = 9;
DOUBLE diams[DIAMS];
void initDiams(){
  diams[0] = 1e-08;
  diams[1] = 1e-05;
  diams[2] = 1e-04;
  diams[3] = 1e-03;
  diams[4] = 1e-02;
  diams[5] = 1e-01;
  diams[6] = 1;
  diams[7] = 2;
  diams[8] = 5;
}

Modes1DContainer calculateThirdDegConvolution(int n, const Modes1DContainer& u){
  int k, k_1;
  ComplexScalar first, second, r;
  Modes1DContainer s(n), t(2*n);
  for(k=0; k <= 2*n; k++){
      t[Index1D(k)] = ComplexScalar(0);
      for(k_1 = -n; k_1 <= n; k_1++){
        if(abs(k-k_1) <= n && abs(k_1) <= n){
            first = u[Index1D(k_1)];
            second = u[Index1D(k-k_1)];
            r = first * second;
            t[Index1D(k)] += r;
        }
    }
    if(k != 0) t[Index1D(-k)] = conjugate(t[Index1D(k)]);
  }
  for(k=0; k <= n; k++){
      s[Index1D(k)] = ComplexScalar(0);
      for(k_1 = -2*n; k_1 <= 2*n; k_1++){
        if(abs(k-k_1) <= n && abs(k_1) <= 2*n){
            first = t[Index1D(k_1)];
            second = u[Index1D(k-k_1)];
            r = first * second;
            s[Index1D(k)] += r;
        }
    }
    s[Index1D(-k)] = conjugate(s[Index1D(k)]);
  }
  return s;
}

Modes1DContainer calculateFourthDegConvolution(int n, const Modes1DContainer& u){
  int k, k_1;
  ComplexScalar first, second, r;
  Modes1DContainer s(2*n), t1(2*n), t2(2*n);
  for(k=0; k <= 2*n; k++){
      t1[Index1D(k)] = ComplexScalar(0);
      for(k_1 = -n; k_1 <= n; k_1++){
        if(abs(k-k_1) <= n && abs(k_1) <= n){

            first = u[Index1D(k_1)];
            second = u[Index1D(k-k_1)];
            r = first * second;
            t1[Index1D(k)] += r;
        }
    }
    if(k != 0) t1[Index1D(-k)] = conjugate(t1[Index1D(k)]);
  }
  for(k=0; k <= 2*n; k++){
      t2[Index1D(k)] = ComplexScalar(0);
      for(k_1 = -2*n; k_1 <= 2*n; k_1++){
        if(abs(k-k_1) <= n && abs(k_1) <= 2*n){

            first = t1[Index1D(k_1)];
            second = u[Index1D(k-k_1)];
            r = first * second;
            t2[Index1D(k)] += r;
        }
    }
    if(k != 0) t2[Index1D(-k)] = conjugate(t2[Index1D(k)]);
  }
  for(k=0; k <= n; k++){
      s[Index1D(k)] = ComplexScalar(0);
      for(k_1 = -2*n; k_1 <= 2*n; k_1++){
        if(abs(k-k_1) <= n && abs(k_1) <= 2*n){

            first = t2[Index1D(k_1)];
            second = u[Index1D(k-k_1)];
            r = first * second;
            s[Index1D(k)] += r;
        }
    }
    if(k != 0) s[Index1D(-k)] = conjugate(s[Index1D(k)]);
  }
  return s;
}

Modes1DContainer calculateFifthDegConvolution(int n, const Modes1DContainer& u){
  int k, k_1;
  ComplexScalar first, second, r;
  Modes1DContainer s(2*n), t1(2*n), t2(4*n), t3(6*n);
  for(k=0; k <= 2*n; k++){
      t1[Index1D(k)] = ComplexScalar(0);
      for(k_1 = -n; k_1 <= n; k_1++){
        if(abs(k-k_1) <= n && abs(k_1) <= n){

            first = u[Index1D(k_1)];
            second = u[Index1D(k-k_1)];
            r = first * second;
            t1[Index1D(k)] += r;
        }
    }
    if(k != 0) t1[Index1D(-k)] = conjugate(t1[Index1D(k)]);
  }
  for(k=0; k <= 3*n; k++){
      t2[Index1D(k)] = ComplexScalar(0);
      for(k_1 = -2*n; k_1 <= 2*n; k_1++){
        if(abs(k-k_1) <= n && abs(k_1) <= 2*n){

            first = t1[Index1D(k_1)];
            second = u[Index1D(k-k_1)];
            r = first * second;
            t2[Index1D(k)] += r;
        }
    }
    if(k != 0) t2[Index1D(-k)] = conjugate(t2[Index1D(k)]);
  }
  for(k=0; k <= 2*n; k++){
      t3[Index1D(k)] = ComplexScalar(0);
      for(k_1 = -2*n; k_1 <= 2*n; k_1++){
        if(abs(k-k_1) <= 2*n && abs(k_1) <= 2*n){

            first = t1[Index1D(k_1)];
            second = t1[Index1D(k-k_1)];
            r = first * second;
            t3[Index1D(k)] += r;
        }
    }
    if(k != 0) t3[Index1D(-k)] = conjugate(t3[Index1D(k)]);
  }
  std::cout << "t3:\n" << t3 << "\n";
  for(k=0; k <= n; k++){
      s[Index1D(k)] = ComplexScalar(0);
      for(k_1 = -2*n; k_1 <= 2*n; k_1++){
        if(abs(k-k_1) <= n && abs(k_1) <= 2*n){

            first = t3[Index1D(k_1)];
            second = u[Index1D(k-k_1)];
            r = first * second;
            s[Index1D(k)] += r;
        }
    }
    if(k != 0) s[Index1D(-k)] = conjugate(s[Index1D(k)]);
  }
  return s;
}

Modes1DContainer calculateSeventhDegConvolution(int n, const Modes1DContainer& u){
  int k, k_1;
  ComplexScalar first, second, r;
  Modes1DContainer s(2*n), t1(2*n), t2(4*n), t3(6*n), t4(2*n);
  for(k=0; k <= 2*n; k++){
      t1[Index1D(k)] = ComplexScalar(0);
      for(k_1 = -n; k_1 <= n; k_1++){
        if(abs(k-k_1) <= n && abs(k_1) <= n){

            first = u[Index1D(k_1)];
            second = u[Index1D(k-k_1)];
            r = first * second;
            t1[Index1D(k)] += r;
        }
    }
    if(k != 0) t1[Index1D(-k)] = conjugate(t1[Index1D(k)]);
  }
  for(k=0; k <= 3*n; k++){
      t2[Index1D(k)] = ComplexScalar(0);
      for(k_1 = -2*n; k_1 <= 2*n; k_1++){
        if(abs(k-k_1) <= n && abs(k_1) <= 2*n){

            first = t1[Index1D(k_1)];
            second = u[Index1D(k-k_1)];
            r = first * second;
            t2[Index1D(k)] += r;
        }
    }
    if(k != 0) t2[Index1D(-k)] = conjugate(t2[Index1D(k)]);
  }
  for(k=0; k <= 4*n; k++){
      t3[Index1D(k)] = ComplexScalar(0);
      for(k_1 = -2*n; k_1 <= 2*n; k_1++){
        if(abs(k-k_1) <= 2*n && abs(k_1) <= 2*n){

            first = t1[Index1D(k_1)];
            second = t1[Index1D(k-k_1)];
            r = first * second;
            t3[Index1D(k)] += r;
        }
    }
    if(k != 0) t3[Index1D(-k)] = conjugate(t3[Index1D(k)]);
  }
  std::cout << "t3:\n" << t3 << "\n";
  
  for(k=0; k <= n; k++){
      s[Index1D(k)] = ComplexScalar(0);
      for(k_1 = -4*n; k_1 <= 4*n; k_1++){
        if(abs(k-k_1) <= 3*n && abs(k_1) <= 4*n){

            first = t3[Index1D(k_1)];
            second = t2[Index1D(k-k_1)];
            r = first * second;
            s[Index1D(k)] += r;
        }
    }
    if(k != 0) s[Index1D(-k)] = conjugate(s[Index1D(k)]);
  }
  
  return s;
}

void singleTest(int numberOfTerms, int fftVariant, int n, int m1, int m2, int m3, int m4, int m5){  
  clock_t start, end;
  Modes1DContainer mcEven(n), mcOdd(n);
  capd::auxil::OutputStream out[DIAMS];
  capd::auxil::OutputStream log(std::cout, false, true);
  std::stringstream ss2;
  ss2  << "terms_" << numberOfTerms << "_n_" << n << "_summary.txt";
  log.logfile(ss2.str().c_str(), true); log.log = true;
  int i;
  for(i=0; i < DIAMS; ++i){
    out[i].show = false;
    out[i].flush = true;
  }
  for(i=0; i < DIAMS; ++i){
    std::stringstream ss;
    ss  << "terms_" << numberOfTerms << "_diam_" << diams[i] << "_n_" << n << ".txt";
    std::cout << "\n\nTEST number of terms: " << numberOfTerms << ", n: " << n << ", diam: " << diams[i] << "\n";
    out[i].logfile(ss.str().c_str(), true); out[i].log = true;
    pickRandomSet(n, diams[i], 2, mcEven, mcOdd);
    time_t rawtime;
    time ( &rawtime );
    out[i] << "The current local time is: " << ctime (&rawtime) << "\n";
    out[i] << "TEST number of terms: " << numberOfTerms << ", n: " << n << ", diam: " << diams[i] << "\n";
    out[i] << "u:\n" << mcOdd << "\n";
    int coord;
    Modes1DContainer& mc(mcOdd);
    Modes1DContainer r1(n), r2(n), r3(n), rc(n), rt(n);
    int m_ps, m_padding;
    switch(numberOfTerms){
      case 2 : m_ps = m1; m_padding = m2; break;
      case 3 : m_ps = m1; m_padding = m3; break;
      case 4 : m_ps = m2; m_padding = m4; break;
      case 5 : m_ps = m2; m_padding = m4; break;
      case 6 : m_ps = m3; m_padding = m4; break;
      case 7 : m_ps = m3; m_padding = m5; break;
      default : m_ps = m3; m_padding = m5;
    }    
    Interval directMaxDiam, paddingMaxDiam, ps1MaxDiam, ps2MaxDiam;
    
    FFT1D fft1(n, m1, Interval::pi()), 
          fft2(n, m_ps, Interval::pi()),
          fft_padding(n, m_padding, Interval::pi());
    start = clock();
    if(numberOfTerms == 2){
      rc = fft1.calculateConvolution(mc, mc);
    }    
    if(numberOfTerms == 3){  
      rc = calculateThirdDegConvolution(n, mc);
    }    
    if(numberOfTerms == 4){
      rc = calculateFourthDegConvolution(n, mc);
    }
    if(numberOfTerms == 5){
      rc = calculateFifthDegConvolution(n, mc);
    }
    if(numberOfTerms == 7){
      rc = calculateSeventhDegConvolution(n, mc);
    }
    end = clock();
    std::cout << end-start << "\n";
    out[i] << "direct convolution:\n" << rc << "\n";
    directMaxDiam = rc.maxDiam(coord);
    std::cout << "direct convolution     MAX DIAM: " << rc.maxDiam(coord);
    out[i] << "direct convolution     MAX DIAM: " << rc.maxDiam(coord);
    std::cout << "IS REACHED AT COORD: " << coord << ", AVG DIAM: " << rc.avgDiam() << "\n";
    out[i] << "IS REACHED AT COORD: " << coord << ", AVG DIAM: " << rc.avgDiam() << "\n";
    log << "diam: " << diams[i] << "\n";
    log << "direct: " << directMaxDiam << "\n";
    
    FFT1D::DFT1DGrid hatu = FFT1D::DFT1DGrid(m_padding);
    FFT1D::DFT1DGrid s = FFT1D::DFT1DGrid(m_padding);
    start = clock();
    fft_padding.fastTransform(mc, hatu, capd::jaco::padding, fftVariant);
    if(numberOfTerms >= 2){
      s.multiply(hatu, hatu);
      if(numberOfTerms >= 3){
        s.multiplyOnly(s, hatu);
        if(numberOfTerms >= 4){
          s.multiplyOnly(s, hatu);
          if(numberOfTerms >= 5){
            s.multiplyOnly(s, hatu);
            if(numberOfTerms >= 6){
              s.multiplyOnly(s, hatu);
              if(numberOfTerms >= 7)
                s.multiplyOnly(s, hatu);
            }
          }
        }
      }
    }else{
      throw std::runtime_error("number of terms is too small.\n");
    }
    fft_padding.fastInverseTransform(s, r3, capd::jaco::padding, fftVariant);
    end = clock();
    out[i] << "padding:\n" << r3 << "\n";
    paddingMaxDiam = r3.maxDiam(coord);
    std::cout << "padding                MAX DIAM: " << r3.maxDiam(coord);
    out[i] << "padding                MAX DIAM: " << r3.maxDiam(coord);
    std::cout << "IS REACHED AT COORD: " << coord << ", AVG DIAM: " << r3.avgDiam() << "\n";
    out[i] << "IS REACHED AT COORD: " << coord << ", AVG DIAM: " << r3.avgDiam() << "\n";
    std::cout << ", MAX DIFF WITH DIRECT: " << calculateMaxDiffBetweenDiams(r3, rc, coord);
    out[i] << ", MAX DIFF WITH DIRECT: " << calculateMaxDiffBetweenDiams(r3, rc, coord);
    std::cout << "IS REACHED AT COORD: " << coord << "\n";
    out[i] << "IS REACHED AT COORD: " << coord << "\n";    
    log << "padding (m=" << m_padding << "): " << paddingMaxDiam << "\n";

    hatu =  FFT1D::DFT1DGrid(m_ps);
    s =  FFT1D::DFT1DGrid(m_ps);
    start = clock();
    fft2.fastTransform(mc, hatu, capd::jaco::phaseShiftOddEven, fftVariant);
    if(numberOfTerms >= 2){
      s.multiply(hatu, hatu);
      if(numberOfTerms >= 3){
        s.multiplyOnly(s, hatu);
        if(numberOfTerms >= 4){
          s.multiplyOnly(s, hatu);
          if(numberOfTerms >= 5){
            s.multiplyOnly(s, hatu);
            if(numberOfTerms >= 6){
              s.multiplyOnly(s, hatu);
              if(numberOfTerms >= 7)
                s.multiplyOnly(s, hatu);
            }
          }
        }
      }
    }else{
      throw std::runtime_error("number of terms is too small.\n");
    }
    fft2.fastInverseTransform(s, r1, capd::jaco::phaseShiftOddEven, fftVariant);
    end = clock();    
    out[i] << "phase shift (odd/even):\n" << r1 << "\n";
    ps1MaxDiam = r1.maxDiam(coord);
    std::cout << "phase shift (odd/even) MAX DIAM: " << r1.maxDiam(coord);
    out[i] << "phase shift (odd/even) MAX DIAM: " << r1.maxDiam(coord);
    std::cout << "IS REACHED AT COORD: " << coord << ", AVG DIAM: " << r1.avgDiam();
    out[i] << "IS REACHED AT COORD: " << coord << ", AVG DIAM: " << r1.avgDiam();
    std::cout << ", MAX DIFF WITH DIRECT: " << calculateMaxDiffBetweenDiams(r1, rc, coord);
    out[i] << ", MAX DIFF WITH DIRECT: " << calculateMaxDiffBetweenDiams(r1, rc, coord);
    std::cout << "IS REACHED AT COORD: " << coord << "\n";
    out[i] << "IS REACHED AT COORD: " << coord << "\n";
    log << "ps1 (m=" << m_ps << "): " << ps1MaxDiam << "\n";
    
    hatu = FFT1D::DFT1DGrid(m1);
    s = FFT1D::DFT1DGrid(m1);
    FFT1D::DFT1DGrid hatu2(m1), hatu3(m1);
    start = clock();
    fft1.fastTransform(mc, hatu, capd::jaco::none, fftVariant);
    if(numberOfTerms >= 2){
      s.multiply(hatu, hatu);
      if(numberOfTerms >= 3){
        s.multiplyOnly(s, hatu);
        if(numberOfTerms >= 4){
          s.multiplyOnly(s, hatu);
          if(numberOfTerms >= 5){
            s.multiplyOnly(s, hatu);
            if(numberOfTerms >= 6){
              s.multiplyOnly(s, hatu);
              if(numberOfTerms >= 7)
                s.multiplyOnly(s, hatu);
            }
          }
        }
      }
    }else{
      throw std::runtime_error("number of terms is too small.\n");
    }
    fft1.fastInverseTransform(s, rt, capd::jaco::none, fftVariant);
    end = clock();    
    start = clock();
    fft1.fastTransform(mc, hatu2, capd::jaco::phaseShiftOddEven, fftVariant); //here i changed to phaseShiftOddEven in order to avoid aliasing
    if(numberOfTerms >= 2){
      s.multiply(hatu2, hatu2);
      if(numberOfTerms >= 3){
        s.multiplyOnly(s, hatu2);
        if(numberOfTerms >= 4){
          s.multiplyOnly(s, hatu2);
          if(numberOfTerms >= 5){
            s.multiplyOnly(s, hatu2);
            if(numberOfTerms >= 6){
              s.multiplyOnly(s, hatu2);
              if(numberOfTerms >= 7)
                s.multiplyOnly(s, hatu2);
            }
          }
        }
      }
    }else{
      throw std::runtime_error("number of terms is too small.\n");
    }
    fft1.fastInverseTransform(s, r2, capd::jaco::phaseShiftOddEven, fftVariant); //here i changed to phaseShiftOddEven in order to avoid aliasing
    end = clock();    
    start = clock();
    fft1.fastTransform(mc, hatu3, capd::jaco::phaseShiftRegular, fftVariant); //here i changed to phaseShiftOddEven in order to avoid aliasing
    if(numberOfTerms >= 2){
      s.multiply(hatu3, hatu3);
      if(numberOfTerms >= 3){
        s.multiplyOnly(s, hatu3);
        if(numberOfTerms >= 4){
          s.multiplyOnly(s, hatu3);
          if(numberOfTerms >= 5){
            s.multiplyOnly(s, hatu3);
            if(numberOfTerms >= 6){
              s.multiplyOnly(s, hatu3);
              if(numberOfTerms >= 7)
                s.multiplyOnly(s, hatu3);
            }
          }
        }
      }
    }else{
      throw std::runtime_error("number of terms is too small.\n");
    }
    fft1.fastInverseTransform(s, r3, capd::jaco::phaseShiftRegular, fftVariant); //here i changed to phaseShiftOddEven in order to avoid aliasing
  
    r2 *= 2.;
    r3 *= 1.;
    rt *= 1.;
    r3 += r2;
    r3 += rt;
    r3 *= 0.25;
    end = clock();
    out[i] << "phase shift (regular):\n" << r3 << "\n";
    ps2MaxDiam = r3.maxDiam(coord);
    std::cout << "phase shift (regular)  MAX DIAM: " << r3.maxDiam(coord);
    out[i] << "phase shift (regular)  MAX DIAM: " << r3.maxDiam(coord);
    std::cout << "IS REACHED AT COORD: " << coord << ", AVG DIAM: " << r3.avgDiam() << "\n";
    out[i] << "IS REACHED AT COORD: " << coord << ", AVG DIAM: " << r3.avgDiam() << "\n";
    std::cout << ", MAX DIFF WITH DIRECT: " << calculateMaxDiffBetweenDiams(r3, rc, coord);
    out[i] << ", MAX DIFF WITH DIRECT: " << calculateMaxDiffBetweenDiams(r3, rc, coord);
    std::cout << "IS REACHED AT COORD: " << coord << "\n";
    out[i] << "IS REACHED AT COORD: " << coord << "\n";
    log << "ps2 (m=" << m1 << "): " << ps2MaxDiam << "\n";
    log << "$" << diams[i] << "$&$" << leftBound(directMaxDiam) << "$&$" << leftBound(paddingMaxDiam) << "$&$" << leftBound(ps1MaxDiam) << "$&$" << leftBound(ps2MaxDiam) << "$\\\\\\\hline\n\n";    
  }
}

void thirdDegSingleTest(int fftVariant, DOUBLE diam, int n, int m1, int m2){
  Modes1DContainer mcEven = Modes1DContainer(n), mcOdd = Modes1DContainer(n);
//  pickRandomSet(n, diam, 2, mcEven, mcOdd);
  mcOdd.set(Index1D(1), ComplexScalar(1, 0));
  mcOdd.set(Index1D(2), ComplexScalar(1, 0));
  mcOdd.subspaceType = capd::jaco::even;
  mcOdd.solutionType = capd::jaco::realValued;
  generalDebug2 << "mcOdd:\n" << mcOdd << "\n";
  int coord;
  Modes1DContainer& mc(mcOdd);
  Modes1DContainer r1(n), r2(n), r3(n), rc(n), rt(n);
  FFT1D fft1(n, m1, Interval::pi()), fft2(n, m2, Interval::pi());
  rc = calculateThirdDegConvolution(n, mc);
  generalDebug2 << "direct convolution:\n" << rc << "\n";
  std::cout << "direct convolution     MAX DIAM: " << rc.maxDiam(coord);
  std::cout << " AT COORD: " << coord << ", AVG DIAM: " << rc.avgDiam() << "\n";
  FFT1D::DFT1DGrid hatu(m1), hatu2(m1), hatu3(m1), s(m1);
  fft1.fastTransform(mc, hatu, capd::jaco::phaseShiftOddEven, fftVariant);
  s.multiply(hatu, hatu);
  s.multiplyOnly(s, hatu);
  generalDebug2 << "s:\n" << s << "\n";
  fft1.fastInverseTransform(s, r1, capd::jaco::phaseShiftOddEven, fftVariant);
  generalDebug2 << "phase shift (odd/even):\n" << r1 << "\n";
  std::cout << "phase shift (odd/even) MAX DIAM: " << r1.maxDiam(coord);
  std::cout << " AT COORD: " << coord << ", AVG DIAM: " << r1.avgDiam();
  std::cout << ", MAX DIFF WITH DIRECT: " << calculateMaxDiffBetweenDiams(r1, rc, coord);
  std::cout << " AT COORD: " << coord << "\n";
  fft1.fastTransform(mc, hatu, capd::jaco::none, fftVariant);
  s.multiply(hatu, hatu);
  s.multiplyOnly(s, hatu);
  generalDebug2 << "s:\n" << s << "\n";
  fft1.fastInverseTransform(s, rt, capd::jaco::none, fftVariant);
  
  fft1.fastTransform(mc, hatu2, capd::jaco::phaseShiftOddEven, fftVariant); //here i changed to phaseShiftOddEven in order to avoid aliasing
  s.multiply(hatu2, hatu2);
  s.multiplyOnly(s, hatu2);
  fft1.fastInverseTransform(s, r2, capd::jaco::phaseShiftOddEven, fftVariant); //here i changed to phaseShiftOddEven in order to avoid aliasing

  fft1.fastTransform(mc, hatu3, capd::jaco::phaseShiftRegular, fftVariant); //here i changed to phaseShiftOddEven in order to avoid aliasing
  s.multiply(hatu3, hatu3);
  s.multiplyOnly(s, hatu3);
  fft1.fastInverseTransform(s, r3, capd::jaco::phaseShiftRegular, fftVariant); //here i changed to phaseShiftOddEven in order to avoid aliasing

  r2 *= 2.;
  r3 *= 1.;
  rt *= 1.;
  r3 += r2;
  r3 += rt;
  r3 *= 0.25;
  generalDebug2 << "phase shift (regular):\n" << r3 << "\n";
  std::cout << "phase shift (regular)  MAX DIAM: " << r2.maxDiam(coord);
  std::cout << " AT COORD: " << coord << ", AVG DIAM: " << r2.avgDiam() << "\n";
  std::cout << ", MAX DIFF WITH DIRECT: " << calculateMaxDiffBetweenDiams(r2, rc, coord);
  std::cout << " AT COORD: " << coord << "\n";
  hatu = FFT1D::DFT1DGrid(m2);
  s = FFT1D::DFT1DGrid(m2);
  fft2.fastTransform(mc, hatu, capd::jaco::padding, fftVariant);
  s.multiply(hatu, hatu);
  s.multiplyOnly(s, hatu);
  generalDebug2 << "s:\n" << s << "\n";
  fft2.fastInverseTransform(s, r3, capd::jaco::padding, fftVariant);
  generalDebug2 << "padding:\n" << r3 << "\n";
  std::cout << "padding                MAX DIAM: " << r3.maxDiam(coord);
  std::cout << " AT COORD: " << coord << ", AVG DIAM: " << r3.avgDiam() << "\n";
  std::cout << ", MAX DIFF WITH DIRECT: " << calculateMaxDiffBetweenDiams(r3, rc, coord);
  std::cout << " AT COORD: " << coord << "\n";
}

void fourthDegSingleTest(int fftVariant, DOUBLE diam, int n, int m1, int m2){
  Modes1DContainer mcEven = Modes1DContainer(n), mcOdd = Modes1DContainer(n);
  pickRandomSet(n, diam, 2, mcEven, mcOdd);
//  mcOdd.set(Index1D(2), ComplexScalar(0, -1));
//  mcOdd.subspaceType = capd::jaco::odd;
//  mcOdd.solutionType = capd::jaco::realValued;
  generalDebug2 << "mcOdd:\n" << mcOdd << "\n";
  int coord;
  Modes1DContainer& mc(mcOdd);
  Modes1DContainer r1(n), r2(n), r3(n), rc(n), rt(n), rt2(n);
  FFT1D fft1(n, m1, Interval::pi()), fft2(n, m2, Interval::pi());
  rc = calculateFourthDegConvolution(n, mc);
  generalDebug2 << "direct convolution:\n" << rc << "\n";
  std::cout << "direct convolution     MAX DIAM: " << rc.maxDiam(coord);
  std::cout << " AT COORD: " << coord << ", AVG DIAM: " << rc.avgDiam() << "\n";
  FFT1D::DFT1DGrid hatu(m1), hatu2(m1), hatu3(m1), s(m1);
  fft1.fastTransform(mc, hatu, capd::jaco::phaseShiftOddEven, fftVariant);
  s.multiply(hatu, hatu);
  s.multiplyOnly(s, hatu);
  s.multiplyOnly(s, hatu);

  generalDebug2 << "s:\n" << s << "\n";
  fft1.fastInverseTransform(s, r1, capd::jaco::phaseShiftOddEven, fftVariant);
  generalDebug2 << "phase shift (odd/even):\n" << r1 << "\n";
  std::cout << "phase shift (odd/even) MAX DIAM: " << r1.maxDiam(coord);
  std::cout << " AT COORD: " << coord << ", AVG DIAM: " << r1.avgDiam();
  std::cout << ", MAX DIFF WITH DIRECT: " << calculateMaxDiffBetweenDiams(r1, rc, coord);
  std::cout << " AT COORD: " << coord << "\n";
  fft1.fastTransform(mc, hatu, capd::jaco::none, fftVariant);
  s.multiply(hatu, hatu);
  s.multiplyOnly(s, hatu);
  s.multiplyOnly(s, hatu);

  generalDebug2 << "s:\n" << s << "\n";
  fft1.fastInverseTransform(s, rt, capd::jaco::none, fftVariant);
  fft1.fastTransform(mc, hatu2, capd::jaco::phaseShiftOddEven, fftVariant); //here i changed to phaseShiftOddEven in order to avoid aliasing
  s.multiply(hatu2, hatu2);
  s.multiplyOnly(s, hatu2);
  s.multiplyOnly(s, hatu2);
  fft1.fastInverseTransform(s, r2, capd::jaco::phaseShiftOddEven, fftVariant); //here i changed to phaseShiftOddEven in order to avoid aliasing

  fft1.fastTransform(mc, hatu3, capd::jaco::phaseShiftRegular, fftVariant); //here i changed to phaseShiftOddEven in order to avoid aliasing
  s.multiply(hatu3, hatu3);
  s.multiplyOnly(s, hatu3);
  s.multiplyOnly(s, hatu3);
  fft1.fastInverseTransform(s, r3, capd::jaco::phaseShiftRegular, fftVariant); //here i changed to phaseShiftOddEven in order to avoid aliasing

  r2 *= 2.;
  r3 *= 1.;
  rt *= 1.;
  r3 += r2;
  r3 += rt;
  r3 *= 0.25;
  generalDebug2 << "phase shift (regular):\n" << r3 << "\n";
  std::cout << "phase shift (regular)  MAX DIAM: " << r2.maxDiam(coord);
  std::cout << " AT COORD: " << coord << ", AVG DIAM: " << r2.avgDiam() << "\n";
  std::cout << ", MAX DIFF WITH DIRECT: " << calculateMaxDiffBetweenDiams(r2, rc, coord);
  std::cout << " AT COORD: " << coord << "\n";
  hatu = FFT1D::DFT1DGrid(m2);
  s = FFT1D::DFT1DGrid(m2);
  fft2.fastTransform(mc, hatu, capd::jaco::padding, fftVariant);
  s.multiply(hatu, hatu);
  s.multiplyOnly(s, hatu);
  s.multiplyOnly(s, hatu);

  generalDebug2 << "s:\n" << s << "\n";

  fft2.fastInverseTransform(s, r3, capd::jaco::padding, fftVariant);
  generalDebug2 << "padding:\n" << r3 << "\n";
  std::cout << "padding                MAX DIAM: " << r3.maxDiam(coord);
  std::cout << " AT COORD: " << coord << ", AVG DIAM: " << r3.avgDiam() << "\n";
  std::cout << ", MAX DIFF WITH DIRECT: " << calculateMaxDiffBetweenDiams(r3, rc, coord);
  std::cout << " AT COORD: " << coord << "\n";
}

void fifthDegSingleTest(int fftVariant, DOUBLE diam, int n, int m1, int m2){
  Modes1DContainer mcEven = Modes1DContainer(n), mcOdd = Modes1DContainer(n);
  pickRandomSet(n, diam, 2, mcEven, mcOdd);
//  mcOdd.set(Index1D(2), ComplexScalar(0, -1));
//  mcOdd.subspaceType = capd::jaco::odd;
//  mcOdd.solutionType = capd::jaco::realValued;
  generalDebug2 << "mcOdd:\n" << mcOdd << "\n";
  int coord;
  Modes1DContainer& mc(mcOdd);
  Modes1DContainer r1(n), r2(n), r3(n), rc(n), rt(n);
  FFT1D fft1(n, m1, Interval::pi()), fft2(n, m2, Interval::pi());
  rc = calculateFifthDegConvolution(n, mc);
  generalDebug2 << "direct convolution:\n" << rc << "\n";
  std::cout << "direct convolution     MAX DIAM: " << rc.maxDiam(coord);
  std::cout << " AT COORD: " << coord << ", AVG DIAM: " << rc.avgDiam() << "\n";
  FFT1D::DFT1DGrid hatu(m1), hatu2(m1), s(m1);
  fft1.fastTransform(mc, hatu, capd::jaco::phaseShiftOddEven, fftVariant);
  s.multiply(hatu, hatu);
  s.multiplyOnly(s, hatu);
  s.multiplyOnly(s, hatu);
  s.multiplyOnly(s, hatu);

  generalDebug2 << "s:\n" << s << "\n";
  fft1.fastInverseTransform(s, r1, capd::jaco::phaseShiftOddEven, fftVariant);
  generalDebug2 << "phase shift (odd/even):\n" << r1 << "\n";
  std::cout << "phase shift (odd/even) MAX DIAM: " << r1.maxDiam(coord);
  std::cout << " AT COORD: " << coord << ", AVG DIAM: " << r1.avgDiam();
  std::cout << ", MAX DIFF WITH DIRECT: " << calculateMaxDiffBetweenDiams(r1, rc, coord);
  std::cout << " AT COORD: " << coord << "\n";
  fft1.fastTransform(mc, hatu, capd::jaco::none, fftVariant);
  s.multiply(hatu, hatu);
  s.multiplyOnly(s, hatu);
  s.multiplyOnly(s, hatu);
  s.multiplyOnly(s, hatu);

  generalDebug2 << "s:\n" << s << "\n";
  fft1.fastInverseTransform(s, rt, capd::jaco::none, fftVariant);
  fft1.fastTransform(mc, hatu2, capd::jaco::phaseShiftOddEven, fftVariant);
  s.multiply(hatu2, hatu2);
  s.multiplyOnly(s, hatu2);
  s.multiplyOnly(s, hatu2);
  s.multiplyOnly(s, hatu2);

  fft1.fastInverseTransform(s, r2, capd::jaco::phaseShiftOddEven, fftVariant);
  r2 *= 0.5;
  rt *= 0.5;
  r2 += rt;
  generalDebug2 << "phase shift (regular):\n" << r2 << "\n";
  std::cout << "phase shift (regular)  MAX DIAM: " << r2.maxDiam(coord);
  std::cout << " AT COORD: " << coord << ", AVG DIAM: " << r2.avgDiam() << "\n";
  std::cout << ", MAX DIFF WITH DIRECT: " << calculateMaxDiffBetweenDiams(r2, rc, coord);
  std::cout << " AT COORD: " << coord << "\n";
  hatu = FFT1D::DFT1DGrid(m2);
  s = FFT1D::DFT1DGrid(m2);
  fft2.fastTransform(mc, hatu, capd::jaco::padding, fftVariant);
  s.multiply(hatu, hatu);
  s.multiplyOnly(s, hatu);
  s.multiplyOnly(s, hatu);
  s.multiplyOnly(s, hatu);

  generalDebug2 << "s:\n" << s << "\n";
  fft2.fastInverseTransform(s, r3, capd::jaco::padding, fftVariant);
  generalDebug2 << "padding:\n" << r3 << "\n";
  std::cout << "padding                MAX DIAM: " << r3.maxDiam(coord);
  std::cout << " AT COORD: " << coord << ", AVG DIAM: " << r3.avgDiam() << "\n";
  std::cout << ", MAX DIFF WITH DIRECT: " << calculateMaxDiffBetweenDiams(r3, rc, coord);
  std::cout << " AT COORD: " << coord << "\n";
}

void seventhDegSingleTest(int fftVariant, DOUBLE diam, int n, int m1, int m2){
  Modes1DContainer mcEven = Modes1DContainer(n), mcOdd = Modes1DContainer(n);
  pickRandomSet(n, diam, 2, mcEven, mcOdd);
//  mcOdd.set(Index1D(2), ComplexScalar(0, -1));
//  mcOdd.subspaceType = capd::jaco::odd;
//  mcOdd.solutionType = capd::jaco::realValued;
  generalDebug2 << "mcOdd:\n" << mcOdd << "\n";
  int coord;
  Modes1DContainer& mc(mcOdd);
  Modes1DContainer r1(n), r2(n), r3(n), rc(n), rt(n), rt2(n);
  FFT1D fft1(n, m1, Interval::pi()), fft2(n, m2, Interval::pi());

  FFT1D::DFT1DGrid hatu(m1), hatu2(m1), hatu3(m1), s(m1);
  fft1.fastTransform(mc, hatu, capd::jaco::phaseShiftOddEven, fftVariant);
  s.multiply(hatu, hatu);
  s.multiplyOnly(s, hatu);
  s.multiplyOnly(s, hatu);
  s.multiplyOnly(s, hatu);
  s.multiplyOnly(s, hatu);
  s.multiplyOnly(s, hatu);

  generalDebug2 << "s:\n" << s << "\n";
  fft1.fastInverseTransform(s, r1, capd::jaco::phaseShiftOddEven, fftVariant);
  generalDebug2 << "phase shift (odd/even):\n" << r1 << "\n";
  std::cout << "phase shift (odd/even) MAX DIAM: " << r1.maxDiam(coord);
  std::cout << " AT COORD: " << coord << ", AVG DIAM: " << r1.avgDiam();
  std::cout << ", MAX DIFF WITH DIRECT: " << calculateMaxDiffBetweenDiams(r1, rc, coord);
  std::cout << " AT COORD: " << coord << "\n";
  fft1.fastTransform(mc, hatu, capd::jaco::none, fftVariant);
  s.multiply(hatu, hatu);
  s.multiplyOnly(s, hatu);
  s.multiplyOnly(s, hatu);
  s.multiplyOnly(s, hatu);
  s.multiplyOnly(s, hatu);
  s.multiplyOnly(s, hatu);

  generalDebug2 << "s:\n" << s << "\n";
  fft1.fastInverseTransform(s, rt, capd::jaco::none, fftVariant);
  fft1.fastTransform(mc, hatu2, capd::jaco::phaseShiftOddEven, fftVariant);
  s.multiply(hatu2, hatu2);
  s.multiplyOnly(s, hatu2);
  s.multiplyOnly(s, hatu2);
  s.multiplyOnly(s, hatu2);
  s.multiplyOnly(s, hatu2);
  s.multiplyOnly(s, hatu2);
  fft1.fastInverseTransform(s, r2, capd::jaco::phaseShiftOddEven, fftVariant);

  fft1.fastTransform(mc, hatu3, capd::jaco::phaseShiftRegular, fftVariant); //here i changed to phaseShiftOddEven in order to avoid aliasing
  s.multiply(hatu3, hatu3);
  s.multiplyOnly(s, hatu3);
  s.multiplyOnly(s, hatu3);
  s.multiplyOnly(s, hatu3);
  s.multiplyOnly(s, hatu3);
  s.multiplyOnly(s, hatu3);
  fft1.fastInverseTransform(s, r3, capd::jaco::phaseShiftRegular, fftVariant); //here i changed to phaseShiftOddEven in order to avoid aliasing

  r2 *= 2.;
  r3 *= 1.;
  rt *= 1.;
  r3 += r2;
  r3 += rt;
  r3 *= 0.25;
  generalDebug2 << "phase shift (regular):\n" << r3 << "\n";
  std::cout << "phase shift (regular)  MAX DIAM: " << r3.maxDiam(coord);
  std::cout << " AT COORD: " << coord << ", AVG DIAM: " << r3.avgDiam() << "\n";
  std::cout << ", MAX DIFF WITH DIRECT: " << calculateMaxDiffBetweenDiams(r3, rc, coord);
  std::cout << " AT COORD: " << coord << "\n";
  hatu = FFT1D::DFT1DGrid(m2);
  s = FFT1D::DFT1DGrid(m2);
  fft2.fastTransform(mc, hatu, capd::jaco::padding, fftVariant);
  s.multiply(hatu, hatu);
  s.multiplyOnly(s, hatu);
  s.multiplyOnly(s, hatu);
  s.multiplyOnly(s, hatu);
  s.multiplyOnly(s, hatu);
  s.multiplyOnly(s, hatu);

  generalDebug2 << "s:\n" << s << "\n";
  fft2.fastInverseTransform(s, r3, capd::jaco::padding, fftVariant);
  generalDebug2 << "padding:\n" << r3 << "\n";
  std::cout << "padding                MAX DIAM: " << r3.maxDiam(coord);
  std::cout << " AT COORD: " << coord << ", AVG DIAM: " << r3.avgDiam() << "\n";
  std::cout << ", MAX DIFF WITH DIRECT: " << calculateMaxDiffBetweenDiams(r3, rc, coord);
  std::cout << " AT COORD: " << coord << "\n";
}

void perturbationsTest(int fftVariant, DOUBLE diam, int n, int m1, int m2){
  Modes1DContainer mcEven = Modes1DContainer(n), mcOdd = Modes1DContainer(n), mcEven2 = Modes1DContainer(n), mcEven3 = Modes1DContainer(n);
  int i;
  for(i=1; i < 4; ++i){
    mcEven.set(Index1D(i), ComplexScalar(1 + Interval(-diam , diam), 0));
    mcEven2.set(Index1D(i), ComplexScalar(1 + Interval(-diam, diam), 0));
  }

  mcEven.subspaceType = capd::jaco::even;
  mcEven.solutionType = capd::jaco::realValued;
  mcEven2.subspaceType = capd::jaco::even;
  mcEven2.solutionType = capd::jaco::realValued;
  mcEven3.subspaceType = capd::jaco::even;
  mcEven3.solutionType = capd::jaco::realValued;

  for(i=(n+1)/2; i <= n; ++i){
    mcEven.set(Index1D(i), ComplexScalar( Interval(-1e-12, 1e-12), 0));
    mcEven3.set(Index1D(i), ComplexScalar( Interval(-1e-12, 1e-12), 0));
  }
  generalDebug2 << "mcEven:\n" << mcEven << "\n";
  generalDebug2 << "mcEven2:\n" << mcEven2 << "\n";
  generalDebug2 << "mcEven3:\n" << mcEven3 << "\n";
  int coord;
  Modes1DContainer r1(n), r2(n), r3(n), rc(n), rt(n);
  FFT1D fft1(n, m1, Interval::pi()), fft2(n, m2, Interval::pi());

  FFT1D::DFT1DGrid hatu(m1), hatu2(m1), hatu3(m1), s(m1), s2(m1);
  fft1.fastTransform(mcEven, hatu, capd::jaco::phaseShiftOddEven, fftVariant);
  fft1.fastTransform(mcEven2, hatu2, capd::jaco::phaseShiftOddEven, fftVariant);
  fft1.fastTransform(mcEven3, hatu3, capd::jaco::phaseShiftOddEven, fftVariant);
  s.multiply(hatu, hatu);
  s2.multiply(hatu2, hatu2);
  generalDebug2 << "s2:\n" << s2 << "\n";
  s -= s2;
  generalDebug2 << "s:\n" << s << "\n";
  fft1.fastInverseTransform(s, r1, capd::jaco::phaseShiftOddEven, fftVariant);
  generalDebug2 << "phase shift (odd/even):\n" << r1 << "\n";
  std::cout << "phase shift (odd/even) MAX DIAM: " << r1.maxDiam(coord);
  std::cout << " AT COORD: " << coord << ", AVG DIAM: " << r1.avgDiam();
  std::cout << ", MAX DIFF WITH DIRECT: " << calculateMaxDiffBetweenDiams(r1, rc, coord);
  std::cout << " AT COORD: " << coord << "\n";

  hatu = FFT1D::DFT1DGrid(m2), hatu2 = FFT1D::DFT1DGrid(m2), hatu3 = FFT1D::DFT1DGrid(m2);
  s = FFT1D::DFT1DGrid(m2), s2 = FFT1D::DFT1DGrid(m2);
  fft2.fastTransform(mcEven, hatu, capd::jaco::padding, fftVariant);
  fft2.fastTransform(mcEven2, hatu2, capd::jaco::padding, fftVariant);
  fft2.fastTransform(mcEven3, hatu3, capd::jaco::padding, fftVariant);
  s.multiply(hatu2, hatu3);
  generalDebug2 << "hatu:\n" << hatu << "\n";
  generalDebug2 << "hatu2:\n" << hatu2 << "\n";
  generalDebug2 << "hatu3:\n" << hatu3 << "\n";
  s *= 2.;
  s2.multiply(hatu3, hatu3);
  s += s2;
  generalDebug2 << "s:\n" << s << "\n";
  fft2.fastInverseTransform(s, r3, capd::jaco::padding, fftVariant);
  generalDebug2 << "padding:\n" << r3 << "\n";
  std::cout << "padding                MAX DIAM: " << r3.maxDiam(coord);
  std::cout << " AT COORD: " << coord << ", AVG DIAM: " << r3.avgDiam() << "\n";
  std::cout << ", MAX DIFF WITH DIRECT: " << calculateMaxDiffBetweenDiams(r3, rc, coord);
  std::cout << " AT COORD: " << coord << "\n";
}

void test(int numberOfTerms, int i){
  Modes1DContainer mcEven, mcOdd, r1, r2, r3, rt, rc;
  FFT1D::DFT1DGrid hatu, hatu2, s;
  std::cout << "================n=" << n[i] << ", m1=" << m1[i] << ", m2=" << m2[i] << "================\n";
//  std::cout << "\n================self-sorting variant test================\n";
//  singleTest(capd::jaco::selfSorting, diam, n[i], m1[i], m2[i]);

// we test only scrambled variant (FFT-B - FFT-C)  
  std::cout << "\n================scrambled variant test================\n";
  singleTest(numberOfTerms, capd::jaco::scrambled, n[i], m1[i], m2[i], m3[i], m4[i], m5[i]);
//  std::cout << "\n================unscrambled variant test================\n";
//  singleTest(capd::jaco::unscrambled, diam, n[i], m1[i], m2[i]);
}

///experimental, not actually used
void wrappingEffectInFFT(){
  int m = 512,
      n = 10;

  FFT1D fft(n, m, Interval::pi());

  IntervalVector u(m), transformedU(m), v(m), transformedV(m);
  ComplexVector complexU(m);

  FFT1D::DFT1DGrid hatu(m), hatv(m), hatuC(m), hatvC(m), hatuR(m), hatvR(m), s(m), s1(m), s2(m), s3(m), s4(m), ss(m);

  int i;
  srand(time(0));

  Modes1DContainer mu(n, 2*n), mv(n, 2*n), uC(n, 2*n), vC(n, 2*n), uR(n, 2*n), vR(n, 2*n), r(n, 2*n), r1(n, 2*n), middle(n, 2*n),
      radius(n, 2*n), rt(n, 2*n), c(n, 2*n);
  Interval diam = Interval(-1, 1);

  for(i=1; i <= n; ++i){
    ComplexScalar s = ComplexScalar(Interval(rand())/RAND_MAX, Interval(rand())/RAND_MAX);
    uC.set(Index1D(i), s);
    uR.set(Index1D(i), ComplexScalar(diam, diam));
    mu.set(Index1D(i), s + ComplexScalar(diam, diam));
    s = ComplexScalar(Interval(rand())/RAND_MAX, Interval(rand())/RAND_MAX);
    vC.set(Index1D(i), s);
    vR.set(Index1D(i), ComplexScalar(diam, diam));
    mv.set(Index1D(i), s + ComplexScalar(diam, diam));
  }

  for(i=0; i < 1; ++i){
    fft.fastTransform(uC, hatuC);
    fft.fastTransform(vC, hatvC);
    fft.fastTransform(uR, hatuR);
    fft.fastTransform(vR, hatvR);
    fft.fastTransform(mu, hatu);
    fft.fastTransform(mv, hatv);
    s.multiply(hatu, hatv);
    s *= ComplexScalar(0.0001);
    s1.multiply(hatuC, hatvC);
    s2.multiply(hatuC, hatvR);
    s3.multiply(hatuR, hatvC);
    s4.multiply(hatuR, hatvR);
    fft.fastInverseTransform(s, r);
    ss = s1+s2+s3+s4;
    fft.fastInverseTransform(s1, rt);
    r1 = rt;
    fft.fastInverseTransform(s2, rt);
    r1 += rt;
    fft.fastInverseTransform(s3, rt);
    r1 += rt;
    fft.fastInverseTransform(s4, rt);
    r1 += rt;
    r1 *= ComplexScalar(0.0001);
  }
  c = fft.calculateConvolution( mu, mv);
  generalDebug << "conv:\n"<<c<<"\n";
  generalDebug << "r:\n" << r << "\n";
  r.split(middle, radius);
  generalDebug << "r1:\n" << r1 << "\n";
  generalDebug << "middle r:\n" << middle << "\nradius r:\n"<<radius<<"\n";
  r1.split(middle, radius);
  generalDebug << "middle r1:\n" << middle << "\nradius r1:\n"<<radius<<"\n";

}

///experimental, not actually used
void overestimateReduction(int n, int m, DOUBLE diam){
  Modes1DContainer mcEven(m), mcOdd(m), r(m), rRad(m), rRad1(m), rRad2(m), rMid(m), mid(m), rad(m), rad1(m), rad2(m);
  pickRandomSet(n, diam, 2, mcEven, mcOdd);
  mcOdd.split(mid, rad);
  rad1 = rad;
  DOUBLE d = 0.5;
  rad1 *= d;
  rad2 = rad;
  rad2 *= (1-d);
  FFT1D::DFT1DGrid hatu(m), hatMid(m), hatRad(m), hatRad1(m), hatRad2(m), s(m);
  FFT1D fft1(n, m, Interval::pi());
  std::cout << "set:\n" << mcOdd << "\n";
  fft1.fastTransform(mid, hatMid, capd::jaco::none, capd::jaco::selfSorting);
//  fft1.fastTransform(rad1, hatRad1, capd::jaco::none, capd::jaco::selfSorting);
//  fft1.fastTransform(rad2, hatRad2, capd::jaco::none, capd::jaco::selfSorting);
  fft1.fastTransform(rad, hatRad, capd::jaco::none, capd::jaco::selfSorting);
  hatMid *= Interval(1) / m;
  hatRad *= Interval(1) / m;
//  hatRad1 *= Interval(1) / m;
//  hatRad2 *= Interval(1) / m;
  fft1.fastInverseTransform(hatMid, rMid, capd::jaco::none, capd::jaco::selfSorting);
  fft1.fastInverseTransform(hatRad, rRad, capd::jaco::none, capd::jaco::selfSorting);
//  fft1.fastInverseTransform(hatRad1, rRad1, capd::jaco::none, capd::jaco::selfSorting);
//  fft1.fastInverseTransform(hatRad2, rRad2, capd::jaco::none, capd::jaco::selfSorting);
  std::cout << "result (mid):\n" << rMid << "\n";
  std::cout << "result (rad):\n" << rRad << "\n";
//  std::cout << "result (rad1):\n" << rRad1 << "\n";
//  std::cout << "result (rad2):\n" << rRad2 << "\n";
}

void jetTest(int i, DOUBLE diam){
  int n_ = n[i];
  ///begin FOJ initialization for FFT integrator
  FOJ1D::dim = Modes1DJetContainer::modes2arraySizeStatic(n_);
  FOJ1D buffer;
  FOJ1D::initialized = 1;
  FOJ1D::buffer = &buffer;
  ///end FOJ initialization
  int m = m1[i];
  IntervalMatrix mat(2 * m, 2 * m), mat2(2 * m, 2 * m);
  Modes1DJetContainer mcEven(n_), mcOdd(n_), r(n_), mid(n_), rad(n_);
  IntervalVec vecr(2 * m), vecm(2 * m);
  JetFFT1D::DFT1DGrid hat(m), s(m);
  pickRandomSet(n_, diam, 2, mcEven, mcOdd);
  std::cout << "r.subT=" << mcOdd.subspaceType << ", r.solT=" << mcOdd.solutionType << "\n";
  JetFFT1D fft1(n_, m, Interval::pi());
  fft1.fastTransform(mcOdd, hat, capd::jaco::phaseShiftOddEven, capd::jaco::scrambled);
  s.multiply(hat, hat);
  generalDebug << "s:\n" << s << "\n";
  fft1.fastInverseTransform(s, r, capd::jaco::phaseShiftOddEven, capd::jaco::scrambled);
  generalDebug << "r:\n" << r << "\n";
  r.monodromyMatrix2(mat);
  generalDebug << "mat:\n" << mat << "\n";
  mcOdd.split(vecm, vecr);
  generalDebug << "vecr:\n" << vecr << "\n";
  vecr = mat * vecr;
  generalDebug << "vecr:\n" << vecr << "\n";
  rad = fft1.calculateConvolution(mcOdd, mcOdd);
  generalDebug << "rad:\n" << rad << "\n";
  rad.monodromyMatrix2(mat2);
  generalDebug << "mat:\n" << mat2 << "\n";
  generalDebug << "diff:\n" << mat - mat2 << "\n";
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

void fft2dTest(){
  int m = 16,
      n = 5,
      d = Modes2DContainer::modes2arraySizeStatic(n);

  FFT2D fft2d(n, m, Interval::pi());

  FFT2D::DFT2DGridType s(m), ru(m), rv(m);
  ComplexVector temp(m), tempR(m), pad(m);

  Modes2DContainer u(n, d), v(n, d), r(n, d), c(n, d);
  Index2D index1, index2, ind;
  index1[0] = 1;
  index1[1] = 1;
  u.set(index1, ComplexScalar(0, -1));
  v.set(index1, u[index1]);
  index2[0] = 1;
  index2[1] = 2;
  u.set(index2, ComplexScalar(0, -1));
  v.set(index2, u[index2]);

  int i,j;
  Interval diam = Interval(-1e-15, 1e-15);
  IndexRange range;
  range.setRange(0, capd::jaco::strong, n, capd::jaco::weak);

  for(ind = fft2d.firstModeIndex(range); !ind.limitReached(range); ind.inc(range, true)){
    u.set(ind, ComplexScalar(double(rand())/RAND_MAX + diam, double(rand())/RAND_MAX + diam));
    v.set(ind, u[ind]);
  }
  clock_t start, end;
  start = clock();
  fft2d.transform(u, ru);
  fft2d.transform(v, rv);

  s.multiply(ru, rv);
  fft2d.inverseExtendedTransform(s, r);
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
  std::cout << "conv t=" << end-start << "\n";
  for(i=0; i <= n; ++i)
    for(j=-n; j <= n; ++j){
      ind[0] = j; ind[1] = i;
      generalDebug << "["<<i<<"]["<<j<<"] "<<r[ind]<<"\n";
    }

}

///experimental, not actually used
void basicEfficiencyTest(){
  clock_t start = clock();
  srand(time(0));
  Interval g = double(rand()) / RAND_MAX,
         m = double(rand()) / RAND_MAX;
  int f = 0;
  std::cout << "g=" << g << ", m=" << m << "\n";
  long int i, j, k;
  for(j=0; j < 10; ++j)
    for(i=0; i < 100000; ++i){
        g = double(rand()) / RAND_MAX;
        m = double(rand()) / RAND_MAX;
      for(k=0; k < 1000; ++k){
        if(f != 0)
          m = g;
      }
    }
  clock_t end = clock();
  std::cout << "time: " << (end - start) << "\n";
}

void gnuplot(int numberOfTerms, int n, int m2){  
  Modes1DContainer mcEven(n), mcOdd(n);
  capd::auxil::OutputStream out[DIAMS];
  capd::auxil::OutputStream log(std::cout, false, true);
  log.logfile("log.txt", true); log.log = true;
  capd::auxil::OutputStream plotlog(std::cout, false, true); plotlog.logfile("plotlog.txt", true); plotlog.log = true;
  capd::auxil::OutputStream plotlog2(std::cout, false, true); plotlog2.logfile("plotlog2.txt", true); plotlog2.log = true;
  capd::auxil::OutputStream plotlog3(std::cout, false, true); plotlog3.logfile("plotlog3.txt", true); plotlog3.log = true;
  capd::auxil::OutputStream plotlog4(std::cout, false, true); plotlog4.logfile("plotlog4.txt", true); plotlog4.log = true;
  int i;
  for(i=0; i < DIAMS; ++i){
    out[i].show = false;
    out[i].flush = true;
  }
  i = 3;
  std::stringstream ss;
  ss << numberOfTerms << "_terms" << diams[i] << "_diam.txt";
  std::cout << "\n\nTEST number of terms: " << numberOfTerms << ", diam: " << diams[i] << "\n";
  out[i].logfile(ss.str().c_str(), true); out[i].log = true;

  double diam = diams[i];
  
  mcOdd[Index1D(1)].im = 1 + Interval(-diam / 2, diam / 2);
  mcOdd[Index1D(-1)].im = -1 + Interval(-diam / 2, diam / 2);
  
  mcOdd[Index1D(2)].im = -1 + Interval(-diam / 2, diam / 2);
  mcOdd[Index1D(-2)].im = 1 + Interval(-diam / 2, diam / 2);
  
  mcOdd[Index1D(3)].im = 1 + Interval(-diam / 2, diam / 2);
  mcOdd[Index1D(-3)].im = -1 + Interval(-diam / 2, diam / 2);
  
  mcOdd[Index1D(4)].im = -1 + Interval(-diam / 2, diam / 2);
  mcOdd[Index1D(-4)].im = 1 + Interval(-diam / 2, diam / 2);
  
  out[i] << "u:\n" << mcOdd << "\n";
  int coord;
  Modes1DContainer& mc(mcOdd);
  Modes1DContainer r1(n), r2(n), r3(n), rc(n), rt(n);
  FFT1D fft2(n, m2, Interval::pi());
  
  if(numberOfTerms == 2){
    rc = fft2.calculateConvolution(mc, mc);
  }    
  if(numberOfTerms == 3){  
    rc = calculateThirdDegConvolution(n, mc);
  }    
  if(numberOfTerms == 4){
    rc = calculateFourthDegConvolution(n, mc);
  }
  if(numberOfTerms == 5){
    rc = calculateFifthDegConvolution(n, mc);
  }
  if(numberOfTerms == 7){
    rc = calculateSeventhDegConvolution(n, mc);
  }
  
   FFT1D::DFT1DGrid hatu =  FFT1D::DFT1DGrid(m2);
  FFT1D::DFT1DGrid s =  FFT1D::DFT1DGrid(m2);
  fft2.fastTransform(rc, hatu, capd::jaco::phaseShiftOddEven, capd::jaco::selfSorting);
  plotlog3 << hatu << "\n";
  
  out[i] << "direct convolution:\n" << rc << "\n";
  std::cout << "direct convolution     MAX DIAM: " << rc.maxDiam(coord);
  out[i] << "direct convolution     MAX DIAM: " << rc.maxDiam(coord);
  std::cout << "IS REACHED AT COORD: " << coord << ", AVG DIAM: " << rc.avgDiam() << "\n";
  out[i] << "IS REACHED AT COORD: " << coord << ", AVG DIAM: " << rc.avgDiam() << "\n";
  log << rc.avgDiam() << "\n";
  
  fft2.fastTransform(mc, hatu, capd::jaco::phaseShiftOddEven, capd::jaco::selfSorting);
//  plotlog << hatu << "\n";
//  int j;
//  for(j=0; j < 9; ++j){
//    plotlog2 << 2*3.14159265358979323846264338*j/9 << " " << hatu[j*m2/9].im.leftBound() << "\n";
//  }
  if(numberOfTerms >= 2){
    s.multiply(hatu, hatu);
    if(numberOfTerms >= 3){
      s.multiplyOnly(s, hatu);
      if(numberOfTerms >= 4){
        s.multiplyOnly(s, hatu);
        if(numberOfTerms >= 5)
          s.multiplyOnly(s, hatu);
      }
    }
  }else{
    throw std::runtime_error("number of terms is too small.\n");
  }
  fft2.fastInverseTransform(s, r3, capd::jaco::phaseShiftOddEven, capd::jaco::selfSorting);
  fft2.fastTransform(r3, hatu, capd::jaco::phaseShiftOddEven, capd::jaco::selfSorting);
  plotlog4 << hatu << "\n";
  
  out[i] << "padding:\n" << r3 << "\n";
  std::cout << "padding                MAX DIAM: " << r3.maxDiam(coord);
  out[i] << "padding                MAX DIAM: " << r3.maxDiam(coord);
  std::cout << "IS REACHED AT COORD: " << coord << ", AVG DIAM: " << r3.avgDiam() << "\n";
  out[i] << "IS REACHED AT COORD: " << coord << ", AVG DIAM: " << r3.avgDiam() << "\n";
  std::cout << ", MAX DIFF WITH DIRECT: " << calculateMaxDiffBetweenDiams(r3, rc, coord);
  out[i] << ", MAX DIFF WITH DIRECT: " << calculateMaxDiffBetweenDiams(r3, rc, coord);
  std::cout << "IS REACHED AT COORD: " << coord << "\n";
  out[i] << "IS REACHED AT COORD: " << coord << "\n";
  log << r3.avgDiam() << "\n";
}

int main(int argc, char * argv[]){

  setLoggers();
  fftDebug.log = true;
  generalDebug2.log = true;
  initN();
  initDiams();
  if(argc != 3){
    throw std::runtime_error("Number of input parameters doesn't match.\nCorrect format is the following:\n"
        "FFTTests numberOfTerms diametersIndex");
  }
  test(atoi(argv[1]), atoi(argv[2]));
  
//  gnuplot( atoi(argv[1]), atoi(argv[2]), atoi(argv[3]) );
  
//  singleTest(capd::jaco::unscrambled, atof(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
  
//    if(argc != 3){
//      throw std::runtime_error("Number of input parameters doesn't match.\n");
//    }
//    jetTest(atoi(argv[1]), atof(argv[2]));
}
