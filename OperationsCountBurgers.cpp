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
 * OperationsCountBurgers.cpp
 *
 *  Created on: Aug 31, 2011
 *      Author: cyranka
 */
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <cstring>
#include "config.h"
//#include "capd/filib/Interval.h"
#include "intervals/Interval.hpp"

#define DOUBLE long double
#if __FILIB__
  typedef capd::filib::Interval<DOUBLE> Interval;
#else
  typedef capd::intervals::Interval<DOUBLE> Interval;
  #define PI Interval::pi()
//  typedef DOUBLE Interval;
//  #define PI 3.1415926535897932384626433832795
#endif

#include "ComplexScalar.h"
#include "FirstOrderJet.h"
#include "Odd.h"
#include "Index.h"
#include "Coefficients.h"
#include "Pair.h"

#include "Odd.h"
#include "dynset/C0Rect2Set.hpp"
#include "dynset/C0Rect2RSet.hpp"

#include "FFT.h"
#include "Equations.h"
#include "FFTDynSys.h"
#include "PolyBd.h"

typedef capd::vectalg::Matrix<Interval,_D,_D> IntervalMatrix;
typedef capd::vectalg::Vector<Interval, _D> IntervalVector;
typedef capd::vectalg::EuclNorm<IntervalVector, IntervalMatrix> EuclNorm;
typedef capd::jaco::ComplexScalar<Interval> ComplexScalar;
typedef capd::jaco::ComplexDerivativePair<ComplexScalar> ComplexDerivativePair;

//Integrators 1D
//Basic Fad Taylor 1D
typedef capd::jaco::Index1D Index1D;
typedef capd::jaco::MaximumNorm<Index1D> MaximumNorm;
//typedef capd::jaco::Odd<Interval, Index1D, MaximumNorm> OddSubspace;
typedef capd::jaco::RealPolynomialBound<ComplexScalar, Index1D, 0> RealPolynomialBound;
//typedef capd::jaco::RealPolynomialBound<ComplexScalar, Index1D, 0> RealPolynomialBound;
typedef capd::jaco::Burgers<RealPolynomialBound> Burgers;
typedef capd::jaco::KS<RealPolynomialBound> KS;

typedef capd::jaco::MaximumNorm<Index1D> Norm1D;
//typedef capd::jaco::Real<Interval, Index1D, Norm1D> Real1D;
typedef capd::jaco::Odd<Interval, Index1D, Norm1D> Real1D;
typedef capd::dynset::C0Rect2RSet<IntervalMatrix> C0Set;
//FFT1D

typedef capd::jaco::FFT1D<ComplexScalar, ComplexScalar, 0, 0, RealPolynomialBound> FFT1D;

typedef capd::jaco::DPDE2<Burgers, FFT1D, 0> BurgersDPDE; ///set equation here

typedef capd::jaco::FirstOrderJet<ComplexDerivativePair, 0> FOJ1D;

typedef capd::jaco::ComplexPolyBdJetOptimized<Interval, FOJ1D, Index1D, 0> ModesContainer;

typedef capd::jaco::FFT1D<FOJ1D, ComplexScalar, 0, 0, ModesContainer> JetFFT1D;

typedef capd::jaco::DPDE2<Burgers, JetFFT1D, 0> JetBurgersDPDE; ///set equation here

typedef capd::jaco::FFTTaylorDynSys<BurgersDPDE, JetBurgersDPDE, 0> FFTDynSys;

typedef FFTDynSys::ModesContainerType ModesContainer1D;

capd::auxil::OutputStream fftApproachOperationsCount(std::cout, false, true);
capd::auxil::OutputStream directApproachOperationsCount(std::cout, false, true);

int ns[100];
int ms[100];
int orders[100];

/**HERE YOU SET WHICH VALUES TO USE IN TESTS*/
const int starti = 3,     
          endi = 19,
          startio = 0,
          endio = 7;

///this is to calculate the number of operations for the padding technique
void initTablesPadding(){
  ns[2] = 14; ms[2] = 45;  
  ns[3] = 19; ms[3] = 60;  
  ns[4] = 24; ms[4] = 75;  
  ns[5] = 29; ms[5] = 90;  
  ns[6] = 34; ms[6] = 108;  
  ns[7] = 39; ms[7] = 120;  
  ns[8] = 44; ms[8] = 135;  
  ns[9] = 49; ms[9] = 150;    
  ns[10] = 54; ms[10] = 180;
  ns[11] = 59; ms[11] = 180;
  ns[12] = 64; ms[12] = 200;
  ns[13] = 69; ms[13] = 216;
  ns[14] = 74; ms[14] = 225;
  ns[15] = 79; ms[15] = 240;
  ns[16] = 84; ms[16] = 256;
  ns[17] = 89; ms[17] = 270;
  ns[18] = 94; ms[18] = 288;
  ns[19] = 99; ms[19] = 300;
  
  
  orders[0] = 3;
  orders[1] = 4;
  orders[2] = 5;
  orders[3] = 6;
  orders[4] = 7;
  orders[5] = 8;
  orders[6] = 9;
  orders[7] = 10;
  orders[8] = 15;
  orders[9] = 12;
}

void initTables(){
  
  ns[2] = 14;ms[2] = 30;  
  ns[3] = 19;ms[3] = 40;  
  ns[4] = 24;ms[4] = 50;  
  ns[5] = 29;ms[5] = 60;  
  ns[6] = 34;ms[6] = 72;  
  ns[7] = 39;ms[7] = 80;  
  ns[8] = 44;ms[8] = 90;  
  ns[9] = 49;ms[9] = 100;    
  ns[10] = 54; ms[10] = 120;
  ns[11] = 59; ms[11] = 120;
  ns[12] = 64; ms[12] = 144;
  ns[13] = 69; ms[13] = 144;
  ns[14] = 74; ms[14] = 150;
  ns[15] = 79; ms[15] = 160;
  ns[16] = 84; ms[16] = 180;
  ns[17] = 89; ms[17] = 180;
  ns[18] = 94; ms[18] = 192;
  ns[19] = 99; ms[19] = 200;
  ns[20] = 109; ms[20] = 225;
  ns[21] = 119; ms[21] = 240;
  ns[22] = 129; ms[22] = 288;
  ns[23] = 139; ms[23] = 288;
  ns[24] = 149; ms[24] = 300;
/*  ns[25] = 159; ms[25] = 320;
  ns[26] = 169; ms[26] = 360;
  ns[27] = 179; ms[27] = 360;
  ns[28] = 189; ms[28] = 400;
  ns[29] = 199; ms[29] = 400;
  ns[30] = 299; ms[30] = 600;
  ns[31] = 399; ms[31] = 800;
  ns[32] = 499; ms[32] = 1000;
  ns[33] = ; ms[33] = ;
  ns[34] = ; ms[34] = ;
  ns[35] = ; ms[35] = ;
  ns[36] = ; ms[36] = ;
  ns[37] = ; ms[37] = ;
  ns[38] = ; ms[38] = ;
  ns[39] = ; ms[39] = ;*/
  orders[0] = 3;
  orders[1] = 4;
  orders[2] = 5;
  orders[3] = 6;
  orders[4] = 7;
  orders[5] = 8;
  orders[6] = 9;
  orders[7] = 10;
  orders[8] = 15;
  orders[9] = 12;
}


int main(){
  fftApproachOperationsCount.logfile("fftOpCnt_Burgers.txt", true);
  fftApproachOperationsCount.log = true;
  directApproachOperationsCount.logfile("directOpCnt_Burgers.txt", true);
  directApproachOperationsCount.log = true;  
  setLoggers();
  initTablesPadding();
  
  int i, j;  
  for(j = startio; j <= endio; ++j){
    for(i = starti; i <= endi; ++i){    
      int n = ns[i], //change FOJ1D stack dimension
           m = ms[i],
           order = orders[j];
       std::cout << "n=" << n << ", m=" << m << ", order=" << order << "\n";
       Interval nu(0.1), //this was chosen empirically such that set will decrease by 10% at the time 1
                step(0.0001);
    
       ///begin FOJ initialization for FFT integrator
       capd::jaco::DPDEContainer container;
       ///2.set here the subspace of the initial condition e.g. setToRealValuedOdd means that the initial condition is real valued odd
       container.setToRealValued();
       FOJ1D::initialize(ModesContainer1D::modes2arraySizeStatic(n), container);
       ///end FOJ initialization   
    
       FFTDynSys dynsys( n, m, step, order, PI, nu);
    
       Index1D idx;
       RealPolynomialBound u_0(n), enclosure(n);
       
       double diam = 1e-08;
       Interval r = Interval(-diam/2, diam/2);
       
       int i;
       for(i=1; i <= n; ++i){
         u_0.set(Index1D(i), ComplexScalar(r, r));
       } 
       
       dynsys.setJetDynSysAlgorithmType(capd::jaco::FFT);        
       clearSums();
       C0Set set(u_0, 1);             
       for(i=0; i < 1; ++i){    
         set.move(dynsys);
       }    
       fftApproachOperationsCount << order << " " << n << " " << " " << RRmultiplicationsSum / 1e+06 << " " << RRadditionsSum / 1e+06 << "\n";
           
       dynsys.setJetDynSysAlgorithmType(capd::jaco::direct);
       clearSums();
       set = C0Set(u_0, 1);             
       for(i=0; i < 1; ++i){    
         set.move(dynsys);
       }    
       directApproachOperationsCount << order << " " << n << " " << " " << RRmultiplicationsSum / 1e+06 << " " << RRadditionsSum / 1e+06 << "\n";
    } 
    fftApproachOperationsCount << "\n";
    directApproachOperationsCount << "\n";
  }
}
