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

#include <time.h>
#include <math.h>
#include <stdlib.h>

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
  
#include "capd/vectalg/Matrix.h"
#include "ComplexScalar.h"
#include "Pair.h"
#include "FirstOrderJet.h"
#include "Odd.h"
#include "Index.h"
#include "Coefficients.h"
#include "PolyBd.h"
  
typedef capd::vectalg::Matrix<Interval,_D,_D> IntervalMatrix;
typedef capd::vectalg::Vector<Interval, _D> IntervalVector;
typedef capd::jaco::ComplexScalar<Interval> ComplexScalar;
typedef capd::jaco::ComplexDerivativePair<ComplexScalar> ComplexDerivativePair;

typedef capd::jaco::FirstOrderJet<ComplexDerivativePair, 60> FOJ1D;

int main(){
  setLoggers();
  int n = 20;
      //m = 64,
      //d = 20,
      //order = 6;
  Interval nu(0.127),
           step(0.000128);
  ///begin FOJ initialization for FFT integrator
  FOJ1D::dim = n;
  capd::jaco::DPDEContainer container;
  container.setToRealValuedOdd();
  FOJ1D::initialCondition = &container;
  FOJ1D buffer;
  FOJ1D::initialized = 1;
  FOJ1D::buffer = &buffer;
  ///end FOJ initialization
  FOJ1D foj1, foj2, fojr, fojr2, fojr3;
  //foj1.val = ComplexScalar(0, 1); odd
  //foj2.val = ComplexScalar(0, 1);
  foj1.val = ComplexScalar(0, 0.5); 
  foj2.val = ComplexScalar(0, 0.5);
  fojr2 = FOJ1D(0);
  fojr3 = FOJ1D(0);
  Interval r;
  int i;
  for(i=1; i < n; i++){
    foj1.setVariationalPartToId(i);
    r = Interval(rand())/RAND_MAX;
    foj1.ksi[i] *= ComplexScalar(r, 0);
    foj2.setVariationalPartToId(i);
    foj2.ksi[i] *= ComplexScalar(r, 0);
    fojr2.setVariationalPartToId(i);
    fojr2.ksi[i] *= ComplexScalar(r, 0);
  }
  FOJ1D::optimizationControl = capd::jaco::local;
  fojr = FOJ1D(0);
  fojr += 4 * foj2;
  fojr3 = conjugate(fojr);
  generalDebug << "fojr:\n"<<fojr << "\n";
  generalDebug << "fojr3:\n"<<fojr3<<"\n";
  FOJ1D::optimizationControl = capd::jaco::global;
  FOJ1D::initialCondition->solutionType = capd::jaco::complexValued;
  FOJ1D::initialCondition->solutionType2 = capd::jaco::other;
  FOJ1D::initialCondition->subspaceType = capd::jaco::non;
  fojr = FOJ1D(0);
  generalDebug2 << "foj1:\n" << foj1 << "\n";
  generalDebug2 << "foj2:\n" << foj1 * foj2 << "\n";
  fojr += 4 * foj2;
  fojr3 = conjugate(fojr);
  generalDebug2 << "fojr:\n" << fojr << "\n";
  generalDebug2 << "fojr3:\n"<<fojr3<<"\n";
}
