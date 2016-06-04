/* DBCPTests.cpp
 *
 *  Created on: May, 12, 2016
 *      Author: Jacek Cyranka
 *
 *  This file contains some test cases for DBCP Heteroclinic connection proof
 *
 */

#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <cstring>

#include <sstream>
#include <stdio.h>
#include <string.h>

//#include "capd/filib/Interval.h"
#include "config.h"
//#include "capd/intervals/Interval.hpp"


bool COUNT_OPERATIONS = false;
#include "intervals/Interval.hpp" // for operations count
#define DOUBLE double

#include "alglib/capd2alglib.h"

#if __FILIB__
  typedef capd::filib::Interval<DOUBLE> Interval;
#else
  typedef capd::intervals::Interval<DOUBLE> Interval;
  #define PI Interval::pi()
//  typedef DOUBLE Interval;
//  #define PI 3.1415926535897932384626433832795
#endif


// MOCK FUNCTIONS NEEDED FOR NONRIGOROUS INTEGRATOR TO WORK
void setLeftBound( Interval& i, DOUBLE b){
  i.setLeftBound(b);
}
void setRightBound( Interval& i, DOUBLE b){
  i.setRightBound(b);
}

// MOCK FUNCTIONS NEEDED FOR NONRIGOROUS INTEGRATOR TO WORK
DOUBLE diam( DOUBLE i ){
  return 0.;
}
DOUBLE power( DOUBLE i, DOUBLE e ){
  return pow(i, e);
}
void setLeftBound( DOUBLE& i, DOUBLE b){
  i = b;
}
void setRightBound( DOUBLE& i, DOUBLE b){
  i = b;
}


long long unsigned TOTAL_ITER;

#include "capd/dynsys/BasicFadTaylor.h"
#include "capd/dynsys/FadTaylor.h"

#include "ComplexScalar.h"
#include "FirstOrderJet.h"
#include "Odd.h"
#include "Even.h"
#include "Index.h"
#include "Coefficients.h"
#include "Pair.h"

#include "capd/dynset/C0Rect2Set.hpp"
#include "capd/dynset/C0Rect2RSet.hpp"

#include "FFT.h"
#include "Equations.h"
#include "FFTDynSys.h"
#include "PolyBd.h"

#include "DPDEInclusionCW.h"
#include "InclRect2Set.hpp"

#include "PolyBdInputReader.h"

#include "FixedPoint.h"

#define _D 0

typedef capd::vectalg::Matrix< DOUBLE, _D, _D > DoubleMatrix;
typedef capd::vectalg::Vector< DOUBLE, _D > DoubleVector;

typedef capd::vectalg::Matrix<Interval,_D,_D> IntervalMatrix;
typedef capd::vectalg::Vector<Interval, _D> IntervalVector;
typedef capd::vectalg::Vector<bool, _D> BoolVector;
typedef capd::vectalg::MaxNorm<IntervalVector, IntervalMatrix> MaxNorm;
typedef capd::vectalg::EuclNorm<IntervalVector, IntervalMatrix> EuclNorm;
typedef capd::jaco::ComplexScalar<Interval> ComplexScalar;
typedef capd::jaco::ComplexDerivativePair<ComplexScalar> ComplexDerivativePair;


//Integrators 1D
//Basic Fad Taylor 1D
typedef capd::jaco::Index1D Index1D;
typedef capd::jaco::MaximumNorm<Index1D> MaximumNorm;
//typedef capd::jaco::Odd<Interval, Index1D, MaximumNorm> OddSubspace;
typedef capd::jaco::Even<Interval, Index1D, MaximumNorm, IntervalMatrix> EvenSubspace;
typedef capd::jaco::RealPolynomialBound<ComplexScalar, Index1D, 0, EvenSubspace> RealPolynomialBound;
//typedef capd::jaco::RealPolynomialBound<ComplexScalar, Index1D, 0> RealPolynomialBound;
typedef capd::jaco::Burgers<RealPolynomialBound> Burgers;
typedef capd::jaco::KS<RealPolynomialBound> KS;
typedef capd::jaco::GL<RealPolynomialBound> GL;
typedef capd::jaco::SH<RealPolynomialBound> SH;
typedef capd::jaco::CH<RealPolynomialBound> CH;
typedef capd::jaco::DBCP<RealPolynomialBound> DBCP;


typedef capd::dynset::C0Rect2RSet<IntervalMatrix> C0Set;
//FFT1D

typedef capd::jaco::FFT1D<ComplexScalar, ComplexScalar, 0, 0, RealPolynomialBound> FFT1D;

typedef capd::jaco::DPDE3<DBCP, FFT1D, 0> SHDPDE; ///set equation here

typedef capd::jaco::FirstOrderJet<ComplexDerivativePair, 0> FOJ1D;

typedef capd::jaco::ComplexPolyBdJetOptimized<Interval, FOJ1D, Index1D, 0, EvenSubspace> ModesContainer;

typedef capd::jaco::FFT1D<FOJ1D, ComplexScalar, 0, 0, ModesContainer> JetFFT1D;

typedef capd::jaco::DPDE3<DBCP, JetFFT1D, 0> JetSHDPDE; ///set equation here

typedef capd::jaco::FFTTaylorDynSys<SHDPDE, JetSHDPDE, 0> FFTDynSys;

typedef FFTDynSys::ModesContainerType ModesContainer1D;

typedef capd::jaco::DPDEInclusionCW<SHDPDE, FFTDynSys> DPDEInclusionCW3;
typedef capd::diffIncl2::InclRect2Set<IntervalMatrix, RealPolynomialBound> InclRect2Set;


typedef capd::jaco::Box<IntervalVector, IntervalMatrix, RealPolynomialBound> Box;
typedef capd::jaco::BoxFinder<SHDPDE, JetSHDPDE, DOUBLE, 0> BoxFinder;

//Needed for partial derivatives -- the cone condition computations
typedef capd::jaco::DPDE2<DBCP, FFT1D, 0> PartialDerivativePDE; ///set equation here

//and nonrigorous integrator


typedef capd::jaco::ComplexScalar<DOUBLE> ComplexScalar_;
typedef capd::jaco::Even<DOUBLE, Index1D, MaximumNorm, DoubleMatrix> EvenSubspace_;
typedef capd::jaco::RealPolynomialBound<ComplexScalar_, Index1D, 0, EvenSubspace_> RealPolynomialBound_;
typedef capd::jaco::FFT1D<ComplexScalar_, ComplexScalar_, 0, 0, RealPolynomialBound_> FFT1D_;
typedef capd::jaco::DBCP<RealPolynomialBound_> DBCP_;
typedef capd::jaco::DPDE3<DBCP_, FFT1D_, 0> SHDPDE_;
typedef capd::jaco::ComplexDerivativePair<ComplexScalar_> ComplexDerivativePair_;
typedef capd::jaco::FirstOrderJet<ComplexDerivativePair_, 0> FOJ1D_;
typedef capd::jaco::ComplexPolyBdJetOptimized<DOUBLE, FOJ1D_, Index1D, 0, EvenSubspace_> ModesContainer_;
typedef capd::jaco::FFT1D<FOJ1D_, ComplexScalar_, 0, 0, ModesContainer_> JetFFT1D_;
typedef capd::jaco::DPDE3<DBCP_, JetFFT1D_, 0> JetSHDPDE_;



//loads Polynomial Bounds from the specified file
int loadDataFromFile(char* fileName, int* m, int* M, int* dftPts, int* dftPts2, int* order, double* step, double* nu, double* sigma, double* piOverL ){

  FILE* file;
  int r;

  if(!(file = fopen(fileName, "r"))) {
    std::cerr << "The specified file does not exists.\n";
    throw std::runtime_error("The specified file does not exists.\n");
  }

  double test;
  char string[100000];

  if(fscanf(file, "m=%d\n", m) <= 0) {
    std::cerr << "Input data file format error. Check if it has proper format.\n";
    throw std::runtime_error("Input data file format error. Check if it has proper format.\n");
  }
  std::cout  << "m=" << *m << "\n";

  if(fscanf(file, "M=%d\n", M) <= 0) {
    std::cerr << "Input data file format error. Check if it has proper format.\n";
    throw std::runtime_error("Input data file format error. Check if it has proper format.\n");
  }
  std::cout  << "M=" << *M << "\n";

  if(fscanf(file, "dftPts=%d\n", dftPts) <= 0) {
    std::cerr << "Input data file format error. Check if it has proper format.\n";
    throw std::runtime_error("Input data file format error. Check if it has proper format.\n");
  }
  std::cout  << "dftPts=" << *dftPts << "\n";

  if(fscanf(file, "dftPts2=%d\n", dftPts2) <= 0) {
    std::cerr << "Input data file format error. Check if it has proper format.\n";
    throw std::runtime_error("Input data file format error. Check if it has proper format.\n");
  }
  std::cout  << "dftPts2=" << *dftPts2 << "\n";

  if(fscanf(file, "order=%d\n", order) <= 0) {
    std::cerr << "Input data file format error. Check if it has proper format.\n";
    throw std::runtime_error("Input data file format error. Check if it has proper format.\n");
  }
  std::cout  << "order=" << *order << "\n";


  if(fscanf(file, "step=%le\n", &test) <= 0) {
    std::cerr << "Input data file format error. Check if it has proper format.\n";
    throw std::runtime_error("Input data file format error. Check if it has proper format.\n");
  }
  *step = test;
  std::cout  << "step=" << *step << "\n";

  if(fscanf(file, "nu=%le\n", &test) <= 0) {
    std::cerr << "Input data file format error. Check if it has proper format.\n";
    throw std::runtime_error("Input data file format error. Check if it has proper format.\n");
  }
  *nu = test;
  std::cout  << "nu=" << *nu << "\n";

  if(fscanf(file, "sigma=%le\n", &test) <= 0) {
    std::cerr << "Input data file format error. Check if it has proper format.\n";
    throw std::runtime_error("Input data file format error. Check if it has proper format.\n");
  }
  *sigma = test;
  std::cout  << "sigma=" << *sigma << "\n";

  if(fscanf(file, "piOverL=%le\n", &test) <= 0) {
    std::cerr << "Input data file format error. Check if it has proper format.\n";
    throw std::runtime_error("Input data file format error. Check if it has proper format.\n");
  }
  *piOverL = test;
  std::cout  << "piOverL=" << *piOverL << "\n";

  fclose(file);
  //end of reading data
  return r;
}


int loadApproximateFixedPoint(char* fileName, IntervalVector& fixedPoint){

  FILE* file;

  if(!(file = fopen(fileName, "r"))) {
    std::cerr << "The specified file does not exists.\n";
    throw std::runtime_error("The specified file does not exists.\n");
  }
  DOUBLE nr;
  int i = 0;
  while(fscanf(file, "%le\n", &nr ) > 0 && i < fixedPoint.size() ){

    fixedPoint[i++] = Interval(nr, nr);

  }

  std::cout << "fixedPoint=\n" << fixedPoint << "\n";

  fclose(file);
}



std::ostream& modes(std::ostream& out, const RealPolynomialBound& pb){
  out << "m=" << pb.n << "\n";
  out << "M=" << pb.N << "\n";
  out << (capd::jaco::DPDEContainer&)pb ;
  for(Index1D index = pb.firstModeIndex(pb.irProjection); !index.limitReached(pb.irProjection); index.inc(pb.irProjection)){
    out << "(" <<  pb[index].re << ", " << pb[index].im << ")\n";
  }
  out << "far_tail=\n";
  out << pb.farTail.getC() << "\n";
  out << pb.farTail.getS() << "\n";
  return out;
}



int main(int argc, char * argv[]){
  setLoggers();

  int m, /*change FOJ1D stack dimension*/ M, dftPts, dftPts2, order;
  double step_, nu_, sigma_, piOverL_;
  Interval nu, step;

  std::ofstream currentSet;
  loadDataFromFile("config.in", &m, &M, &dftPts, &dftPts2, &order, &step_, &nu_, &sigma_, &piOverL_ );


  step = step_;
  nu = nu_;

  IntervalVector zeroCenteredBox(m + 1), eigenvaluesRe(m + 1), fixedPoint(m + 1);
  DoubleVector dfPoint( m + 1 );

  int i,j;
  for(i = 0; i < m + 1; i++)
    dfPoint[i] = rightBound( fixedPoint[i] );

  ///begin FOJ initialization for FFT integrator
  capd::jaco::DPDEContainer container;
  ///2.set here the subspace of the initial condition e.g. setToRealValuedOdd means that the initial condition is real valued odd
  container.setToRealValuedEven();
  FOJ1D::initialize(ModesContainer1D::modes2arraySizeStatic(m), container);
  FOJ1D_::initialize(ModesContainer_::modes2arraySizeStatic(m), container);
  ///end FOJ initialization

  RealPolynomialBound tail(m, M, container);

  JetSHDPDE jetDPDE(m, dftPts, nu, PI, order);
  SHDPDE dpde(m, M, dftPts, dftPts2, nu, PI, order);
  PartialDerivativePDE pdpde(m, M, dftPts, dftPts2, nu_, PI, order);

  jetDPDE.m_r = 0;
  jetDPDE.m_N_coeff = 1.;
  dpde.m_r = 0;
  dpde.m_N_coeff = 1.;
  pdpde.m_r = 0;
  pdpde.m_N_coeff = 1.;

  BoxFinder finder(dpde, jetDPDE);

  Box Qbox(m + 1);

  IntervalVector V = Qbox;

  //PARTIAL DERVIATIVE TESTS IN FINITE DIMENSIONS

  //TEST 1
  std::cout << "PARTIAL DERIVATIVE TEST 1\n";

  V[0] = 1.;
  V[1] = 2.;
  V[2] = 3.;
  V[3] = 4.;
  V[4] = 5.;
  V[5] = 6.;
  V[6] = 7.;
  V[7] = 8.;
  V[8] = 9.;
  V[9] = 10.;
  V[10] = -1.;
  V[11] = -2.;
  V[12] = -3.;
  V[13] = -4.;
  V[14] = -5.;
  V[15] = -6.;
  V[16] = -7.;
  V[17] = -8.;
  V[18] = -9.;
  V[19] = -10.;

  tail = V;
  std::cout << "V:\n" << V << "\n";
  std::cout << "tail:\n" << tail << "\n";
  IntervalMatrix jacobian = jetDPDE.jacobian( V );
  IntervalMatrix jacobianFFT = jetDPDE.jacobianFFT( V );
  std::cout << "jacobian:\n" << jacobian << "\n";
  std::cout << "jacobianFFT:\n" << jacobianFFT << "\n";
  RealPolynomialBound r( tail ), rFFT( tail ), rr(tail), tmcN(tail);

  //cleans the tail part of 'tail'
  tail.cleanTail();
  pdpde.CalculateNonlinearTermDirectly( tail, r , capd::jaco::full, false);
  IntervalMatrix expJacobian(Qbox.m_n, Qbox.m_n);
  //std::cout << "r:\n" << r << "\n";
  for(int i = 0; i < Qbox.m_n; i++){
    for(int j = 0; j < Qbox.m_n; j++){
      expJacobian[i][j] = 3 * (  r[ Index1D( abs(i-j) ) ]  +  r[ Index1D( abs(i+j) ) ] ).re ;
    }
    expJacobian[i][0] = 0.;
    expJacobian[i][i] += dpde.lambda_k( i );
  }
  expJacobian[0][0] = 1.;
  std::cout << "expJacobian:\n" << expJacobian << "\n";
  //COMPARISION OF PARTIAL DERIVATIVES MATRICES OBTAINED USING TWO DIFFERENT METHODS
  for(int i = 0; i < Qbox.m_n; i++){
    for(int j = 0; j < Qbox.m_n; j++){
      if( ! expJacobian[i][j].subset(jacobian[i][j]) ){
        std::cout << "PARTIAL DERIVATIVE TEST 1 is not satisfied (in jacobian matrix).\n " << expJacobian[i][j] << " is not a subset of " << jacobian[i][j] << "\n";
        exit(1);
      }

      if( ! expJacobian[i][j].subset(jacobianFFT[i][j]) ){
        std::cout << "PARTIAL DERIVATIVE TEST 1 is not satisfied (in jacobianFFT matrix).\n " << expJacobian[i][j] << " is not a subset of " << jacobianFFT[i][j] << "\n";
        exit(1);
      }
    }
  }

  //TEST 2
  std::cout << "PARTIAL DERIVATIVE TEST 2\n";

  V[0] = 1.;
  V[1] = 2.;
  V[2] = 1.;
  V[3] = 2.;
  V[4] = 1.;
  V[5] = 2.;
  V[6] = 1.;
  V[7] = 2.;
  V[8] = 1.;
  V[9] = 2.;
  V[10] = -1.;
  V[11] = -2.;
  V[12] = -1.;
  V[13] = -2.;
  V[14] = -1.;
  V[15] = -2.;
  V[16] = -1.;
  V[17] = -2.;
  V[18] = -1.;
  V[19] = -2.;

  double rad = 1e-6;
  V += Interval(-rad, rad);

  tail = V;
  std::cout << "V:\n" << V << "\n";
  std::cout << "tail:\n" << tail << "\n";
  jacobian = jetDPDE.jacobian( V );
  jacobianFFT = jetDPDE.jacobianFFT( V );
  std::cout << "jacobian:\n" << jacobian << "\n";
  std::cout << "jacobianFFT:\n" << jacobianFFT << "\n";

  r = tail ;
  //cleans the tail part of 'tail'
  tail.cleanTail();
  pdpde.CalculateNonlinearTermDirectly( tail, r , capd::jaco::full, false);

  //std::cout << "r:\n" << r << "\n";
  for(int i = 0; i < Qbox.m_n; i++){
    for(int j = 0; j < Qbox.m_n; j++){
      expJacobian[i][j] = 3 * (   r[ Index1D( abs(i-j) ) ]  +  r[ Index1D( abs(i+j) ) ] ).re ;
    }
    expJacobian[i][0] = 0.;
    expJacobian[i][i] += dpde.lambda_k( i );
  }
  expJacobian[0][0] = 1.;

  std::cout << "expJacobian:\n" << expJacobian << "\n";
  //COMPARISION OF PARTIAL DERIVATIVES MATRICES OBTAINED USING TWO DIFFERENT METHODS
  for(int i = 0; i < Qbox.m_n; i++){
    for(int j = 0; j < Qbox.m_n; j++){
      if( ! jacobian[i][j].subset(expJacobian[i][j] + 10 * Interval(-rad, rad)) ){
        std::cout << "PARTIAL DERIVATIVE TEST 2 is not satisfied (expJacobian not subset of jacobain matrix).\n " << jacobian[i][j] << " is not a subset of " << expJacobian[i][j] << "\n";
        exit(1);
      }

      if( ! expJacobian[i][j].subset(jacobianFFT[i][j]) ){
        std::cout << "PARTIAL DERIVATIVE TEST 2 is not satisfied (jacobianFFT matrix).\n " << expJacobian[i][j] << " is not a subset of " << jacobianFFT[i][j] << "\n";
        exit(1);
      }
    }
  }




  //VECTOR FIELD TEST IN FINITE DIMENSIONS
  std::cout << "VECTOR FIELD TEST 1\n";

  V[0] = 1.;
  V[1] = 2.;
  V[2] = 3.;
  V[3] = 4.;
  V[4] = 5.;
  V[5] = 6.;
  V[6] = 7.;
  V[7] = 8.;
  V[8] = 9.;
  V[9] = 10.;
  V[10] = -1.;
  V[11] = -2.;
  V[12] = -3.;
  V[13] = -4.;
  V[14] = -5.;
  V[15] = -6.;
  V[16] = -7.;
  V[17] = -8.;
  V[18] = -9.;
  V[19] = -10.;
  V[20] = 1.;
  V[21] = 2.;
  V[22] = 1.;
  V[23] = 2.;
  V[24] = 1.;
  V[25] = 2.;
  V[26] = 1.;
  V[27] = 2.;
  V[28] = 1.;
  V[29] = 2.;
  V[30] = 0.0001;
  V[31] = 0.0002;
  V[32] = 0.0001;
  V[33] = 0.0002;
  V[34] = 0.0001;
  V[35] = 0.0002;
  V[36] = 0.0001;
  V[37] = 0.0002;
  V[38] = 0.0001;
  V[39] = 0.0002;

  rad = 1e-10;
  V += Interval(-rad, rad);

  tail = V;
  tail.infiniteDimensional = false;

  dpde(tail, r);

  std::cout << "r:\n" << r << "\n";

  dpde.useFFT = true;


  dpde(tail, rFFT);
  std::cout << "rFFT:\n" << rFFT << "\n";

  for(int i = 0; i <= m; i++){
    if(! r[i].subset( rFFT[i] )){
      std::cout << "VECTOR FIELD TEST IN FINITE DIMENSION NOT SATISFIED.\n On " << i << " coordinate, " << r[i] << " is not a subset of " << rFFT[i] << "\n";
      exit(1);
    }
  }


  //VECTOR FIELD TESTS IN INFINITE DIMENSIONS

  std::cout << "VECTOR FIELD INFINITE DIMENSION TEST 1\n";

  tail[0] = 1.;
  tail[1] = 2.;
  tail[2] = 3.;
  tail[3] = 4.;
  tail[4] = 5.;
  tail[5] = 6.;
  tail[6] = 7.;
  tail[7] = 8.;
  tail[8] = 9.;
  tail[9] = 10.;
  tail[10] = -1.;
  tail[11] = -2.;
  tail[12] = -3.;
  tail[13] = -4.;
  tail[14] = -5.;
  tail[15] = -6.;
  tail[16] = -7.;
  tail[17] = -8.;
  tail[18] = -9.;
  tail[19] = -10.;
  tail[20] = 1.;
  tail[21] = 2.;
  tail[22] = 1.;
  tail[23] = 2.;
  tail[24] = 1.;
  tail[25] = 2.;
  tail[26] = 1.;
  tail[27] = 2.;
  tail[28] = 1.;
  tail[29] = 2.;
  tail[30] = 0.0001;
  tail[31] = 0.0002;
  tail[32] = 0.0001;
  tail[33] = 0.0002;
  tail[34] = 0.0001;
  tail[35] = 0.0002;
  tail[36] = 0.0001;
  tail[37] = 0.0002;
  tail[38] = 0.0001;
  tail[39] = 0.0002;
  tail[40] = -0.0001;
  tail[41] = -0.0002;
  tail[42] = -0.0001;
  tail[43] = -0.0002;
  tail[44] = -0.0001;
  tail[45] = -0.0002;
  tail[46] = -0.0001;
  tail[47] = -0.0002;
  tail[48] = -0.0001;
  tail[49] = -0.0002;
  tail[50] = 0.000001;
  tail[51] = 0.000002;
  tail[52] = 0.000001;
  tail[53] = 0.000002;
  tail[54] = 0.000001;
  tail[55] = 0.000002;
  tail[56] = 0.000001;
  tail[57] = 0.000002;
  tail[58] = 0.000001;
  tail[59] = 0.000002;
  tail[60] = -0.000001;
  tail[61] = -0.000002;
  tail[62] = -0.000001;
  tail[63] = -0.000002;
  tail[64] = -0.000001;
  tail[65] = -0.000002;
  tail[66] = -0.000001;
  tail[67] = -0.000002;
  tail[68] = -0.000001;
  tail[69] = -0.000002;
  tail[70] = 0.00000001;
  tail[71] = 0.00000002;
  tail[72] = 0.00000001;
  tail[73] = 0.00000002;
  tail[74] = 0.00000001;
  tail[75] = 0.00000002;
  tail[76] = 0.00000001;
  tail[77] = 0.00000002;
  tail[78] = 0.00000001;
  tail[79] = 0.00000002;
  tail[80] = -0.00000001;
  tail[81] = -0.00000002;
  tail[82] = -0.00000001;
  tail[83] = -0.00000002;
  tail[84] = -0.00000001;
  tail[85] = -0.00000002;
  tail[86] = -0.00000001;
  tail[87] = -0.00000002;
  tail[88] = -0.00000001;
  tail[89] = -0.00000002;
  tail[90] = 0.00000001;
  tail[91] = 0.00000002;
  tail[92] = 0.00000001;
  tail[93] = 0.00000002;
  tail[94] = 0.00000001;
  tail[95] = 0.00000002;
  tail[96] = 0.00000001;
  tail[97] = 0.00000002;
  tail[98] = 0.00000001;
  tail[99] = 0.00000002;
  tail[100] = -0.000000001;
  tail[101] = -0.000000002;
  tail[102] = -0.000000001;
  tail[103] = -0.000000002;
  tail[104] = -0.000000001;
  tail[105] = -0.000000002;
  tail[106] = -0.000000001;
  tail[107] = -0.000000002;
  tail[108] = -0.000000001;
  tail[109] = -0.000000002;

  rad = 1e-10;
  for(int i = 0; i < 110; i++)
    tail[i] += Interval(-rad, rad);

  tail.infiniteDimensional = true;
  r.infiniteDimensional = true;
  setS(tail, 9 );

  dpde.useFFT = false;
  dpde(tail, r);
  const RealPolynomialBound& rconst= r;


  for(int k = M + 1; k < 2 * M; k++){
    if( !r.redundantMode( Index1D(k) ).re.subset( power(k, s(r)) * rconst[Index1D(k)].re ) ){
      std::cout << "VECTOR FIELD TEST 1 IN INFINITE DIMENSIONS is not satisfied for k=" << k << " and re value\n" ;
      std::cout << r.redundantMode( Index1D(k) ).re << "is not subset of " << power(k, s(r)) * rconst[Index1D(k)].re  << "\n";
      exit(1);
    }
    if( !r.redundantMode( Index1D(k) ).im.subset( power(k, s(r)) * rconst[Index1D(k)].im ) ){
      std::cout << "VECTOR FIELD TEST 1 IN INFINITE DIMENSIONS is not satisfied for k=" << k << " and im value\n" ;
      std::cout << "" << r.redundantMode( Index1D(k) ).im << " " << power(k, s(r)) * rconst[Index1D(k)].im << "\n";
      exit(1);
    }
  }


  std::cout << "VECTOR FIELD INFINITE DIMENSION TEST 2\n";

  tail[0] = 1.;
  tail[1] = 2.;
  tail[2] = 3.;
  tail[3] = 4.;
  tail[4] = 5.;
  tail[5] = 6.;
  tail[6] = 7.;
  tail[7] = 8.;
  tail[8] = 9.;
  tail[9] = 10.;
  tail[10] = -1.;
  tail[11] = -2.;
  tail[12] = -3.;
  tail[13] = -4.;
  tail[14] = -5.;
  tail[15] = -6.;
  tail[16] = -7.;
  tail[17] = -8.;
  tail[18] = -9.;
  tail[19] = -10.;
  tail[20] = 1.;
  tail[21] = 2.;
  tail[22] = 1.;
  tail[23] = 2.;
  tail[24] = 1.;
  tail[25] = 2.;
  tail[26] = 1.;
  tail[27] = 2.;
  tail[28] = 1.;
  tail[29] = 2.;
  tail[30] = 0.0001;
  tail[31] = 0.0002;
  tail[32] = 0.0001;
  tail[33] = 0.0002;
  tail[34] = 0.0001;
  tail[35] = 0.0002;
  tail[36] = 0.0001;
  tail[37] = 0.0002;
  tail[38] = 0.0001;
  tail[39] = 0.0002;
  tail[40] = -0.0001;
  tail[41] = -0.0002;
  tail[42] = -0.0001;
  tail[43] = -0.0002;
  tail[44] = -0.0001;
  tail[45] = -0.0002;
  tail[46] = -0.0001;
  tail[47] = -0.0002;
  tail[48] = -0.0001;
  tail[49] = -0.0002;
  tail[50] = 0.000001;
  tail[51] = 0.000002;
  tail[52] = 0.000001;
  tail[53] = 0.000002;
  tail[54] = 0.000001;
  tail[55] = 0.000002;
  tail[56] = 0.000001;
  tail[57] = 0.000002;
  tail[58] = 0.000001;
  tail[59] = 0.000002;
  tail[60] = -0.000001;
  tail[61] = -0.000002;
  tail[62] = -0.000001;
  tail[63] = -0.000002;
  tail[64] = -0.000001;
  tail[65] = -0.000002;
  tail[66] = -0.000001;
  tail[67] = -0.000002;
  tail[68] = -0.000001;
  tail[69] = -0.000002;
  tail[70] = 0.00000001;
  tail[71] = 0.00000002;
  tail[72] = 0.00000001;
  tail[73] = 0.00000002;
  tail[74] = 0.00000001;
  tail[75] = 0.00000002;
  tail[76] = 0.00000001;
  tail[77] = 0.00000002;
  tail[78] = 0.00000001;
  tail[79] = 0.00000002;
  tail[80] = -0.000001;
  tail[81] = -0.000002;
  tail[82] = -0.000001;
  tail[83] = -0.000002;
  tail[84] = -0.000001;
  tail[85] = -0.000002;
  tail[86] = -0.000001;
  tail[87] = -0.000002;
  tail[88] = -0.000001;
  tail[89] = -0.000002;
  tail[90] = 0.0001;
  tail[91] = 0.0002;
  tail[92] = 0.0001;
  tail[93] = 0.0002;
  tail[94] = 0.0001;
  tail[95] = 0.0002;
  tail[96] = 0.0001;
  tail[97] = 0.0002;
  tail[98] = 0.0001;
  tail[99] = 0.0002;
  tail[100] = -1;
  tail[101] = -2;
  tail[102] = -1;
  tail[103] = -2;
  tail[104] = -1;
  tail[105] = -2;
  tail[106] = -1;
  tail[107] = -2;
  tail[108] = -1;
  tail[109] = -2;

  rad = 1e-10;
  for(int i = 0; i < 110; i++)
    tail[i] += Interval(-rad, rad);

  tail.infiniteDimensional = true;
  r.infiniteDimensional = true;
  setS(tail, 6 );

  dpde.useFFT = false;
  dpde(tail, r);


  std::cout << "r:\n" << r << "\n";
  for(int k = M + 1; k < 2 * M; k++){

    if( !r.redundantMode( Index1D(k) ).re.subset( power(k, s(r)) * ((const RealPolynomialBound&)r)[Index1D(k)].re ) ){
      std::cout << "VECTOR FIELD TEST 1 IN INFINITE DIMENSIONS is not satisfied for k=" << k << " and re value\n" ;
      std::cout << r.redundantMode( Index1D(k) ).re << "is not subset of " << power(k, s(r)) * rconst[Index1D(k)].re  << "\n";
      exit(1);
    }
    if( !r.redundantMode( Index1D(k) ).im.subset( power(k, s(r)) * ((const RealPolynomialBound&)r)[Index1D(k)].im ) ){
      std::cout << "VECTOR FIELD TEST 1 IN INFINITE DIMENSIONS is not satisfied for k=" << k << " and im value\n" ;
      std::cout << "" << r.redundantMode( Index1D(k) ).im << " " << power(k, s(r)) * rconst[Index1D(k)].im << "\n";
      exit(1);
    }
  }


  //ITEGRATOR TEST (CONVERGENCE TO A FIXED POINT)

  //loading approximate fixed point position
  loadApproximateFixedPoint( "fixedPoint.in", fixedPoint );

  tail.cleanFinitePart();
  tail.cleanTail();

  tail = fixedPoint;

  std::cout << "ITEGRATOR TEST (CONVERGENCE TO A FIXED POINT)\n Integrating the set: " << tail << "\n";

  //setting back the correct parameters


  const int STEPS = 100;
  DPDEInclusionCW3 diffIncl(m, dftPts, M, dftPts2, PI, nu, order, step, MaxNorm());
  diffIncl.getDynamicalSystem().setJetDynSysAlgorithmType(capd::jaco::direct);

  InclRect2Set set( tail );

  tail.infiniteDimensional = true;
  r.infiniteDimensional = true;
  r.cleanFinitePart();
  r.cleanTail();

  rr.infiniteDimensional = true;
  rr.cleanFinitePart();
  rr.cleanTail();

  dpde.m_r = 2;
  dpde.m_N_coeff = nu;
  dpde(tail, r);

  for(i=0; i < STEPS; ++i){
    set.move(diffIncl);
    std::cout << (IntervalVector)set << "\n";
    if(i % 100 == 0){
      std::cout << set.getPerturbationParams() << "\n";
    }
  }


  dpde(set.getPerturbationParams(), rr);

  std::cout << "RESULT OF INTEGRATION TEST:\n rr.sumOfNorms()=" << rr.sumOfNorms() << ", r.sumOfNorms()=" << r.sumOfNorms() << "\n";
  if( rr.sumOfNorms().rightBound() > r.sumOfNorms().rightBound() ){
    std::cout << "ITEGRATOR TEST (CONVERGENCE TO A FIXED POINT) is not satisfied, the euclidean norm of F(later set) is larger than that of F(smaller set)\n";
    exit(1);
  }


  //THE CONE CONDITION TEST


  FOJ1D::destroy();
}
