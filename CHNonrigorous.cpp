/*
 * CHNonrigorous.cpp
 *
 *  Created on: Sep 26, 2013
 *      Author: cyranka
 */


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

//FFT1D
typedef capd::jaco::FFT1D<ComplexScalar, ComplexScalar, 0, 0, RealPolynomialBound> FFT1D;
typedef capd::jaco::DPDE3<DBCP, FFT1D, 0> SHDPDE; ///set equation here
typedef capd::jaco::FFTBasicDynSys<SHDPDE, 0> FFTDynSys;

//JetVersion (needed to calculate jacobian)
typedef capd::jaco::ComplexDerivativePair<ComplexScalar> ComplexDerivativePair;
typedef capd::jaco::FirstOrderJet<ComplexDerivativePair, 0> FOJ1D;
typedef capd::jaco::ComplexPolyBdJetOptimized<Interval, FOJ1D, Index1D, 0, EvenSubspace> ModesContainer;
typedef capd::jaco::FFT1D<FOJ1D, ComplexScalar, 0, 0, ModesContainer> JetFFT1D;
typedef capd::jaco::DPDE3<DBCP, JetFFT1D, 0> JetDBCP; ///set equation here


template< class MatrixT, class DoubleMatrixT >
class Diagonalizator {
public:
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef DoubleMatrixT DoubleMatrixType;
  typedef typename DoubleMatrixType::RowVectorType DoubleVectorType;
  typedef typename VectorType::ScalarType ScalarType;

  static void schurDecompose(const DoubleMatrixType& A, DoubleMatrixType& T, DoubleMatrixType& S){
     if(A.numberOfRows() != A.numberOfColumns())
        throw std::invalid_argument("schurDecompose works only for square matrices");
      ap::real_2d_array a;
      int n = A.numberOfRows();
      a.setbounds(0,n-1, 0, n-1) ;
      int i,j,k;
      for(i =0; i < n; i++)
        for(j=0; j< n; j++){
          a(i,j) = A[i][j];
        }
      ap::real_2d_array t;
      ap::real_2d_array s;
      ap::real_1d_array tau;
      ap::real_2d_array h;
      ap::real_2d_array h1;
      ap::real_2d_array s1;
      ap::real_2d_array s2;
      ap::real_2d_array s3;
      h1.setbounds(1, n, 1, n);
      s1.setbounds(1, n, 1, n);
      s2.setbounds(1, n, 1, n);
      s3.setbounds(1, n, 1, n);

      rmatrixhessenberg(a, n, tau);

      rmatrixhessenbergunpackq(a, n, tau, s);

      rmatrixhessenbergunpackh(a, n, h);

      //changing indexing of h and s from [0...n-1] to [1...n]
      for(i=0; i<n; i++)
        for(j=0; j<n; j++){
          h1(i+1, j+1)=h(i,j);
          s1(i+1, j+1)=s(i,j);
        }

      if(! upperhessenbergschurdecomposition(h1, n, s2))
          throw std::runtime_error("algorithm for schur decomposition did not converge!");

      ///multiplying orthogonal s1 by orthogonal s2
      for(i=1; i<=n; i++)
        for(j=1; j<=n; j++){
          s3(i, j)=0;
          for(k=1; k<=n; k++){
            s3(i, j)+=s1(i, k)*s2(k, j);
          }
        }

      for(int i =0; i < n; i++){
        for(int j=0; j< n; j++){
          T[j][i]=h1(j+1,i+1);
          S[j][i]=s3(j+1,i+1);
        }
      }
      return true;
  }

  ///calculates schurMatrix of 'der' matrix
  static bool schurMatrix(DoubleMatrixType& der, DoubleMatrixType& T, DoubleMatrixType& S, std::ostream& out){
    bool r=schurDecompose(der, T, S);
    return r;
  }

  static void translateComplexChOfCoordIntoRealOne(const DoubleMatrixType& rVec, const DoubleMatrixType& iVec){
    //we scan imaginary matrix iVec looking for blocks 2x2 corresponding to i*w and -i*w part from eigenvectors v+i*w, v-i*w
    int i, j;
    int m_dim = rVec.numberOfRows();
    DoubleMatrixType r(m_dim, m_dim);
    for(i = 0; i < m_dim - 1; i += 2)
      for(j = 0; j < m_dim; j++) {
        if(iVec[i][j] != 0 || iVec[i + 1][j] != 0) { //block with eigenvectors in form v+iw, v-iw is detected
          //instead of block [v+iw,v-iw] we take block [v,w]
          r[i][j] = rVec[i][j]; //v_1
          r[i + 1][j] = rVec[i + 1][j]; //v_2
          r[i][j + 1] = iVec[i][j]; //w_1
          r[i + 1][j + 1] = iVec[i + 1][j]; //w_2
          j++;
        } else { //we rewrite real eigenvectors
          r[i][j] = rVec[i][j];
          r[i + 1][j] = rVec[i + 1][j];
        }
      }
    //in case of odd dimension we check last row
    if(m_dim % 2 == 1){
      for(j = 0; j < m_dim; j++) {
        if(iVec[i][j] != 0){
          r[i][j + 1] = iVec[i][j];
          r[i][j] = rVec[i][j];
          j++;
        }else{
          r[i][j]=rVec[i][j];
        }
      }
    }
    return r;
  }

  static void krawczykRefineMatrix(DoubleMatrixType& A, const DoubleMatrixType& approxAinv, MatrixType& Ainv){
    int i;
    VectorType x;
    int dim=Ainv.numberOfColumns();
    MatrixType AinvT=capd::vectalg::transpose(Ainv),
               candidate;
    for(i=0; i<dim; i++){
      x=AinvT[i];
      if(!krawczykRefineColumn( A, approxAinv, x, i ))
        return false;
      AinvT[i]=x;
    }
    candidate=capd::vectalg::transpose(AinvT);
    intersection(Ainv, candidate, Ainv);
    return true;
  }

  static void changeOfCoordinates(const MatrixType& jacobian, MatrixType& Q, MatrixType& Qinv, DoubleVectorType& rEigen, DoubleVectorType& iEigen, std::ostream& out){
    int d = jacobian.numberOfRows();
    DoubleMatrixType dT(d, d), dS(d, d), dSinv(d, d), rVec(d, d), iVec(d, d), rVecInv(d, d), dJac(d, d);
    rEigen = DoubleVectorType(d),
    iEigen = DoubleVectorType(d);

    MatrixType eigenvectors, eigenvectorsInv, S, Sinv;

    Q = MatrixType(d, d);
    Qinv = MatrixType(d, d);

    for(int i = 0; i < d; i++)
      for(int j = 0; j < d; j++){
        dJac[i][j] = jacobian[i][j].mid().rightBound();
      }

    schurMatrix(dJac, dT, dS, out);
    dS=capd::vectalg::transpose(dS);
    capd::alglib::computeEigenvaluesAndEigenvectors(dT, rEigen, iEigen, rVec, iVec);
    //we have a complex change of coordinates in two matrices rVec the real part and iVec the imaginary part
    //we unite them into one matrix which generates block diagonal form of T
    rVec = translateComplexChOfCoordIntoRealOne(rVec, -iVec);
    rVec = capd::vectalg::transpose(rVec);

    eigenvectors = MatrixType(rVec);
    rVecInv = capd::matrixAlgorithms::inverseMatrix(rVec);
    eigenvectorsInv = MatrixType(rVecInv);

    S = MatrixType(dS);
    //matrices S together with eigenvectors is our change of coordinates, inverses are calculated rigorously
    Q = eigenvectors * S;
    dSinv = capd::matrixAlgorithms::inverseMatrix(dS);
    Sinv = MatrixType(dSinv);
    std::cout<<"Refining inverse matrices using the Krawczyk operator.\n";
    //inflating matrix Sinv
    Sinv += ScalarType(-1, 1);
    //refining it using the Krawczyk operator
    if(!krawczykRefineMatrix(dS, dSinv, Sinv)){
      out << "Procedure of rigorously inverting a matrix using the Krawczyk operator failed.\nTry to rerun program with different"
          << " projection size m.\n";
      std::cerr << "Procedure of rigorously inverting a matrix using the Krawczyk operator failed.\nTry to rerun program with different"
                << " projection size m.\n";
      throw std::runtime_error("Procedure of rigorously inverting a matrix using the Krawczyk operator failed.\nTry to rerun program with different projection size m.\n");
    }

    //the same with rVec matrix
    eigenvectorsInv+=ScalarType(-1, 1);

    if(!krawczykRefineMatrix(rVec, rVecInv, eigenvectorsInv)){
      out << "Procedure of rigorously inverting a matrix using the Krawczyk operator failed.\nTry to rerun program with different"
          << " projection size m.\n";
      std::cerr << "Procedure of rigorously inverting a matrix using the Krawczyk operator failed.\nTry to rerun program with different"
                << " projection size m.\n";
      throw std::runtime_error("Procedure of rigorously inverting a matrix using the Krawczyk operator failed.\nTry to rerun program with different projection size m.\n");
    }
    Qinv = Sinv * eigenvectorsInv;

    //out << "Q=" << Q << "\nQinv=" << Qinv << "\n";
    //MatrixType diagonalization = Q * jacobian * Qinv;
    //out << "diagonalization=" << diagonalization << "\n";
  }

};


///partial specialization designed for nonrigorous computations
template<class DoubleMatrixT>
class Diagonalizator<DoubleMatrixT, DoubleMatrixT> {
public:
  typedef DoubleMatrixT DoubleMatrixType;
  typedef typename DoubleMatrixType::RowVectorType DoubleVectorType;

  static bool schurDecompose(const DoubleMatrixType& A, DoubleMatrixType& T, DoubleMatrixType& S){

     if(A.numberOfRows() != A.numberOfColumns())
        throw std::invalid_argument("schurDecompose works only for square matrices");
      ap::real_2d_array a;
      int n = A.numberOfRows();
      a.setbounds(0,n-1, 0, n-1) ;
      int i,j,k;
      for(i =0; i < n; i++)
        for(j=0; j< n; j++){
          a(i,j) = A[i][j];
        }
      ap::real_2d_array t;
      ap::real_2d_array s;
      ap::real_1d_array tau;
      ap::real_2d_array h;
      ap::real_2d_array h1;
      ap::real_2d_array s1;
      ap::real_2d_array s2;
      ap::real_2d_array s3;
      h1.setbounds(1, n, 1, n);
      s1.setbounds(1, n, 1, n);
      s2.setbounds(1, n, 1, n);
      s3.setbounds(1, n, 1, n);

      rmatrixhessenberg(a, n, tau);

      rmatrixhessenbergunpackq(a, n, tau, s);

      rmatrixhessenbergunpackh(a, n, h);

      //changing indexing of h and s from [0...n-1] to [1...n]
      for(i=0; i<n; i++)
        for(j=0; j<n; j++){
          h1(i+1, j+1)=h(i,j);
          s1(i+1, j+1)=s(i,j);
        }

      if(! upperhessenbergschurdecomposition(h1, n, s2))
          throw std::runtime_error("algorithm for schur decomposition did not converge!");

      ///multiplying orthogonal s1 by orthogonal s2
      for(i=1; i<=n; i++)
        for(j=1; j<=n; j++){
          s3(i, j)=0;
          for(k=1; k<=n; k++){
            s3(i, j)+=s1(i, k)*s2(k, j);
          }
        }

      for(int i =0; i < n; i++){
        for(int j=0; j< n; j++){
          T[j][i]=h1(j+1,i+1);
          S[j][i]=s3(j+1,i+1);
        }
      }
      return true;
  }

  ///calculates schurMatrix of 'der' matrix
  static bool schurMatrix(const DoubleMatrixType& der, DoubleMatrixType& T, DoubleMatrixType& S, std::ostream& out){
    bool r=schurDecompose(der, T, S);
    return r;
  }

  static DoubleMatrixType translateComplexChOfCoordIntoRealOne(const DoubleMatrixType& rVec, const DoubleMatrixType& iVec){
    //we scan imaginary matrix iVec looking for blocks 2x2 corresponding to i*w and -i*w part from eigenvectors v+i*w, v-i*w
    int m_dim = rVec.numberOfColumns();
    int i, j;
    DoubleMatrixType r(m_dim, m_dim);
    for(i = 0; i < m_dim - 1; i += 2)
      for(j = 0; j < m_dim; j++) {
        if(iVec[i][j] != 0 || iVec[i + 1][j] != 0) { //block with eigenvectors in form v+iw, v-iw is detected
          //instead of block [v+iw,v-iw] we take block [v,w]
          r[i][j] = rVec[i][j]; //v_1
          r[i + 1][j] = rVec[i + 1][j]; //v_2
          r[i][j + 1] = iVec[i][j]; //w_1
          r[i + 1][j + 1] = iVec[i + 1][j]; //w_2
          j++;
        } else { //we rewrite real eigenvectors
          r[i][j] = rVec[i][j];
          r[i + 1][j] = rVec[i + 1][j];
        }
      }
    //in case of odd dimension we check last row
    if(m_dim % 2 == 1){
      for(j = 0; j < m_dim; j++) {
        if(iVec[i][j] != 0){
          r[i][j + 1] = iVec[i][j];
          r[i][j] = rVec[i][j];
          j++;
        }else{
          r[i][j]=rVec[i][j];
        }
      }
    }
    return r;
  }

  static void changeOfCoordinates(const DoubleMatrixType& jacobian, DoubleMatrixType& Q, DoubleMatrixType& Qinv, DoubleVectorType& rEigen, DoubleVectorType& iEigen, std::ostream& out){
    int d = jacobian.numberOfRows();
    DoubleMatrixType dT(d, d), dS(d, d), dSinv(d, d), rVec(d, d), iVec(d, d), rVecInv(d, d), dJac(d, d);
    rEigen = DoubleVectorType(d),
    iEigen = DoubleVectorType(d);

    Q = DoubleMatrixType(d, d);
    Qinv = DoubleMatrixType(d, d);

    schurMatrix(jacobian, dT, dS, out);
    dS=capd::vectalg::transpose(dS);
    capd::alglib::computeEigenvaluesAndEigenvectors(dT, rEigen, iEigen, rVec, iVec);
    //we have a complex change of coordinates in two matrices rVec the real part and iVec the imaginary part
    //we unite them into one matrix which generates block diagonal form of T
    rVec = translateComplexChOfCoordIntoRealOne(rVec, -iVec);
    rVecInv = capd::matrixAlgorithms::inverseMatrix(rVec);

    //matrices S together with eigenvectors is our change of coordinates, inverses are calculated rigorously
    Q = rVecInv * dS;
    dSinv = capd::matrixAlgorithms::inverseMatrix(dS);
    Qinv = dSinv * rVec;

  }

};



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
  int n = 12, //change FOJ1D stack dimension
      m = 25,
      order = 7;
  Interval nu,
           step(0.0001);

  ///begin FOJ initialization for FFT integrator
  capd::jaco::DPDEContainer container;
  ///2.set here the subspace of the initial condition e.g. setToRealValuedOdd means that the initial condition is real valued odd
  container.setToRealValuedEven();
  FOJ1D::initialize(ModesContainer::modes2arraySizeStatic(n), container);
  ///end FOJ initialization

  clock_t start, end;

  Index1D idx;
  RealPolynomialBound u_0(n), enclosure(n);
  idx[0] = 0;
  //u_0[idx].re = 1e-10;
  idx[0] = 1;
  u_0[idx].re =  1e-10;

  int STEPS;
  capd::auxil::OutputStream log(std::cout, false, true);
  std::stringstream ss;
  if(testNumber == 0){
    //a test case
    ///3.choose the algorithm type, see enum AlgorithmType in FFTDynSys.h file
    //nu = 6.085;
    //nu = 3.546241427; //scaling 2 * pi
    nu = 3.56;
    //nu = 4;
    //nu = 14.18496571; //scaling pi
    //nu = 56.73986284; //scaling pi/2
    STEPS = 100000;
    //n = 17; //n has to be the same as above
    step = 0.00005;
    ss << "test1_CHproj_";
    srand(time(0));

  /*u_0[Index1D(2)] = 0.330288;
    u_0[Index1D(6)] = -0.0161974;
    u_0[Index1D(10)] = 0.000751484;
    u_0[Index1D(14)] = -3.57331e-05; */

  /*u_0[Index1D(3)] = 0.250095;
    u_0[Index1D(9)] = -0.0030652;
    u_0[Index1D(15)] = 3.72671e-05;*/

    //a case for which the solution goes really close to one attr fixed point, but ends at another fixed point
    for(int i = 1; i < 5; i++){
      //only one mode is forced , should be i = 2, i < 4
      if(i==3){
        if( rand() % 2 == 0 )
          u_0[ Index1D(i) ].re = 0 ;
        else
          u_0[ Index1D(i) ].re = 0 ;
      }else{
        if( rand() % 2 == 0 )
          u_0[ Index1D(i) ].re = 1e-10 ;
        else
          u_0[ Index1D(i) ].re = -1e-10 ;
      }
    }


  }
  if(testNumber == 1){
    //the left panel from DBCP paper p. 3691
    nu = 8;
    STEPS = 100000;
    n = 22;
    step = 0.00001;

    //n at least 17
/*    u_0[Index1D(0)].re = 1e-10;
    u_0[Index1D(1)].re = 0.381048;
    u_0[Index1D(2)].re = 3.02336e-09;
    u_0[Index1D(3)].re = 0.190014;
    u_0[Index1D(4)].re = 2.50279e-09;
    u_0[Index1D(5)].re = 0.0508882;
    u_0[Index1D(6)].re = 1.04236e-09;
    u_0[Index1D(7)].re = 0.0141308;
    u_0[Index1D(8)].re = 3.86516e-10;
    u_0[Index1D(9)].re = 0.0039876;
    u_0[Index1D(10)].re = 1.36241e-10;
    u_0[Index1D(11)].re = 0.00112793;
    u_0[Index1D(12)].re = 4.62271e-11;
    u_0[Index1D(13)].re = 0.000319338;
    u_0[Index1D(14)].re = 1.52424e-11;
    u_0[Index1D(15)].re = 9.03778e-05;
    u_0[Index1D(16)].re = 4.9e-12;
    u_0[Index1D(17)].re = 2.54792e-05;*/

    ss << "test2_CHproj_";
  }
  if(testNumber == 2){
    //the right panel from DBCP paper p. 3691
    //nu = 26.37; //(without scaling the eigenvalues
    nu = 6.59; //(with scaling)
    STEPS = 400000;
    n = 15;
    step = 0.00001;
    ss << "test3_CHproj_";
    srand(time(0));
    rand() % 2;
    for(int i = 1; i < 5; i++){
      if( rand() % 2 == 0 )
        u_0[ Index1D(i) ].re = 1e-10;
      else
        u_0[ Index1D(i) ].re = -1e-10;
    }
  }

  JetDBCP jetDBCP( n, m, nu, PI, order);
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
    if(i % 100 == 0){
      std::cout << dynsys.step << "\n";
      std::cout << u << "\n";
    }
  }

  log << "\nset at the end: " << u << "\n";
  calculateDiams(u, max);
  log << "max diameter of the set at the end: " << max << "\n";

  IntervalMatrix jacobian, Q, Qinv, diag;
  IntervalVector rEigen, iEigen;
  jacobian = jetDBCP.jacobian( u );
  ///diagonalize jacobian calculated at the last point
  Diagonalizator<IntervalMatrix, IntervalMatrix>::changeOfCoordinates( jacobian, Q, Qinv, rEigen, iEigen, std::cout);

  log << "jacobian:\n" << jacobian << "\n";
  diag = Q * jacobian * Qinv;
  for(int i = 0; i < diag.numberOfRows(); i++)
    for(int j = 0; j < diag.numberOfColumns(); j++)
      if( (diag[i][j] > 0. ? diag[i][j] : -diag[i][j]) < 1e-12)
        diag[i][j] = 0.;
  log << "diagonalized:\n" << diag << "\n";

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
