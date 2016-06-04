/* DBCPModelHetConProof.cpp
 *
 *  Created on: Jun 23, 2015
 *      Author: Jacek Cyranka
 *
 *  This file contains an algorithm for proving heteroclinic connections in diblock copolymer model
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
  while( fscanf(file, "%le\n", &nr ) > 0 && i < fixedPoint.size() ){

    fixedPoint[i++] = Interval(nr, nr);

  }

  std::cout << "fixedPoint=\n" << fixedPoint << "\n";

  fclose(file);
}


std::ostream& printModes( const char* string, std::ostream& out, const RealPolynomialBound& pb ){
  out << string << "\n";
  out << "m=" << pb.n << "\n";
  out << "M=" << pb.N << "\n";
  out << "DPDEContainer=\n";
  out << (capd::jaco::DPDEContainer&)pb ;
  for(Index1D index = pb.firstModeIndex(pb.irProjectionPlusFiniteTail); !index.limitReached(pb.irProjectionPlusFiniteTail); index.inc(pb.irProjectionPlusFiniteTail)){
    out << "(" <<  pb[index].re << "," << pb[index].im << ")\n";
  }
  out << pb.farTail.getC() << "\n";
  out << pb.farTail.getS() << "\n";
  return out;
}





void integrate(int approach){

  int m, /*change FOJ1D stack dimension*/ M, dftPts, dftPts2, order;
  double step_, nu_, sigma_, piOverL_;
  Interval nu, step;

  std::ofstream currentSet;
  loadDataFromFile("config.in", &m, &M, &dftPts, &dftPts2, &order, &step_, &nu_, &sigma_, &piOverL_ );

  step = step_;
  nu = nu_;

  ///1.initialize FOJ, the dimension of the first order jet and the initial condition type defined by a DPDEContainer instance
  ///should be provided, see FOJ1D::initialize description

  clock_t start, end;

  int MAXSTEPS = 200000;;
  capd::auxil::OutputStream log(std::cout, false, true);
  std::stringstream ss;

  ///begin FOJ initialization for FFT integrator
  capd::jaco::DPDEContainer container;
  ///2.set here the subspace of the initial condition e.g. setToRealValuedOdd means that the initial condition is real valued odd
  container.setToRealValuedEven();
  FOJ1D::initialize(ModesContainer1D::modes2arraySizeStatic(m), container);
  ///end FOJ initialization

  RealPolynomialBound u_0(m, M, container), basin(m, M, container), enclosure(m);

  PolyBdInputReader<RealPolynomialBound> inputReader("manifold.in");
  u_0 =  inputReader.polyBd;

  u_0.changeM(100, false);

  std::cout << "Input polynomial bound=\n" << u_0 << "\n";

  PolyBdInputReader<RealPolynomialBound> basinReader("basin.in");
  basin = basinReader.polyBd;


  std::cout << "basin of attraction=\n" << basin << "\n";


  ss << "DBCPModelHetConProof_";
  DPDEInclusionCW3 diffIncl(m, dftPts, M, dftPts2, PI, nu, order, step, MaxNorm());
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
  log << "Taylor method order=" << order << ", constant time step=" << step << "\n";
  log << "Galerkin projection dimension m=" << m << ", M_{FFT}=" << dftPts << "\n";
  log << "the whole infinite dimensional system is being integrated (the Lohner algorithm for differential inclusions is used).\n";
  if(approach == 0){
    diffIncl.getDynamicalSystem().setJetDynSysAlgorithmType(capd::jaco::direct);
    log << "Using the direct approach\n";
  }else{
    if(approach == 1){
      diffIncl.getDynamicalSystem().setJetDynSysAlgorithmType(capd::jaco::FFT);
      log << "Using the FFT approach\n";
    }else{
      diffIncl.getDynamicalSystem().setJetDynSysAlgorithmType(capd::jaco::FFTButFirstOrderDirect);
      log << "Using the FFT approach, but the first order normalized derivative is calculated directly to avoid a blowup\n";
    }
  }

  log << "initial condition (infinite dimensional):\n" << u_0 << "\nitegration started, the output below is the Galerkin projection of the set at each timestep\n";
  int i;
  //  //end init mo
  InclRect2Set set(u_0);

  Interval max;

  start = clock();
  for(i=0; i < MAXSTEPS; ++i){
//    std::cout << "step #" << i << " ";
    set.move(diffIncl);
    log << (IntervalVector)set << "\n";
    if(i % 100 == 0){
      log << "The set at time step #" << i << " ( h = " << step << " )\n";
      log << set.getPerturbationParams() << "\n";
      currentSet.open("current_set.out", std::ofstream::out);
      std::stringstream ss;
      ss << "current set step #" << i << " of the integration process";
      printModes(ss.str().c_str(), currentSet, set.getPerturbationParams() );
      currentSet.close();

      if( set.getPerturbationParams().subset(basin) ){
        log << "\nTHE CURRENT SET AT THE TIME STEP i=" << i << " IS WITHIN THE BASIN OF ATTRACTION.\n END OF INTEGRATION IN TIME PROCESS.\n";
        break;
      }
    }
  }
  if( i == MAXSTEPS ){
    log << "\nTHE MAXIMAL NUMBER OF TIME STEPS REACHED BEFORE THE SET WAS INTEGRATED INTO THE BASIN OF ATTRACTION.\n";
    log << "\nset at the end (infinite dimensional): " << set.getPerturbationParams() << "\n";
  }
  time ( &rawtime );
  log << "End of the integration current local time is: " << ctime (&rawtime) << "\n";
  calculateDiams((IntervalVector)set, max);
  log << "max diameter of the set at the end: " << max << "\n";
  end = clock();
  log << "time: "<<(end-start)/CLOCKS_PER_SEC<<"\n";
  log << "RRmultiplicationsSum=" << RRmultiplicationsSum << "\n";
  log << "RRadditionsSum=" << RRadditionsSum << "\n";

}


//calculates the logarithmic norm for the provided self-consistent bounds
//todo: NOW IT IS DONE IN FINITE DIMENSION
Interval calculateLogNorm(Box& Qbox, RealPolynomialBound& tail, IntervalVector& eigenvaluesRe, JetSHDPDE& jetDPDE, std::ostream& out){

  std::ofstream currentSet;

  IntervalVector V = Qbox.wrapAffine();

  out << "V:\n" << V << "\n";

  jetDPDE.useFFT = false;

  IntervalMatrix jacobian = jetDPDE.jacobian( V );
  int n = jacobian.numberOfRows();


  IntervalMatrix dFtilde = Qbox.m_Q * jacobian * Qbox.m_Qinv;

  std::cout << "dFtilde:\n" << dFtilde << "\n";
  std::cout << "mdisc:\n" << Qbox.m_disc << "\n";
  std::cout << "eigenvaluesRe:\n" << eigenvaluesRe << "\n";


  int m, /*change FOJ1D stack dimension*/ M, dftPts, dftPts2, order;
  double step_, nu_, sigma_, piOverL_;

  loadDataFromFile("config.in", &m, &M, &dftPts, &dftPts2, &order, &step_, &nu_, &sigma_, &piOverL_ );

  PartialDerivativePDE pdpde(m, M, dftPts, dftPts2, nu_, PI, order);
  pdpde.m_N_coeff = 1.;

  tail = Qbox.wrapAffine();
  tail.copyTailPartFrom( Qbox.getTail() );

  RealPolynomialBound r( tail );
  const RealPolynomialBound& rc(r);

  pdpde.CalculateNonlinearTermDirectly( tail, r , capd::jaco::full, false);

  IntervalMatrix expJacobian(M, M);

  for(int i = 0; i < M; i++){
    for(int j = 0; j < M; j++){
      expJacobian[i][j] = i * i * jetDPDE.piOver_l * jetDPDE.piOver_l * jetDPDE.m_N_coeff.re * 3 * ( r[ Index1D( abs(i-j) ) ]  +  r[ Index1D( abs(i+j) ) ] ).re ;
    }
  }


  IntervalVector sums( M ), conecond( M ), sumst( Qbox.m_n );
  Interval sumi, psum1, psum2;


  int i = 38; //manually checked where is the smallest term on diagonal
  sumi = 0.; psum1 = 0.; psum2 = 0.;
  //the 1st component of the cone condition sum
  for(int j=1; j < M + i; j++){ //this loop has to start from 1 , as we do not differentiate with resp to a_0 -- it is constant
    if( i != j ){
      if(j < Qbox.m_n)
        sumi += iabs(dFtilde[i][j] + dFtilde[j][i]);
      else{
        for(int k=0; k < Qbox.m_n; k++){
          psum1 += Qbox.m_Q[i][k] * ( r[ Index1D( abs(k-j) ) ]  +  rc[ Index1D( abs(k+j) ) ] ).re * k * k * jetDPDE.piOver_l * jetDPDE.piOver_l * jetDPDE.m_N_coeff.re * 3;
        }

        for(int k=0; k < Qbox.m_n; k++){
          psum2 += Qbox.m_Qinv[k][i] * ( r[ Index1D( abs(j - k) ) ]  +  rc[ Index1D( abs(k + j) ) ] ).re * j * j * jetDPDE.piOver_l * jetDPDE.piOver_l * jetDPDE.m_N_coeff.re * 3;
        }
        sumi += iabs( psum1 + psum2 );
      }
    }
  }

  Interval logNorm = sumi  + dFtilde[i][i];

  return logNorm;

}


void proveStableFixedPoint(){
  int m, /*change FOJ1D stack dimension*/ M, dftPts, dftPts2, order;
  double step_, nu_, sigma_, piOverL_;
  Interval nu, step;

  std::ofstream currentSet;
  loadDataFromFile("config.in", &m, &M, &dftPts, &dftPts2, &order, &step_, &nu_, &sigma_, &piOverL_ );

  step = step_;
  nu = nu_;

  IntervalVector fixedPoint(m + 1), zeroCenteredBox(m + 1), eigenvaluesRe(m + 1);
  DoubleVector dfPoint( m + 1 );

  loadApproximateFixedPoint( "fixedPoint.in", fixedPoint );
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

  std::cout << "ModesContainer1D::modes2arraySizeStatic(m)=" << ModesContainer1D::modes2arraySizeStatic(m) << "\n";
  std::cout << "ModesContainer1D::modes2arraySizeStatic(m)=" << ModesContainer_::modes2arraySizeStatic(m) << "\n";

  JetSHDPDE jetDPDE(m, dftPts, nu, PI, order);

  IntervalMatrix jacobian = jetDPDE.jacobian( fixedPoint );

  int n = jacobian.numberOfRows();
  IntervalMatrix T(n, n);

  DoubleMatrix dJacobian( n, n ), E(n, n) ;

  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      dJacobian[i][j] = rightBound( jacobian[i][j] );

  SHDPDE dpde(m, M, dftPts, dftPts2, nu, PI, order);

  BoxFinder finder(dpde, jetDPDE);

  Box Qbox(m + 1);

  IntervalVector FevaluatedAtFP(m + 1);


  DOUBLE small = 1e-10;
  zeroCenteredBox = Interval(-small, small);

  RealPolynomialBound tail(m, M, container), tr(m, M, container);
  setS ( tail, 6 );

  SHDPDE_ basicDpde(m, M, dftPts, dftPts2, nu_, 3.1415926535897932384626433832795, order);
  JetSHDPDE_ basicJetDpde(m, dftPts, nu_, 3.1415926535897932384626433832795, order);

  RealPolynomialBound_ u( m ), f( m );
  DoubleVector delta( m + 1 ), dF( m + 1 );
  u = dfPoint;
  DoubleMatrix uJac = basicJetDpde.jacobian(u), uJacInv;

  basicDpde(u, f);
  dF = f;
  uJacInv = capd::matrixAlgorithms::inverseMatrix( uJac );
  delta = uJacInv * dF;
  dfPoint = dfPoint - delta;
  u = dfPoint;
  basicDpde(u, f);


  finder.validateBox(dfPoint, dJacobian,
      zeroCenteredBox, ///a zero centred box (box around fixed point - fixed point)
      Qbox, tail, E, eigenvaluesRe, T,
      FevaluatedAtFP, std::cout);

  std::cout << "FevaluatedAtFP: " << FevaluatedAtFP << "\n";

  finder.inflateTrappingRegion(Qbox, eigenvaluesRe, T, dfPoint, FevaluatedAtFP, std::cout , 23);

  Interval logNorm = calculateLogNorm(Qbox, Qbox.getTail(), eigenvaluesRe, jetDPDE, std::cout);
  if(logNorm < 0){
    std::cout << "logarithmic norm is NEGATIVE, logNorm =" << logNorm << "\n";
  }else{
    std::cout << "logarithmic norm is POSITIVE, logNorm =" << logNorm << "\n";
  }

  //the found trapping region is in Q coordinates (diagonalizing the coordinate system at the fixed point)
  //we find a set given in canonical coordinates, which is contained in the trapping region in Q coordinates.
  IntervalVector Qwrap = Qbox.wrapAffine();
  std::cout << "Qwrap:\n" << Qwrap << "\n";

  const DOUBLE FACTOR = 1000;
  IntervalVector v = ( Qbox.m_Qinv * Qbox.m_v ) / FACTOR;

  //check if the box is inside the trapping region
  IntervalVector t = Qbox.m_Q * v;

  for(int i=0; i < t.size(); i++)
    if( ! Qbox.m_v[i].contains( t[i] ) )
      std::cout << "t[" << i << "] is not a subset of v[" << i << "] (" << t[i] << ", and " << Qbox.m_v[i] << "\n";

  tr = Qwrap;
  tr.copyTailPartFrom( Qbox.getTail() );
  std::cout << "m_v:\n" << Qbox.m_v;
  std::cout << "tail:\n" << tr;

  std::ofstream basin;
  basin.open("basin.in", std::ofstream::out);
  printModes("Proved_basin_of_attraction_of_the_stable_fixed_point:", basin, tr );
  basin.close();

  std::ofstream qinv;
  qinv.open("Qinv.in", std::ofstream::out);

  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++){
      qinv << Qbox.m_Qinv[i][j] << " ";
    }

  qinv.close();

}


//verifies the cone condition in the provided self-consistent bounds
//todo: NOW IT IS DONE IN FINITE DIMENSION
bool verifyConeCondition(Box& Qbox, RealPolynomialBound& tail, IntervalVector& eigenvaluesRe, JetSHDPDE& jetDPDE, std::ostream& out){

  std::ofstream currentSet;

  IntervalVector V = Qbox;

  out << "V:\n" << V << "\n";

  jetDPDE.useFFT = false;

  IntervalMatrix jacobian = jetDPDE.jacobian( V );
  int n = jacobian.numberOfRows();

  RealPolynomialBound tr(tail);
  tr = Qbox.wrapAffine();
  tr.copyTailPartFrom( Qbox.getTail() );

  IntervalVector Q( Qbox.m_n );

  for(int i = 0; i < Qbox.m_n; i++){
    if( eigenvaluesRe[i] < 0 ){
      Q[i] = -1;
    }else{
      Q[i] = 1;
    }
  }

  IntervalVector condition( Qbox.m_n );

  Interval sum = 0.;
  for(int i = 0; i < Qbox.m_n; i++){
    condition[i] = 2 * ( jacobian[i][i] > 0 ? jacobian[i][i] : - jacobian[i][i]);

    sum = 0.;
    for(int j = 0; j < Qbox.m_n; j++){
      if(i != j){
        sum += Q[i] * jacobian[i][j] + Q[j] * jacobian[j][i];
      }
    }

    condition[i] -= ( sum > 0 ? sum : -sum );

  }

  int m, /*change FOJ1D stack dimension*/ M, dftPts, dftPts2, order;
  double step_, nu_, sigma_, piOverL_;
  loadDataFromFile("config.in", &m, &M, &dftPts, &dftPts2, &order, &step_, &nu_, &sigma_, &piOverL_ );

  PartialDerivativePDE pdpde(m, M, dftPts, dftPts2, nu_, PI, order);
  pdpde.m_N_coeff = 1.;

  RealPolynomialBound r( tail );
  const RealPolynomialBound& rc(r);

  pdpde.CalculateNonlinearTermDirectly( tr, r , capd::jaco::full, false);
  r[0] = 0.;

  IntervalMatrix expJacobian(M, M);

  out << "r:\n" << r << "\n";

  for(int i = 0; i < M; i++){
    for(int j = 0; j < M; j++){
      expJacobian[i][j] = i * i * jetDPDE.piOver_l * jetDPDE.piOver_l * jetDPDE.m_N_coeff.re * 3 * ( r[ Index1D( abs(i-j) ) ]  +  r[ Index1D( abs(i+j) ) ] ).re ;
    }
  }

  IntervalVector sums( M ), conecond( M ), sumst( Qbox.m_n );
  Interval sumi, psum1, psum2;
  for(int i=0; i < M; i++){
    sumi = 0.;

    for(int j=1; j < M + i; j++){ //this loop has to start from 1 , as we do not differentiate with resp to a_0 -- it is constant
      if( i != j ){
        //the 1st component of the cone condition sum \partial F_i / \partial x_j
        if(i < Qbox.m_n)
          psum1 = Q[i] * ( r[ Index1D( abs(i - j) ) ]  +  rc[ Index1D( abs(i + j) ) ] ).re;
        else //all eigenvalues are stable for j > m, so Q[j] = -1
          psum1 = -1. * ( r[ Index1D( abs(i - j) ) ]  +  rc[ Index1D( abs(i + j) ) ] ).re;

        psum1 *= i * i * jetDPDE.piOver_l * jetDPDE.piOver_l * jetDPDE.m_N_coeff.re * 3;

        //the 2nd component of the cone condition sum \partial F_j / \partial x_i
        if(j < Qbox.m_n)
          psum2 = Q[j] * ( r[ Index1D( abs(j - i) ) ]  +  rc[ Index1D( abs(i + j) ) ] ).re;
        else //all eigenvalues are stable for j > m, so Q[j] = -1
          psum2 = -1. * ( r[ Index1D( abs(j - i) ) ]  +  rc[ Index1D( abs(i + j) ) ] ).re;

        psum2 *= j * j * jetDPDE.piOver_l * jetDPDE.piOver_l * jetDPDE.m_N_coeff.re * 3;
      }
      sumi += iabs(psum1 + psum2);
    }
    // j > M + i part (all terms are in the tail)
    //add the infinite series contribution for the 1st component sum \partial F_i / \partial x_j
    sumi += 2 * C(r) * ( power( M, -s(r) + 1 ) / (s(r) - 1) ) * i * i * jetDPDE.piOver_l * jetDPDE.piOver_l * jetDPDE.m_N_coeff.re * 3;
    //add the infinite series contribution for the 2nd component sum \partial F_j / \partial x_i -- have to decrease the power by two due to the j^2 factor
    //additional factor (1 + i/M) comes from estimating (j+i)^2 \leq j^2( 1 + i/M )
    sumi += 2 * C(r) * ( power( M, -(s(r) - 2) + 1 ) / ((s(r) - 2) - 1) ) * jetDPDE.piOver_l * jetDPDE.piOver_l * jetDPDE.m_N_coeff.re * 3 * ( 1 + i / M );

    sums[i] = sumi;
  }


  //VERIFICATION FOR i >= M

  //2 * suma modow z r + nu*C/k^{s-2} mniejsze niz 0?

  for(int i=0; i < Qbox.m_n; i++){
    sumi = 0.;
    for(int j=0; j < Qbox.m_n; j++){
      if(i != j)
        sumi += Q[j] * jacobian[i][j];
    }
    sumst[i] = sumi;
  }

  Interval v;
  for(int i=0; i < M; i++){
    v = ( pdpde.lambda(i) + rc[ Index1D( abs(i+i) ) ].re * i * i * jetDPDE.piOver_l * jetDPDE.piOver_l * jetDPDE.m_N_coeff.re * 3 );
    v = iabs( v );
    conecond[i] = 2 * v ;

    conecond[i] -= sums[i];
  }

  //verify the cone condition
  //for each i we check if conecond[i] > epsilon

  Interval epsilon = 1e-05;
  bool satisfied = true;
  for(int i=0; i < M; i++){
    if( leftBound( conecond[i] ) < epsilon )
      satisfied = false;
  }

  out << "conecond full vector: \n" << conecond << "\n";

  Interval min = 10000000000.;
  for(int i=0; i < M; i++){
    if( leftBound(conecond[i]) < min )
      min = leftBound(conecond[i]);
  }

  out << "CONE CONDITION EPSILON=" << min << "\n";

  out.flush();
  //out << "sums: " << sums << "\n";
  //out << "sumst: " << sumst << "\n";

  return satisfied;

}




void proveUnstableManifold(int direction){
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

  std::cout << "ModesContainer1D::modes2arraySizeStatic(m)=" << ModesContainer1D::modes2arraySizeStatic(m) << "\n";
  std::cout << "ModesContainer1D::modes2arraySizeStatic(m)=" << ModesContainer_::modes2arraySizeStatic(m) << "\n";

  JetSHDPDE jetDPDE(m, dftPts, nu, PI, order);

  IntervalMatrix jacobian = jetDPDE.jacobian( fixedPoint );
  int n = jacobian.numberOfRows();
  IntervalMatrix T(n, n),
                 Q(n, n);

  DoubleMatrix dJacobian( n, n ), E(n, n) ;

  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      dJacobian[i][j] = rightBound( jacobian[i][j] );

  SHDPDE dpde(m, M, dftPts, dftPts2, nu, PI, order);

  BoxFinder finder(dpde, jetDPDE);

  Box Qbox(m + 1);

  IntervalVector FevaluatedAtFP(m + 1);

  const DOUBLE SMALL = 1e-12;

  const DOUBLE MANIFOLD = 0.75e-01;

  zeroCenteredBox = Interval( -SMALL, SMALL );

  zeroCenteredBox[direction] = Interval( -MANIFOLD, MANIFOLD );

  RealPolynomialBound tail(m, M, container), tr(m, M, container);
  setS ( tail, 6 );

  BoolVector directions(n);

  for(i = 0; i < n; i++){
    if(dJacobian[i][i] > 0){
      directions[i] = 1;
      Q[i][i] = 1;
    }else{
      if(dJacobian[i][i] < 0){
        directions[i] = -1;
        Q[i][i] = -1;
      }else{
        std::cerr << "eigenvalue 0 -- fixed point zero is not hyperbolic ??? \n";
        exit(1);
      }
    }
  }

  std::ofstream boxDebug("box_debug.txt");

  finder.validateBox(dfPoint, dJacobian,
      zeroCenteredBox, ///a zero centred box (box around fixed point - fixed point)
      Qbox, tail, E, eigenvaluesRe, T,
      FevaluatedAtFP, boxDebug);

  std::cout << "VERIFYING CONE CONDITION!\n";

  bool ccsatisfy = verifyConeCondition( Qbox, tail, eigenvaluesRe, jetDPDE, boxDebug );

  tr = Qbox;

  tr.copyTailPartFrom( tail );

  //detaching the piece of the boundary from the set in the unstable direction
  tr.set(Index1D(direction), ComplexScalar( rightBound(tr[ Index1D(direction) ].re), 0. ) );
  tr.set(Index1D(0), ComplexScalar(0.,0.) );

  std::cout << "the computed self-consistent bounds for the unstable (zero) fixed point:\n";
  std::cout << tr << "\n";

  std::ofstream manifold;
  manifold.open("manifold.in", std::ofstream::out);

  printModes("isolating_neighbourhood_about_zero", manifold, tr);

  manifold.close();

  if(ccsatisfy){
    std::cout << "CONE CONDITION SATISFIED\n";
  }else{
    std::cerr << "error - cone condition is not satisfied!\n";
    exit(1);
  }

  boxDebug.close();

}




int main(int argc, char * argv[]){
  setLoggers();

  if(argc != 3){
      std::cerr << "Number of arguments wrong.\nUsage: ./DBCPModelHetConProof [integrate|fixedpoint|manifold] approach_number\n" <<
                   "\napproach_number is the approach type, 0 - the direct approach, 1 - the FFT approach, 2 - the FFT approach, but the first normalized derivative is calculated directly, in order to avoid blow-ups.\n";
    }else{
      if( strcmp(argv[1], "integrate") == 0 )
        integrate( atoi( argv[2] ) );
      else
        if( strcmp(argv[1], "fixedpoint") == 0 )
          proveStableFixedPoint();
        else
          if( strcmp(argv[1], "manifold") == 0 ){
            proveUnstableManifold(3);
          }
    }

  FOJ1D::destroy();

}
