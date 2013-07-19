/*
 * PolyBd2DTests.cpp
 *
 *  Created on: Jun 27, 2013
 *      Author: cyranka
 */

#include <time.h>

#include "capd/filib/Interval.h"
#include "capd/intervals/Interval.hpp"
#include "config.h"
#include <iomanip>

///here 'long double' has to be used instead of 'double' because of some rounding problems
///(result from FFT is correct , but direct evaluation result was not subset of FFT result due to rounding error issues)
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
#include <assert.h>

typedef capd::vectalg::Vector<Interval, 0> IntervalVector;
typedef capd::jaco::ComplexScalar<Interval> ComplexScalar;
typedef capd::jaco::Index2D Index;
typedef capd::jaco::MaximumNorm<Index> MaximumNorm;
typedef capd::jaco::Odd<Interval, Index, MaximumNorm > OddSubspace;
typedef capd::jaco::Even<Interval, Index, MaximumNorm > EvenSubspace;
typedef capd::jaco::RealPolynomialBound<ComplexScalar, Index, 0> RealPolynomialBound;// there is MaximumNormFixed
//typedef capd::jaco::RealPolynomialBound<ComplexScalar, Index, 0, EvenSubspace> RealPolynomialBound;

typedef capd::jaco::Burgers<RealPolynomialBound> Burgers;

//typedef capd::jaco::FFT2DOneComponent<ComplexScalar, ComplexScalar, MaximumNorm, 0, 0, RealPolynomialBound> FFT2D;
typedef capd::jaco::FFT2DOneComponent<ComplexScalar, ComplexScalar, MaximumNorm, 0, 0, RealPolynomialBound> FFT2D;

typedef capd::jaco::DPDE2<Burgers, FFT2D, 0> DPDE2;

typedef FFT2D::DFTGridType DFT2DGridType;
typedef FFT2D::ModesContainerType Modes2DContainer;


enum testPolyBds{odd, even};

struct TestCase{
  int m, M, w, dft_m, dft_M, s;
  Interval C, diam;
  TestCase(int m_, int M_, int w_, int dft_m_, int dft_M_, int s_, Interval C_, Interval diam_) : m(m_), M(M_), w(w_), dft_m(dft_m_), dft_M(dft_M_), s(s_), C(C_), diam(diam_){}
  TestCase(){}
};








void printTC(capd::auxil::OutputStream& out, TestCase tc){
  out << "Current test case:\n" << "m=" << tc.m << ", M=" << tc.M << ", w=" << tc.w << ", dft_m=" << tc.dft_m << ", dft_M=" << tc.dft_M << ", s=" << tc.s << "\n";
}

TestCase testCases[100];
int numberOfTestCases;

void setTestCases(){
  numberOfTestCases = 10;
  testCases[0] = TestCase(5, 49, 2, 16, 150, 5, 1e+02, 1e-10 * Interval(-1., 1.));
  testCases[1] = TestCase(10, 100, 2, 32, 320, 5, 1e+02, 1e-10 * Interval(-1., 1.));
  testCases[2] = TestCase(5, 49 , 2, 16, 150, 7, 1e+3, 1e-10 * Interval(-1., 1.));
  testCases[3] = TestCase(10, 100, 2, 32, 320, 7, 1e+3, 1e-10 * Interval(-1., 1.));
  testCases[4] = TestCase(5, 49, 3, 16, 150, 5, 1e+3, 1e-10 * Interval(-1., 1.));
  testCases[5] = TestCase(10, 100, 3, 32, 320, 5, 1e+3, 1e-10 * Interval(-1., 1.));
  testCases[6] = TestCase(5, 24, 2, 16, 75, 5, 10, 1e-10 * Interval(-1., 1.));
  testCases[7] = TestCase(5, 24, 2, 16, 75, 5, 10, 1e-6 * Interval(-1., 1.));
  testCases[8] = TestCase(5, 49, 2, 16, 150, 5, 10, 1e-6 * Interval(-1., 1.));
  testCases[9] = TestCase(5, 49, 4, 16, 150, 5, 1e+3, 1e-10 * Interval(-1., 1.));
}

capd::auxil::OutputStream testDebug(std::cout, false, true);
capd::auxil::OutputStream testData(std::cout, false, true);

RealPolynomialBound constructPolyBdForResult(TestCase tc, int subcase){
  capd::jaco::DPDEContainer container;
  container.setToRealValuedEven();//works for assumption that functions are either odd or even

  RealPolynomialBound polybd(tc.m, tc.M, container, tc.w);
  return polybd;
}

DFT2DGridType constructDFT(TestCase tc){
  DFT2DGridType dft(tc.M);
  return dft;
}



/**Given C>0, s>0 Constructs a polynomial bound with random values satisfying <= C/|k|^s and <= 1
 */
RealPolynomialBound constructRandomPolyBd(TestCase tc, int subcase){
  capd::jaco::DPDEContainer container;
  container.setToRealValuedEven();

  RealPolynomialBound polybd(tc.m, tc.M, container, tc.w);
  srand(time(0));
  for(int i = -tc.M; i <= tc.M; ++i){
    for(int j = -tc.M; j <= tc.M; ++j){
      Index ind = Index(i, j);
      if(ind.upperHalfspace() && polybd.irProjectionPlusFiniteTail.withinRange(ind) && !ind.isZero()){
        Interval c = (double(rand())/RAND_MAX);
        Interval sign = rand() % 2 == 0 ? 1. : -1.;
        c = tc.C / power(ind.maxNorm(), tc.s) > 1. ? 1. : c * tc.C / power(ind.maxNorm(), tc.s);
        c *= sign;
        c += tc.diam;

        ComplexScalar cs;
        if(subcase == even)
          cs = ComplexScalar(c, 0.);
        if(subcase == odd)
          cs = ComplexScalar(0., c);
        polybd[ind] = cs;
      }
    }
  }
  setC(polybd, tc.C);
  setS(polybd, tc.s);
  return polybd;
}


FFT2D constructFft(TestCase tc){
  FFT2D fft(tc.M, tc.dft_M, Interval::pi());
  return fft;
}

DPDE2 constructDPDE(TestCase tc){
  DPDE2 dpde(tc.m, tc.M, tc.dft_m, tc.dft_M, 1, Interval::pi(), true, 1, tc.w);
  return dpde;
}


void basicTest(){
  int m = 5,
      M = 49,
      dft_m = 16,
      dft_M = 150,
      w = 2;

  TestCase tc(m, M, w, dft_m, dft_M, 6, 100, 1e-10 * Interval(-1., 1.));

  RealPolynomialBound polybd = constructRandomPolyBd(tc, even),
                      polybd2 = constructRandomPolyBd(tc, even),
                      r = constructPolyBdForResult(tc, even),
                      r_opt = constructPolyBdForResult(tc, even),
                      rPS = constructPolyBdForResult(tc, even),
                      c = constructPolyBdForResult(tc, even);

  FFT2D fft = constructFft(tc);
  
  DFT2DGridType grid = constructDFT(tc),
                s = constructDFT(tc),
                sPS = constructDFT(tc),
                grid2 = constructDFT(tc),
                gridPS = constructDFT(tc),
                grid2PS = constructDFT(tc);
  DPDE2 dpde2 = constructDPDE(tc);
  std::cout << "dpde2.w=" << dpde2.w << "\n";

  typename DPDE2::IndexRangeType range1, range2, range3, range4;
  range1 = polybd.irFull;
  range2.setRange(tc.w * tc.M, capd::jaco::strong, (tc.w) * tc.M + 10, capd::jaco::weak);
  range3.setRange(tc.w * tc.M, capd::jaco::weak, -1, capd::jaco::weak);
  range4.setRange(tc.w * tc.M, capd::jaco::weak, tc.w * tc.w * tc.w * tc.M, capd::jaco::weak);

  typedef DPDE2::IndexType IndexType;
  Interval hs = IndexType::harmonicSumK_1<Interval, typename DPDE2::IndexRangeType, capd::jaco::EuclideanNorm<IndexType> >(range3, tc.s);

  Interval sum(0);
  for(DPDE2::IndexType i = dpde2.firstModeIndex(range4); !i.limitReached(range4); i.inc(range4)){
    sum += 1. / power(i.squareEuclNorm(), tc.s / 2.);
  }
  std::cout << "hs=" << hs << ", sum=" << sum << "\n";
  assert(sum <= hs);

  time_t start = clock();
  //dpde2.addBound(polybd, polybd2, r, false);
  time_t end = clock();
  //std::cout << "convoluiton time=" << (end-start) / CLOCKS_PER_SEC << "\n";

  start = clock();
  dpde2.calculateConvolution(polybd, polybd2, r_opt);
  end = clock();
  std::cout << "optimized convoluiton time=" << (end-start) / CLOCKS_PER_SEC << "\n";

  //testDebug << "polyBd: " << polybd << ", polyBd2: " << polybd2 << "\n";
  testDebug << "(optimizedFFT) r:\n" << r_opt << "\n\n";

  std::cout << "r_opt subset r: " << r_opt.subset(r) << "\n";

}

enum NormTypes{maximumNorm, euclideanNorm};

/**
 * This test is testing the w_1 constant from the paper
 */
void constantsTest(TestCase tc, int normType){
  RealPolynomialBound polybd = constructPolyBdForResult(tc, even);
  DPDE2 dpde2 = constructDPDE(tc);

  typename DPDE2::IndexRangeType range1, range2, range3, range4;
  range1 = polybd.irProjectionPlusFiniteTail;
  range2.setRange(tc.w * tc.M, capd::jaco::strong, (tc.w) * tc.M + 10, capd::jaco::weak);


  Interval w1;
  if(normType == maximumNorm)
    w1 = 1. - 1. / tc.w;
  if(normType == euclideanNorm)
    w1 = 1. - sqrt(2.) / tc.w;
  Interval norm1, norm2;

  for(DPDE2::IndexType k = dpde2.firstModeIndex(range2); !k.limitReached(range2); k.inc(range2)){
    for(DPDE2::IndexType k1 = dpde2.firstModeIndex(range1); !k1.limitReached(range1); k1.inc(range1)){
      if(normType == maximumNorm){
        norm1 = (k-k1).maxNorm();
        norm2 = w1 * k.maxNorm();
      }
      if(normType == euclideanNorm){
        norm1 = (k-k1).euclNorm<Interval>();
        norm2 = w1 * k.euclNorm<Interval>();
      }
      assert(norm1 >= norm2);
      generalDebug << k[0] << " " << k[1] << " " << k1[0] << " " << k1[1] << " " << norm1 << " " << norm2 << "\n";
    }
  }

}


void optimizedFFTIncludedInFFTTest(TestCase tc, int subspace){
  RealPolynomialBound polybd = constructRandomPolyBd(tc, subspace),
                      polybd2 = constructRandomPolyBd(tc, subspace),
                      r = constructPolyBdForResult(tc, subspace),
                      r_opt = constructPolyBdForResult(tc, subspace),
                      rPS = constructPolyBdForResult(tc, subspace),
                      c = constructPolyBdForResult(tc, subspace);

  FFT2D fft = constructFft(tc);

  DFT2DGridType grid = constructDFT(tc),
                s = constructDFT(tc),
                sPS = constructDFT(tc),
                grid2 = constructDFT(tc),
                gridPS = constructDFT(tc),
                grid2PS = constructDFT(tc);
  DPDE2 dpde2 = constructDPDE(tc);

  time_t start = clock();
  dpde2.calculateConvolution(polybd, polybd2, r, false);
  time_t end = clock();

  start = clock();
  dpde2.calculateConvolution(polybd, polybd2, r_opt, true);
  end = clock();

  testData.logfile("lastTestData.txt", true);
  testData << "optimizedFFTIncludedInFFTTest\n";
  testData << "polybd:\n" << polybd << "\n" << "polybd2:\n" << polybd2 << "\n" << "r_opt:\n" << r_opt << "\nr:\n" << r << "\n";

  assert(r_opt.subset(r));

  testDebug << "assert(r_opt.subset(r)) passed succesfully\n";

}


void ConvolutionWithTheSameFunctionTest(TestCase tc, int subspace){
  RealPolynomialBound polybd = constructRandomPolyBd(tc, subspace),
                      r1 = constructPolyBdForResult(tc, subspace),
                      r2 = constructPolyBdForResult(tc, subspace),
                      rPS = constructPolyBdForResult(tc, subspace),
                      c = constructPolyBdForResult(tc, subspace);

  FFT2D fft = constructFft(tc);

  DFT2DGridType grid = constructDFT(tc),
                s = constructDFT(tc),
                sPS = constructDFT(tc),
                grid2 = constructDFT(tc),
                gridPS = constructDFT(tc),
                grid2PS = constructDFT(tc);
  DPDE2 dpde2 = constructDPDE(tc);

  time_t start = clock();
  dpde2.calculateConvolution(polybd, r1, true);
  time_t end = clock();

  start = clock();
  dpde2.calculateConvolution(polybd, polybd, r2, true);
  end = clock();
  DPDE2::IndexRangeType range1  = r2.irProjectionPlusFiniteTail;

  Interval rad(-1e-16, 1e-16); //there is difference sometimes at 1e-7 , but this is considered negligible in this test
  for(DPDE2::IndexType i = dpde2.firstModeIndex(range1); !i.limitReached(range1); i.inc(range1)){
    r2[i] += ComplexScalar(rad, rad);
  }

  testData.logfile("lastTestData.txt", true);
  testData << "ConvolutionWithTheSameFunctionTest\n";
  testData << "polybd:\n" << polybd << "\n" << "r1:\n" << r1 << "\nr2:\n" << r2 << "\n";


  setC(r1, 0);
  setC(r2, 0);

  assert(r1.subset(r2));

  testDebug << "assert(r1.subset(r2)) passed succesfully\n";

}


void boundsTest(TestCase tc, int subspace){
  RealPolynomialBound polybd = constructRandomPolyBd(tc, subspace),
                      polybd2 = constructRandomPolyBd(tc, subspace),
                      r1 = constructPolyBdForResult(tc, subspace),
                      r2 = constructPolyBdForResult(tc, subspace),
                      r3 = constructPolyBdForResult(tc, subspace),
                      r4 = constructPolyBdForResult(tc, subspace),
                      rPS = constructPolyBdForResult(tc, subspace),
                      c = constructPolyBdForResult(tc, subspace);

  FFT2D fft = constructFft(tc);

  DFT2DGridType grid = constructDFT(tc),
                s = constructDFT(tc),
                sPS = constructDFT(tc),
                grid2 = constructDFT(tc),
                gridPS = constructDFT(tc),
                grid2PS = constructDFT(tc);
  DPDE2 dpde2 = constructDPDE(tc);

  dpde2.addBound(polybd, r1);

  dpde2.addBound(polybd, polybd2, r2);


  DPDE2::IndexRangeType range1, range2, range3, range4;
  range1 = polybd.irFull;
  range2.setRange(tc.w * tc.M, capd::jaco::strong, (tc.w) * tc.M + 5, capd::jaco::weak);
  range3.setRange(tc.w * tc.M, capd::jaco::weak, -1, capd::jaco::weak);
  range4.setRange(tc.w * tc.M, capd::jaco::weak, tc.w * tc.w * tc.M, capd::jaco::weak);

  const RealPolynomialBound& cr1 = polybd,
                             cr2 = polybd2;

  typedef DPDE2::IndexType IndexType;
  Interval hs = IndexType::harmonicSumK_1<Interval, typename DPDE2::IndexRangeType, capd::jaco::EuclideanNorm<IndexType> >(range3, tc.s);

  Interval sum(0);
  for(DPDE2::IndexType i = dpde2.firstModeIndex(range4); !i.limitReached(range4); i.inc(range4)){
    sum += 1. / power(DPDE2::VNormType::norm(i), tc.s);
  }
  std::cout << "hs=" << hs << ", sum=" << sum << "\n";
  assert(sum <= hs);

  int counter = 0;

  for(DPDE2::IndexType i = dpde2.firstModeIndex(range1); !i.limitReached(range1); i.inc(range1)){
    range2.k = i;
    if(range2.k.upperHalfspace()){
      for(DPDE2::IndexType k_1 = dpde2.firstModeIndex(range2); !k_1.limitReached(range2); k_1.inc(range2)){
        assert((range2.k - k_1).maxNorm() > tc.w * tc.M || k_1.maxNorm() > tc.w * tc.M);
        r3[range2.k] += cr1[range2.k - k_1] * cr1[k_1];
        r4[range2.k] += cr1[range2.k - k_1] * cr2[k_1];
        counter++;
      }
    }
  }
  r3.projectOntoSubspace();
  r4.projectOntoSubspace();


  testData.logfile("lastTestData.txt", true);
  testData << "boundsTest\n";
  testData << "polybd:\n" << polybd << "\n" << "r1:\n" << r1 << "\nr2:\n" << r2 << "\nr3:\n" << r3 << "\nr4:\n" << r4 << "\n";

  testDebug << "summed " << counter << " values in this case\n";
  std::cout << "summed " << counter << " values in this case\n";

  assert(r3.subset(r1));
  testDebug << "assert(r3.subset(r1)) passed succesfully\n";

  assert(r4.subset(r2));
  testDebug << "assert(r4.subset(r2)) passed succesfully\n";

}

/**
 * TODO: this test is not passing - check it
 */
void FFTIsReasonableTest(TestCase tc, int subspace){

  RealPolynomialBound polybd = constructRandomPolyBd(tc, subspace),
                      polybd2 = constructRandomPolyBd(tc, subspace),
                      rFFT = constructPolyBdForResult(tc, subspace),
                      rDirect = constructPolyBdForResult(tc, subspace),
                      rPS = constructPolyBdForResult(tc, subspace),
                      c = constructPolyBdForResult(tc, subspace);

  FFT2D fft = constructFft(tc);

  DFT2DGridType grid = constructDFT(tc),
                s = constructDFT(tc),
                sPS = constructDFT(tc),
                grid2 = constructDFT(tc),
                gridPS = constructDFT(tc),
                grid2PS = constructDFT(tc);
  DPDE2 dpde2 = constructDPDE(tc);

  DPDE2::IndexRangeType range1, range2;
  range1 = polybd.irProjectionPlusFiniteTail;
  range2 = polybd.irFull;

  setC(polybd, 0.); //this can be tested only for 0, because otherwise the result is only true for odd or even cases

  const RealPolynomialBound& cr1 = polybd,
                             cr2 = polybd2;


  for(DPDE2::IndexType i = dpde2.firstModeIndex(range1); !i.limitReached(range1); i.inc(range1)){
    range2.k = i;
    if(range2.k.upperHalfspace()){

      for(DPDE2::IndexType k_1 = dpde2.firstWithinRange(range2); !k_1.limitReached(range2); k_1.inc(range2)){
        rDirect[range2.k] += cr1[range2.k - k_1] * cr1[k_1];
      }
    }
  }


  dpde2.calculateConvolution(polybd, polybd, rFFT);

  testData.logfile("lastTestData.txt", true);
  testData << "FFTIsReasonableTest\n";
  testData << "polybd:\n" << polybd << "\nrDirect:\n" << rDirect << "\nrFFT:\n" << rFFT << "\n";

  assert(rDirect.subset(rFFT));
  testDebug << "assert(rDirect.subset(rFFT)) passed succesfully\n";

}



int main(){
  setLoggers();
  testDebug.logfile("PolyBd2DTests_output.txt", true);
  testDebug.log = true;
  testData.log = true;

  setTestCases();
  time_t rawtime;
  time ( &rawtime );
  testDebug << "The current local time is: " << ctime (&rawtime) << "\n";
  testDebug << "This file contains output from the tests performed by PolyBd2DTests program.\n";
  /*for (int var = 0; var < numberOfTestCases; ++var) {
    printTC(testDebug, testCases[var]);
    FFTIsReasonableTest(testCases[var], even);
    //boundsTest(testCases[var], even);
    optimizedFFTIncludedInFFTTest(testCases[var], even);
    ConvolutionWithTheSameFunctionTest(testCases[var], even);
  }*/

  int m = 5,
      M = 49,
      dft_m = 16,
      dft_M = 150,
      w = 2;

  TestCase tc(m, M, w, dft_m, dft_M, 6, 100, 1e-10 * Interval(-1., 1.));
  boundsTest(tc, even);

}
