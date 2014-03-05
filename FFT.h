/*
 * FFT.h
 *
 *  Created on: Oct 12, 2011
 *      Author: cyranka
 */

#ifndef FFT_H_
#define FFT_H_

#include "capd/vectalg/Vector.hpp"
#include "capd/vectalg/Matrix.hpp"
#include "Index.h"
#include "DFTGrid.h"
#include "ComplexPolyBdJetOptimized.h"

namespace capd{
namespace jaco{

enum PermutationStep { withoutUnscrambling, withUnscrambling};
enum SmallTransformsOmegas { sin60, sin90, sin270, sin72, sin36};

/**
 * In all FFT algorithms implemented here, there is a convention that N is the maximal absolute value of each index
 * coefficient, it means that it is varied from -N/2 upto N/2. Whereas M is the number of discrete points in the interval
 * [0, 2\pi], that is the DFT is calculated at 2\pij/M, where j=0, ..., M-1.
 *
 * The aliasing error can be avoided in various ways. Either PhaseShift (by half or quarter grid cell) or Padding techniques
 * can be used. Refer the paper for details.
 */

enum FFTVariant{unscrambled, scrambled, selfSorting};
enum AliasingRemoval{none, phaseShiftOddEven, phaseShiftRegular, padding};
enum Shift{quarterGridCell, halfGridCell};

enum TransformType{realValuedTransform, complexValuedTransform};//if transform is real-valued, or complex-valued
/**Real odd boundary conditions*/

/**Calculates one-dimensional FFT.
 * First template parameter (ScalarT) is a type that is transformed (complex number, Taylor coefficients),
 * second template parameter (ComplexScalarT) is a class representing a complex number (needed in order to construct
 * the complex matrices W and hatW - representing the discrete transform and inverse transform).
 *
 * Above types are provided separately in order to allow calculating FFT of types that are not necessarily ComplexScalars
 * (for example Jets).
 *
 * ModesContainerT is a type of modes container, input is not a vector, but an abstract ModesContainerType. This is due to
 * the fact in a vector modes can be stored in various ways (negatively indexed modes for example, real and complex parts...)
 *
 */
template< class ScalarT, class ComplexScalarT, int N, int M,
    class ModesContainerT = capd::jaco::ComplexPolyBdJetOptimized<typename ComplexScalarT::ScalarType, ScalarT, capd::jaco::Index1D, (N == 0 ? 0 : N + 1)> >
class FFT1D : public ModesContainerT::SubspaceType{
public:
  typedef ComplexScalarT ComplexScalarType;
  typedef ScalarT ScalarType;
  typedef typename ComplexScalarT::ScalarType IntervalType;
  typedef ModesContainerT ModesContainerType;
  typedef capd::vectalg::Vector<ComplexScalarT, M> ComplexVectorType;
  typedef capd::vectalg::Matrix<ComplexScalarT, 2, 2> QuadraticMatrixType;
  typedef capd::jaco::DFT1DGrid<IntervalType, ScalarType, M> DFT1DGrid;
  typedef DFT1DGrid DFTGridType;
  typedef typename DFT1DGrid::VectorType VectorType;
  typedef typename ModesContainerType::IndexType IndexType;
  typedef typename ModesContainerType::SubspaceType SubspaceType;
  typedef typename SubspaceType::IndexRangeType IndexRangeType;

  int n;
  int m;
  typename capd::vectalg::Vector<int, M> factors; ///< factors of the number N
  typename capd::vectalg::Vector<int, M> factorsInv; ///< inversed order factors
  int p; ///< number of factors
  QuadraticMatrixType W2;
  ComplexVectorType omegas; ///array with complex numbers exp(i * j * k)
  ComplexVectorType omegasInv; ///array with complex numbers exp(-i * j * k)
  ComplexVectorType omegasExt; ///exp(i * j * k) and exp(-i * j * k) together in one array (positive and negative indices)
  ComplexVectorType smallTransformsOmegas;
  ComplexVectorType smallTransformsOmegasInv;
  ComplexVectorType phaseShiftHalfGridCell;
  ComplexVectorType phaseShiftQuarterGridCell;
  ComplexVectorType phaseShiftHalfGridCellInv;
  ComplexVectorType phaseShiftQuarterGridCellInv;
  ComplexScalarType sqrt5D4;
  VectorType temp; ///< temporary results
  VectorType tempR; ///< temporary results
  VectorType pad; ///< temporary results
  DFTGridType s; ///< for temporary results
  typename capd::vectalg::Vector<ScalarType, 15> auxiliary;

  IntervalType pi;
  bool fTmp;

  static bool isAPowerOfTwo(int a, int& p){
    int d;
    p = 0;
    if(a == 1)
      return false;
    while((a / 2) > 0){
      d = (a % 2) && (a % 3) && (a % 5);
      if(d != 0 )
        return false;
      if(a % 5 == 0){
        a /= 5;
        p ++;
        continue;
      }
      if(a % 4 == 0){
        a /= 4;
        p ++;
        continue;
      }
      if(a % 3 == 0){
        a /= 3;
        p ++;
        continue;
      }
      if(a % 2 == 0){
        a /= 2;
        p ++;
        continue;
      }

    }
    return true;
  }

  void fulfillFactors(int a, int p){
    int i;
    for(i=0; i < p; ++i){
      if(a % 5 == 0){
        factors[i + 1] = 5;
        factorsInv[p - i] = 5;
        a /= 5;
        continue;
      }
      if(a % 4 == 0){
        factors[i + 1] = 4;
        factorsInv[p - i] = 4;
        a /= 4;
        continue;
      }
      if(a % 3 == 0){
        factors[i + 1] = 3;
        factorsInv[p - i] = 3;
        a /= 3;
        continue;
      }
      if(a % 2 == 0){
        factors[i + 1] = 2;
        factorsInv[p - i] = 2;
        a /= 2;
        continue;
      }

    }
  }

  FFT1D(){}

  FFT1D(int n_, int m_, IntervalType pi_) : SubspaceType(n_, 2*n_), n(n_), m(m_), factors(m), factorsInv(m), W2(2, 2), omegas(m),
      omegasInv(m), smallTransformsOmegas(10), smallTransformsOmegasInv(10), phaseShiftHalfGridCell(m), phaseShiftQuarterGridCell(m),
      phaseShiftHalfGridCellInv(m), phaseShiftQuarterGridCellInv(m), temp(m), tempR(m), pad(m), s(m),
      auxiliary(15), pi(pi_){
    if(!isAPowerOfTwo(m, p)){
      std::cerr << "Number of points for the discrete transform is not a power of two.\n";
      throw std::runtime_error("Number of points for the discrete transform is not a power of two.\n");
    }
    //important  warning about the aliasing error possibility

    if(! (n + 1 <= m/2)){
      std::cerr << "ATTENTION! ALIASING ERROR WILL BE PRESENT! INCREASE m\n n="<<n<<", m="<<m<<"\n";
    }
    if(n + 1 <= m/2 && 3*n >= m){
      std::cerr << "Warning! Possible aliasing error, PHASE SHIFT ALIASING REMOVAL TECHNIQUE HAS TO BE USED\n n="<<n<<", m="<<m<<"\n";
    }

    W2[0][0] = 1.;
    W2[0][1] = 1.;
    W2[1][0] = 1.;
    W2[1][1] = -1.;

    //fulfill the omega array
    int k, z; IntervalType r;
    fulfillFactors(m, p);
    sqrt5D4 = sqrt(IntervalType(5)) / 4;
    smallTransformsOmegas[sin60] = sin(pi/3);
    smallTransformsOmegasInv[sin60] = -sin(pi/3);

    smallTransformsOmegas[sin90] = 1.;
    smallTransformsOmegasInv[sin90] = -1.;

    smallTransformsOmegas[sin270] = -1.;
    smallTransformsOmegasInv[sin270] = 1.;

    smallTransformsOmegas[sin72] = sin(2*pi/5);
    smallTransformsOmegasInv[sin72] = -sin(2*pi/5);

    smallTransformsOmegas[sin36] = sin(pi/5);
    smallTransformsOmegasInv[sin36] = -sin(pi/5);

    int i;
    IntervalType Delta = pi / m;
    phaseShiftHalfGridCell[0] = 1;
    phaseShiftHalfGridCellInv[0] = 1;
    for(i = 1; i <= n; ++i){
      phaseShiftHalfGridCell[i] = ComplexScalarType(cos(i * Delta), sin(i * Delta));
      phaseShiftHalfGridCellInv[i] = conjugate(phaseShiftHalfGridCell[i]);
    }
    for(i = -n; i <= -1; ++i){
      phaseShiftHalfGridCell[m + i] = ComplexScalarType(cos(i * Delta), sin(i * Delta));
      phaseShiftHalfGridCellInv[m + i] = conjugate(phaseShiftHalfGridCell[m + i]);
    }
    Delta = pi / (2 * m);
    phaseShiftQuarterGridCell[0] = 1;
    phaseShiftQuarterGridCellInv[0] = 1;
    for(i = 1; i <= n; ++i){
      phaseShiftQuarterGridCell[i] = ComplexScalarType(cos(i * Delta), sin(i * Delta));
      phaseShiftQuarterGridCellInv[i] = conjugate(phaseShiftQuarterGridCell[i]);
    }
    for(i = -n; i <= -1; ++i){
      phaseShiftQuarterGridCell[m + i] = ComplexScalarType(cos(i * Delta), sin(i * Delta));
      phaseShiftQuarterGridCellInv[m + i] = conjugate(phaseShiftQuarterGridCell[m + i]);
    }

    for(k=0; k < m; k++){
      r =  pi * (IntervalType(2*k) / IntervalType(m));
      if((12 * k) % m == 0){ ///this means that (2k/m) pi = (z/6) pi, i.e. at least one of the values is rational
        z = (12 * k) / m;
        switch(z){
          case 0 : omegas[k].re = 1; omegas[k].im = 0; break;
          case 1 : omegas[k].re = cos(r); omegas[k].im = 0.5; break;
          case 2 : omegas[k].re = 0.5; omegas[k].im = sin(r); break;
          case 3 : omegas[k].re = 0; omegas[k].im = 1; break;
          case 4 : omegas[k].re = -0.5; omegas[k].im = sin(r); break;
          case 5 : omegas[k].re = cos(r); omegas[k].im = 0.5; break;
          case 6 : omegas[k].re = -1; omegas[k].im = 0; break;
          case 7 : omegas[k].re = cos(r); omegas[k].im = -0.5; break;
          case 8 : omegas[k].re = -0.5; omegas[k].im = sin(r); break;
          case 9 : omegas[k].re = 0; omegas[k].im = -1; break;
          case 10 : omegas[k].re = 0.5; omegas[k].im = sin(r); break;
          case 11 : omegas[k].re = cos(r); omegas[k].im = -0.5; break;
          default : throw std::runtime_error("fatal error in FFT1D constructor.\n");
        }
        omegasInv[k].re = omegas[k].re;
        omegasInv[k].im = -omegas[k].im;
      }else{
        omegas[k].re = cos(r);
        omegas[k].im = sin(r);
        omegasInv[k].re = cos(r);
        omegasInv[k].im = -sin(r);
     }
    }

  }


  /**Returns the smallest number of discrete points in FFT that avoids aliasing error.
   */
  static int numberOfPointsAvoidingAliasing(int n_){
    int r = 3 * n_ + 1, t;
    while(!isAPowerOfTwo(r, t))
      r++;
    return r;
  }

  ///W_2
  inline const void smallTransform2(const ScalarType& first, const ScalarType& second, ScalarType& rFirst, ScalarType& rSecond,
                                    bool calculateFirst = true, bool calculateSecond = true){
    auxiliary[1] = first;
    if(calculateFirst) rFirst = first + second;
    if(calculateSecond) rSecond = auxiliary[1] - second;
  }

  ///W_3
  inline const void smallTransform3(const ScalarType& first, const ScalarType& second, const ScalarType& third,
                                    ScalarType& rFirst, ScalarType& rSecond, ScalarType& rThird,
                                    const ComplexVectorType& smallTransformsOmegas,
                                    bool calculateFirst = true, bool calculateSecond = true, bool calculateThird = true){
    auxiliary[1] = second + third;
    auxiliary[2] = first - 0.5 * auxiliary[1];
    auxiliary[3] = smallTransformsOmegas[sin60] * (second - third);
    if(calculateFirst) rFirst = first + auxiliary[1];
    if(calculateSecond) rSecond = auxiliary[2] + ComplexScalarType::i() * auxiliary[3];
    if(calculateThird) rThird = auxiliary[2] - ComplexScalarType::i() * auxiliary[3];
  }

  ///W_4
  inline const void smallTransform4(ScalarType& first, ScalarType& second, ScalarType& third, ScalarType& fourth,
                                    ScalarType& rFirst, ScalarType& rSecond, ScalarType& rThird, ScalarType& rFourth,
                                    const ComplexVectorType& smallTransformsOmegas,
                                    bool calculateFirst = true, bool calculateSecond = true, bool calculateThird = true,
                                    bool calculateFourth = true){
    auxiliary[1] = first + third;
    auxiliary[2] = second + fourth;
    auxiliary[3] = first - third;
    auxiliary[4] = second - fourth;
    if(calculateFirst) rFirst = auxiliary[1] + auxiliary[2];
    if(calculateSecond) rSecond = auxiliary[3] + ComplexScalarType::i() * smallTransformsOmegas[sin90] * auxiliary[4];
    if(calculateThird) rThird = auxiliary[1] - auxiliary[2];
    if(calculateFourth) rFourth = auxiliary[3] + ComplexScalarType::i() * smallTransformsOmegas[sin270] * auxiliary[4];
  }

  ///W_5
  inline const void smallTransform5(const ScalarType& first, const ScalarType& second, const ScalarType& third, const ScalarType& fourth, const ScalarType& fifth,
                                    ScalarType& rFirst, ScalarType& rSecond, ScalarType& rThird, ScalarType& rFourth, ScalarType& rFifth,
                                    const ComplexVectorType& smallTransformsOmegas,
                                    bool calculateFirst = true, bool calculateSecond = true, bool calculateThird = true,
                                    bool calculateFourth = true, bool calculateFifth = true){
    auxiliary[1] = second + fifth;
    auxiliary[2] = third + fourth;
    auxiliary[3] = second - fifth;
    auxiliary[4] = third - fourth;
    auxiliary[5] = auxiliary[1] + auxiliary[2];
    auxiliary[6] = sqrt5D4 * (auxiliary[1] - auxiliary[2]);
    auxiliary[7] = first - 0.25 * auxiliary[5];
    auxiliary[8] = auxiliary[7] + auxiliary[6];
    auxiliary[9] = auxiliary[7] - auxiliary[6];
    ///ATTENTION: the following operation has been splitted into two lines, because there is a buffer used to add/substract
    ///FirstOrderJets, and this buffer works only up to three operations per line.
    auxiliary[10] = (smallTransformsOmegas[sin72] * auxiliary[3]);
    auxiliary[10] = auxiliary[10] + (smallTransformsOmegas[sin36] * auxiliary[4]);
    ///ATTENTION: the following operation has been splitted into two lines, because there is a buffer used to add/substract
    ///FirstOrderJets, and this buffer works only up to three operations per line.
    auxiliary[11] = (smallTransformsOmegas[sin36] * auxiliary[3]);
    auxiliary[11] = auxiliary[11] - (smallTransformsOmegas[sin72] * auxiliary[4]);
    if(calculateFirst) rFirst = first + auxiliary[5];
    if(calculateSecond) rSecond = auxiliary[8] + ComplexScalarType::i() * auxiliary[10] ;
    if(calculateThird) rThird = auxiliary[9] + ComplexScalarType::i() * auxiliary[11];
    if(calculateFourth) rFourth = auxiliary[9] - ComplexScalarType::i() * auxiliary[11];
    if(calculateFifth) rFifth = auxiliary[8] - ComplexScalarType::i() * auxiliary[10];
  }

  ///``Self-sorting'' variant
  inline const void FFTSS(const VectorType& u, const ComplexVectorType& omegas, const ComplexVectorType& smallTransformsOmegas,
                          typename capd::vectalg::Vector<int, M>& factors, VectorType& temp, VectorType& r){
    int j = 0,
        n_j, l_j, m_j, l, k,
        firstIndex, secondIndex, thirdIndex, fourthIndex, fifthIndex,
        rFirstIndex, rSecondIndex, rThirdIndex, rFourthIndex, rFifthIndex,
        x;
    r = u;
    for(j = 1, l_j = 1; j <= p; ++j){
      n_j = factors[j];
      m_j = m / (l_j*n_j);
      for(k=0; k < m_j; k++){
        for(l=0; l < l_j; l++){
          firstIndex = l + k*l_j;
          //s = k / n_j;
          //t = k % n_j;
          //rFirstIndex = t*m/n_j + s*l_j + l;
          rFirstIndex = n_j*k*l_j + l;
          secondIndex = l + k*l_j + m/n_j;
          //s = (k + m_j)/ n_j;
          //t = (k + m_j) % n_j;
          //rSecondIndex = t*m/n_j + s*l_j + l;
          rSecondIndex = n_j*k*l_j + l_j + l;
          if(n_j == 2){
            smallTransform2(r[firstIndex], r[secondIndex], temp[rFirstIndex], temp[rSecondIndex]);
          }else{
            thirdIndex = l + k*l_j + 2*m/n_j;
            //s = (k + 2*m_j)/ n_j;
            //t = (k + 2*m_j) % n_j;
            //rThirdIndex = t*m/n_j + s*l_j + l;
            rThirdIndex = n_j*k*l_j + 2*l_j + l;
            if(n_j == 3){
              smallTransform3(r[firstIndex], r[secondIndex], r[thirdIndex],
                              temp[rFirstIndex], temp[rSecondIndex], temp[rThirdIndex], smallTransformsOmegas);
            }else{
              fourthIndex = l + k*l_j + 3*m/n_j;
              //s = (k + 3*m_j)/ n_j;
              //t = (k + 3*m_j) % n_j;
              //rFourthIndex = t*m/n_j + s*l_j + l;
              rFourthIndex = n_j*k*l_j + 3*l_j + l;
              if(n_j == 4){
                smallTransform4(r[firstIndex], r[secondIndex], r[thirdIndex], r[fourthIndex],
                                temp[rFirstIndex], temp[rSecondIndex], temp[rThirdIndex], temp[rFourthIndex], smallTransformsOmegas);
              }
              if(n_j == 5){
                fifthIndex = l + k*l_j + 4*m/n_j;
                //s = (k + 4*m_j)/ n_j;
                //t = (k + 4*m_j) % n_j;
                //rFourthIndex = t*m/n_j + s*l_j + l;
                rFifthIndex = n_j*k*l_j + 4*l_j + l;
                smallTransform5(r[firstIndex], r[secondIndex], r[thirdIndex], r[fourthIndex], r[fifthIndex],
                                temp[rFirstIndex], temp[rSecondIndex], temp[rThirdIndex], temp[rFourthIndex], temp[rFifthIndex],
                                smallTransformsOmegas);
                temp[rFifthIndex] = omegas[4 * k * l_j] * temp[rFifthIndex];
              }
              temp[rFourthIndex] = omegas[3 * k * l_j] * temp[rFourthIndex];
            }
            temp[rThirdIndex] = omegas[2 * k * l_j] * temp[rThirdIndex];
          }
          temp[rSecondIndex] = omegas[k * l_j] * temp[rSecondIndex];
        }
      }
      r = temp;
      l_j *= n_j;
      fftDebug << "j="<<j<<", factors[j]="<<factors[j]<<" l_j+1="<<l_j<<", m_j="<<m_j<<", r:\n";
      for(x=0; x < r.size(); ++x)
        fftDebug<<r[x]<<"\n";
      fftDebug<<"\n";
    }
  }

  /** The transform procedure body, without normalizing.
   *
   * Note: when withoutUnscrambling variant of FFT is used (the output data are not unscrambled - the order in which modes appear on
   * the output is not at all relevant) the permutation step (multiplication by permutation matrix - the last 'if') can be skipped.
   * 
   * In case the permutation step is skipped, the inverse transform must be executed without performing the permutation step at all,
   * IMPORTANT: FFTv2 variant has to be used in order for calculations to match.
   *
   * @param u vector that is being transformed
   * @param omegas vector consisting of values of exp(i * j * k), depends if a transform or inverse transform is being performed
   * @param temp temporary buffer vector, used when permuting elements of the vector u
   * @param r a reference, when the result is returned
   * @return
   */
  inline const void FFT(const VectorType& u, const ComplexVectorType& omegas, const ComplexVectorType& smallTransformsOmegas,
                        typename capd::vectalg::Vector<int, M>& factors, VectorType& temp, VectorType& r, int unscrambl = withoutUnscrambling,
                        int transformType = realValuedTransform){
    int j = 0,
        n_j = 2,
        l_j, m_j, l, k, firstIndex, secondIndex, thirdIndex, fourthIndex, fifthIndex, s, t, x;
    r = u;
    ///scrambled data are calculated here
    for(j=1, l_j = 1; j <= p; ++j){
      n_j = factors[j];
      m_j = m / (l_j*n_j);
      for(l=0; l < l_j; ++l){
        for(k=0; k <= (unscrambl == withoutUnscrambling && transformType == realValuedTransform ? m_j/2 : m_j-1); ++k){ //this is an optimization, but works ONLY for real valued series
        ///TODO: check why we cannot optimize this way when unscrambling is used (it looks like the theorem is not fulfilled)
//      for(k=0; k < m_j; ++k){
          firstIndex = k + l*m_j*n_j;
          secondIndex = k + m_j + l*m_j*n_j;
          if(n_j == 2 ){
            smallTransform2(r[firstIndex], r[secondIndex], r[firstIndex], r[secondIndex]);
            r[secondIndex] = omegas[l_j * k] * r[secondIndex];
          }
          if(n_j == 3){
            thirdIndex = k + 2*m_j + l*m_j*n_j;
            smallTransform3(r[firstIndex], r[secondIndex], r[thirdIndex], r[firstIndex], r[secondIndex], r[thirdIndex], smallTransformsOmegas);
            r[secondIndex] = omegas[l_j * k] * r[secondIndex];
            r[thirdIndex] = omegas[2 * l_j * k] * r[thirdIndex];
          }
          if(n_j == 4){
            thirdIndex = k + 2*m_j + l*m_j*n_j;
            fourthIndex = k + 3*m_j + l*m_j*n_j;
            smallTransform4(r[firstIndex], r[secondIndex], r[thirdIndex], r[fourthIndex],
                            r[firstIndex], r[secondIndex], r[thirdIndex], r[fourthIndex], smallTransformsOmegas);
            r[secondIndex] = omegas[l_j * k] * r[secondIndex];
            r[thirdIndex] = omegas[2 * l_j * k] * r[thirdIndex];
            r[fourthIndex] = omegas[3 * l_j * k] * r[fourthIndex];
          }
          if(n_j == 5){
            thirdIndex = k + 2*m_j + l*m_j*n_j;
            fourthIndex = k + 3*m_j + l*m_j*n_j;
            fifthIndex = k + 4*m_j + l*m_j*n_j;
            smallTransform5(r[firstIndex], r[secondIndex], r[thirdIndex], r[fourthIndex], r[fifthIndex],
                            r[firstIndex], r[secondIndex], r[thirdIndex], r[fourthIndex], r[fifthIndex],
                            smallTransformsOmegas);
            r[secondIndex] = omegas[l_j * k] * r[secondIndex];
            r[thirdIndex] = omegas[2 * l_j * k] * r[thirdIndex];
            r[fourthIndex] = omegas[3 * l_j * k] * r[fourthIndex];
            r[fifthIndex] = omegas[4 * l_j * k] * r[fifthIndex];
          }
        }
        if(unscrambl == withoutUnscrambling && transformType == realValuedTransform){
          
          for(k = m_j/2 + 1; k < m_j; ++k){
            for(s=0; s < n_j; ++s){
              r[k + s*m_j + l*m_j*n_j] = conjugate(r[m_j - k + s*m_j + l*m_j*n_j]);
            }
          }
        }
      }
      fftDebug << "j="<<j<<", factors[j]="<<factors[j]<<" l_j+1="<<l_j<<", m_j="<<m_j<<", r:\n";
      
      if(unscrambl == withoutUnscrambling && transformType == realValuedTransform){
        for(l=0; l < m / m_j; ++l){
          r[l * m_j].setImaginaryPartToZero();
          if(m_j % 2 == 0) r[l * m_j + m_j/2].setImaginaryPartToZero();
        }
      }
      l_j *= n_j;
      for(x=0; x < r.size(); ++x)
        fftDebug<<r[x]<<"\n";
      fftDebug<<"\n";
    }

    ///unscrambling the data
    if(unscrambl == withUnscrambling){
      if(p > 1){
        for(j=p-1, l_j = m/(factors[p]); j >= 1; --j){
          n_j = factors[j];
          l_j /= n_j;
          m_j = m / (l_j*n_j);
          temp = r;
          for(l=0; l < l_j; ++l){
            for(k=1; k < m_j*n_j - 1; ++k){
                t = k / n_j;
                s = k % n_j;
                firstIndex = k + l*m_j*n_j;
                secondIndex = m_j*s + t + l*m_j*n_j;
                if(firstIndex != secondIndex){
                  r[firstIndex] = temp[secondIndex];
                }
            }
          }

        }
      }
      fftDebug << "unscrambled:\n";
      for(x=0; x < r.size(); ++x)
        fftDebug<<r[x]<<", diam="<<diam(r[x].value().re)<<", "<<diam(r[x].value().im)<<"\n";
      fftDebug<<"\n";
    }
  }


  /**This is the variant of the FFT algorithm, which first multiplies by permutation matrices and then by the sparse matrices with
   * ''omegas''.
   * 
   * Note: when withoutUnscrambling variant of FFT is used (the output data are not unscrambled - the order in which modes appear on
   * the output is not at all relevant) the permutation step (multiplication by permutation matrix - the first 'if') can be skipped. 
   */
  inline const void FFTv2(const VectorType& u, const ComplexVectorType& omegas, const ComplexVectorType& smallTransformsOmegas,
                          typename capd::vectalg::Vector<int, M>& factors, VectorType& temp, VectorType& r, int unscrambl = withoutUnscrambling,
                          int transformType = realValuedTransform){
    int j = 0,
        n_j = 2,
        l_j, m_j, l, k, firstIndex, secondIndex, thirdIndex, fourthIndex, fifthIndex, s, t, x;
    r = u;
    
    //it seems that this function is implemented only for realValued transform
    if(transformType != realValuedTransform)
      throw std::runtime_error("The function FFTv2 is not implemented for complex input (this function is optimized for real-valued transform)");
    
    if(unscrambl == withUnscrambling){
      if(p > 1){
        for(j=p, l_j = m; j > 1; --j){ //this is different than in v1, because in this case P_1=Id
          n_j = factors[j];
          l_j /= n_j;
          m_j = m / (l_j*n_j);
          temp = r;
          for(k=0; k < m_j; ++k){
            for(l=1; l < l_j*n_j-1; ++l){
              t = l / n_j;
              s = l % n_j;
              firstIndex  = l         + k*l_j*n_j;
              secondIndex = s*l_j + t + k*l_j*n_j;
              if(firstIndex != secondIndex){
                r[secondIndex] = temp[firstIndex];//this is very important that we use the indices in this way
                                                  //because this refers to multiplication by P_{n_j}^{l_j}
                                                  //and this is different that P_{m_j}^{n_j} (from the first variant)
              }
            }
          }
        }
      }
    }

    for(j=1, l_j=1; j <= p; ++j){
      n_j = factors[j];
      m_j = m / (l_j*n_j);
      for(k=0; k < m_j; ++k){
        for(l=0; l < l_j; ++l){
          firstIndex = l + k*n_j*l_j;
          secondIndex = l + l_j + k*n_j*l_j;
          ///booleans firstIndex <= sth, secondIndex <= sth, etc... are here, because we want to calculate only half, and obtain
          ///the rest by the conjugacy condition.
          ///TODO: works only for real valued solutions
          if(n_j == 2){
            r[secondIndex] = omegas[m_j * l] * r[secondIndex];
            smallTransform2(r[firstIndex], r[secondIndex], r[firstIndex], r[secondIndex],
                            (firstIndex <= n_j*l_j / 2 + k*n_j*l_j), (secondIndex <= n_j*l_j / 2 + k*n_j*l_j));
          }
          if(n_j == 3){
            thirdIndex = l + 2*l_j + k*n_j*l_j;
            r[secondIndex] = omegas[m_j * l] * r[secondIndex];
            r[thirdIndex] = omegas[2 * m_j * l] * r[thirdIndex];
            smallTransform3(r[firstIndex], r[secondIndex], r[thirdIndex], r[firstIndex], r[secondIndex], r[thirdIndex], smallTransformsOmegasInv,
                            (firstIndex <= n_j*l_j / 2 + k*n_j*l_j), (secondIndex <= n_j*l_j / 2 + k*n_j*l_j),
                            (thirdIndex <= n_j*l_j / 2 + k*n_j*l_j));
          }
          if(n_j == 4){
            thirdIndex = l + 2*l_j+k*n_j*l_j;
            fourthIndex = l + 3*l_j+k*n_j*l_j;
            r[secondIndex] = omegas[m_j * l] * r[secondIndex];
            r[thirdIndex] = omegas[2 * m_j * l] * r[thirdIndex];
            r[fourthIndex] = omegas[3 * m_j * l] * r[fourthIndex];
            smallTransform4(r[firstIndex], r[secondIndex], r[thirdIndex], r[fourthIndex],
                            r[firstIndex], r[secondIndex], r[thirdIndex], r[fourthIndex], smallTransformsOmegasInv,
                            (firstIndex <= n_j*l_j / 2 + k*n_j*l_j), (secondIndex <= n_j*l_j / 2 + k*n_j*l_j),
                            (thirdIndex <= n_j*l_j / 2 + k*n_j*l_j), (fourthIndex <= n_j*l_j / 2 + k*n_j*l_j));
          }
          if(n_j == 5){
            thirdIndex = l + 2*l_j+k*n_j*l_j;
            fourthIndex = l + 3*l_j+k*n_j*l_j;
            fifthIndex = l + 4*l_j+k*n_j*l_j;
            r[secondIndex] = omegas[m_j * l] * r[secondIndex];
            r[thirdIndex] = omegas[2 * m_j * l] * r[thirdIndex];
            r[fourthIndex] = omegas[3 * m_j * l] * r[fourthIndex];
            r[fifthIndex] = omegas[4 * m_j * l] * r[fifthIndex];
            smallTransform5(r[firstIndex], r[secondIndex], r[thirdIndex], r[fourthIndex], r[fifthIndex],
                            r[firstIndex], r[secondIndex], r[thirdIndex], r[fourthIndex], r[fifthIndex], smallTransformsOmegasInv,
                            (firstIndex <= n_j*l_j / 2 + k*n_j*l_j), (secondIndex <= n_j*l_j / 2 + k*n_j*l_j),
                            (thirdIndex <= n_j*l_j / 2 + k*n_j*l_j), (fourthIndex <= n_j*l_j / 2 + k*n_j*l_j),
                            (fifthIndex <= n_j*l_j / 2 + k*n_j*l_j));
          }
        }
        //here the rest of values are calculated by conjugating the already calcualated values (works only for real-valued transform)
        for(l = 1; l < (n_j*l_j + 1)/ 2; l++){
          r[l + n_j*l_j / 2 + k*n_j*l_j] = conjugate(r[(n_j*l_j + 1)/ 2 + k*n_j*l_j - l]);
        }

      }
      fftDebug << "j="<<j<<", factors[j]="<<factors[j]<<", l_j+1="<<l_j<<", m_j="<<m_j<<", r:\n";
      l_j *= n_j;
      for(x=0; x < r.size(); ++x)
        fftDebug<<r[x]<<"\n";
      fftDebug<<"\n";
    }
  }


  /** Calculates the result of the matrix-vector multiplication
   * /f[
   *   \hat{u}=W_M \cdot u
   * /f]
   * /f[
   *   W_M(k,l)=e^{ikl2\pi/M}
   * /f]
   */
  inline void transform(const VectorType& u, VectorType& r, int variant, int transformType = realValuedTransform){
    switch(variant){
      case unscrambled : FFT(u, omegas, smallTransformsOmegas, factors, temp, r, withUnscrambling, transformType); break;
      case scrambled : FFT(u, omegas, smallTransformsOmegas, factors, temp, r, withoutUnscrambling, transformType); break;
      case selfSorting : FFTSS(u, omegas, smallTransformsOmegas, factors, temp, r); break;
      default: FFTSS(u, omegas, smallTransformsOmegas, factors, temp, r);
    }
  }

  /** Calculates the result of the matrix-vector multiplication (the inverse transform)
   * /f[
   *   u=\hat{W}_M \cdot \hat{u}
   * /f]
   * /f[
   *   \hat{W}_M(k,l)=e^{ikl2\pi/M}
   * /f]
   */
  inline void inverseTransform(const VectorType& u, VectorType& r, int variant, int transformType = realValuedTransform){
    switch(variant){
      case unscrambled : FFT(u, omegasInv, smallTransformsOmegasInv, factorsInv, temp, r, withUnscrambling, transformType); break;
      case scrambled : FFTv2(u, omegasInv, smallTransformsOmegasInv, factorsInv, temp, r, withoutUnscrambling, transformType); break;
      case selfSorting : FFTSS(u, omegasInv, smallTransformsOmegasInv, factorsInv, temp, r); break;
      default: FFTSS(u, omegasInv, smallTransformsOmegasInv, factorsInv, temp, r);
    }
  }


  /** Returns v = [u_0, u_1, ..., u_{N/2-1}, 0, ..., 0] or v = [0, u_0, u_1, ..., u_{N/2-1}, 0, ..., 0]. Translates
   * a ModesContainerVector to a Vector, which is then being transformed.
   *
   * @param u is the vector which is being padded
   * @param v is the padded u
   * @param positiveIndexed if u are the components of v indexed by the positive indices (true) or the negative indices (false)
   * @param zeroMode if the zero indexed mode is taken into account
   */
  inline void padVector(const ModesContainerType& u, VectorType& v, bool upperHalfspace, bool zeroMode){
    int i;
    if(upperHalfspace){
      for(i = 0; i <= n; ++i){
        v[i] = u[Index1D(i)];
      }
      for(i = -n; i <= -1; ++i){
        v[m + i] = u[Index1D(i)];
      }
    }else{
      for(i=1; i <= n; ++i)
        v[i] = u[Index1D(-i)];
      if(zeroMode)
        v[0] = u[Index1D(0)];
      else
        v[0] = ComplexScalarType(0); ///<this line is important, otherwise the zero-indexed element may be added twice
    }
  }

  /**Performs a phase shift of the input vector of size m;
   */
  inline void phaseShiftVector(VectorType& v, int phaseShift, int direction){
    int i;
    for(i = 0; i < m; ++i){
      if(phaseShift == phaseShiftOddEven){
        if(direction == 1)
          v[i] *= phaseShiftQuarterGridCell[i];
        else
          v[i] *= phaseShiftQuarterGridCellInv[i];
      }
      if(phaseShift == phaseShiftRegular){
        if(direction == 1)
          v[i] *= phaseShiftHalfGridCell[i];
        else
          v[i] *= phaseShiftHalfGridCellInv[i];
      }
    }
  }

  
  
  /**THIS FUNCTION IS A COMPLEX VERSION OF FFT, THIS IS IMPLEMENTED FOR 2D TRANSFORM.
   * 
   * default variant here is the self-sorting (because this transform is called for complex data and optimizations cannot be performed)
   */
  inline void extendedTransform(const ModesContainerType& modes, DFTGridType& r, int aliasingRemoval = padding, int variant = selfSorting){
    
    padVector(modes, pad, 1, 1);
     if(aliasingRemoval == padding){
       if(m <= 3 * n){
         std::cerr << "Relation m <= 3 * n is not satisfied, m=" << m << ", 3*n=" << 3*n << "\n";
         throw std::runtime_error("Relation m <= 3 * n is not satisfied\n");
       }
     }else{
       if(aliasingRemoval == phaseShiftOddEven)
         phaseShiftVector(pad, phaseShiftOddEven, 1);
       if(aliasingRemoval == phaseShiftRegular)
         phaseShiftVector(pad, phaseShiftRegular, 1);
     }
     fftDebug << "input:\n";
     int x;
     for(x=0; x < pad.size(); ++x)
       fftDebug<<pad[x]<<", diam="<<diam(pad[x].value().re)<<", "<<diam(pad[x].value().im)<<"\n";
     transform(pad, r, variant, complexValuedTransform);
     r.setSubspaceType(modes);
     r.projectOntoSubspace();
    
  }

  /**assumes that the input modes 'modes' represents a real-valued function
   * 
   */
  inline void fastTransform(const ModesContainerType& modes, DFTGridType& r, int aliasingRemoval = phaseShiftOddEven, int variant = scrambled){
//    padVector(modes, pad, 1, 0);
//    transform(pad, tempR);
//    int i;
//    for(i=0; i < m; ++i){
//      r[i] = 2 * tempR[i] + modes[Index1D(0)];
//    }
//    r.setImaginaryPartToZero();
    padVector(modes, pad, 1, 1);
    if(aliasingRemoval == padding){
      if(m <= 3 * n){
        std::cerr << "Relation m <= 3 * n is not satisfied, m=" << m << ", 3*n=" << 3*n << "\n";
        throw std::runtime_error("Relation m <= 3 * n is not satisfied\n");
      }
    }else{
      if(aliasingRemoval == phaseShiftOddEven)
        phaseShiftVector(pad, phaseShiftOddEven, 1);
      if(aliasingRemoval == phaseShiftRegular)
        phaseShiftVector(pad, phaseShiftRegular, 1);
    }
    fftDebug << "input:\n";
    int x;
    for(x=0; x < pad.size(); ++x)
      fftDebug<<pad[x]<<", diam="<<diam(pad[x].value().re)<<", "<<diam(pad[x].value().im)<<"\n";
    transform(pad, r, variant);
    r.setSubspaceType(modes);
    r.projectOntoSubspace();
  }

  /**used in 2D FFT transform. 'extended' means that the optimizations which can be performed in real-valued case are not performed here
   */
  inline void inverseExtendedTransform(const DFTGridType& s, ModesContainerType& r, int aliasingRemoval = padding, int variant = selfSorting){
    if(!ModesContainerType::storesLowerHalfspaceIndependently){
      std::cerr << " Error in FFT1D inverseExtendedTransform function. It is allowed only for containers that are storing the lower " <<
          "halfspace of modes independently (they are not dependent by the conjugating relation a_k=\\overline{a_{-k}}\n";
      throw std::runtime_error(" Error in FFT1D inverseExtendedTransform function. It is allowed only for containers that are storing the lower halfspace of modes independently (they are not dependent by the conjugating relation a_k=\\overline{a_{-k}}\n");
    }
    fftDebug << "input:\n";
    int x;
    for(x=0; x < s.size(); ++x)
      fftDebug<<s[x]<<", diam="<<diam(s[x].value().re)<<", "<<diam(s[x].value().im)<<"\n";
    
    inverseTransform(s, tempR, variant, complexValuedTransform);
    if(aliasingRemoval == padding){
      if(m <= 3 * n){
        std::cerr << "Relation m <= 3 * n is not satisfied, m=" << m << ", 3*n=" << 3*n << "\n";
        throw std::runtime_error("Relation m <= 3 * n is not satisfied\n");
      }
    }else{
      if(aliasingRemoval == phaseShiftOddEven)
        phaseShiftVector(tempR, phaseShiftOddEven, -1);
      if(aliasingRemoval == phaseShiftRegular)
        phaseShiftVector(tempR, phaseShiftRegular, -1);
    }
    int i;
    for(i = 0; i <= n; ++i){
      r[Index1D(i)] = tempR[i];
    }
    if(ModesContainerType::storesLowerHalfspaceIndependently){
      for(i = -n; i < 0; ++i){
        r[Index1D(i)] = tempR[m + i];
      }
    }
    r.setSubspaceType(s);
    r.projectOntoSubspace();    
  }
  
  inline void fastInverseTransform(const DFTGridType& s, ModesContainerType& r, int aliasingRemoval = phaseShiftOddEven, int variant = scrambled){
//    inverseTransform((const VectorType&)s, tempR);
//    int i;
//    for(i=0; i <= n; ++i){
//      r[IndexType(i)] = tempR[i];
//      if(i > 0 && ModesContainerType::storesLowerHalfspaceIndependently){
//        r[IndexType(-i)] = conjugate(tempR[i]);
//      }
//    }
//    VectorType& v = (VectorType&)s;
    fftDebug << "input:\n";
    int x;
    for(x=0; x < s.size(); ++x)
      fftDebug<<s[x]<<", diam="<<diam(s[x].value().re)<<", "<<diam(s[x].value().im)<<"\n";
    
    inverseTransform(s, tempR, variant);
    if(aliasingRemoval == padding){
      if(m <= 3 * n){
        std::cerr << "Relation m <= 3 * n is not satisfied, m=" << m << ", 3*n=" << 3*n << "\n";
        throw std::runtime_error("Relation m <= 3 * n is not satisfied\n");
      }
    }else{
      if(aliasingRemoval == phaseShiftOddEven)
        phaseShiftVector(tempR, phaseShiftOddEven, -1);
      if(aliasingRemoval == phaseShiftRegular)
        phaseShiftVector(tempR, phaseShiftRegular, -1);
    }
    int i;
    for(i = 0; i <= n; ++i){
      r[Index1D(i)] = tempR[i];
    }
    if(ModesContainerType::storesLowerHalfspaceIndependently){
      for(i = -n; i < 0; ++i){
        r[Index1D(i)] = tempR[m + i];
      }
    }
    r.setSubspaceType(s);
    r.projectOntoSubspace();

//    for(i = -n; i <= -1; ++i){
//      intersection(tempR[-i].value().re, tempR[m + i].value().re, tempR[m + i].value().re);
//      tempR[-i].value().re = tempR[m + i].value().re;
//      r[Index1D(i)] = tempR[m + i];
//    }
//    generalDebug << "tempR:\n";
//    for(j = 0; j < tempR.size(); ++j)
//      generalDebug << tempR[j] <<", diam="<<diam(tempR[j].value().re)<< "\n";
  }

  /**
   * for test purpose only. Calculates naively a convolution (using O(N^2) operations).
   */
  const ModesContainerType calculateConvolution(const ModesContainerType& u, const ModesContainerType& v) const {
    int k, k_1;
    ScalarType first, second, r;
    //ModesContainerType s(n, n, true); //this was before, was working with Polynomial Bounds?
    ModesContainerType s(n);
    s.multiply(u, v);
    for(k=0; k <= n; k++){
        s[IndexType(k)] = ComplexScalarType(0);
        for(k_1 = -n; k_1 <= n; k_1++){
          if(abs(k-k_1) <= n && abs(k_1) <= n){

              first = (k_1==0 ? ScalarType(0) : u[IndexType(k_1)]);
              second = (k-k_1==0 ? ScalarType(0) : v[IndexType(k-k_1)]);
              r = first * second;
              s[IndexType(k)] += r;
          }
        }
        s[IndexType(-k)] = conjugate(s[IndexType(k)]);
    }
    return s;
  }

};


/*
 * 17.04.2013 tested to work good
 */
template< class ScalarT, class ComplexScalarT, int N, int M,
  class ModesContainerT = capd::jaco::ComplexPolyBdJetOptimized<typename ComplexScalarT::ScalarType, ScalarT, capd::jaco::Index2D, 2*2*2*N*N> >
class FFT2DOneComponent : public ModesContainerT::SubspaceType{
public:
  typedef ScalarT ScalarType;
  typedef ComplexScalarT ComplexScalarType;
  typedef typename ComplexScalarType::ScalarType IntervalType;
  typedef FFT1D<ScalarType, ComplexScalarType, N, M> FFT1DType;
  typedef typename FFT1DType::VectorType Vector1DType;
  typedef ModesContainerT ModesContainer2DType;
  typedef typename FFT1DType::ModesContainerType ModesContainer1DType;
  typedef ModesContainer2DType ModesContainerType;
  typedef capd::jaco::DFT2DGrid<IntervalType, ScalarType, M> DFT2DGridType;
  typedef DFT2DGridType DFTGridType;
  typedef typename DFT2DGridType::DFT1DGridType DFT1DGridType;
  typedef capd::jaco::Index2D Index2DType;
  typedef Index2DType IndexType;
  typedef typename FFT1DType::IndexType Index1DType;

  typedef capd::jaco::DFT1DGrid<IntervalType, ModesContainer1DType, M> Grid1DOfModesContainer1DType; ///<is that not too complicated?
  typedef capd::jaco::ComplexPolyBdJetOptimized<IntervalType, DFT1DGridType, Index1DType, 2*N > ModesContainer1DOfGrid1DType;
  typedef typename ModesContainerType::SubspaceType SubspaceType;
  typedef typename SubspaceType::IndexRangeType IndexRangeType;


  int n;
  int m;
  FFT1DType fft1d;
  DFT1DGridType rProjection;
  ModesContainer1DType mc1d;
  Grid1DOfModesContainer1DType gridOfModes;
  ModesContainer1DOfGrid1DType modesContainerOfGrid;
  DFT2DGridType s; ///< for temporary results
  ModesContainer2DType t; ///< for temporary results

  FFT2DOneComponent(){}
  
  FFT2DOneComponent(int n_, int m_, IntervalType pi_) : SubspaceType(n_, 2*n_), n(n_), m(m_), fft1d(n_, m_, pi_), rProjection(m_), mc1d(n_),
      gridOfModes(m_, false), modesContainerOfGrid(n_), s(m_), t(n_){
    int i;
    for(i=0; i < m; ++i)
      gridOfModes[i] = ModesContainer1DType(n);    
    for(i=-n; i <= n; ++i){     
      modesContainerOfGrid[Index1DType(i)] = DFT1DGridType(m);      
    }
    
  }

  /**Takes 1D projection of 2D discrete set of modes, and returns 1D ModesContainer.
   *
   * @param component this is the component of modes for which the transform is being calculated (either 0 or 1)
   */
  inline void takeProjection(const ModesContainer2DType& mc2d, ModesContainer1DType& mc1d, int fixedSecondComponent, int component) const{
    int k_1;
    Index2DType index;
    index.l = component;
    index[1] = fixedSecondComponent;
    
    for(k_1 = -n; k_1 <= n; ++k_1){
      index[0] = k_1;
      mc1d[Index1D(k_1)] = mc2d[index];
    }
    
  }

  /**Sets a fixed 1D projection of 2D discrete set of modes with provided array of values.
   * @param direction this is the component of modes for which the transform is being calculated (either 0 or 1)
   */
  inline void setProjection(const ModesContainer1DType& mc1d, ModesContainer2DType& mc2d, int fixedSecondComponent) const{
    int k_1;
    Index2DType index;
    index[1] = fixedSecondComponent;
    for(k_1 = -n; k_1 <= n; ++k_1){
      index[0] = k_1;
      index.l = 0;
      mc2d.set(index, mc1d[Index1D(k_1)]);
      index.l = 1;
      mc2d.set(index, mc1d[Index1D(k_1)]);
    }
  }


  /**calculates the scalar product
   * /f[
   *  \sum_{k_1\in\text{a Projection}}{a_{k_1}\cdot b_{k-k_1}}
   * /f]
   */
  inline void scalarProduct(const DFT2DGridType& a, const DFT2DGridType& b, ModesContainer2DType& r){
    s.multiply(a, b);

    IndexType index;
    IndexRangeType ir;
    ir.setRange(0, capd::jaco::strong, n, capd::jaco::weak);
    
    inverseExtendedTransform(s, r);

  }

  /**Write what are differences with extendedTransform. (it calculates only one coordinate)
   *
   */
  inline void transform(const ModesContainer2DType& modes, DFT2DGridType& r, int aliasingRemoval = padding){
    int k_2, j_1;

    ///first, calculate N independent FFTs, one for each fixed second index component
    for(k_2 = 0; k_2 <= n; ++k_2){
      takeProjection(modes, mc1d, k_2, 0);           
      
      //IMPORTANT: here we have to use extended version, because a 1D projection of modes {a_k} (obtained by fixing the second component)
      //from a 2D set of modes may not satisfy the condition a_{-k}=\overline{a_k}.
      fft1d.extendedTransform(mc1d, rProjection, aliasingRemoval, selfSorting);
      
      //We want the result to be saved in Grid1D, representing a DFT values calculated for all of the discrete points, and for all
      //possible values of the modes index second component (which was fixed).

      for(j_1=0; j_1 < m; ++j_1){
        gridOfModes[j_1].set(Index1D(k_2), rProjection[j_1]);
      }
      
    }         

    ///second, calculate M independent FFTs, one for each j_1 discrete point index that were calculated in the previous step
    for(j_1=0; j_1 < m; ++j_1){
      //IMPORTANT: we know that result is going to be real, therefore we use ''fast version'' which avoids some calculations.
      fft1d.fastTransform(gridOfModes[j_1], rProjection, aliasingRemoval);

      r[j_1]=rProjection; ///r doesn't have more components therefore r[j_1] not r[0][j_1] or r[1][j_1]
    }

    r.setSubspaceType(modes);
    r.projectOntoSubspace();

  }

  inline void fastTransform(const ModesContainer2DType& modes, DFT2DGridType& r, int aliasingRemoval = padding){
    transform(modes, r, aliasingRemoval);
  }
  
  /**
   *
   * @param s two dimensional DFT grid
   * @param r values of modes after inverse transform (direction parameter determines which direction [in 2D case either 0 or 1])
   *          is being calculated)
   */
  inline void inverseTransform(const DFT2DGridType& s, ModesContainer2DType& r, int aliasingRemoval = padding){
    int j_1, k_2;

    for(j_1=0; j_1 < m; ++j_1){
      //IMPORTANT: we know that data is real, therefore we use ''fast version'' which avoids some calculations.
      fft1d.fastInverseTransform(s[j_1], mc1d, aliasingRemoval);
//      fft1d.inverseExtendedTransform(s[j_1], mc1d);
      for(k_2 = -n; k_2 <= n; ++k_2){
        modesContainerOfGrid[Index1D(k_2)][j_1] = mc1d[Index1D(k_2)];
      }
    }
    for(k_2 = 0; k_2 <= n; ++k_2){
      //IMPORTANT: here we cannot use ''fast version'', because we are calculating a 1D modes projection {a_k} (obtained by
      //fixing the second component) from a 2D set of modes may not satisfy the condition a_{-k}=\overline{a_k}.
      fft1d.inverseExtendedTransform(modesContainerOfGrid[Index1D(k_2)], mc1d, aliasingRemoval, selfSorting);
      setProjection(mc1d, r, k_2);
    }
  }

  ///TODO: temporary functions, for debugging
  void printModes(const ModesContainerType& mc) const{
    IndexType index;
    IndexRangeType ir;
    int j;
    ir.setRange(0, capd::jaco::strong, n, capd::jaco::weak);
//    for(j=0; j < index.d(); ++j){
      for(index = firstModeIndex(ir); !index.limitReached(ir); index.inc(ir)) {
        generalDebug << index << " " << mc[index] << "\n";
      }
 //   }
  }

  inline void fastInverseTransform(const DFT2DGridType& s, ModesContainer2DType& r, int aliasingRemoval = padding){
    inverseTransform(s, r, aliasingRemoval);
    r.setSubspaceType(s);
    r.projectOntoSubspace();
  }


};





/**Calculates 2D FFT for M which is a power of two. This is a variant dedicated for use in the Taylor integrators, it
 * uses specific data representation. First template parameter (ScalarT) is a type that is transformed (complex number,
 * Taylor coefficients), second template parameter (ComplexScalarT) is a class representing a complex number (needed in
 * order to construct the complex matrices W and hatW - representing the discrete transform and inverse transform.
 *
 */
template< class ScalarT, class ComplexScalarT, int N, int M,
  class ModesContainerT = capd::jaco::ComplexPolyBdJetOptimized<typename ComplexScalarT::ScalarType, ScalarT, capd::jaco::Index2DTwoComponents, 2*2*2*N*N> >
class FFT2D : public ModesContainerT::SubspaceType{
public:
  typedef ScalarT ScalarType;
  typedef ComplexScalarT ComplexScalarType;
  typedef typename ComplexScalarType::ScalarType IntervalType;
  typedef FFT1D<ScalarType, ComplexScalarType, N, M> FFT1DType;
  typedef typename FFT1DType::VectorType Vector1DType;
  typedef ModesContainerT ModesContainer2DType;
  typedef typename FFT1DType::ModesContainerType ModesContainer1DType;
  typedef ModesContainer2DType ModesContainerType;
  typedef capd::jaco::DFT2DGrid<IntervalType, ScalarType, M> DFT2DGridType;
  typedef capd::jaco::ComponentGrid<DFT2DGridType, 2> DFT2DComponentGridType;
  typedef DFT2DComponentGridType DFTGridType;
  typedef typename DFT2DGridType::DFT1DGridType DFT1DGridType;
  typedef capd::jaco::Index2DTwoComponents Index2DType;
  typedef Index2DType IndexType;
  typedef typename FFT1DType::IndexType Index1DType;

  typedef capd::jaco::DFT1DGrid<IntervalType, ModesContainer1DType, M> Grid1DOfModesContainer1DType; ///<is that not too complicated?
  typedef capd::jaco::ComplexPolyBdJetOptimized<IntervalType, DFT1DGridType, Index1DType, 2*N > ModesContainer1DOfGrid1DType;
  typedef typename ModesContainerType::SubspaceType SubspaceType;
  typedef typename SubspaceType::IndexRangeType IndexRangeType;


  int n;
  int m;
  FFT1DType fft1d;
  DFT1DGridType rProjection;
  ModesContainer1DType mc1d;
  ModesContainer2DType projected, grad1, grad2;
  Grid1DOfModesContainer1DType gridOfModes;
  ModesContainer1DOfGrid1DType modesContainerOfGrid;
  capd::vectalg::Matrix<DFT2DGridType, 0, 0> s; ///< for temporary results
  capd::vectalg::Matrix<ModesContainer2DType, 0, 0> t; ///< for temporary results

  FFT2D(){}
  
  FFT2D(int n_, int m_, IntervalType pi_) : n(n_), m(m_), SubspaceType(n_, 2 * n_), fft1d(n_, m_, pi_), rProjection(m_), mc1d(n_),
      gridOfModes(m_, false), modesContainerOfGrid(n_), s(2, 2, false), t(2, 2, false), projected(n_), grad1(n_), grad2(n_){
    int i, j;
    for(i=0; i < m; ++i)
      gridOfModes[i] = ModesContainer1DType(n);
    for(i=-n; i <= n; ++i)
      modesContainerOfGrid[Index1D(i)] = DFT1DGridType(m);
    for(i=0; i < 2; i++)
      for(j=0; j < 2; j++){
        s[i][j] = DFT2DGridType(m);
        t[i][j] = ModesContainer2DType( n );
      }
  }


  /**Takes 1D projection of 2D discrete set of modes, and returns 1D ModesContainer.
   *
   * @param direction this is the component of modes for which the transform is being calculated (either 0 or 1)
   */
  inline void takeProjection(const ModesContainer2DType& mc2d, ModesContainer1DType& mc1d, int fixedSecondComponent, int direction) const{
    int k_1;
    Index2DType index;
    index.l = direction;
    index[1] = fixedSecondComponent;
    for(k_1 = -n; k_1 <= n; ++k_1){
      index[0] = k_1;
      mc1d[Index1D(k_1)] = mc2d[index];
    }
  }

  /**Sets a fixed 1D projection of 2D discrete set of modes with provided array of values.
   * @param direction this is the component of modes for which the transform is being calculated (either 0 or 1)
   */
  inline void setProjection(const ModesContainer1DType& mc1d, ModesContainer2DType& mc2d, int fixedSecondComponent) const{
    int k_1;
    Index2DType index;
    index[1] = fixedSecondComponent;
    for(k_1 = -n; k_1 <= n; ++k_1){
      index[0] = k_1;
      index.l = 0;
      mc2d.set(index, mc1d[Index1D(k_1)]);
      index.l = 1;
      mc2d.set(index, mc1d[Index1D(k_1)]);
    }
  }

  inline void setProjection(const ModesContainer1DType& mc1d, ModesContainer2DType& mc2d, int fixedSecondComponent, int component) const{
    int k_1;
    Index2DType index;
    index[1] = fixedSecondComponent;
    for(k_1 = -n; k_1 <= n; ++k_1){
      index[0] = k_1;
      index.l = component;
      mc2d.set(index, mc1d[Index1D(k_1)]);
    }
  }


/*  inline void scalarProduct(const DFTGridType& gridJ, const DFTGridType& grid1ImJ, const DFTGridType& grid2ImJ, ModesContainer2DType& r){
    s[0][0].multiply(gridJ[0], grid1ImJ[0]);
    s[1][0].multiply(gridJ[1], grid1ImJ[1]);
    s[0][1].multiply(gridJ[0], grid2ImJ[0]);
    s[1][1].multiply(gridJ[1], grid2ImJ[1]);

    IndexType index, index_projected;
    IndexRangeType ir;
    ir.setRange(0, capd::jaco::strong, n, capd::jaco::weak);

    inverseExtendedTransform(s[0][0], t[0][0]);
    inverseExtendedTransform(s[1][0], t[1][0]);
    inverseExtendedTransform(s[0][1], t[0][1]);
    inverseExtendedTransform(s[1][1], t[1][1]);
    //attention: all t's above have one component

    for(index = this->firstModeIndex(ir); !index.limitReached(ir); index.inc(ir)) {
      //t[0][index.l].set(index, index[0]*t[0][index.l][index]);
      //t[1][index.l].set(index, index[1]*t[1][index.l][index]);
      index_projected = index;
      index_projected.l = 0;
      r[index] = t[0][index.l][index_projected] + t[1][index.l][index_projected];
    }
  } */

  /**
   * calculates the scalar product
   * (U\cdot\nabla)V
   *
   * grid stores the Fourier transform of U (two component vector, as U has two components)
   *
   * grid.gradient[0] stores the Fourier transform of \nabla V_1 (also two component vector)
   *
   * grid.gradient[1] stores the Fourier transform of \nabla V_2 (also two component vector)
   */
  inline void scalarProduct(const DFTGridType& grid, ModesContainer2DType& r){
    s[0][0].multiply(grid[0], grid.gradient[0][0]);
    s[1][0].multiply(grid[1], grid.gradient[0][1]);
    s[0][1].multiply(grid[0], grid.gradient[1][0]);
    s[1][1].multiply(grid[1], grid.gradient[1][1]);

    IndexType index, index_projected;
    IndexRangeType ir;
    ir.setRange(0, capd::jaco::strong, n, capd::jaco::weak);

    inverseExtendedTransform(s[0][0], t[0][0]);
    inverseExtendedTransform(s[1][0], t[1][0]);
    inverseExtendedTransform(s[0][1], t[0][1]);
    inverseExtendedTransform(s[1][1], t[1][1]);
    //attention: all t's above have one component

    for(index = this->firstModeIndex(ir); !index.limitReached(ir); index.inc(ir)) {
      //t[0][index.l].set(index, index[0]*t[0][index.l][index]);
      //t[1][index.l].set(index, index[1]*t[1][index.l][index]);
      index_projected = index;
      index_projected.l = 0;
      r[index] = t[0][index.l][index_projected] + t[1][index.l][index_projected];
    }
  }


  inline void scalarProduct(const DFTGridType& grid, DFTGridType& r){
    s[0][0].multiply(grid[0], grid.gradient[0][0]);
    s[1][0].multiply(grid[1], grid.gradient[0][1]);
    s[0][1].multiply(grid[0], grid.gradient[1][0]);
    s[1][1].multiply(grid[1], grid.gradient[1][1]);

    r[0] = s[0][0] + s[1][0];
    r[1] = s[0][1] + s[1][1];
  }

  inline void scalarProduct(const DFTGridType& grid1, const DFTGridType& grid2, DFTGridType& r){
    s[0][0].multiply(grid1[0], grid2.gradient[0][0]);
    s[1][0].multiply(grid1[1], grid2.gradient[0][1]);
    s[0][1].multiply(grid1[0], grid2.gradient[1][0]);
    s[1][1].multiply(grid1[1], grid2.gradient[1][1]);

    r[0] = s[0][0] + s[1][0];
    r[1] = s[0][1] + s[1][1];
  }


  /**
   * calculates the scalar product
   * (U\cdot\nabla)V
   *
   * grid1 stores the Fourier transform of U (two component vector, as U has two components)
   *
   * grid2.gradient[0] stores the Fourier transform of \nabla V_1 (also two component vector)
   *
   * grid2.gradient[1] stores the Fourier transform of \nabla V_2 (also two component vector)
   */
  inline void scalarProduct(const DFTGridType& grid1, const DFTGridType& grid2, ModesContainer2DType& r){
    s[0][0].multiply(grid1[0], grid2.gradient[0][0]);
    s[1][0].multiply(grid1[1], grid2.gradient[0][1]);
    s[0][1].multiply(grid1[0], grid2.gradient[1][0]);
    s[1][1].multiply(grid1[1], grid2.gradient[1][1]);

    IndexType index, index_projected;
    IndexRangeType ir;
    ir.setRange(0, capd::jaco::strong, n, capd::jaco::weak);

    inverseExtendedTransform(s[0][0], t[0][0]);
    inverseExtendedTransform(s[1][0], t[1][0]);
    inverseExtendedTransform(s[0][1], t[0][1]);
    inverseExtendedTransform(s[1][1], t[1][1]);
    //attention: all t's above have one component

    for(index = this->firstModeIndex(ir); !index.limitReached(ir); index.inc(ir)) {
      //t[0][index.l].set(index, index[0]*t[0][index.l][index]);
      //t[1][index.l].set(index, index[1]*t[1][index.l][index]);
      index_projected = index;
      index_projected.l = 0;
      r[index] = t[0][index.l][index_projected] + t[1][index.l][index_projected];
    }


  }


  /**
   * More basic version of the transform function. The component which is being transformed is provided as the input.
   */
  inline void extendedTransform(const ModesContainer2DType& modes, DFT2DGridType& r, int component){
    int k_2, j_1;
    ///first, calculate N independent FFTs, one for each fixed second index component
    for(k_2 = 0; k_2 <= n; ++k_2){
      takeProjection(modes, mc1d, k_2, component);
      //IMPORTANT: here we have to use extended version, because a 1D projection of modes {a_k} (obtained by fixing the second component)
      //from a 2D set of modes may not satisfy the condition a_{-k}=\overline{a_k}.
      fft1d.extendedTransform(mc1d, rProjection);

      //We want the result to be saved in Grid1D, representing a DFT values calculated for all of the discrete points, and for all
      //possible values of the modes index second component (which was fixed).
      for(j_1=0; j_1 < m; ++j_1){
        gridOfModes[j_1].set(Index1D(k_2), rProjection[j_1]);// we have to set the conjugates also (unnecessary to calculate them)
      }
    }

    ///second, calculate M independent FFTs, one for each j_1 discrete point index that were calculated in the previous step
    for(j_1=0; j_1 < m; ++j_1){
      //IMPORTANT: we know that result is going to be real, therefore we use ''fast version'' which avoids some calculations.
      fft1d.fastTransform(gridOfModes[j_1], rProjection);

      r[j_1]=rProjection;
    }
  }

  inline void extendedTransform(const ModesContainer2DType& modes, DFTGridType& r){
    extendedTransform(modes, r[0], 0);
    extendedTransform(modes, r[1], 1);
  }

  inline void CalculateGradients(const ModesContainer2DType& u, ModesContainer2DType& grad1, ModesContainer2DType& grad2){
    IndexRangeType ir;
    ir.setRange(0, capd::jaco::strong, n, capd::jaco::weak);
    for(IndexType ind = this->firstModeIndex(ir), ind_2; !ind.limitReached(ir); ind.inc(ir)){
      if(ind.l == 0){
        grad1.set(ind,  ind[0] * (ComplexScalarType::i() * u[ind])); //partial u_1 / partial x_1
        ind_2 = ind;
        ind_2.l = 1;
        grad1.set(ind_2,  ind[1] * (ComplexScalarType::i() * u[ind])); //partial u_1 / partial x_2
      }else{
        grad2.set(ind,   ind[1] * (ComplexScalarType::i() * u[ind])); //partial u_1 / partial x_1
        ind_2 = ind;
        ind_2.l = 0;
        grad2.set(ind_2, ind[0] * (ComplexScalarType::i() * u[ind])); //partial u_1 / partial x_2
      }
    }
  }

  inline void project(const ModesContainer2DType& in, ModesContainer2DType& projected){
    IndexRangeType range;
    range.setRange(0, capd::jaco::strong, n, capd::jaco::weak);
    IndexType i, i2;

    for(i = this->firstModeIndex(range), i.l = 0; !i.limitReached(range); i.inc(range, true)){
      //calculate scalar products
      ScalarType scpr = i[0] * in[i] + i[1] * in[i.nextComponent()];
      IntervalType norm = i.squareEuclNorm();
      projected[i] = in[i] - (i[0] / norm) * scpr;
      projected[i.nextComponent()] = in[i.nextComponent()] - (i[1] / norm) * scpr;

    }
  }

  /* This function is the extended version from the one-dimensional transform. It calculates transform of the
   * function components stored in modes, and all gradients (these are stored in r.gradient).
   */
  inline void fastTransform(const ModesContainer2DType& modes, DFTGridType& r){
    CalculateGradients(modes, grad1, grad2);

    project(modes, projected);
    extendedTransform(projected, r[0], 0);
    extendedTransform(projected, r[1], 1);

    extendedTransform(grad1, r.gradient[0][0], 0);
    extendedTransform(grad1, r.gradient[0][1], 1);

    extendedTransform(grad2, r.gradient[1][0], 0);
    extendedTransform(grad2, r.gradient[1][1], 1);

  }

  /**
   *
   * @param s two dimensional DFT grid
   * @param r values of modes after inverse transform (direction parameter determines which direction [in 2D case either 0 or 1])
   *          is being calculated)
   */
  inline void inverseExtendedTransform(const DFT2DGridType& s, ModesContainer2DType& r){
    int j_1, k_2;

    for(j_1=0; j_1 < m; ++j_1){
      //IMPORTANT: we know that data is real, therefore we use ''fast version'' which avoids some calculations.
      fft1d.fastInverseTransform(s[j_1], mc1d);
//      fft1d.inverseExtendedTransform(s[j_1], mc1d);
      for(k_2 = -n; k_2 <= n; ++k_2){
        modesContainerOfGrid[Index1D(k_2)][j_1] = mc1d[Index1D(k_2)];
      }
    }
    for(k_2 = 0; k_2 <= n; ++k_2){
      //IMPORTANT: here we cannot use ''fast version'', because we are calculating a 1D modes projection {a_k} (obtained by
      //fixing the second component) from a 2D set of modes may not satisfy the condition a_{-k}=\overline{a_k}.
      fft1d.inverseExtendedTransform(modesContainerOfGrid[Index1D(k_2)], mc1d);
      setProjection(mc1d, r, k_2);
    }
  }


  /**
   *
   * @param s two dimensional , two component DFT grid
   * @param r values of modes after inverse transform (component parameter determines which direction [in 2D case either 0 or 1])
   *          is being calculated)
   */
  inline void inverseExtendedTransform(const DFTGridType& s, ModesContainerType& r, int component){
    int j_1, k_2;

    for(j_1=0; j_1 < m; ++j_1){
      //IMPORTANT: we know that data is real, therefore we use ''fast version'' which avoids some calculations.
      fft1d.fastInverseTransform(s[component][j_1], mc1d);
//      fft1d.inverseExtendedTransform(s[j_1], mc1d);
      for(k_2 = -n; k_2 <= n; ++k_2){
        modesContainerOfGrid[Index1D(k_2)][j_1] = mc1d[Index1D(k_2)];
      }
    }
    for(k_2 = 0; k_2 <= n; ++k_2){
      //IMPORTANT: here we cannot use ''fast version'', because we are calculating a 1D modes projection {a_k} (obtained by
      //fixing the second component) from a 2D set of modes may not satisfy the condition a_{-k}=\overline{a_k}.
      fft1d.inverseExtendedTransform(modesContainerOfGrid[Index1D(k_2)], mc1d);
      setProjection(mc1d, r, k_2, component);
    }
  }


  void printModes(const ModesContainerType& mc) const{
    IndexType index;
    IndexRangeType ir;
    int j;
    ir.setRange(0, capd::jaco::strong, n, capd::jaco::weak);
//    for(j=0; j < index.d(); ++j){
      for(index = firstModeIndex(ir); !index.limitReached(ir); index.inc(ir)) {
        generalDebug << index << " " << mc[index] << "\n";
      }
 //   }
  }


};


}
}

#endif /* FFT_H_ */
