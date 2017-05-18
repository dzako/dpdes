/*
 * Equations.h
 *
 *  Created on: Dec 6, 2011
 *      Author: cyranka
 */

#ifndef EQUATIONS_H_
#define EQUATIONS_H_

#include "PolyBd.h"
#include "FFT.h"
#include <typeinfo>

namespace capd {
namespace jaco {

enum Range {
  full, finitePart, finitePartPlusTail, redundantRange
};

template<class PolyBdT>
class DPDEDefinition {
public:
  typedef PolyBdT PolyBdType;
  typedef typename PolyBdType::RealType RealType;
  typedef typename PolyBdType::ComplexScalarType ComplexType;
  typedef typename PolyBdType::IndexType IndexType;

  bool uproject;

  virtual RealType ni(int k) const = 0;
  virtual RealType ni(const IndexType& index) const = 0;
  virtual bool isDissipative(int k) = 0;
  virtual RealType lambda(int i) const = 0;
  virtual RealType lambda_k(int k) const = 0;
  virtual RealType lambda_k(const IndexType& k) const = 0;
  virtual RealType V(int K) const = 0;
  virtual RealType V(const IndexType& k) const = 0;
  virtual int maximumPoint(const RealType& h, int r, int k) const = 0;
};

///the Burgers PDE
template<class PolyBdT>
class Burgers: public PolyBdT::SubspaceType, public capd::jaco::DPDEDefinition<
    PolyBdT> {
public:
  typedef PolyBdT PolyBdType;
  typedef typename PolyBdType::SubspaceType SubspaceType;
  typedef typename PolyBdType::ComplexScalarType ComplexType;
  typedef typename PolyBdType::RealType RealType;
  typedef typename PolyBdType::IndexType IndexType;

  RealType nu;
  int m_p;
  int m_d;
  int m_r;
  int m_sufficientlyLarge;
  int m_factorMaxM;
  ComplexType m_N_coeff;
  double S_DISSIPATIVE;
  RealType piOver_l;

  Burgers(RealType nu_) :
      nu(nu_), S_DISSIPATIVE(-1), piOver_l(1.) {
    m_p = 2.;
    m_d = 1.;
    m_r = 1.;
    m_sufficientlyLarge = m_d + m_p + 1;
    m_N_coeff = -0.5;
    m_factorMaxM = 3;
  }

  Burgers(int m, int M, RealType nu_) :
      SubspaceType(m, M), nu(nu_), S_DISSIPATIVE(-1), piOver_l(1.) {
    m_p = 2.;
    m_d = 1.;
    m_r = 1.;
    m_sufficientlyLarge = m_d + m_p + 1;
    m_N_coeff = -0.5;
    m_factorMaxM = 3;
  }

  RealType ni(int k) const {
    return nu;
  }

  RealType ni(const IndexType& index) const {
    return nu;
  }

  const ComplexType& Ncoeff() const {
    return m_N_coeff;
  }

  bool isDissipative(int k) {
    if (lambda(k) < S_DISSIPATIVE) {
      return true;
    }
    return false;
  }

  /**i is position in an array of a mode.
   * Returns real part of the eigenvalue \lambda
   */
  RealType lambda(int i) const {
    IndexType j = array2modeIndex(i);
    //return -j.squareEuclNorm() * ni(j);
    return -j.squareEuclNormAlpha( __ALPHA__ ) * ni(j);  ///TODO: changed to flat torus domain
  }

  ///!!!!!! IMPORTANT HERE CHANGED EIGENVALUES TO BILAPLACIAN !!!!!

  ///k is the index of mode, lambda_k from the paper.
  RealType lambda_k(int k) const {
    std::cout << "FATAL ERROR  !\n";
    return -k * k * ni(k);
  }

  RealType lambda_k(const IndexType& k) const {
    //return -k.squareEuclNorm() * ni(k);
    return -k.squareEuclNormAlpha( __ALPHA__ ) * ni(k);  ///TODO: changed to flat torus domain
  }

  ///function returning V(K)=\{ \inf{v(|k|} | |k|>=K} \}, see definition of dissipative PDE in the paper.
  ///Eigenvalues of a dPDE satisfies \lambda_k=-\nu(|k|)|k|^p .
  RealType V(int K) const {
    return leftBound(ni(K));
  }

  RealType V(const IndexType& k) const {
    return leftBound(ni(k));
  }

  /**returns smallest integer K larger than k, such that f(i) is monotonously non increasing function for i>K.
   * Where f(i) = e^{h\lambda_k} k^r
   */
  int maximumPoint(const RealType& h, int r, int k) const {
    return rightBound(ceil(rightBound(sqrt(RealType(r / (2. * nu * h))))));
  }
  using SubspaceType::array2modeIndex;
};



///the Burgers PDE
template<class PolyBdT>
class BurgersBilaplacian: public PolyBdT::SubspaceType, public capd::jaco::DPDEDefinition<
    PolyBdT> {
public:
  typedef PolyBdT PolyBdType;
  typedef typename PolyBdType::SubspaceType SubspaceType;
  typedef typename PolyBdType::ComplexScalarType ComplexType;
  typedef typename PolyBdType::RealType RealType;
  typedef typename PolyBdType::IndexType IndexType;

  RealType nu;
  int m_p;
  int m_d;
  int m_r;
  int m_sufficientlyLarge;
  int m_factorMaxM;
  ComplexType m_N_coeff;
  double S_DISSIPATIVE;
  RealType piOver_l;

  BurgersBilaplacian(RealType nu_) :
      nu(nu_), S_DISSIPATIVE(-1), piOver_l(1.) {
    m_p = 2.;
    m_d = 1.;
    m_r = 1.;
    m_sufficientlyLarge = m_d + m_p + 1;
    m_N_coeff = -0.5;
    m_factorMaxM = 3;
  }

  BurgersBilaplacian(int m, int M, RealType nu_) :
      SubspaceType(m, M), nu(nu_), S_DISSIPATIVE(-1), piOver_l(1.) {
    m_p = 2.;
    m_d = 1.;
    m_r = 1.;
    m_sufficientlyLarge = m_d + m_p + 1;
    m_N_coeff = -0.5;
    m_factorMaxM = 3;
  }

  RealType ni(int k) const {
    return nu;
  }

  RealType ni(const IndexType& index) const {
    return nu;
  }

  const ComplexType& Ncoeff() const {
    return m_N_coeff;
  }

  bool isDissipative(int k) {
    if (lambda(k) < S_DISSIPATIVE) {
      return true;
    }
    return false;
  }

  /**i is position in an array of a mode.
   * Returns real part of the eigenvalue \lambda
   */
  RealType lambda(int i) const {
    IndexType j = array2modeIndex(i);
    //return - j.squareEuclNorm() * j.squareEuclNorm() * j.squareEuclNorm() * j.squareEuclNorm() * j.squareEuclNorm() * j.squareEuclNorm() * ni(j);
    return - j.squareEuclNorm() * j.squareEuclNorm() * ni(j);
  }

  ///k is the index of mode, lambda_k from the paper.
  RealType lambda_k(int k) const {
    //return - (k * k) * (k * k) * (k * k) * (k * k) * (k * k) * (k * k) * ni(k);
    return - (k * k) * (k * k) * ni(k);
  }

  RealType lambda_k(const IndexType& k) const {
    //return - k.squareEuclNorm() * k.squareEuclNorm() * k.squareEuclNorm() * k.squareEuclNorm() * k.squareEuclNorm() * k.squareEuclNorm() * ni(k);
    return - k.squareEuclNorm() * k.squareEuclNorm() * ni(k);
  }

  ///function returning V(K)=\{ \inf{v(|k|} | |k|>=K} \}, see definition of dissipative PDE in the paper.
  ///Eigenvalues of a dPDE satisfies \lambda_k=-\nu(|k|)|k|^p .
  RealType V(int K) const {
    return leftBound(ni(K));
  }

  RealType V(const IndexType& k) const {
    return leftBound(ni(k));
  }

  /**returns smallest integer K larger than k, such that f(i) is monotonously non increasing function for i>K.
   * Where f(i) = e^{h\lambda_k} k^r
   */
  int maximumPoint(const RealType& h, int r, int k) const {
    return rightBound(ceil(rightBound(sqrt(RealType(r / (2. * nu * h))))));
  }
  using SubspaceType::array2modeIndex;
};


///the Euler PDE
template<class PolyBdT>
class Euler : public PolyBdT::SubspaceType, public capd::jaco::DPDEDefinition<
    PolyBdT> {
public:
  typedef PolyBdT PolyBdType;
  typedef typename PolyBdType::SubspaceType SubspaceType;
  typedef typename PolyBdType::ComplexScalarType ComplexType;
  typedef typename PolyBdType::RealType RealType;
  typedef typename PolyBdType::IndexType IndexType;

  RealType nu;
  int m_p;
  int m_d;
  int m_r;
  int m_sufficientlyLarge;
  int m_factorMaxM;
  ComplexType m_N_coeff;
  double S_DISSIPATIVE;
  RealType piOver_l;

  Euler(RealType nu_) :
      nu(nu_), S_DISSIPATIVE(-1), piOver_l(1.) {
    m_p = 2.;
    m_d = 1.;
    m_r = 1.;
    m_sufficientlyLarge = m_d + m_p + 1;
    m_N_coeff = -1.;
    m_factorMaxM = 3;
  }

  Euler(int m, int M, RealType nu_) :
      SubspaceType(m, M), nu(nu_), S_DISSIPATIVE(-1), piOver_l(1.) {
    m_p = 2.;
    m_d = 1.;
    m_r = 1.;
    m_sufficientlyLarge = m_d + m_p + 1;
    m_N_coeff = -1.;
    m_factorMaxM = 3;
  }

  RealType ni(int k) const {
    return nu;
  }

  RealType ni(const IndexType& index) const {
    return nu;
  }

  const ComplexType& Ncoeff() const {
    return m_N_coeff;
  }

  bool isDissipative(int k) {
    if (lambda(k) < S_DISSIPATIVE) {
      return true;
    }
    return false;
  }

//  /**i is position in an array of a mode.
//   * Returns real part of the eigenvalue \lambda
//   */
//  RealType lambda(int i) const {
//    return 0.;
//  }
//
//
//  ///k is the index of mode, lambda_k from the paper.
//  RealType lambda_k(int k) const {
//    return 0.;
//  }
//
//  RealType lambda_k(const IndexType& k) const {
//    return 0.;
//  }
//
//  ///function returning V(K)=\{ \inf{v(|k|} | |k|>=K} \}, see definition of dissipative PDE in the paper.
//  ///Eigenvalues of a dPDE satisfies \lambda_k=-\nu(|k|)|k|^p .
//  RealType V(int K) const {
//    throw std::runtime_error("not implemented, Euler PDE is NOT dissipative");
//  }
//
//  RealType V(const IndexType& k) const {
//    throw std::runtime_error("not implemented, Euler PDE is NOT dissipative");
//  }
//
//  /**returns smallest integer K larger than k, such that f(i) is monotonously non increasing function for i>K.
//   * Where f(i) = e^{h\lambda_k} k^r
//   */
//  int maximumPoint(const RealType& h, int r, int k) const {
//    throw std::runtime_error("not implemented, Euler PDE is NOT dissipative");
//  }
//
//
//  /**i is position in an array of a mode.
//     * Returns real part of the eigenvalue \lambda
//     */
//    RealType lambda(int i) const {
//      IndexType j = array2modeIndex(i);
//      return -j.squareEuclNorm() * ni(j);
//    }


    RealType lambda(int i) const {
      IndexType j = array2modeIndex(i);
      return - j.squareEuclNormAlpha( __ALPHA__ ) * ni(j);
      //return 0;
    }


    ///k is the index of mode, lambda_k from the paper.
    RealType lambda_k(int k) const {
      return - k * k * ni(k);
      //return 0;
    }

    RealType lambda_k(const IndexType& k) const {
      return - k.squareEuclNormAlpha( __ALPHA__ ) * ni(k);
      //return 0;
    }

    ///function returning V(K)=\{ \inf{v(|k|} | |k|>=K} \}, see definition of dissipative PDE in the paper.
    ///Eigenvalues of a dPDE satisfies \lambda_k=-\nu(|k|)|k|^p .
    RealType V(int K) const {
      return leftBound(ni(K));
    }

    RealType V(const IndexType& k) const {
      return leftBound(ni(k));
    }

    /**returns smallest integer K larger than k, such that f(i) is monotonously non increasing function for i>K.
     * Where f(i) = e^{h\lambda_k} k^r
     */
    int maximumPoint(const RealType& h, int r, int k) const {
      return rightBound(ceil(rightBound(sqrt(RealType(r / (2. * nu * h))))));
    }


  using SubspaceType::array2modeIndex;
};



///the alpha-Euler PDE
template<class PolyBdT>
class AlphaEuler : public PolyBdT::SubspaceType, public capd::jaco::DPDEDefinition<
    PolyBdT> {
public:
  typedef PolyBdT PolyBdType;
  typedef typename PolyBdType::SubspaceType SubspaceType;
  typedef typename PolyBdType::ComplexScalarType ComplexType;
  typedef typename PolyBdType::RealType RealType;
  typedef typename PolyBdType::IndexType IndexType;

  RealType nu;
  int m_p;
  int m_d;
  int m_r;
  int m_sufficientlyLarge;
  int m_factorMaxM;
  ComplexType m_N_coeff;
  double S_DISSIPATIVE;
  RealType piOver_l;

  AlphaEuler(RealType nu_) :
      nu(nu_), S_DISSIPATIVE(-1), piOver_l(1.) {
    m_p = 2.;
    m_d = 1.;
    m_r = 1.;
    m_sufficientlyLarge = m_d + m_p + 1;
    m_N_coeff = -1.;
    m_factorMaxM = 3;
  }

  AlphaEuler(int m, int M, RealType nu_) :
      SubspaceType(m, M), nu(nu_), S_DISSIPATIVE(-1), piOver_l(1.) {
    m_p = 2.;
    m_d = 1.;
    m_r = 1.;
    m_sufficientlyLarge = m_d + m_p + 1;
    m_N_coeff = -1.;
    m_factorMaxM = 3;
  }

  RealType ni(int k) const {
    return nu;
  }

  RealType ni(const IndexType& index) const {
    return nu;
  }

  const ComplexType& Ncoeff() const {
    return m_N_coeff;
  }

  bool isDissipative(int k) {
    if (lambda(k) < S_DISSIPATIVE) {
      return true;
    }
    return false;
  }

  RealType lambda(int i) const {
    IndexType j = array2modeIndex(i);
    return - ni(j);
    //return 0;
  }


  ///k is the index of mode, lambda_k from the paper.
  RealType lambda_k(int k) const {
    return - ni(k);
    //return 0;
  }

  RealType lambda_k(const IndexType& k) const {
    return - ni(k);
    //return 0;
  }

  ///function returning V(K)=\{ \inf{v(|k|} | |k|>=K} \}, see definition of dissipative PDE in the paper.
  ///Eigenvalues of a dPDE satisfies \lambda_k=-\nu(|k|)|k|^p .
  RealType V(int K) const {
    return leftBound(ni(K));
  }

  RealType V(const IndexType& k) const {
    return leftBound(ni(k));
  }

  /**returns smallest integer K larger than k, such that f(i) is monotonously non increasing function for i>K.
   * Where f(i) = e^{h\lambda_k} k^r
   */
  int maximumPoint(const RealType& h, int r, int k) const {
    return rightBound(ceil(rightBound(sqrt(RealType(r / (2. * nu * h))))));
  }


  using SubspaceType::array2modeIndex;
};



///the real Ginzburg-Landau PDE
template<class PolyBdT>
class GL: public PolyBdT::SubspaceType, public capd::jaco::DPDEDefinition<
    PolyBdT> {
public:
  typedef PolyBdT PolyBdType;
  typedef typename PolyBdType::SubspaceType SubspaceType;
  typedef typename PolyBdType::ComplexScalarType ComplexType;
  typedef typename PolyBdType::RealType RealType;
  typedef typename PolyBdType::IndexType IndexType;

  RealType nu;
  int m_p;
  int m_d;
  int m_r;
  int m_sufficientlyLarge;
  int m_factorMaxM;
  ComplexType m_N_coeff;
  double S_DISSIPATIVE;
  RealType piOver_l;

  GL(RealType nu_) :
      nu(nu_), S_DISSIPATIVE(-0.01), piOver_l(1.) {
    m_p = 2.;
    m_d = 1.;
    m_r = 0.;
    m_sufficientlyLarge = m_d + m_p + 1;
    m_N_coeff = -1.;
    m_factorMaxM = 3;
  }

  GL(int m, int M, RealType nu_) :
      SubspaceType(m, M), nu(nu_), S_DISSIPATIVE(-0.01), piOver_l(1.) {
    m_p = 2.;
    m_d = 1.;
    m_r = 0.;
    m_sufficientlyLarge = m_d + m_p + 1;
    m_N_coeff = -1.;
    m_factorMaxM = 3;
  }

  RealType ni(int k) const {
    return nu;
  }

  RealType ni(const IndexType& k) const {
    return nu;
  }

  const ComplexType& Ncoeff() const {
    return m_N_coeff;
  }

  bool isDissipative(int k) {
    if (lambda(k) < S_DISSIPATIVE) {
      return true;
    }
    return false;
  }

  ///i is position in an array of a mode.
  RealType lambda(int i) const {
    IndexType j = array2modeIndex(i);
    return -nu * j.squareEuclNorm() + 1;
  }

  ///k is the index of mode, lambda_k from the paper.
  RealType lambda_k(int k) const {
    return -k * k + 1.;
  }

  RealType lambda_k(const IndexType& k) const {
    return -nu * k.squareEuclNorm() + 1.;
  }

//    ///i is position in the near tail array of a mode.
//    ComplexType lambdaTail(int i) const{ int j=array2modeTail(i); return -nu * j * j + 1.; }

  ComplexType getNCoeff() const {
    return Ncoeff;
  }

  ///function returning V(K)=\{ \inf{v(|k|} | |k|>=K} \}, see definition of dissipative PDE in the paper.
  ///Eigenvalues of a dPDE satisfies \lambda_k=-\nu(|k|)|k|^p .
  RealType V(int K) const {
    return leftBound(ni(K));
  }

  RealType V(const IndexType& k) const {
    return leftBound(ni(k));
  }

  /**returns smallest integer K larger than k, such that f(i) is monotonously non increasing function for all i>K.
   * Where f(i) = e^{h\lambda_k} k^r
   */
  int maximumPoint(const RealType& h, int r, int k) const {
    return rightBound(ceil(sqrt(RealType(r / (2. * nu * h)))));
  }
  using SubspaceType::array2modeIndex;
};

///the Kuramoto-Shivasinsky equation
template<class PolyBdT>
class KS: public PolyBdT::SubspaceType, public capd::jaco::DPDEDefinition<
    PolyBdT> {
public:
  typedef PolyBdT PolyBdType;
  typedef typename PolyBdType::SubspaceType SubspaceType;
  typedef typename PolyBdType::ComplexScalarType ComplexType;
  typedef typename PolyBdType::RealType RealType;
  typedef typename PolyBdType::IndexType IndexType;

  RealType nu;
  int m_p;
  int m_d;
  int m_r;
  int m_sufficientlyLarge;
  int m_factorMaxM;
  ComplexType m_N_coeff;
  double S_DISSIPATIVE;
  RealType piOver_l;

  KS(RealType nu_) :
      nu(nu_), S_DISSIPATIVE(-0.01), piOver_l(1.) {
    m_p = 4.;
    m_d = 1.;
    m_r = 1.;
    m_sufficientlyLarge = m_d + m_p + 1;
    m_N_coeff = 1.;
    m_factorMaxM = 3;
  }

  KS(int m, int M, RealType nu_) :
      SubspaceType(m, M), nu(nu_), piOver_l(1.) {
    m_p = 4.;
    m_d = 1.;
    m_r = 1.;
    m_sufficientlyLarge = m_d + m_p + 1;
    m_N_coeff = 1.;
    S_DISSIPATIVE = -0.01;
    m_factorMaxM = 3;
  }

  RealType ni(int k) const {
    IndexType j = array2modeIndex(k);
    return (j.isZero() ? 0 : nu - 1. / j.squareEuclNorm());
  }

  RealType ni(const IndexType& k) const {
    return (k.isZero() ? 0 : nu - 1. / k.squareEuclNorm());
  }

  const ComplexType& Ncoeff() const {
    return m_N_coeff;
  }

  bool isDissipative(int k) {
    if (lambda(k) < S_DISSIPATIVE) {
      return true;
    }
    return false;
  }

  ///i is position in an array of a mode.
  RealType lambda(int i) const {
    IndexType j = array2modeIndex(i);
    return -j.squareEuclNorm() * (j.squareEuclNorm() * nu - 1.);
  }

  ///k is the index of mode, lambda_k from the paper.
  RealType lambda_k(int k) const {
    return -k * k * (k * k * nu - 1.);
  }

  RealType lambda_k(const IndexType& k) const {
    return -k.squareEuclNorm() * (k.squareEuclNorm() * nu - 1.);
  }

  ComplexType getNCoeff() const {
    return Ncoeff;
  }

  ///function returning V(K)=\{ \inf{v(|k|} | |k|>=K} \}, see definition of dissipative PDE in the paper.
  ///Eigenvalues of a dPDE satisfies \lambda_k=-\nu(|k|)|k|^p .
  RealType V(int K) const {
    return leftBound(ni(K));
  }

  RealType V(const IndexType& k) const {
    return leftBound(ni(k));
  }

  /**returns smallest integer K larger than k, such that f(i) is monotonously non increasing function for all i>K.
   * Where f(i) = e^{h\lambda_k} k^r
   */
  int maximumPoint(const RealType& h, int r, int k) const {
    if (-4 * rightBound(h) * nu * k * k * k * k + 2 * rightBound(h) * k * k + r
        <= 0 && k * k >= 1 / (4 * nu)) {
      return k;
    } else {
      int c = k;
      while (!(-4 * rightBound(h) * nu * c * c * c * c
          + 2 * rightBound(h) * c * c + r <= 0 && c * c >= 1 / (4 * nu))) {
        c++;
      }
      return c;
    }
  }
  using SubspaceType::array2modeIndex;
};

///the Swift-Hohenberg dPDE
template<class PolyBdT>
class SH: public PolyBdT::SubspaceType, public capd::jaco::DPDEDefinition<
    PolyBdT> {
public:
  typedef PolyBdT PolyBdType;
  typedef typename PolyBdType::SubspaceType SubspaceType;
  typedef typename PolyBdType::ComplexScalarType ComplexType;
  typedef typename PolyBdType::RealType RealType;
  typedef typename PolyBdType::IndexType IndexType;

  RealType nu;
  int m_p;
  int m_d;
  int m_r;
  int m_sufficientlyLarge;
  int m_factorMaxM;
  ComplexType m_N_coeff;
  double S_DISSIPATIVE;
  RealType piOver_l;

  SH(RealType nu_) :
      nu(nu_), S_DISSIPATIVE(-0.01), piOver_l(1.) {
    m_p = 4.;
    m_d = 1.;
    m_r = 0.;
    m_sufficientlyLarge = m_d + m_p + 1;
    m_N_coeff = -1.;
    m_factorMaxM = 10;
  }

  SH(int m, int M, RealType nu_) :
      SubspaceType(m, M), nu(nu_), S_DISSIPATIVE(-0.01), piOver_l(1.) {
    m_p = 4.;
    m_d = 1.;
    m_r = 0.;
    m_sufficientlyLarge = m_d + m_p + 1;
    m_N_coeff = -1.;
    m_factorMaxM = 10;
  }

  RealType ni(int k) const {
    IndexType j = array2modeIndex(k);
    double sqN = j.squareEuclNorm();
    return (j.isZero() ? 0 : 1 - 2. / sqN + (1 - nu) / (sqN * sqN));
  }

  RealType ni(const IndexType& k) const {
    double sqN = k.squareEuclNorm();
    return (k.isZero() ? 0 : 1 - 2. / sqN + (1 - nu) / (sqN * sqN));
  }

  const ComplexType& Ncoeff() const {
    return m_N_coeff;
  }

  bool isDissipative(int k) {
    if (lambda(k) < S_DISSIPATIVE) {
      return true;
    }
    return false;
  }

  ///i is position in an array of a mode.
  RealType lambda(int i) const {
    IndexType j = array2modeIndex(i);
    double sqN = j.squareEuclNorm();
    return -sqN * (sqN - 2.) - 1. + nu;
  }

  ///k is the index of mode, lambda_k from the paper.
  RealType lambda_k(int k) const {
    return -k * k * (k * k - 2.) - 1. + nu;
  }

  RealType lambda_k(const IndexType& k) const {
    double sqN = k.squareEuclNorm();
    return -sqN * (sqN - 2.) - 1. + nu;
  }

  ComplexType getNCoeff() const {
    return Ncoeff;
  }

  ///Returns V(K)=\{ \inf{v(|k|} | |k|>=K} \}, see definition of dissipative PDE in the paper.
  ///Eigenvalues of a dPDE satisfies \lambda_k=-\nu(|k|)|k|^p .
  RealType V(int K) const {
    return leftBound(ni(K));
  }

  RealType V(const IndexType& k) const {
    return leftBound(ni(k));
  }

  /**returns smallest integer K larger than k, such that f(i) is monotonously non increasing function for all i>K.
   * Where f(i) = e^{h\lambda_k} k^r
   */
  int maximumPoint(const RealType& h, int r, int k) const {
    if (-4 * rightBound(h) * k * k * k * k + 4 * rightBound(h) * k * k + r <= 0
        && k * k >= 0.5) {
      return k;
    } else {
      int c = k;
      while (!(-4 * rightBound(h) * c * c * c * c + 4 * rightBound(h) * c * c
          + r <= 0 && c * c >= 0.5)) {
        c++;
      }
      return c;
    }
  }
  using SubspaceType::array2modeIndex;
};

/////the Cahn-Hillard dPDE
//template<class PolyBdT>
//class CH: public PolyBdT::SubspaceType, public capd::jaco::DPDEDefinition<
//    PolyBdT> {
//public:
//  typedef PolyBdT PolyBdType;
//  typedef typename PolyBdType::SubspaceType SubspaceType;
//  typedef typename PolyBdType::ComplexScalarType ComplexType;
//  typedef typename PolyBdType::RealType RealType;
//  typedef typename PolyBdType::IndexType IndexType;
//
//  RealType nu;
//  int m_p;
//  int m_d;
//  int m_r;
//  int m_sufficientlyLarge;
//  int m_factorMaxM;
//  ComplexType m_N_coeff;
//  double S_DISSIPATIVE;
//  RealType piOver_l; //TODO: temporary - should be defined by boundary conditions
//
//  CH(RealType nu_) :
//      nu(nu_), S_DISSIPATIVE(-0.01) {
//    m_p = 4.;
//    m_d = 1.;
//    m_r = 2.;
//    m_sufficientlyLarge = m_d + m_p + 1;
//    m_N_coeff = -nu;
//    piOver_l = 1.;
//    m_factorMaxM = 10;
//  }
//
//  CH(int m, int M, RealType nu_) :
//      SubspaceType(m, M), nu(nu_), S_DISSIPATIVE(-0.01) {
//    m_p = 4.;
//    m_d = 1.;
//    m_r = 2.;
//    m_sufficientlyLarge = m_d + m_p + 1;
//    m_N_coeff = -nu;
//    piOver_l = 1.;
//    m_factorMaxM = 10;
//  }
//
//  RealType ni(int k) const {
//    IndexType j = array2modeIndex(k);
//    RealType sqN = j.squareEuclNorm() * piOver_l * piOver_l * piOver_l
//        * piOver_l;
//    return (j.isZero() ? 0 : 1 - nu / sqN);
//  }
//
//  RealType ni(const IndexType& k) const {
//    RealType sqN = k.squareEuclNorm() * piOver_l * piOver_l * piOver_l
//        * piOver_l;
//    return (k.isZero() ? 0 : 1 - nu / sqN);
//  }
//
//  const ComplexType& Ncoeff() const {
//    return m_N_coeff;
//  }
//
//  bool isDissipative(int k) {
//    if (lambda(k) < S_DISSIPATIVE) {
//      return true;
//    }
//    return false;
//  }
//
//  ///i is position in an array of a mode.
//  RealType lambda(int i) const {
//    IndexType j = array2modeIndex(i);
//    RealType sqN = j.squareEuclNorm() * piOver_l * piOver_l;
//    return -sqN * (sqN - nu);
//  }
//
//  ///k is the index of mode, lambda_k from the paper.
//  RealType lambda_k(int k) const {
//    return -k * k * piOver_l * piOver_l * (k * k * piOver_l * piOver_l - nu);
//  }
//
//  RealType lambda_k(const IndexType& k) const {
//    RealType sqN = k.squareEuclNorm() * piOver_l * piOver_l;
//    return -sqN * (sqN - nu);
//  }
//
//  ComplexType getNCoeff() const {
//    return Ncoeff;
//  }
//
//  //TODO: this functions below are implemented provisionally
//  ///Returns V(K)=\{ \inf{v(|k|} | |k|>=K} \}, see definition of dissipative PDE in the paper.
//  ///Eigenvalues of a dPDE satisfies \lambda_k=-\nu(|k|)|k|^p .
//  RealType V(int K) const {
//    return leftBound(ni(K));
//  }
//
//  RealType V(const IndexType& k) const {
//    return leftBound(ni(k));
//  }
//
//  /**returns smallest integer K larger than k, such that f(i) is monotonously non increasing function for all i>K.
//   * Where f(i) = e^{h\lambda_k} k^r
//   */
//  //TODO: different for this equation
//  int maximumPoint(const RealType& h, int r, int k) const {
//    if (-4 * rightBound(h) * k * k * k * k + 4 * rightBound(h) * k * k + r <= 0
//        && k * k >= 0.5) {
//      return k;
//    } else {
//      int c = k;
//      while (!(-4 * rightBound(h) * c * c * c * c + 4 * rightBound(h) * c * c
//          + r <= 0 && c * c >= 0.5)) {
//        c++;
//      }
//      return c;
//    }
//  }
//  using SubspaceType::array2modeIndex;
//};



///the FitzHugh-Nagumo PDE
#define _EPSILON_  0.01
#define _A_ 0.1
#define _L_CONST_ _A_
#define _GAMMA_ 5.
#define _NU_ 0.049087385212340519350978
//template<class PolyBdT>
//class FN: public PolyBdT::SubspaceType, public capd::jaco::DPDEDefinition<
//    PolyBdT> {
//public:
//  typedef PolyBdT PolyBdType;
//  typedef typename PolyBdType::SubspaceType SubspaceType;
//  typedef typename PolyBdType::ComplexScalarType ComplexType;
//  typedef typename PolyBdType::RealType RealType;
//  typedef typename PolyBdType::IndexType IndexType;
//
//  RealType nu;
//  int m_p;
//  int m_d;
//  int m_r;
//  int m_sufficientlyLarge;
//  int m_factorMaxM;
//  int mhalf;
//  ComplexType m_N_coeff;
//  double S_DISSIPATIVE;
//  RealType piOver_l;
//
//  FN(RealType nu_) :
//      nu(nu_), S_DISSIPATIVE(-0.01), piOver_l(1.) {
//    m_p = 2.;
//    m_d = 1.;
//    m_r = 0.;
//    m_sufficientlyLarge = m_d + m_p + 1;
//    m_N_coeff = -1.;
//    m_factorMaxM = 10;
//    uproject = true;
//  }
//
//  FN(int m, int M, RealType nu_) :
//      SubspaceType(m, M), nu(nu_), S_DISSIPATIVE(-0.01), piOver_l(1.) {
//
//    if( m % 2 == 0 || M % 2 == 0 ){
//      std::cerr << "For FitzHugh-Nagumo m (" << m << ") and M (" << M <<") has to be ODD (there are 0...m modes for u, and 0...m modes for v).\n";
//      exit(1);
//    }
//    mhalf = m / 2;
//    m_p = 2.;
//    m_d = 1.;
//    m_r = 0.;
//    m_sufficientlyLarge = m_d + m_p + 1;
//    m_N_coeff = -1.;
//    m_factorMaxM = 10;
//    uproject = true;
//  }
//
//  RealType ni(int k) const {
//    IndexType j = array2modeIndex(k);
//
//    double sqN = j.squareEuclNorm();
//    return - nu + _L_CONST_ / sqN;
//
//  }
//
//  RealType ni(const IndexType& k) const {
//    double sqN = k.squareEuclNorm();
//
//    return - nu + _L_CONST_ / sqN;
//  }
//
//  const ComplexType& Ncoeff() const {
//    return m_N_coeff;
//  }
//
//  bool isDissipative(int k) {
//    if (lambda(k) < S_DISSIPATIVE) {
//      return true;
//    }
//    return false;
//  }
//
//  //lambda functions return a nonzero value only for first part of the vector (indexed 0...m)
//  //storing solutions first component (u), for the second part (component v) return 0.
//
//  ///i is position in an array of a mode.
//  RealType lambda(int i) const {
//    IndexType j = array2modeIndex(i);
//
//    if( j[0] <= mhalf ){
//      double sqN = j.squareEuclNorm();
//      return - nu * sqN + _L_CONST_;
//    }else{
//      return 0;
//    }
//  }
//
//  ///k is the index of mode, lambda_k from the paper.
//  RealType lambda_k(int k) const {
//    if ( k <= mhalf ){
//      return - nu * k * k + _L_CONST_;
//    }else{
//      return 0;
//    }
//  }
//
//  RealType lambda_k(const IndexType& k) const {
//    if( k[0] <= mhalf ){
//      double sqN = k.squareEuclNorm();
//      return - nu * sqN + _L_CONST_;
//    }else{
//      return 0;
//    }
//  }
//
//  ComplexType getNCoeff() const {
//    return Ncoeff;
//  }
//
//  ///Returns V(K)=\{ \inf{v(|k|} | |k|>=K} \}, see definition of dissipative PDE in the paper.
//  ///Eigenvalues of a dPDE satisfies \lambda_k=-\nu(|k|)|k|^p .
//  RealType V(int K) const {
//    return leftBound(ni(K));
//  }
//
//  RealType V(const IndexType& k) const {
//    return leftBound(ni(k));
//  }
//
//  /**TODO: Old version !!!
//   */
//  int maximumPoint(const RealType& h, int r, int k) const {
//    if (-4 * rightBound(h) * k * k * k * k + 4 * rightBound(h) * k * k + r <= 0
//        && k * k >= 0.5) {
//      return k;
//    } else {
//      int c = k;
//      while (!(-4 * rightBound(h) * c * c * c * c + 4 * rightBound(h) * c * c + r <= 0 && c * c >= 0.5)) {
//        c++;
//      }
//      return c;
//    }
//  }
//  using SubspaceType::array2modeIndex;
//  using DPDEDefinition<PolyBdT>::uproject;
//};



///the Diblock Copolymer dPDE
template<class PolyBdT>
class DBCP: public PolyBdT::SubspaceType, public capd::jaco::DPDEDefinition<
    PolyBdT> {
public:
  typedef PolyBdT PolyBdType;
  typedef typename PolyBdType::SubspaceType SubspaceType;
  typedef typename PolyBdType::ComplexScalarType ComplexType;
  typedef typename PolyBdType::RealType RealType;
  typedef typename PolyBdType::IndexType IndexType;

  RealType nu;
  RealType sigma;
  int m_p;
  int m_d;
  int m_r;
  int m_sufficientlyLarge;
  int m_factorMaxM;
  ComplexType m_N_coeff;
  double S_DISSIPATIVE;
  RealType piOver_l;

  DBCP(RealType nu_) :
      nu(nu_), S_DISSIPATIVE(-0.01) {
    m_p = 4.;
    m_d = 1.;
    m_r = 2.;
    m_factorMaxM = 10;
    m_sufficientlyLarge = m_d + m_p + 1;
    m_N_coeff = nu_;
    piOver_l = 0.5;
    RealType small_interval;
    setLeftBound(small_interval, -1e-16);
    setRightBound(small_interval, 1e-16);
    sigma = 0.40528473456935108  + small_interval;
    std::cout << "sigma=" << sigma << "\n";
  }

  DBCP(int m, int M, RealType nu_) :
      SubspaceType(m, M), nu(nu_), S_DISSIPATIVE(-0.01) {
    m_p = 4.;
    m_d = 1.;
    m_r = 2.;
    m_factorMaxM = 10;
    m_sufficientlyLarge = m_d + m_p + 1;
    m_N_coeff = nu_;
    piOver_l = 0.5;
    RealType small_interval;
    setLeftBound(small_interval, -1e-16);
    setRightBound(small_interval, 1e-16);
    sigma = 0.40528473456935108 + small_interval;
  }

  RealType ni(int k) const {
    IndexType j = array2modeIndex(k);
    RealType sqN = j.squareEuclNorm();
    return (
        j.isZero() ? 0 :
            piOver_l * piOver_l * piOver_l * piOver_l
                - piOver_l * piOver_l * nu / sqN + (nu * sigma) / (sqN * sqN));
  }

  RealType ni(const IndexType& k) const {
    RealType sqN = k.squareEuclNorm();
    return (
        k.isZero() ? 0 :
            piOver_l * piOver_l * piOver_l * piOver_l
                - piOver_l * piOver_l * nu / sqN + (nu * sigma) / (sqN * sqN));
  }

  const ComplexType& Ncoeff() const {
    return m_N_coeff;
  }

  bool isDissipative(int k) {
    if (lambda(k) < S_DISSIPATIVE) {
      return true;
    }
    return false;
  }

  ///i is position in an array of a mode.
  RealType lambda(int i) const {
    IndexType j = array2modeIndex(i);
    RealType sqN = j.squareEuclNorm() * piOver_l * piOver_l;
    return -sqN * (sqN - nu) - nu * sigma;
  }

  ///k is the index of mode, lambda_k from the paper.
  RealType lambda_k(int k) const {
    return -k * k * piOver_l * piOver_l * (k * k * piOver_l * piOver_l - nu)
        - nu * sigma;
  }

  RealType lambda_k(const IndexType& k) const {
    RealType sqN = k.squareEuclNorm() * piOver_l * piOver_l;
    return -sqN * (sqN - nu) - nu * sigma;
  }

  ComplexType getNCoeff() const {
    return Ncoeff;
  }


  ///Returns V(K)=\{ \inf{v(|k|} | |k|>=K} \}, see definition of dissipative PDE in the paper.
  ///Eigenvalues of a dPDE satisfies \lambda_k=-\nu(|k|)|k|^p .
  RealType V(int K) const {
    return leftBound(ni(K));
  }

  RealType V(const IndexType& k) const {
    return leftBound(ni(k));
  }

  /**returns smallest integer K larger than k, such that f(i) is monotonously non increasing function for all i>K.
   * Where f(i) = e^{h\lambda_k} k^r
   */
  int maximumPoint(const RealType& h, int r, int k) const {
    if (-4 * rightBound(h) * k * k * k * k * piOver_l * piOver_l * piOver_l
        * piOver_l + 2 * nu * rightBound(h) * k * k * piOver_l * piOver_l + r
        <= 0 && k * k >= nu / (4 * piOver_l * piOver_l)) {
      return k;
    } else {
      int c = k;
      while (!(-4 * rightBound(h) * c * c * c * c * piOver_l * piOver_l
          * piOver_l * piOver_l
          + 2 * nu * rightBound(h) * c * c * piOver_l * piOver_l + r <= 0
          && c * c >= nu / (4 * piOver_l * piOver_l))) {
        c++;
      }
      return c;
    }
  }
  using SubspaceType::array2modeIndex;
};

/** This is a class representing a nonlinear PDE (this is the input to the integrator). It is the base for other classes representing
 * pdes consisting higher order nonlinear term.
 *
 * Elliptic part is a linear operator, and nonlinear part is the second degree polynomial and its derivatives
 * u^2
 * (u^2)_x
 * (u^2)_{xx}
 * etc...
 *
 * Elliptic part dominates the nonlinear part, in sense the maximal order of the derivative in the elliptic operator is larger than the
 * order of the derivative of u^2 in the nonlinear part.
 *
 * EquationT class defines what specific equation this is. Like for example what are eigenvalues and what is constant in front
 * of the nonlinear term, and what is the order of the polynomial in the nonlinearity.
 * In particular the value m_r in EquationT defines the order of the derivative of u^2 in the nonlinear part (>= 0).
 *
 * FFTT is a FFT class.
 *
 * ORD is the maximal order.
 *
 * When VNormT is set to EuclideanNorm the program intentionally will not compile (because of efficiency reasons - the
 * function norm is not returning int anymore, it must return interval)  - but the case of EculideanNorm was tested too
 */
template<class EquationT, class FFTT, int ORD,
    class VNormT = capd::jaco::MaximumNorm<typename FFTT::IndexType> >
class DPDE2: public EquationT {
public:
  typedef EquationT EquationType;
  typedef FFTT FFTType;
  typedef VNormT VNormType;
  typedef typename FFTType::DFTGridType DFTGridType;
  typedef typename FFTType::ModesContainerType ModesContainerType;
  typedef capd::vectalg::Container<DFTGridType, ORD> GridsContainerType;
  typedef capd::vectalg::Container<ModesContainerType, ORD> ModesContainerContainerType;
  typedef typename EquationType::PolyBdType PolyBdType;
  typedef PolyBdType ParamType;
  typedef typename FFTType::IndexType IndexType;
  typedef typename FFTType::IndexRangeType IndexRangeType;
  typedef typename FFTType::ScalarType ScalarT; ///< this is actual scalar, ScalarType typename is reserved for enclosure functions
  typedef typename EquationType::ComplexType ComplexScalarType;
  typedef typename EquationType::RealType RealType;
  typedef typename PolyBdType::MatrixType MatrixType;
  typedef RealType ScalarType; ///< ScalarType has to be this, because enclosure functions are using this
  typedef typename PolyBdType::RealContainerType VectorType;

  int m, M, w; ///<w is the constant >= 2 from the paper , which appears in the G set definition, i.e. G(wM)
  int dftPts1; ///<number of discrete points for FFT, such that aliasing is omitted
  int& dftPts; ///<an alias for dftPts1
  FFTType fft1; ///<used for calculating only the finite dimensional part (takes m elements and transforms them)
  int dftPts2; ///<number of discrete points for FFT, such that aliasing is omitted
  FFTType fft2; ///<used for calculating the infinite dimensional part (takes M elements into transform)
  int dftPts3; ///<number of discrete points for FFT, such that aliasing is omitted
  FFTType fft3; ///<used for calculating the infinite dimensional part (takes 2*M elements into transform)

  IndexRangeType ir_m, ir_M, ir_finiteTail, ir_overFft3; ///<this is the range of indices of modes that are not included in the fft3 transform (|k|>2M)

  PolyBdType m_Nt, m_bt, m_gt, tpb; ///<auxiliary variables, remember to initialize them!
  DFTGridType tg1, tg2, tg3, tg1_2, tg2_2, tg3_2, tg1_3, tg2_3, tg3_3, sdft; ///<auxiliary variables, remember to initialize them!
  ModesContainerType tmc, tmc2, tmc3; ///<auxiliary variables, remember to initialize them!

  //below are temporary variables for overestimation optimized FFT
  ModesContainerType delta_u, delta_v, mu, mv, r1;

  DFTGridType r_delta_u, r_abs_mu, r_delta_v, r_abs_mv, r_mu, r_mv, s1, s2, s3,
      s4;

  DFTGridType empty, ///<auxiliary variable
      rhs; ///<auxiliary variable
  bool useFFT;
  ModesContainerType td, tl;
  VectorType yc, ///<this is y_c ([\delta] midpoint)
      forcing; ///<this is an optional forcing
  bool initializedHigherDFT;
  RealType pi;

  std::vector<double> inflates;

  RealType w_1;
  int w_2, w_3; ///< the constants w1, w2, w3 from the paper, depending on the norms

  void initializeW(const int w, RealType& w1, int& w2, int& w3) {
    if (typeid(VNormType)
        == typeid(typename FFTType::IndexRangeType::NormType)) {
      w_1 = 1. - 1. / w;
      w_2 = 1;
      w_3 = 1;
    } else {
      if (typeid(VNormType) == typeid(capd::jaco::EuclideanNorm<IndexType>)
          && typeid(typename FFTType::IndexRangeType::NormType)
              == typeid(capd::jaco::MaximumNorm<IndexType>)) {
        IndexType i;
        w_1 = 1. - sqrt(RealType(i.d())) / w;
        w_2 = 1;
        w_3 = 1;
      } else {
        std::cerr
            << "The is no code implemented for the provided norms ( |.|_v and |.|_n)";
        throw std::runtime_error(
            "The is no code implemented for the provided norms ( |.|_v and |.|_n)");
      }
    }
  }

  /**The constructor for FINITE dimensional integrator - only the projection is taken into account
   */
  DPDE2(int m_, int dftPts1_, RealType nu_, RealType pi_, int order) :
      EquationType(m_, m_, nu_), m(m_), dftPts1(dftPts1_), dftPts(dftPts1),
      fft1(m, dftPts1, pi_), ir_m(IndexType::zero()), ir_M(IndexType::zero()), tpb(m_, m_),
          tg1(dftPts1), tg1_2(dftPts1), sdft(dftPts1), tmc(m), tmc2(m), tmc3(m),
          empty(dftPts1), rhs(dftPts1), useFFT(false), td(m), tl(m),
          yc(PolyBdType::modes2realArraySizeStatic(m_)),
          forcing(PolyBdType::modes2realArraySizeStatic(m_)),
          initializedHigherDFT(false), pi(pi_),
          inflates(PolyBdType::modes2realArraySizeStatic(m_), 1){
    ir_m.setRange(0, strong, m, weak);
    initializeW(w, w_1, w_2, w_3);
  }

  /**The constructor for INFINITE dimensional integrator (differential inclusion, with tails etc...)
   *
   * w is the constant >= 2 from the paper , which appears in the G set definition, i.e. G(wM)
   *
   * dftPts3 is the number of points used by FFT for calculating the polynomial bounds convolutions , this number
   * >= (w + 1)*dftPts2,
   *
   * when dftPts3 == w*dftPts2 terms such that a_k\in G(m) and b_{k-k_1}\nin G(w*M) are missing in the FFT result
   */
  DPDE2(int m_, int M_, int dftPts1_, int dftPts2_, RealType nu_, RealType pi_,
      int order, bool initializeHigherDFT = true, int w_ = 2) :
      EquationType(m_, M_, nu_), m(m_), M(M_), w(w_), dftPts1(dftPts1_), dftPts(
          dftPts1), fft1(m, dftPts1, pi_), dftPts2(dftPts2_), dftPts3(
          (w_ + 1) * dftPts2_), ir_m(IndexType::zero()), ir_M(
          IndexType::zero()), m_Nt(m, M), m_bt(m, M), m_gt(m, M), tpb(m, M), tg1(
          dftPts1), tg2(dftPts2), tg3(dftPts3), tg1_2(dftPts1), tg2_2(dftPts2), tg3_2(
          dftPts3), tg1_3(dftPts1), tg2_3(dftPts2), tg3_3(dftPts3), sdft(
          dftPts1), tmc(m), tmc2(m), tmc3(m), delta_u(m, M, w), delta_v(m, M, w), mu(m,
          M, w), mv(m, M, w), r1(m, M, w), r_delta_u(dftPts3), r_abs_mu(
          dftPts3), r_delta_v(dftPts3), r_abs_mv(dftPts3), r_mu(dftPts3), r_mv(
          dftPts3), s1(dftPts3), s2(dftPts3), s3(dftPts3), s4(dftPts3), empty(
          dftPts1), rhs(dftPts1), useFFT(false), td(m), tl(m), yc(
          PolyBdType::modes2realArraySizeStatic(m_)), forcing(
          PolyBdType::modes2realArraySizeStatic(m_)), initializedHigherDFT(
          initializeHigherDFT), pi(pi_), inflates(PolyBdType::modes2realArraySizeStatic(M_), 1) {
    if (initializeHigherDFT) {
      fft2 = FFTType(M, dftPts2, pi);
      fft3 = FFTType((w + 1) * M, dftPts3, pi);
    }
    initializeW(w, w_1, w_2, w_3);
    ir_m.setRange(0, strong, m, weak);
    updateRange();
  }

  void updateRange() {
    ir_M.setRange(0, strong, M, weak);
    ir_finiteTail.setRange(m, strong, M, weak);
    ir_overFft3.setRange((w + 1) * M, strong, -1, strong);
  }

  int getFinitePart() const {
    return capd::jaco::finitePart;
  }

  virtual DPDE2& getVectorField() {
    return *this;
  }

  inline void scalarProduct(int j, int imj, const DFTGridType& grid1,
      const DFTGridType& grid2, DFTGridType& r) {
    r.multiply(grid1, grid2);
  }

  inline void calculateGrids(int i, const ModesContainerContainerType& modes,
      GridsContainerType& grids) {
    int j;
    for (j = 0; j <= i; ++j) {
      fft1.fastTransform(modes[j], grids[j]);
    }
  }

  /**Procedure for calculating i-th normalized derivative. USES FFT.
   *
   * @param i order for which the calculations are done
   * @param grids an array containing L_2 coefficients of the right-hand side of the PDE
   * @param modes an array containing l_2 coefficients of the right-hand side of the PDE
   * @param rhsSeries l_2 coefficients of the result
   * @param rhsFunctionSpace L_2 coefficients of the result
   * @param calculateRhsFunctionSpace if L_2 coefficients of the calculated rhs should be calculated
   */
  inline void rightHandSide(int i, const GridsContainerType& grids,
      const ModesContainerContainerType& modes, ModesContainerType& rhsSeries,
      DFTGridType& rhsFunctionSpace, bool calculateRhsFunctionSpace = true) {
    int j;
    //if the solution is real valued and even/odd then there are a lot of zeros in the jets, and the additions/multiplication by
    //zeros should be avoided
    if (ScalarT::initialConditionIsRealValued()) ///optimizing thing
      ScalarT::switchToRealValuedL2(); //if solution is real valued then imaginary part of jets is zero
    for (j = 0; j <= (i - 1) / 2; ++j) {
      ///remark: there is already multiplication by k inside this procedure
      ///< calculates \sum_{k_1\in a Projection}{(k, a_k)\cdot a_{k-k_1}}
      scalarProduct(j, i - j, grids[j], grids[i - j], sdft);
      if (j == 0)
        rhs = sdft; //this is needed in order to s have the same DPDEContainer
      else
        rhs += sdft;
    }
    if (i > 0) {
      rhs *= 2;
    }
    if (i % 2 == 0) {
      scalarProduct(i / 2, i / 2, grids[i / 2], grids[i / 2], sdft);
      if (i == 0)
        rhs = sdft;
      else
        rhs += sdft;
    }

    //we switch back to complex valued, because the FFT is complex, regardless the solution is real valued
    ScalarT::switchToComplexValued();
    fft1.fastInverseTransform(rhs, rhsSeries);
    generalDebug2 << "rhsSeries:\n" << rhsSeries << "\n";
    //here we cannot switch to real valued, because multiplication by i (switches zeros re to/from im)
    ComplexScalarType constant = Ncoeff();
    for (j = 0; j < this->m_r; j++)
      constant *= ComplexScalarType::i();
    rhsSeries *= constant;
    if (ScalarT::initialConditionIsRealValued())
      ScalarT::switchToRealValued();
    multiplyComponents(rhsSeries);
    IndexType index;
    for (index = firstModeIndex(ir_m); !index.limitReached(ir_m);
        index.inc(ir_m)) {
      rhsSeries.set(index,
          rhsSeries[index] + lambda_k(index) * modes[i][index]);
    }
    if (i == 0) { //vector field is calculated, and thus y_c and forcing have to be added
      rhsSeries += yc;
      rhsSeries += forcing;
    }
    rhsSeries *= RealType(1) / RealType(i + 1);
    //we switch back to complex valued, because the FFT is complex
    ScalarT::switchToComplexValued();
    if (calculateRhsFunctionSpace)
      fft1.fastTransform(rhsSeries, rhsFunctionSpace);
  }

  /**Procedure for calculating i-th normalized derivative. NOT USING FFT.
   *
   * @param i order for which the calculations are done
   * @param grids an array containing L_2 coefficients of the right-hand side of the PDE
   * @param modes an array containing l_2 coefficients of the right-hand side of the PDE
   * @param rhsSeries l_2 coefficients of the result
   * @param rhsFunctionSpace L_2 coefficients of the result
   * @param calculateRhsFunctionSpace if L_2 coefficients of the calculated rhs should be calculated
   */
  inline void rightHandSide(int i, const ModesContainerContainerType& modes,
      ModesContainerType& rhsSeries) {
    int j;
    for (j = 0; j <= i; ++j) {
      CalculateNonlinearTermDirectly(modes[j], modes[i - j], tmc);
      if (j == 0)
        rhsSeries = tmc;
      else
        rhsSeries += tmc;
    }
    IndexType index;
    for (index = firstModeIndex(ir_m); !index.limitReached(ir_m);
        index.inc(ir_m)) {
      rhsSeries.set(index,
         rhsSeries[index] + lambda_k(index) * modes[i][index]);
    }
    if (i == 0) { //vector field is calculated, and thus y_c and forcing have to be added
      rhsSeries += yc;
      rhsSeries += forcing;
    }
    rhsSeries *= RealType(1) / RealType(i + 1);
  }

  ///Calculates g_k^{+-} = (T0_k^{+-}-b_k^{+-})e^{-\lambda_k h}+b_k^{+-}, version used by enclosure algorithms.
  ///Returns either Re{a_k} or Im{a_k}, depending on which value is at index i in the modes storage.
  inline RealType g_encl(const RealType& h, int i, const RealType& x,
      const RealType& b) {
    RealType g(0);
    RealType ex = exp(h * lambda(i));
    g.setLeftBound(
        (((x.leftBound() - b.leftBound()) * ex) + b.leftBound()).leftBound());
    g.setRightBound(
        (((x.rightBound() - b.rightBound()) * ex) + b.rightBound()).rightBound());
    return g;
  }

  /** Multiplies each component in the provided container by k (k-th component is multiplied by k).
   *
   * Stores the result in r.
   */
  inline void multiplyComponents(const ModesContainerType& pb,
      ModesContainerType& r) const {
    r = pb;
    IndexType i;
    ScalarType piOver_lPower = power(this->piOver_l, this->m_r);

    if (this->m_r > 0) {
      const IndexRangeType& range((pb.infiniteDimensional ? ir_M : ir_m));
      for (i = firstModeIndex(range); !i.limitReached(range); i.inc(range)) {
        r[i] *= ComplexScalarType(power(i[0], this->m_r) * piOver_lPower);
      }
      pb[IndexType::zero()] *= ComplexScalarType(0);
    }
  }

  /** Multiplies each component in the provided container by k (k-th component is multiplied by k). Multiplications are in place
   * (are saved in pb)
   *
   */
  inline void multiplyComponents(ModesContainerType& pb) const {
    IndexType i;
    ScalarType piOver_lPower = power(this->piOver_l, this->m_r);

    const IndexRangeType& range((pb.infiniteDimensional ? ir_M : ir_m));
    if (this->m_r > 0) {
      for (i = firstModeIndex(range); !i.limitReached(range); i.inc(range)) {
        pb[i] *= ComplexScalarType(power(i[0], this->m_r) * piOver_lPower);
      }
      pb[IndexType::zero()] *= ComplexScalarType(0);
    }
  }

  /**Multiplies the far tail by k.
   */
  inline void multiplyComponentsFarTail(ModesContainerType& pb) const {
    ScalarType piOver_lPower = power(this->piOver_l, this->m_r);
    if (pb.infiniteDimensional) {
      setS(pb, s(pb) - this->m_r); //substracts m_r (constant depending on the order of the nonlinearity derivative) from s(pb)
      setC(pb, C(pb) * piOver_lPower);
    }
  }

  /**Bound of infinite convolution.
   */
  inline ComplexScalarType bound(const ModesContainerType& pb,
      const IndexType& k, ModesContainerType& out) const {
    ///TODO: zamienic if d()==1 na IndexType::harmonicSum()
    ComplexScalarType r;
    if (k.d() == 1) {
      RealType C_ = C(pb);
      int s_ = s(pb);
      RealType unit;
      setLeftBound(unit, -1.);
      setRightBound(unit, 1.);
      RealType t = 2. * C_ * C_ * (1. / (2. * s_ - 1.))
          * power(1. / ScalarType((M + k[0]) * (M)), s_ - 0.5) * unit;

      if (!out.baseReZero)
        r.re = t;
      if (!out.baseImZero)
        r.im = t;
      if (!out.baseReZero && !out.baseImZero)
        r *= 2.;
      //if modes which are bounded are complex we multiply the bound by 2, because when zero-centered complex modes are
      //multiplied the diam of the result is twice as the diam of the input.

    } else {
      if (k.d() == 2) {
        ///TODO: hard-coded constants - should depend on the norm used, this is for |.|_n=|.|_\infty
        ///case |k|_n = |k| = |k|_\infty
        int w_3 = 1;

        RealType C_ = C(pb);
        int s_ = s(pb);

        IndexRangeType range;
        range.setK_1Range(w_3 * w * M, weak, -1, weak);
        RealType unit;
        setLeftBound(unit, -1.);
        setRightBound(unit, 1.);

        RealType divisor = k.maxNorm() == 0 ? 1 : power(VNormType::norm(k), s_);

        RealType t = C_ * C_ * (power(2, s_) / divisor)
            * IndexType::template harmonicSumK_1<RealType, IndexRangeType,
                VNormType>(range, s_) * unit;

        if (!out.baseReZero)
          r.re = t;
        if (!out.baseImZero)
          r.im = t;
        if (!out.baseReZero && !out.baseImZero)
          r *= 2.;
        //if modes which are bounded are complex we multiply the bound by 2, because when zero-centered complex modes are
        //multiplied the diam of the result is twice as the diam of the input.

      } else {
        std::cerr << "Bound can be obtained only in 1D and 2D case.\n";
        throw std::runtime_error(
            "Bound can be obtained only in 1D and 2D case.\n");
      }
    }
    out[k] += r;
    return r;
  }

  /**Bound of one dimensional infinite convolution.
   */
  inline ComplexScalarType bound(const ModesContainerType& pb1,
      const ModesContainerType& pb2, const IndexType& k,
      ModesContainerType& out) const {
    ComplexScalarType r;
    if (k.d() == 1) {
      RealType C1 = C(pb1), C2 = C(pb2);
      int s1 = s(pb1), s2 = s(pb2);

      ScalarType unit;
      setLeftBound(unit, -1.);
      setRightBound(unit, 1.);

      RealType t = 2. * C1 * C2 * sqrt(1. / ((2. * s1 - 1.) * (2. * s2 - 1.)))
          * (power(1. / ScalarType(M + k[0]), s1 - 0.5)
              * power(1. / ScalarType(M), s2 - 0.5)
              + power(1. / ScalarType(M + k[0]), s2 - 0.5)
                  * power(1. / ScalarType(M), s1 - 0.5)) * unit;

      if (!out.baseReZero)
        r.re = t;
      if (!out.baseImZero)
        r.im = t;
      if ((!out.baseReZero && !out.baseImZero))
        r *= 2.;
    } else {
      if (k.d() == 2) {
        ///TODO: hard-coded constants - should depend on the norm used, this is for |.|_n=|.|_\infty
        ///case |k|_n = |k| = |k|_\infty

        RealType Ca = C(pb1), Cb = C(pb2);
        int sa = s(pb1), sb = s(pb2), s_ = (sa > sb ? sb : sa), //min(sa, sb)
            sm = (sa > sb ? sa : sb); //max(sa, sb)

        IndexRangeType range;
        range.setK_1Range(w_3 * w * M, weak, -1, weak);
        RealType unit;
        setLeftBound(unit, -1.);
        setRightBound(unit, 1.);

        RealType divisor = k.maxNorm() == 0 ? 1 : power(VNormType::norm(k), s_);
        RealType t = Ca * Cb / power(w_2 * w * M, sm - s_)
            * (power(2, s_) / divisor)
            * IndexType::template harmonicSumK_1<RealType, IndexRangeType,
                VNormType>(range, s_) * unit;

        //std::cout << "harmonic sum=" << IndexType::template harmonicSumK_1<RealType, IndexRangeType>(range, s_) << "\n";

        if (!out.baseReZero)
          r.re = t;
        if (!out.baseImZero)
          r.im = t;
        if (!out.baseReZero && !out.baseImZero)
          r *= 2.;
        //if modes which are bounded are complex we multiply the bound by 2, because when zero-centered complex modes are
        //multiplied the diam of the result is twice as the diam of the input.

      } else {
        std::cerr << "Bound can be obtained only in 1D and 2D case.\n";
        throw std::runtime_error(
            "Bound can be obtained only in 1D and 2D case.\n");
      }
    }
    out[k] += r;
    return r;
  }

  /**Adds a small bound enclosing infinite dimensional sums. This is called only for infinite dimensional PolyBds.
   */
  inline void addBound(const ModesContainerType& pb, ModesContainerType& out,
      int range = full, bool doubleTheRange = false) const {
    IndexType i;
    IndexRangeType r = out.irFull;
    switch (range) {
    case full:
      r = out.irFull;
      break;
    case finitePart:
      r = out.irProjection;
      break;
    case finitePartPlusTail:
      r = out.irProjectionPlusFiniteTail;
      break;
    case redundantRange:
      r = out.irRedundantRange;
      break;
    }
    if (doubleTheRange && (range == finitePart || range == finitePartPlusTail))
      r.doubleK_1Range();
    for (i = firstModeIndex(r); !i.limitReached(r); i.inc(r)) {
      bound(pb, i, out);
    }
  }

  /**Adds a small bound enclosing infinite dimensional sums. This is called only for infinite dimensional PolyBds.
   */
  inline void addBound(const ModesContainerType& pb1,
      const ModesContainerType& pb2, ModesContainerType& out,
      int range = full) const {
    IndexType i;
    ScalarT t(0.);
    IndexRangeType r = out.irFull;
    switch (range) {
    case full:
      r = out.irFull;
      break;
    case finitePart:
      r = out.irProjection;
      break;
    case finitePartPlusTail:
      r = out.irProjectionPlusFiniteTail;
      break;
    case redundantRange:
      r = out.irRedundantRange;
      break;
    }
    for (i = firstModeIndex(out.irFull); !i.limitReached(out.irFull);
        i.inc(out.irFull)) {
      bound(pb1, pb2, i, out);
    }
  }

  /**THIS FUNCTIONS ARE UNNECESSARY WHEN FFT IS USED WITH M_{FFT} > wM + M
   *
   * When the FFT algorithm is used to calculate the convolutions for k>M, then a remainder should be added, which include terms
   * that the FFT doesn't calculate (because include a term |\cdot| > 2M).
   * The terms are as follows
   * \f[
   *  \sum_{2M < k_1 \leq M+k}{\frac{C}{|k_1|^s}|a_{k-k_1}|} + \sum_{-M \leq k_1 < k-2M}{\frac{C}{|k-k_1|^s}|a_{k_1}|}.
   * \f]

   inline void addFFTRemainderRedundantRange(const ModesContainerType& in, ModesContainerType& out) const{
   */

  /**Calculates the polynomial bound for the convolution of two polynomial bounds.
   *
   *
   * @param r has to include calculated redundant modes
   */
  inline void estimatePolynomialBoundForTheInfinitePart(
      const ModesContainerType& pb, ModesContainerType& r) const {
    setS(r, s(pb));
    RealType C_ = C(pb);
    int s_ = s(pb);

    RealType bound(0), max(0), t;
    IndexType index;
    if (index.d() == 1) {
      for (index = firstModeIndex(r.irRedundantRange);
          !index.limitReached(r.irRedundantRange);
          index.inc(r.irRedundantRange)) {
        r.set(index, power(index[0], s_) * r.redundantMode(index));
        if ((t = r.redundantMode(index).normMax()) > max)
          max = t;
      }
      //tailDebug << "estimate1\n" << "max=" << max << "\n";
      //bound for k<=2M is already stored in redundant modes, now we need to calculate bound for k>2M
      RealType A = pb.sumOfNorms();
      //tailDebug << "A: " << A << "\n";
      RealType twoToS = power(2., s_);
      //multiplication by extra 2. , because we are calculating the maximum norm, and |a\cdot b|<= 2|a||b|
      bound =
          2. * C_
              * (C_ * twoToS / ((s_ - 1.) * power(M, s_ - 1))
                  + 0.5 * twoToS * twoToS * C_ / (power(2 * M + 1, s_))
                  + twoToS * A);
      //tailDebug << "estimate1\n" << "bound=" << bound << "\n";
      if (bound > max)
        max = bound;

      ///!add the bound for the infinite part
      max += 2. * C_ * C_ / ((s_ - 1.) * power(M, s_ - 1));
      setC(r, max);
    } else {
      if (index.d() == 2) {
        RealType A = pb.sumOfNorms(), powerW_1 = power(w_1, s_);
        IndexRangeType range;
        range.setK_1Range(w_3 * M, weak, -1, weak);

        for (index = firstModeIndex(r.irRedundantRange);
            !index.limitReached(r.irRedundantRange);
            index.inc(r.irRedundantRange)) {
          r.set(index,
              power(VNormType::norm(index), s_) * r.redundantMode(index));

          if ((t = r.redundantMode(index).normMax()) > max)
            max = t;
        }

        bound = (2. * C_ * A) / powerW_1
            + C_ * C_ * power(2, s_)
                * IndexType::template harmonicSumK_1<RealType, IndexRangeType,
                    VNormType>(range, s_);

        if (bound > max)
          max = bound;

        setC(r, max);
      } else {
        std::cerr
            << "EstimatePolynomialBoundForTheInfinitePart function (DPDE class) is implemented only for 1D and 2D cases.\n";
        throw std::runtime_error(
            "EstimatePolynomialBoundForTheInfinitePart function (DPDE class) is implemented only for 1D and 2D cases.\n");
      }
    }
  }

  /**Calculates the polynomial bound for the convolution of two polynomial bounds.
   *
   * TODO: estimates here should be changed to IndexType::harmonicSumK_1 and IndexType::harmonicSumKmK_1 to support other dimensions
   * then the ifs like if(index.d() == 1) are not needed.
   *
   * @param r has to include calculated redundant modes
   */
  inline void estimatePolynomialBoundForTheInfinitePart(
      const ModesContainerType& pb1, const ModesContainerType& pb2,
      ModesContainerType& r) const {

    setS(r, min(s(pb1), s(pb2)));
    RealType C1 = C(pb1), C2 = C(pb2);
    int s1 = s(pb1), s2 = s(pb2);

    RealType bound(0), max(0), t;

    IndexType index;
    if (index.d() == 1) {
      for (index = firstModeIndex(r.irRedundantRange);
          !index.limitReached(r.irRedundantRange);
          index.inc(r.irRedundantRange)) {
        r.set(index, power(index[0], s(r)) * r.redundantMode(index));
        if ((t = r.redundantMode(index).normMax()) > max)
          max = t;
      }
      //tailDebug << "estimate2\n" << "max=" << max << "\n";
      //bound for k<=2M is already stored in redundant modes, now we need to calculate bound for k>2M
      RealType A1 = pb1.sumOfNorms(), A2 = pb2.sumOfNorms();
      if (s1 == s2) { //simpler case
        RealType twoToS = power(2., s1);
        //multiplication by extra 2. , because we are calculating the maximum norm, and |a\cdot b|<= 2|a||b|
        bound = C1
            * (C2 * twoToS / ((s1 - 1.) * power(M, s1 - 1))
                + 0.5 * twoToS * twoToS * C2 / (power(2 * M + 1, s1)))
            + C2 * twoToS * A1 + C1 * twoToS * A2;
        //tailDebug << "A1: " << A1 << "\n" << "A2: " << A2 << "\n";
      } else { //pb1.m_s != pb2.m_s
        RealType twoToS1 = power(2., s1), twoToS2 = power(2., s2);
        //multiplication by extra 2. , because we are calculating the maximum norm, and |a\cdot b|<= 2|a||b|
        if (s2 < s1) {
          bound = C2
              * (C1 * twoToS2 / ((s1 - 1.) * power(M, s1 - 1))
                  + C1 * twoToS1
                      / ((s2 - 1.) * power(M, s2 - 1)
                          * power(2 * M + 1, s1 - s2))
                  + twoToS1 * twoToS2 * C1 / power(2 * M + 1, s1))
              + 2. * C1 * A2 * twoToS1 + 2. * C2 * A1 * twoToS2;
        } else {
          bound = C1
              * (C2 * twoToS1 / ((s2 - 1.) * power(M, s2 - 1))
                  + C2 * twoToS2
                      / ((s1 - 1.) * power(M, s1 - 1)
                          * power(2 * M + 1, s2 - s1))
                  + twoToS1 * twoToS2 * C2 / power(2 * M + 1, s2))
              + 2. * C1 * A2 * twoToS1 + 2. * C2 * A1 * twoToS2;
        }
        //tailDebug << "A1: " << A1 << "\n" << "A2: " << A2 << "\n";
      }
      //tailDebug << "estimate2\n" << "bound=" << bound << "\n";
      if (bound > max)
        max = bound;

      ///!add the bound for the infinite part
      max += C1 * C2 / ((s1 - 1.) * power(M, s1 - 1))
          + C1 * C2 / ((s2 - 1.) * power(M, s2 - 1));

      setC(r, max);
    } else {
      if (index.d() == 2) {
        int s_ = (s1 > s2 ? s2 : s1), //min(sa, sb)
            sm = (s1 > s2 ? s1 : s2), k = w * M + 1; //max(sa, sb); //TODO: check this (2d convolutions)
        RealType A = pb1.sumOfNorms(), B = pb2.sumOfNorms();
        IndexRangeType range;
        range.setK_1Range(w_3 * M, weak, -1, weak);

        for (index = firstModeIndex(r.irRedundantRange);
            !index.limitReached(r.irRedundantRange);
            index.inc(r.irRedundantRange)) {
          r.set(index,
              power(VNormType::norm(index), s_) * r.redundantMode(index));

          if ((t = r.redundantMode(index).normMax()) > max)
            max = t;
        }

        bound = (C2 * A) / (power(w_1, s2) * power(k, s2 - s_))
            + (C1 * B) / (power(w_1, s1) * power(k, s1 - s_))
            + (C1 * C2 * power(2, s_)
                * IndexType::template harmonicSumK_1<RealType, IndexRangeType,
                    VNormType>(range, s_)) / power(w_2 * M, sm - s_);

        if (bound > max)
          max = bound;

        setC(r, max);

      } else {
        std::cerr
            << "EstimatePolynomialBoundForTheInfinitePart function (DPDE class) is implemented only for 1D and 2D cases.\n";
        throw std::runtime_error(
            "EstimatePolynomialBoundForTheInfinitePart function (DPDE class) is implemented only for 1D and 2D cases.\n");
      }
    }
  }

  /**Calculates the polynomial bound for the convolution of two polynomial bounds.
   * This functions returns a bound for |k| > 2*M (this is not real value, but is easily calculated and then used for heuristics).
   *
   * TODO: estimates here should be changed to IndexType::harmonicSumK_1 and IndexType::harmonicSumKmK_1 to support other dimensions
   * then the ifs like if(index.d() == 1) are not needed.
   */
  inline virtual RealType estimatePolynomialBoundForTheInfinitePart(
      const PolyBdType& pb) const {

    RealType C = pb.farTail.m_c;
    int s = pb.farTail.m_s;
    IndexType index;
    if (index.d() == 1) {
      RealType bound(0);
      //bound for k<=2M is already stored in redundant modes, now we need to calculate bound for k>2M
      RealType A = pb.sumOfNorms();
      RealType twoToS = power(2., s);
      //multiplication by extra 2. , because we are calculating the maximum norm, and |a\cdot b|<= 2|a||b|
      bound = 2. * C
          * (C * twoToS / ((s - 1.) * power(M, s - 1))
              + 0.5 * twoToS * twoToS * C / (power(2 * M + 1, s)) + twoToS * A);
      return bound;
    } else {
      std::cerr
          << "EstimatePolynomialBoundForTheInfinitePart function (DPDE class) is implemented only for one dimension.\n";
      throw std::runtime_error(
          "EstimatePolynomialBoundForTheInfinitePart function (DPDE class) is implemented only for one dimension.\n");
    }
  }

  /**Nonlinear part of the vector field.
   */
  inline virtual void N(const ModesContainerType& in, ModesContainerType& out,
      int range = full) {

    if (useFFT) {
      CalculateNonlinearTermUsingFFT(in, out);
    } else {
      CalculateNonlinearTermDirectly(in, out, range);
    }
  }

  /**Linear part of the vector field.
   */
  inline virtual void L(const ModesContainerType& in,
      ModesContainerType& out) const {
    IndexType i;
    const IndexRangeType& range(ir_m);

    for (i = this->firstModeIndex(range); !i.limitReached(range);
        i.inc(range)) {
      out[i] = this->lambda_k(i) * in[i];
    }
    ///here we dont care about the infinite part, it is handled in the moveParamsFunction
  }

  /**
   * Whole vector field.
   */
  inline void operator()(const ModesContainerType& in, ModesContainerType& out,
      int range = full) {
    L(in, out);

    N(in, tmc, range);
    out += tmc;
    out += yc; ///add y_c
    out += forcing; ///add forcing

    //assign to the output the calculated farTail
    out.farTail = tmc.farTail;
  }

  /**This function performs FFT transform, which is overestimations optimized. In order to reduce overestimations more fft's are required, therefore more parameters are required.
   *
   * TODO: optimize this - where temporary variables should be? mu, delta_u etc...
   */
  inline void transformOptimized(const ModesContainerType& u,
      const ModesContainerType& v, ModesContainerType& r, int aliasingRemoval =
          padding) {

    u.split(mu, delta_u);
    v.split(mv, delta_v);

    fft3.fastTransform(mu, r_mu, aliasingRemoval);
    fft3.fastTransform(delta_u, r_delta_u, aliasingRemoval);

    mu.useAbsValues = true;
    fft3.fastTransform(mu, r_abs_mu, aliasingRemoval);

    fft3.fastTransform(mv, r_mv, aliasingRemoval);
    fft3.fastTransform(delta_v, r_delta_v, aliasingRemoval);
    mv.useAbsValues = true;
    fft3.fastTransform(mv, r_abs_mv, aliasingRemoval);

    s1.multiply(r_mu, r_mv);

    s2.multiply(r_abs_mu, r_delta_v);
    s3.multiply(r_delta_u, r_abs_mv);
    s4.multiply(r_delta_u, r_delta_v);

    fft3.fastInverseTransform(s1, r1, aliasingRemoval);
    r = r1;

    fft3.fastInverseTransform(s2, r1, aliasingRemoval);
    r1 *= ScalarType(-1., 1.);
    r += r1;
    fft3.fastInverseTransform(s3, r1, aliasingRemoval);
    r1 *= ScalarType(-1., 1.);
    r += r1;
    fft3.fastInverseTransform(s4, r1, aliasingRemoval);
    r1 *= ScalarType(-1., 1.);
    r += r1;
  }

  /**This function performs FFT transform, which is overestimations optimized. In order to reduce overestimations more fft's are required, therefore more parameters are required.
   *
   * TODO: optimize this - where temporary variables should be? mu, delta_u etc...
   */
  inline void transformOptimized(const ModesContainerType& u,
      ModesContainerType& r, int aliasingRemoval = padding) {

    u.split(mu, delta_u);

    fft3.fastTransform(mu, r_mu, aliasingRemoval);
    fft3.fastTransform(delta_u, r_delta_u, aliasingRemoval);
    mu.useAbsValues = true;
    fft3.fastTransform(mu, r_abs_mu, aliasingRemoval);

    s1.multiply(r_mu, r_mu);
    s2.multiply(r_abs_mu, r_delta_u);
    s4.multiply(r_delta_u, r_delta_u);

    fft3.fastInverseTransform(s1, r1, aliasingRemoval);
    r = r1;

    fft3.fastInverseTransform(s2, r1, aliasingRemoval);
    r1 *= ScalarType(-2., 2.);
    r += r1;
    fft3.fastInverseTransform(s4, r1, aliasingRemoval);
    r1 *= ScalarType(-1., 1.);
    r += r1;
  }

  inline void calculateConvolution(const ModesContainerType& pb,
      ModesContainerType& out, bool useOptimizedFFT = true) {
    if (!pb.infiniteDimensional) {
      fft1.fastTransform(pb, tg1);
      tg1_2.multiply(tg1, tg1);
      fft1.fastInverseTransform(tg1_2, out);
    } else {
      if (useOptimizedFFT)
        transformOptimized(pb, out, padding);
      else {
        fft3.fastTransform(pb, tg3, capd::jaco::none);
        tg3_2.multiply(tg3, tg3);
        fft3.fastInverseTransform(tg3_2, out, capd::jaco::none);
      }
      addBound(pb, out);
      estimatePolynomialBoundForTheInfinitePart(pb, out);
    }
  }

  inline void calculateConvolution(const ModesContainerType& pb1,
      const ModesContainerType& pb2, ModesContainerType& out,
      bool useOptimizedFFT = true) {
    if (!pb1.infiniteDimensional and !pb2.infiniteDimensional) {
      fft1.fastTransform(pb1, tg1);
      fft1.fastTransform(pb2, tg1_2);
      tg1_3.multiply(tg1, tg1_2);
      fft1.fastInverseTransform(tg1_3, out);
    } else {
      if (useOptimizedFFT)
        transformOptimized(pb1, pb2, out, padding);
      else {
        fft3.fastTransform(pb1, tg3);
        fft3.fastTransform(pb2, tg3_2);
        tg3_3.multiply(tg3, tg3_2);
        fft3.fastInverseTransform(tg3_3, out);
      }
      addBound(pb1, pb2, out);
      estimatePolynomialBoundForTheInfinitePart(pb1, pb2, out);
    }
  }

  inline void CalculateNonlinearTermUsingFFT(const ModesContainerType& in,
      ModesContainerType& out) {
    if (!in.infiniteDimensional) {
      fft1.fastTransform(in, tg1);
      tg1_2.multiply(tg1, tg1);
      fft1.fastInverseTransform(tg1_2, out);
    } else {
      fft3.fastTransform(in, tg3, capd::jaco::none);
      tg3_2.multiply(tg3, tg3);
      fft3.fastInverseTransform(tg3_2, out, capd::jaco::none);
      addBound(in, out);
      estimatePolynomialBoundForTheInfinitePart(in, out);
    }
    multiplyComponents(out);
    if (in.infiniteDimensional)
      multiplyComponentsFarTail(out);
    ComplexScalarType constant = Ncoeff();
    for (int j = 0; j < this->m_r; j++)
      constant *= ComplexScalarType::i();
    out *= constant;
  }

  inline void CalculateNonlinearTermDirectly(const ModesContainerType& in,
      ModesContainerType& out, int range = full, bool mComponents = true) {
    ((DPDEContainer&) out).multiply(in, in);
    if (!in.infiniteDimensional) {
      IndexType k, k_1;
      ScalarT sum;
      IndexRangeType irK_1(IndexType::zero());
      irK_1.setRange(0, weak, m, weak);
      IndexRangeType irK(IndexType::zero());
      irK.setRange(0, weak, m, weak);
      IndexType first = this->firstModeIndex(irK);

      for (k = first; !k.limitReached(irK); k.inc(irK)) {
        irK_1.k = k;
        sum = ComplexScalarType(0.);
        for (k_1 = k / 2, k_1.inc(irK_1, true); !k_1.limitReached(irK_1);
            k_1.inc(irK_1, true)) {
//        for(k_1=this->firstWithinRange(irK_1), k_1.l=k.l; !k_1.limitReached(irK_1); k_1.inc(irK_1, true)){
          sum += in[k_1] * in[k - k_1];
        }
        sum *= 2.; //this line if for (k_1 = k/2... only!
        if (k.isDivisibleBy2()) {
          sum += in[k / 2] * in[k / 2];
        }

        out[k] = sum;
      }
    } else { //PolyBd is infinite dimensional
      if (!out.infiniteDimensional)
        out = in;
      IndexType k, k_1;
      ScalarT sum;
      IndexRangeType irK_1(IndexType::zero());
      irK_1.setK_1orKmK_1Range(0, weak, M, weak);
      IndexRangeType irK(IndexType::zero());
      irK = out.irFull;
      switch (range) {
      case full:
        irK = out.irFull;
        break;
      case finitePart:
        irK = out.irProjection;
        break;
      case finitePartPlusTail:
        irK = out.irProjectionPlusFiniteTail;
        break;
      case redundantRange:
        irK = out.irRedundantRange;
        break;
      }
      for (k = firstModeIndex(irK); !k.limitReached(irK); k.inc(irK)) {
        irK_1.k = k;
        sum = ComplexScalarType(0.);
        for (k_1 = k / 2, k_1.inc(irK_1, true); !k_1.limitReached(irK_1);
            k_1.inc(irK_1, true)) {
//        for(k_1=this->firstWithinRange(irK_1), k_1.l=k.l; !k_1.limitReached(irK_1); k_1.inc(irK_1, true)){
          if (in.isRealValuedOdd())
            sum.value().re -= in[k_1].value().im * in[k - k_1].value().im;
          else if (in.isRealValuedEven())
            sum.value().re += in[k_1].value().re * in[k - k_1].value().re;
          else
            sum += in[k_1] * in[k - k_1];
        }
        sum *= 2.; //this line if for (k_1 = k/2... only!
        if (k.isDivisibleBy2()) {
          if (in.isRealValuedOdd())
            sum.value().re -= in[k / 2].value().im * in[k / 2].value().im;
          else if (in.isRealValuedEven())
            sum.value().re += in[k / 2].value().re * in[k / 2].value().re;
          else
            sum += in[k / 2] * in[k / 2];
        }
        out[k] = sum;
      }

      if (range != redundantRange)
        addBound(in, out);
      estimatePolynomialBoundForTheInfinitePart(in, out);
      if(mComponents)
        multiplyComponentsFarTail(out);
    }
    if(mComponents){
      if (range != redundantRange)
        multiplyComponents(out);
    }

    ComplexScalarType constant( 1. );
    if(mComponents)
      constant *= Ncoeff();

    for (int j = 0; j < this->m_r; j++)
      constant *= ComplexScalarType::i();
    out *= constant; //multiplication by a coefficient in front of the convolution sum
  }

  inline void CalculateNonlinearTermDirectly(const ModesContainerType& in1,
      const ModesContainerType& in2, ModesContainerType& out,
      int range = full) {
    ((DPDEContainer&) out).multiply(in1, in2);
    if (!in1.infiniteDimensional || !in2.infiniteDimensional) {
      IndexType k, k_1;
      ScalarT sum;
      IndexRangeType irK_1(IndexType::zero());
      irK_1.setRange(0, weak, m, weak);
      IndexRangeType irK(IndexType::zero());
      irK.setRange(0, weak, m, weak);
      IndexType first = this->firstModeIndex(irK);

      for (k = first; !k.limitReached(irK); k.inc(irK)) {
        irK_1.k = k;
        sum = ComplexScalarType(0.);
        for (k_1 = this->firstWithinRange(irK_1), k_1.l = k.l;
            !k_1.limitReached(irK_1); k_1.inc(irK_1, true)) {
          sum += in1[k_1] * in2[k - k_1];
        }
        out[k] = sum;
      }
    } else { //PolyBd is infinite dimensional, we calculate the in
      std::cerr
          << "DPDE2.CalculateNonlinearTermDirectly - Forbidden call\nThe function for infinite-dimensional PolyBds is not implemented.\n";
      throw std::runtime_error(
          "DPDE2.CalculateNonlinearTermDirectly - Forbidden call\nThe function for infinite-dimensional PolyBds is not implemented.\n");
      if (range != redundantRange)
        addBound(in1, in2, out);
      estimatePolynomialBoundForTheInfinitePart(in1, in2, out);
      multiplyComponentsFarTail(out);
    }
    if (range != redundantRange)
      multiplyComponents(out);

    ComplexScalarType constant = Ncoeff();
    for (int j = 0; j < this->m_r; j++)
      constant *= ComplexScalarType::i();
    out *= constant;
  }

  inline void CalculatePerturbationsDirectly(const PolyBdType& in,
      VectorType& out) {
    IndexType k, k_1;
    ComplexScalarType sum, a;
    RealType sm;
    IndexRangeType irK_1(IndexType::zero());
    IndexRangeType irK(IndexType::zero());
    irK.setRange(0, strong, m, weak);
    irK_1.setK_1orKmK_1Range(m, strong, M, weak);
    IndexType first = firstModeIndex(irK);
    IndexRangeType infiniteSumRange(k);
    infiniteSumRange.setRange(M, strong, -1, strong);
    for (k = first; !k.limitReached(irK); k.inc(irK)) {
      irK_1.k = k;
      sum = ComplexScalarType(0.);
      //we take into account only terms from the sum that contains at least one mode from the tail
      for (k_1 = k / 2, k_1.inc(irK_1); !k_1.limitReached(irK_1);
          k_1.inc(irK_1)) {
        sum += in[k_1] * in[k - k_1];
      }
      sum *= 2.;
      //now we estimate the infinite series of terms such that both modes are from the far tail
      //(of order C^2/(|M+1|^2s) [-1, 1])
      sm = 4 * C(in) * C(in)
          * sqrt(
              IndexType::template harmonicSumK_1<RealType, IndexRangeType,
                  VNormType>(infiniteSumRange, 2 * s(in))
                  * IndexType::template harmonicSumKmk_1<RealType,
                      IndexRangeType, VNormType>(infiniteSumRange, 2 * s(in)))
          * RealType(-1., 1.);
      if (in.baseReZero || in.baseImZero)
        sum.re += sm;

      tpb[k] = sum;
    }
    multiplyComponents(tpb);
    ComplexScalarType constant = Ncoeff();
    for (int j = 0; j < this->m_r; j++)
      constant *= ComplexScalarType::i();
    tpb *= constant;
    copyFinitePart(tpb, out);
  }

  void perturbations(const PolyBdType& in, VectorType& out) {
    if (useFFT) {
      std::cerr << "DPDE2.perturbations() using FFT is not implemented.\n";
      throw std::runtime_error(
          "DPDE2.perturbations() using FFT is not implemented.\n");
    } else {
      CalculatePerturbationsDirectly(in, out);
    }
  }

  inline MatrixType jacobian(const VectorType& in) {

    ScalarT::switchToLocalOptimization();
    bool f = useFFT;
    if (f)
      useFFT = false; ///calculating the jacobian using FFT is unnecessary (this is calculated fast)
    tmc = in;
    L(tmc, tmc2);
    N(tmc, tmc3);
    tmc2 += tmc3;
    MatrixType m(in.size(), in.size());
    tmc2.monodromyMatrix(m);

    if (f)
      useFFT = true;
    ScalarT::switchToGlobalOptimization();
    return m;
  }

  inline MatrixType jacobianFFT(const VectorType& in) {

    ScalarT::switchToLocalOptimization();
    bool f = useFFT;
    if (!f)
      useFFT = true; ///calculating the jacobian using FFT is unnecessary (this is calculated fast)
    tmc = in;
    L(tmc, tmc2);
    N(tmc, tmc3);
    tmc2 += tmc3;
    MatrixType m(in.size(), in.size());
    tmc2.monodromyMatrix(m);

    if (!f)
      useFFT = false;
    ScalarT::switchToGlobalOptimization();
    return m;
  }

  ///===TAIL MANAGER===
  inline ComplexScalarType b(const IndexType& index,
      const ModesContainerType& N) const {
    return N[index] / -lambda_k(index);
  }
  inline void bt(const ModesContainerType& N, ModesContainerType& bt) const {
    IndexType i;
    bt = N;
    for (i = firstModeIndex(ir_finiteTail); !i.limitReached(ir_finiteTail);
        i.inc(ir_finiteTail)) {
      bt[i] = b(i, N);
    }
    setC(bt, C(N) / V(M));
    setS(bt, s(N) + m_p);
  }

  inline ComplexScalarType g(const RealType& h, const IndexType& index,
      const ModesContainerType& x_0, const ModesContainerType& N) const {
    RealType g_re, g_im;
    ComplexScalarType bk = b(index, N);
    RealType ex = exp(h * lambda_k(index));
    g_re.setLeftBound(
        (((x_0[index].re.leftBound() - bk.re.leftBound()) * ex)
            + bk.re.leftBound()).leftBound());
    g_re.setRightBound(
        (((x_0[index].re.rightBound() - bk.re.rightBound()) * ex)
            + bk.re.rightBound()).rightBound());

    g_im.setLeftBound(
        (((x_0[index].im.leftBound() - bk.im.leftBound()) * ex)
            + bk.im.leftBound()).leftBound());
    g_im.setRightBound(
        (((x_0[index].im.rightBound() - bk.im.rightBound()) * ex)
            + bk.im.rightBound()).rightBound());
    return ComplexScalarType(g_re, g_im);
  }

  inline void gt(const RealType& h, const ModesContainerType& x_0,
      const ModesContainerType& N, ModesContainerType& gt) const {
    gt = x_0; //just to set the proper dimensions etc...
    IndexType i;
    for (i = firstModeIndex(ir_finiteTail); !i.limitReached(ir_finiteTail);
        i.inc(ir_finiteTail)) {
      gt[i] = g(h, i, x_0, N);
    }
    int sb = s(N) + m_p;
    RealType Cb = C(N) / V(M);
    if (i.d() == 1) {
      if (sb > s(x_0)) {
        int r = sb - s(x_0);
        int maximum = maximumPoint(h, r, M);
        if (maximum < M)
          maximum = M;
        setC(gt, C(x_0) * exp(h * lambda_k(maximum)) * power(maximum, r) + Cb);
      } else {
        setC(gt, C(x_0) * exp(h * lambda_k(M)) * power(M, sb - s(x_0)) + Cb);
      }
      setS(gt, sb);
    } else {
      std::cerr << "Function gt is implemented only for one dimension.\n";
      throw std::runtime_error(
          "Function gt is implemented only for one dimension.\n");
    }
  }

  inline void changeS(ModesContainerType& pb, int newS) const {
    int oldS = s(pb);
    setS(pb, newS);
    setC(pb, C(pb) * power(M, s(pb) - oldS));
  }

  inline RealType estimateCN(const PolyBdType& W_2) {
    return estimatePolynomialBoundForTheInfinitePart(W_2);
  }

  inline RealType estimateCG(const RealType& h, const PolyBdType& x_0,
      const ModesContainerType& W_2) {
    RealType CN = estimateCN(W_2), Cb = CN / V(M);
    int sb = s(W_2) + m_p - 1;
    IndexType i;
    if (i.d() == 1) {
      if (sb >= s(x_0)) {
        int r = sb - s(x_0);
        int maximum = maximumPoint(h, r, M);
        if (maximum < M)
          maximum = M;
        return C(x_0) * exp(h * lambda_k(maximum)) * power(maximum, r) + Cb;
      } else {
        return C(x_0) * exp(h * lambda_k(M)) * power(M, sb - s(x_0)) + Cb;
      }
    } else {
      std::cerr << "DPDE2.estimateCG is implemented only for one dimension.\n";
      throw std::runtime_error(
          "DPDE2.estimateCG is implemented only for one dimension.\n");
    }
  }

  inline RealType investigate(const RealType& step, const PolyBdType& x_0,
      const PolyBdType& W_2) {
    ScalarType Cg = estimateCG(step, x_0, W_2);
    int sg = s(W_2) + m_p - 1.;
    if (C(W_2) != 0 && Cg != 0)
      return power(C(W_2) / Cg, RealType(1. / (s(W_2) - sg))).rightBound();
    else
      return -1;
  }

  int potentialM(RealType basis, bool decrease) const {
    if (!decrease)
      return static_cast<int>(leftBound(basis) * __INCREASE_L__ + __L_CONST__);
    else
      return static_cast<int>(leftBound(basis) * __DECREASE_L__
          + __L_CONST_DECREASE__);
  }

  inline int findS(const RealType& step, const ModesContainerType& x_0,
      ModesContainerType& W_2, bool& changedLastStep) {

    //23.10.2013 commented out, doesn't work ok for DPDE3 - keeps s(W_2) constant
    //16.06.2016 fixed initial s  for tail validation to s(T)-1
    changeS(W_2, s(W_2) - this->m_p + this->m_r);

    /*RealType dL = investigate(step, x_0, W_2);
     int currentM = potentialM(dL, false);
     int maxM = M; ///here we do not allow M to grow, it is forced to keep its initial value
     tailDebug << "findS start.\ncurrent dL=" << dL << "\n";
     while(currentM > maxM && s(W_2) > m_sufficientlyLarge) {
     changeS(W_2, s(W_2) - 1.);
     dL = investigate(step, x_0, W_2);
     currentM = potentialM(dL, false);
     tailDebug << "examinating s=" << s(W_2) << "\ncurrent dL=" << dL << "\n";
     }*/
    tailDebug << "findS finish, final s=" << s(W_2) << "\n";
    /*
     #if !__CONSTANT_M__
     //setting potentially good M
     if(currentM > M) {
     if(suitableM(currentM, false)) {
     changeM(currentM, T0, T, changedLastStep);
     }
     }

     if(currentM < M) {
     if(suitableM(currentM, true)) {
     changeM(currentM, T0, T, changedLastStep);
     }
     }
     #endif
     */
    return s(W_2);
  }

  int changeM(int newM, PolyBdType& T0, PolyBdType& T, bool& mWasChanged,
      bool guard = false) {
    tailDebug << "changing M\n";
    mWasChanged = true;
    int val;
    val = T.changeM(newM, guard);
    T0.changeM(val, guard);
    m_Nt.changeM(val, guard);
    m_bt.changeM(val, guard);
    m_gt.changeM(val, guard);
    this->M = val;
    updateRange();
    return val;
  }

  bool updateM(PolyBdType& T, PolyBdType& T0, int& L, RealType& dL,
      int previousL, RealType previousdL, bool& validated,
      bool& farTailValidate, bool& changedLastStep, bool& mWasChanged) {
    int newM = -1;
    if (truncate(leftBound(dL), 5) > truncate(leftBound(previousdL), 5)) {
      int maxM = this->m_factorMaxM * this->m;
      if (L < maxM)
        newM = potentialM(L, false);
      else {
        if (L > 1000) {
          std::cerr << "Blow up! Validate tail function.\n";
          throw capd::jaco::SetBlowUpException<ScalarType>(dL, "validate tail");
        }
        if (L < 2 * this->M) //if L is not much larger than current M we increase M by little
          newM = potentialM(this->M,
              false)/*static_cast<int>(this->M*__INCREASE_L__+__L_CONST__)*/;
        else
          //otherwise if L is much larger than current M we increase M by large, to make it larger than L
          newM = potentialM(L, false);
      }
    } else {
      if (truncate(leftBound(dL)) == truncate(leftBound(previousdL))) {
        if (L < this->M) {
          newM = potentialM(L,
              true)/*static_cast<int>(L*__INCREASE_L__+__L_CONST__)*/;
        } else {
          if (L > this->M) {
            newM = potentialM(L, false);
          }
        }
      }
    }

    tailDebug << "newM=" << newM << "\n";
    tailDebug << "previousdL=" << previousdL << " dL=" << dL << "\n";

    if (!changedLastStep && newM > 0 && previousdL > 0) {
      if (newM < potentialM(this->m, false))
        newM = potentialM(this->m, false);
      tailDebug << "current M=" << this->M << "\n";
      this->changeM(newM, T0, T, mWasChanged);
      std::cout << "M was updated, current M: " << newM << "\n";
      tailDebug << "M was updated, current M: " << newM << "\n";
      tailDebug << "validate of the far tail failed, increasing M, newM="
          << newM << ", close tail size=" << T.dimension() << "\n";
      tailDebug << "tail NOT validated.\n";

      validated = false;
      farTailValidate = false;
      return true;
    } else {
      return false;
    }
  }


  bool updateT(const RealType& h, const PolyBdType& x_0, PolyBdType& W_2,
      const PolyBdType& gt_b, const RealType& newC, int& L, RealType& dL,
      bool& validated, bool& farTailValidate) {
    bool r = false;
    if (!(newC <= C(W_2))) {
      r = true;
      setCLarger(W_2, __D_2__ * newC);


      //31/05/2016 was here before, I think should be outside the loop

      //if (Cg != 0 && C(W_2) != 0) {
      //  pow = power(C(W_2) / Cg, RealType(1. / (s(W_2) - s(gt_b))));
      //  L = static_cast<int>(pow.rightBound()) + 1;
      //  dL = pow.rightBound();
      //} else {
      //  L = -1;
      //  dL = -1;
      //}

      validated = false;
      farTailValidate = false;
    }

    RealType Cg = estimateCG(h, x_0, W_2), pow;
    if (Cg != 0 && C(W_2) != 0) {
      pow = power(C(W_2) / Cg, RealType(1. / (s(W_2) - s(gt_b))));
      L = static_cast<int>(pow.rightBound()) + 1;
      dL = pow.rightBound();
    } else {
      L = -1;
      dL = -1;
    }

    return r;
  }

  template<typename DoubleT>
  DoubleT truncate(DoubleT d, int decimalPlaces = 2) const {
    int i;
    DoubleT factor = 1;
    for (i = 0; i < decimalPlaces; ++i)
      factor *= 10;
    return ::floor(d * factor) / factor;
  }

  bool validateT(const RealType& step, ModesContainerType& x_0,
      ModesContainerType& W_2, bool& validate, bool& farTailValidate,
      RealType& guessedC, int& previousL, RealType& previousdL,
      bool& changedLastStep, bool& firstStep, bool& recalculateN) {
    bool validated = true;
    int first = W_2.finiteTailBegin(), last = W_2.finiteTailEnd();
    int i;
    IndexType index, tindex;
    bool inflate = false;

    bool mWasChanged = false;

    RealType bk, g;
    tailDebug << "current W_2:\n" << W_2 << "\n";
    tailDebug << "current T0:\n" << x_0 << "\n";

    ///===NEAR TAIL VALIDATION===
    {
      N(W_2, m_Nt);
      tailDebug << "Calculated n (Tail):\n" << m_Nt;
      this->bt(m_Nt, m_bt);
      gt(step, x_0, m_Nt, m_gt);

      for (i = first; i <= last; ++i) {
        inflates[i] = 1.;
        bk = m_bt[i];
        g = m_gt[i];
        tailDebug << "Validating i=" << i << " in tail.\n" << "bkRe=" << bk
            << "\n" << "gRe=" << g << "\n";
        inflate = false;
        if (x_0[i].rightBound() < bk.rightBound()) {
          tailDebug << "W_2[i].rightBound() " << W_2[i].rightBound()
              << "<g.rightBound() " << g.rightBound() << " odp="
              << (W_2[i].rightBound() < g.rightBound()) << "\n";
          if (W_2[i].rightBound() < g.rightBound()) {
            validate = false;
            validated = false;
            tailDebug << "validate failed, setting T(i, 1).rightBound="
                << (1 - __D_G__) * g.rightBound() + __D_G__ * bk.rightBound()
                << "\n";
            W_2.setRightBound(i,
                (1 - __D_G__) * g.rightBound() + __D_G__ * bk.rightBound());
            inflate = true;
            recalculateN = true;
          }
        }
        if (x_0[i].leftBound() > bk.leftBound()) {
          tailDebug << "W_2.leftBound() " << W_2[i].leftBound()
              << "<g.leftBound() " << g.leftBound() << " odp="
              << (W_2[i].leftBound() > g.leftBound()) << "\n";
          if (W_2[i].leftBound() > g.leftBound()) {
            validate = false;
            validated = false;
            tailDebug << "validate failed, setting T(i, 1).leftBound="
                << (1 - __D_G__) * g.leftBound() + __D_G__ * bk.leftBound()
                << "\n";
            W_2.setLeftBound(i,
                (1 - __D_G__) * g.leftBound() + __D_G__ * bk.leftBound());
            inflate = true;
            recalculateN = true;
          }
        }
        //20.10.2013 - a special technique of reducing number of validateTail iterations
        if (inflate) {
          index = array2modeIndex(i);
          W_2[i] = capd::jaco::inflate(W_2[i], __INFLATE_C__);
          for (int j = -__INFLATE_RADIUS__; j <= __INFLATE_RADIUS__; j++) {
            for (int k = 0; k < index.d(); k++) {
              tindex = index;
              tindex[k] = index[k] + j;
              if (W_2.irFiniteTail.withinRange(tindex) && j != 0) {
                //inflates[this->mode2array(tindex, true)] *= 1. + (__INFLATE_C__ - 1.) / power(ScalarType(abs(j)), s(m_Nt) - 1).rightBound();
                //variant below works OK for RADIUS = 2
                inflates[this->mode2array(tindex, true)] *= 1.
                    + (__INFLATE_C__ - 1.) / ScalarType(abs(j)).rightBound();
              }
            }
          }
        }
      }
      for (i = first; i <= last; ++i) {
        tailDebug << "inflates[" << i << "]=" << inflates[i] << "\n";
        if (inflates[i] != 0) {
          RealType inflated = capd::jaco::inflate(W_2[i], inflates[i]);
          W_2[i] = inflated;
        }
      }
    }

    ///===FAR TAIL VALIDATION===
    RealType Cb = C(m_bt);
    RealType Cg = C(m_gt);
    int L = -1;
    RealType L2 = -1;
    RealType gk, max;
    RealType dL = -1;

    //recalculateL(T, gt_b, Cg, L, dL);
    if (C(x_0) != 0. && Cb != 0. && s(m_bt) != s(x_0))
      L2 = leftBound(power(Cb / C(x_0), RealType(1. / (s(m_bt) - s(x_0))))); //recalculates L_2 using current values

    if (L2 != L2) {
      tailDebug << "Critical error.\nL2 escaped the range. L2=" << L2 << "\n";
      std::cerr << "Critical error.\nL2 escaped the range. L2=" << L2 << "\n";
      throw std::runtime_error("Critical error.\nL2 escaped the range.\n");
    }

    tailDebug << "Validating far tail\n" << "C(N)=" << C(m_Nt) << "\n"
        << "C(b)=" << Cb << "\n" << "C(g)=" << Cg << "\n" << "C(T0)=" << C(x_0)
        << "\n";
    tailDebug << "C(T)=" << C(W_2) << "\n" << "s(N)=" << s(m_Nt) << "\n"
        << "s(b)=" << s(m_bt) << "\n" << "s(g)=" << s(m_gt) << "\n" << "s(T)="
        << s(W_2) << "\n";
    tailDebug << "s(T0)=" << s(x_0) << "\n" << "L=" << L << "\n" << "dL=" << dL
        << "\n" << "truncated dL=" << truncate(leftBound(previousdL), 2)
        << "\n";
    tailDebug << "previousL=" << previousL << "\n";
//    if(L>=M) tailDebug<<"T(L,1)="<<W_2(L, 1)<<"\n";
    tailDebug << "L2=" << L2 << "\n" << "previousdL=" << previousdL << "\n"
        << "changedLastStep=" << changedLastStep << "\n";

    RealType newC;
    farTailValidate = true;
    if (s(x_0) > s(m_bt)) {
      if (!(C(x_0) > Cb * power(M, s(x_0) - s(m_bt)))) { //C(T0)<=C(b)(M)^(s(T0)-s(b))
        newC = Cg * power(M, s(W_2) - s(m_gt));
        tailDebug << "case s(T0)>s(b), T0_{M}<=b_{M}\n" << "newC="
            << __D_2__ * newC << "\n";
        if (updateT(step, x_0, W_2, m_gt, newC, L, dL, validated,
            farTailValidate)) {
#if !__CONSTANT_M__
          updateM(W_2, x_0, L, dL, previousL, previousdL, validated,
              farTailValidate, changedLastStep, mWasChanged);
#else
          tailDebug<<"previousdL="<<previousdL<<" dL="<<dL<<"\n";
#endif
        }
      }
      if (!(C(x_0) <= Cb * power(M, s(x_0) - s(m_bt)))) { //C(T0)>C(b)(M)^(s(T0)-s(b))
        //b_{M}<T(0)_{M}
        if (L2 < M && L2 > 0) {
          std::cerr << "Case s(T0)>s(b), T0_{M}>b_{M}\n";
          std::cerr << "Unexpected L2 " << L2 << "<=M " << M << "!\n";
          throw std::runtime_error("Unexpected L2<=M!\n");
        }
        //we check if T_i<T(0)_i for i=M,...,L_2
        //first case T0 is decreasing faster than T, it is enough to check for i=M
        bool t;
        if (s(x_0) >= s(W_2)) {
          newC = C(x_0) * power(M, s(W_2) - s(x_0));
          t = updateT(step, x_0, W_2, m_gt, newC, L, dL, validated,
              farTailValidate);
        }
        //otherwise we check for i=L_2, T is decreasing faster than T0
        else {
          if (L2 > 0) {
            newC = C(x_0) * power(L2, s(W_2) - s(x_0));
            t = updateT(step, x_0, W_2, m_gt, newC, L, dL, validated,
                farTailValidate);
          } else {
            std::cerr << "Case s(b)>s(T0), T0_{M}<b_{M}\n";
            std::cerr << "Unexpected values of L2=infty and s(T0)<s(T).\n";
            std::cerr << "with s(T0)<s(T) validation is not possible.\n";
            throw std::runtime_error(
                "Unexpected L2<=M!\nWith s(T0)<s(T) validation is not possible.\n");
          }
        }
        tailDebug << "case s(T0)>s(b), T0_{M}>b_{M}\nfirst step\nnewC="
            << __D_2__ * newC << "\nt=" << t << "\n";
        if (L2 > 0) {
          newC = Cg * power(L2, s(W_2) - s(m_gt));
          tailDebug << "case s(T0)>s(b), T0_{M}>b_{M}\nsecond step\nnewC="
              << __D_2__ * newC << "\n";
          t = (t
              || updateT(step, x_0, W_2, m_gt, newC, L, dL, validated,
                  farTailValidate));
        }
        if (t) {
          recalculateN = false;
#if !__CONSTANT_M__
          updateM(W_2, x_0, L, dL, previousL, previousdL, validated,
              farTailValidate, changedLastStep, mWasChanged);
#else
          tailDebug<<"previousdL="<<previousdL<<" dL="<<dL<<"\n";
#endif
        }
      }
    }
    if (s(m_bt) == s(x_0)) {
      ///we check the condition !(C(T0)<C(b)) instead of C(T0)>=C(b) because both may be the case
      ///(if intervals C(T0) and C(b) overlap
      if (!(C(x_0) < Cb)) {      //C(T0)>=C(b)
        newC = C(x_0) * power(M, s(W_2) - s(x_0));
        tailDebug << "case s(b)==s(T0), C(T0)>=C(b)\n" << "newC="
                    << __D_2__ * newC << "\n";
        updateT(step, x_0, W_2, m_gt, newC, L, dL, validated, farTailValidate);


      }

      if (!(C(x_0) >= Cb)) {      //C(T0)<C(b)
        if (s(m_gt) < s(W_2)) {
          std::cerr << "Case s(b)==s(T0), C(T0)<C(b)\n";
          std::cerr << "Unexpected s(g).\n";
          std::cerr << "With s(g)<s(T) validation is not possible.\n";
          throw std::runtime_error(
              "Unexpected s(T)!\nWith s(g)<s(T) validation is not possible.\n");
        }
        newC = Cg * power(M, s(W_2) - s(m_gt));
        tailDebug << "case s(b)==s(T0), C(T0)<C(b)\n" << "newC="
            << __D_2__ * newC << "\n";
        if (updateT(step, x_0, W_2, m_gt, newC, L, dL, validated,
            farTailValidate)) {
          recalculateN = false;
#if !__CONSTANT_M__
          updateM(W_2, x_0, L, dL, previousL, previousdL, validated,
              farTailValidate, changedLastStep, mWasChanged);
#else
          tailDebug<<"previousdL="<<previousdL<<" dL="<<dL<<"\n";
#endif
        }
      }
    }

    if (s(m_bt) > s(x_0)) {
      if (!(C(x_0) < Cb * power(M, s(x_0) - s(m_bt)))) {
        newC = C(x_0) * power(M, s(W_2) - s(x_0));
        tailDebug << "case s(b)>s(T0), T0_{M}>=b_{M}\n" << "newC="
            << __D_2__ * newC << "\n";
        updateT(step, x_0, W_2, m_gt, newC, L, dL, validated, farTailValidate);
      }
      if (!(C(x_0) >= Cb * power(M, s(x_0) - s(m_bt)))) {
        if (L2 < M && L2 > 0) {      //L2<\infty
          std::cerr << "Case s(b)>s(T0), T0_{M}<b_{M}\n";
          std::cerr << "Unexpected L2 " << L2 << "<=M " << M << "!\n";
          throw std::runtime_error("Unexpected L2<=M!\n");
        }
        bool t;
        if (s(m_gt) >= s(W_2)) {
          newC = Cg * power(M, s(W_2) - s(m_gt));
          t = updateT(step, x_0, W_2, m_gt, newC, L, dL, validated,
              farTailValidate);
        }
        //otherwise we check for i=L_2, T is decreasing faster than g
        else {
          if (L2 > 0) {      //in fact we check if L2<\infty
            newC = Cg * power(L2, s(W_2) - s(m_gt));
            t = updateT(step, x_0, W_2, m_gt, newC, L, dL, validated,
                farTailValidate);
          } else {
            std::cerr << "Case s(b)>s(T0), T0_{M}<b_{M}\n";
            std::cerr << "Unexpected values of L2=infty and s(g)<s(T).\n";
            std::cerr
                << "Equation is probably not dissipative, reason: p<=q.\n";
            throw std::runtime_error(
                "Unexpected L2<=M!Equation is probably not dissipative.\n");
          }
        }
        tailDebug << "case s(b)>s(T0), T0_{M}<b_{M}\nfirst step\nt=" << t
            << "\n" << "newC=" << __D_2__ * newC << "\n";
        if (L2 > 0) {      //L2<\infty
          newC = C(x_0) * power(L2, s(W_2) - s(x_0));
          tailDebug << "case s(b)>s(T0), T0_{M}<b_{M}\nsecond step\n" << "newC="
              << __D_2__ * newC << "\n";
          t = (t
              || updateT(step, x_0, W_2, m_gt, newC, L, dL, validated,
                  farTailValidate));
        }
        if (t) {
          recalculateN = false;
#if !__CONSTANT_M__
          tailDebug << "L=" << L << "\ndL=" << dL << "\nC(T)/Cg="
              << power(C(W_2) / Cg, 1. / (s(W_2) - s(m_gt))) << "\n";
          updateM(W_2, x_0, L, dL, previousL, previousdL, validated,
              farTailValidate, changedLastStep, mWasChanged);
#else
          tailDebug<<"previousdL="<<previousdL<<" dL="<<dL<<"\n";
#endif
        }
      }
    }
    changedLastStep = mWasChanged;
    previousdL = dL;
    previousL = L;
    ///near and far tail validated
    if (validated)
      tailDebug << "tail VALIDATED.\n";
    else
      tailDebug << "tail NOT validated.\n";

    return validated;
  }

  ///===END TAIL MANAGER===

  inline int dimension() const {
    return PolyBdType::modes2realArraySizeStatic(m);
  }

  using EquationType::m_sufficientlyLarge;
  using EquationType::nu;
  using EquationType::array2modeIndex;
  using EquationType::m_p;
  using EquationType::V;
  using EquationType::lambda;
  using EquationType::lambda_k;
  using EquationType::Ncoeff;
  using EquationType::firstModeIndex;
  using EquationType::firstWithinRange;
  using EquationType::maximumPoint;
};

/** This is a class representing a nonlinear PDE (this is the input to the integrator). It is the base for other classes representing
 * pdes consisting higher order nonlinear term.
 *
 * Elliptic part is a linear operator, and nonlinear part is the third degree polynomial and its derivatives
 * u^3
 * (u^3)_x
 * (u^3)_{xx}
 * etc...
 *
 * Elliptic part dominates the nonlinear part, in sense the maximal order of the derivative in the elliptic operator is larger than the
 * order of the derivative of u^3 in the nonlinear part.
 * * In particular the value m_r in EquationT defines the order of the derivative of u^3 in the nonlinear part (>= 0).
 *
 * EquationT class defines what specific equation this is. Like for example what are eigenvalues and what is constant in front
 * of the nonlinear term, and what is the order of the polynomial in the nonlinearity.
 *
 * FFTT is a FFT class.
 *
 * ORD is the maximal order.
 *
 * When VNormT is set to EuclideanNorm the program intentionally will not compile (because of efficiency reasons - the
 * function norm is not returning int anymore, it must return interval)  - but the case of EculideanNorm was tested too
 *
 * This works only for one dimensional case, as some optimizations are performed, exclusively for this case.
 */
template<class EquationT, class FFTT, int ORD,
    class VNormT = capd::jaco::MaximumNorm<typename FFTT::IndexType> >
class DPDE3: public DPDE2<EquationT, FFTT, ORD> {
public:
  typedef DPDE2<EquationT, FFTT, ORD, VNormT> BaseClass;
  typedef VNormT VNormType;
  typedef typename BaseClass::RealType RealType;
  typedef typename BaseClass::GridsContainerType GridsContainerType;
  typedef typename BaseClass::ModesContainerContainerType ModesContainerContainerType;
  typedef typename BaseClass::ModesContainerType ModesContainerType;
  typedef typename BaseClass::DFTGridType DFTGridType;
  typedef typename BaseClass::ScalarT ScalarT;
  typedef typename BaseClass::ComplexScalarType ComplexScalarType;
  typedef typename BaseClass::IndexType IndexType;
  typedef typename BaseClass::IndexRangeType IndexRangeType;
  typedef typename BaseClass::VectorType VectorType;
  typedef typename BaseClass::PolyBdType PolyBdType;
  typedef typename BaseClass::MatrixType MatrixType;

  GridsContainerType convolutions; ///<this is for storing in an array convolutions that are used as auxiliary values for calculating of the third order polynomial
  PolyBdType tpb2;      ///<auxiliary variable
  PolyBdType tpb3;      ///<auxiliary variable
  ModesContainerType tmcN;      ///<auxiliary variable
  DFTGridType tg3_3;      ///<auxiliary variable
  DFTGridType tg3_4;      ///<auxiliary variable
  DFTGridType tg3_5;      ///<auxiliary variable
  ModesContainerType tmc;      ///<auxiliary variable
  ModesContainerType tmc2;      ///<auxiliary variable
  ModesContainerType tmc3;      ///<auxiliary variable
  ModesContainerType pm1;      ///<auxiliary variable
  ModesContainerType pm2;      ///<auxiliary variable


  /**The constructor for FINITE dimensional integrator - only the projection is taken into account
   */
  DPDE3(int m_, int dftPts1_, RealType nu_, RealType pi_, int order) :
      BaseClass(m_, dftPts1_, nu_, pi_, order), convolutions(order + 1), tpb2(m_, m_),
      //tpb3(m_, m_), tmcN(m_), tmc(2 * m_), tmc2(2 * m_), tmc3(m_) {
      tpb3(m_, m_), tmcN(2 * m_), tmc( m_), tmc2( m_), tmc3(m_), pm1(m_), pm2(m_) {
    int i;
    for (i = 0; i <= order; ++i) {
      convolutions[i] = DFTGridType(dftPts1_);
    }
    ir_m.setRange(0, weak, m, weak);
  }

  /**The constructor for INFINITE dimensional integrator (differential inclusion, with tails etc...)
   *
   * w is the constant >= 2 from the paper , which appears in the G set definition, i.e. G(wM)
   *
   * dftPts3 is the number of points used by FFT for calculating the polynomial bounds convolutions , this number
   * >= (w + 1)*dftPts2,
   *
   * when dftPts3 == w*dftPts2 terms such that a_k\in G(m) and b_{k-k_1}\nin G(w*M) are missing in the FFT result
   */
  DPDE3(int m_, int M_, int dftPts1_, int dftPts2_, RealType nu_, RealType pi_,
      int order, bool initializeHigherDFT = true, int w_ = 2) :
      BaseClass(m_, M_, dftPts1_, dftPts2_, nu_, pi_, order,
          initializeHigherDFT, w_), convolutions(order + 1), tpb2(m, M), tpb3(m,
          M), tmcN(2 * m_), tg3_3(2 * dftPts2_), tg3_4(2 * dftPts2_), tg3_5(
          //2 * dftPts2_), tmc(2 * m_), tmc2(2 * m_), tmc3(m_) {
          2 * dftPts2_), tmc( m_), tmc2( m_), tmc3(m_), pm1(m_), pm2(m_)  {
    int i;
    for (i = 0; i <= order; ++i) {
      convolutions[i] = DFTGridType(dftPts1_);
    }

    ir_m.setRange(0, weak, m, weak);
  }

  DPDE3& getVectorField() {
    return *this;
  }

  inline void calculateGrids(int i, const ModesContainerContainerType& modes,
      GridsContainerType& grids) {
    int j;
    for (j = 0; j <= i; ++j) {
      fft1.fastTransform(modes[j], grids[j]);
      calculateConvolution(j, grids);
    }
  }

  inline void calculateConvolution(int i, const GridsContainerType& grids) {
    int j;
    for (j = 0; j <= (i - 1) / 2; ++j) {
      ///remark: there is already multiplication by k inside this procedure
      ///< calculates \sum_{k_1\in a Projection}{(k, a_k)\cdot a_{k-k_1}}
      sdft.multiply(grids[j], grids[i - j]);
      if (j == 0)
        rhs = sdft; //this is needed in order to s have the same DPDEContainer
      else
        rhs += sdft;
    }
    if (i > 0)
      rhs *= 2;
    if (i % 2 == 0) {
      sdft.multiply(grids[i / 2], grids[i / 2]);
      if (i == 0)
        rhs = sdft;
      else
        rhs += sdft;
    }
    convolutions[i] = rhs;
  }

  inline void rightHandSide(int i, const GridsContainerType& grids,
      const ModesContainerContainerType& modes, ModesContainerType& rhsSeries,
      DFTGridType& rhsFunctionSpace, bool calculateRhsFunctionSpace = true) {
    //if the solution is real valued and even/odd then there are a lot of zeros in the jets, and the additions/multiplication by
    //zeros should be avoided
    if (grids[0].isRealValued()) ///optimizing thing
      ScalarT::switchToRealValuedL2(); //if solution is real valued then imaginary part of jets is zero
    calculateConvolution(i, grids);
    int j;
    for (j = 0; j <= i; ++j) {
      sdft.multiplyOnly(convolutions[j], grids[i - j]);
      if (j == 0)
        rhs = sdft;
      else
        rhs += sdft;
    }

    //we switch back to complex valued, because the FFT is complex, regardless if solution is real valued
    if (grids[0].isRealValued()) ///optimizing thing
      ScalarT::switchToComplexValued();
    fft1.fastInverseTransform(rhs, rhsSeries);

    //here we cannot switch to real valued, because multiplication by i (switches zeros re to/from im)
    ComplexScalarType constant = Ncoeff();
    for (int j = 0; j < this->m_r; j++)
      constant *= ComplexScalarType::i();
    rhsSeries *= constant;

    if (rhsSeries.isRealValued())
      ScalarT::switchToRealValued();
    this->multiplyComponents(rhsSeries);

    IndexType index;
    for (index = firstModeIndex(ir_m); !index.limitReached(ir_m);
        index.inc(ir_m)) {
      rhsSeries.set(index,
          rhsSeries[index] + lambda_k(index) * modes[i][index]);
    }
    if (i == 0) { //vector field is calculated, and thus y_c and forcing have to be added
      rhsSeries += yc;
      rhsSeries += forcing;
    }
    rhsSeries *= RealType(1) / RealType(i + 1);
    //we switch back to complex valued, because the FFT is complex
    if (rhsSeries.isRealValued())
      ScalarT::switchToComplexValued();
    if (calculateRhsFunctionSpace)
      fft1.fastTransform(rhsSeries, rhsFunctionSpace);
  }

  ///special function for FN equation (calculates the linear part [-v, epsilon*u]
  inline void FNlinear(const ModesContainerType& modes, ModesContainerType& rhsSeries) {

    const double EPSILON = _EPSILON_;

    IndexType k, l;
    int range = modes.n / 2;

    IndexRangeType irK(IndexType::zero());
    irK.setRange(0, weak, range, weak);
    ScalarT sc;

    IndexType first = this->firstModeIndex(irK);
    for (k = first; !k.limitReached(irK); k.inc(irK)) {
      l = k;
      l[0] += range+1;

      rhsSeries[k] = -1. * modes[l];
      //rhsSeries[l] = EPSILON * (modes[k] - modes[l]);
      rhsSeries[l] = EPSILON * ( modes[k] - _GAMMA_ * modes[l] );

    }

  }

  ///direct calculation of the convolution, if Equation is (FitzHugh-Nagumo)
  ///when calculating nonlinearity, projects the vector onto u part
  inline void rightHandSide(int i, const ModesContainerContainerType& modes,
      ModesContainerType& rhsSeries) {
    int j, k;
    //case of third degree polynomial, two loops have to be used
    for (j = 0; j <= i; ++j) {
      for (k = 0; k <= j; ++k) {

        CalculateConvolutionDirectly(modes[k], modes[j - k], tmc);
        if (k == 0)
          tmc2 = tmc;
        else
          tmc2 += tmc;
      }
      CalculateConvolutionDirectly(tmc2, modes[i - j], tmc3);

      if (j == 0){
        rhsSeries = tmc3;

        //rhsSeries = tmc3;
        //tmc2 *= - RealType(1. + _A_);
        //rhsSeries += tmc2; //!!!!!! FN modification

      }else{
        rhsSeries += tmc3;

        //rhsSeries += tmc3;
        //tmc2 *= - RealType(1. + _A_);
        //rhsSeries += tmc2; //!!!!!! FN modification

      }
    }

    this->multiplyComponents(rhsSeries);

    ComplexScalarType constant = Ncoeff();
    for (int j = 0; j < this->m_r; j++)
      constant *= ComplexScalarType::i();
    rhsSeries *= constant;

    if( this->uproject ){
      FNlinear(modes[i], tmc);
    }


    IndexType index;
    for (index = firstModeIndex(ir_m); !index.limitReached(ir_m); index.inc(ir_m)) {
      //rhsSeries.set(index,  rhsSeries[index] + lambda_k(index) * modes[i][index] + tmc[index] );
      rhsSeries.set(index,  rhsSeries[index] + lambda_k(index) * modes[i][index] );
    }

    if (i == 0) { //vector field is calculated, and thus y_c and forcing have to be added
      rhsSeries += yc;
      rhsSeries += forcing;
    }
    rhsSeries *= RealType(1) / RealType(i + 1);
  }

  /**Works for real-valued solutions only.
   */
  inline void CalculateConvolutionDirectly(const ModesContainerType& in,
      ModesContainerType& out, int range = full, bool doubleTheRange = false) {
    ((DPDEContainer&) out).multiply(in, in);
    if (!in.infiniteDimensional) {
      IndexType k, k_1;
      ScalarT sum;
      IndexRangeType irK_1(IndexType::zero());
      irK_1.setRange(0, weak, in.n, weak);
      IndexRangeType irK(IndexType::zero());
      irK.setRange(0, weak, out.n, weak);
      IndexType first = this->firstModeIndex(irK);

      for (k = first; !k.limitReached(irK); k.inc(irK)) {
        irK_1.k = k;
        sum = ComplexScalarType(0.);
        for (k_1 = k / 2, k_1.inc(irK_1, true); !k_1.limitReached(irK_1);
            k_1.inc(irK_1, true)) {
          //        for(k_1=this->firstWithinRange(irK_1), k_1.l=k.l; !k_1.limitReached(irK_1); k_1.inc(irK_1, true)){
          sum += in[k_1] * in[k - k_1];
        }
        sum *= 2.; //this line if for (k_1 = k/2... only!
        if (k.isDivisibleBy2()) {
          sum += in[k / 2] * in[k / 2];
        }
        out.set(k, sum);
      }
    } else { //PolyBd is infinite dimensional
      if (!out.infiniteDimensional)
        out = in;
      IndexType k, k_1;
      ScalarT sum;
      IndexRangeType irK_1(IndexType::zero());
      irK_1.setK_1orKmK_1Range(0, weak, M, weak);
      IndexRangeType irK(IndexType::zero());
      irK = out.irFull;
      switch (range) {
      case full:
        irK = out.irFull;
        break;
      case finitePart:
        irK = out.irProjection;
        break;
      case finitePartPlusTail:
        irK = out.irProjectionPlusFiniteTail;
        break;
      case redundantRange:
        irK = out.irRedundantRange;
        break;
      }
      if (doubleTheRange
          && (range == finitePart || range == finitePartPlusTail))
        irK.doubleK_1Range();

      for (k = firstModeIndex(irK); !k.limitReached(irK); k.inc(irK)) {
        irK_1.k = k;
        sum = ComplexScalarType(0.);
        for (k_1 = k / 2, k_1.inc(irK_1, true); !k_1.limitReached(irK_1);
            k_1.inc(irK_1, true)) {
          //        for(k_1=this->firstWithinRange(irK_1), k_1.l=k.l; !k_1.limitReached(irK_1); k_1.inc(irK_1, true)){
          if (in.isRealValuedOdd()) {
            sum.value().re -= in[k_1].value().im * in[k - k_1].value().im;
          } else {
            if (in.isRealValuedEven()) {
              sum.value().re += in[k_1].value().re * in[k - k_1].value().re;
            } else {
              sum += in[k_1] * in[k - k_1];
            }
          }
        }
        sum *= 2.; //this line if for (k_1 = k/2... only!
        if (k.isDivisibleBy2()) {
          if (in.isRealValuedOdd()) {
            sum.value().re -= in[k / 2].value().im * in[k / 2].value().im;
          } else {
            if (in.isRealValuedEven()) {
              sum.value().re += in[k / 2].value().re * in[k / 2].value().re;
            } else {
              sum += in[k / 2] * in[k / 2];
            }
          }
        }
        out[k] = sum;
      }
      if (range != redundantRange)
        addBound(in, out, range, doubleTheRange);
      estimatePolynomialBoundForTheInfinitePart(in, out);
    }
  }


  /**
   * Works for real-valued solutions only.
   */
  inline void CalculateConvolutionDirectly(const ModesContainerType& in1,
      const ModesContainerType& in2, ModesContainerType& out,
      int range = full) {
    ((DPDEContainer&) out).multiply(in1, in2);
    if (!in1.infiniteDimensional || !in2.infiniteDimensional) {

      IndexType k, k_1;
      ScalarT sum;
      IndexRangeType irK_1(IndexType::zero());

      ///15.04.16 added limiting the range (calculated only for part of the vector where u is)
      int range1 = in1.n,
          range2 = in2.n,
          rangeout = out.n;
      if( this->uproject ){
        range1 = in1.n / 2;
        range2 = in2.n / 2;
        rangeout = out.n / 2;
      }

      irK_1.setK_1Range(0, weak, range1, weak);
      irK_1.setKmk_1Range(0, weak, range2, weak);

      IndexRangeType irK(IndexType::zero());
      irK.setRange(0, weak, rangeout, weak);

      IndexType first = this->firstModeIndex(irK);
      for (k = first; !k.limitReached(irK); k.inc(irK)) {
        irK_1.k = k;
        sum = ComplexScalarType(0.);
        for (k_1 = this->firstWithinRange(irK_1), k_1.l = k.l;
            !k_1.limitReached(irK_1); k_1.inc(irK_1, true)) {
          sum += in1[k_1] * in2[k - k_1];
        }
        out.set(k, sum);
      }


    } else { //PolyBd is infinite dimensional, we calculate the in
      if (!out.infiniteDimensional) {
        if (in1.infiniteDimensional)
          out = in1;
        else if (in2.infiniteDimensional)
          out = in2;
      }

      IndexType k, k_1;
      ScalarT sum;
      IndexRangeType irK_1(IndexType::zero());
      irK_1.setK_1orKmK_1Range(0, weak, M, weak);
      IndexRangeType irK(IndexType::zero());
      irK = out.irFull;
      switch (range) {
      case full:
        irK = out.irFull;
        break;
      case finitePart:
        irK = out.irProjection;
        break;
      case finitePartPlusTail:
        irK = out.irProjectionPlusFiniteTail;
        break;
      case redundantRange:
        irK = out.irRedundantRange;
        break;
      }
      for (k = firstModeIndex(irK); !k.limitReached(irK); k.inc(irK)) {
        irK_1.k = k;
        sum = ComplexScalarType(0.);
        for (k_1 = this->firstWithinRange(irK_1), k_1.l = k.l;
            !k_1.limitReached(irK_1); k_1.inc(irK_1, true)) {
          if (in1.isRealValuedOdd() && in2.isRealValuedOdd()) {
            sum.value().re -= in1[k_1].value().im * in2[k - k_1].value().im;
          } else {
            if (in1.isRealValuedEven() && in2.isRealValuedEven()) {
              sum.value().re += in1[k_1].value().re * in2[k - k_1].value().re;
            } else {
              sum += in1[k_1] * in2[k - k_1];
            }
          }
        }
        out[k] = sum;
      }

      if (range != redundantRange)
        addBound(in1, in2, out, range);
      estimatePolynomialBoundForTheInfinitePart(in1, in2, out);
    }
  }

  /**Adds a bound to the calculated L_2 coefficients, which is a contribution of the far modes, i.e. a_k: |k|>2M. This
   * assures that the result is rigorous, and all modes from the tail have been taken into account
   * \f[
   *   \sum_{|k|>2M}{\frac{C}{|k|^s}e^{ikx_j}}\subset\left[-2\sum_{k>2M}{\frac{C}{|k|^s}},2\sum_{k>2M}{\frac{C}{|k|^s}}\right]
   * \f]
   */
  inline void addL2Bound(const ModesContainerType& in, DFTGridType& out) const {
    RealType harmonicSum = IndexType::template harmonicSumK_1<RealType,
        IndexRangeType, VNormType>(ir_overFft3, s(in)), t(
        -rightBound(2 * C(in) * harmonicSum),
        rightBound(2 * C(in) * harmonicSum));
    ComplexScalarType r;
    if (!out.baseReZero)
      r.re = t;
    if (!out.baseImZero)
      r.im = t;
    out += r;
  }


  inline void CalculateNonlinearTermUsingFFT(const ModesContainerType& in,
      ModesContainerType& out) {
    if (!in.infiniteDimensional) {
      fft1.fastTransform(in, tg1);
      tg1_2.multiply(tg1, tg1, tg1);
      fft1.fastInverseTransform(tg1_2, out);
    } else {
      std::cerr
          << "GLDPDE.CalculateNonlinearTermUsingFFT is not implemented.\n";
      throw std::runtime_error(
          "GLDPDE.CalculateNonlinearTermUsingFFT is not implemented.\n");
    }
    this->multiplyComponents(out);
    if (in.infiniteDimensional)
      this->multiplyComponentsFarTail(out);
    ComplexScalarType constant = Ncoeff();
    for (int j = 0; j < this->m_r; j++)
      constant *= ComplexScalarType::i();
    out *= constant;
  }

  /**Calculates the polynomial bound for the convolution of two polynomial bounds.
   * Variant for u^3, i.e. given bound C(u)/|k|^s(u) returns the bound C(u^3)/|k|^s(u^3).
   * This functions returns a bound for |k| > 2*M (this is not real value, but is easily calculated and then used for heuristics).
   *
   * TODO: estimates here should be changed to IndexType::harmonicSumK_1 and IndexType::harmonicSumKmK_1 to support other dimensions
   * then the ifs like if(index.d() == 1) are not needed.
   */
  inline virtual RealType estimatePolynomialBoundForTheInfinitePart(
      const PolyBdType& u) const {

    RealType C = u.farTail.m_c;
    int s = u.farTail.m_s;
    IndexType index;
    if (index.d() == 1) {
      RealType bound(0), boundSquare(0);
      //bound for k<=2M is already stored in redundant modes, now we need to calculate bound for k>2M
      RealType A = u.sumOfNorms();
      RealType twoToS = power(2., s);
      //multiplication by extra 2. , because we are calculating the maximum norm, and |a\cdot b|<= 2|a||b|
      //recursive calculation, first calculates bound for b = u^2, and then for b*u
      boundSquare = 2. * C
          * (C * twoToS / ((s - 1.) * power(M, s - 1))
              + 0.5 * twoToS * twoToS * C / (power(2 * M + 1, s)) + twoToS * A);

      //Estimating for b*u
      RealType A1 = A, A2 = A, C1 = boundSquare, C2 = C;
      int s1 = s, s2 = s;

      if (s1 == s2) { //simpler case
        RealType twoToS = power(2., s1);
        //multiplication by extra 2. , because we are calculating the maximum norm, and |a\cdot b|<= 2|a||b|
        bound = C1
            * (C2 * twoToS / ((s1 - 1.) * power(M, s1 - 1))
                + 0.5 * twoToS * twoToS * C2 / (power(2 * M + 1, s1)))
            + C2 * twoToS * A1 + C1 * twoToS * A2;
        //tailDebug << "A1: " << A1 << "\n" << "A2: " << A2 << "\n";
      } else { //pb1.m_s != pb2.m_s
        std::cerr
            << "error in estimatePolynomialBoundForTheInfinitePart(u) function. Numbers s1 and s2 should be equal\n";
        throw std::runtime_error(
            "error in estimatePolynomialBoundForTheInfinitePart(u) function. Numbers s1 and s2 should be equal\n");
      }
      ///!add the bound for the infinite part
      bound += C1 * C2 / ((s1 - 1.) * power(M, s1 - 1))
          + C1 * C2 / ((s2 - 1.) * power(M, s2 - 1));

      return bound;
    } else {
      std::cerr
          << "EstimatePolynomialBoundForTheInfinitePart function (DPDE class) is implemented only for one dimension.\n";
      throw std::runtime_error(
          "EstimatePolynomialBoundForTheInfinitePart function (DPDE class) is implemented only for one dimension.\n");
    }
  }

  inline virtual void N(const ModesContainerType& in, ModesContainerType& out,
      int range = full) {

    if ( useFFT ) {

      CalculateNonlinearTermUsingFFT(in, out);
    } else {

      ///IMPORTANT: 11.05.2016 found bug -- size of tmcN has to be double the size of in (in order to consider all possible terms in the first convolution)
      ///for calculating the final result.

      if( !in.infiniteDimensional && tmcN.n < 2 * in.n ){
    //    std::cout << "RECALCULATING dimension of tmcN, previous dim=" << tmcN.n << "\n";
        tmcN = ModesContainerType(2 * in.n);
      }
      CalculateConvolutionDirectly(in, tmcN);
      CalculateConvolutionDirectly(in, tmcN, out, range);

      if (range != redundantRange)
        this->multiplyComponents(out);
      if (in.infiniteDimensional)
        this->multiplyComponentsFarTail(out);
      ComplexScalarType constant = Ncoeff();

      for (int j = 0; j < this->m_r; j++)
        constant *= ComplexScalarType::i();
      out *= constant; //multiplication by a coefficient in front of the convolution sum

    }
  }

  inline virtual void L(const ModesContainerType& in,
      ModesContainerType& out) const {
    IndexType i;
    const IndexRangeType& range(ir_m);

    for (i = this->firstModeIndex(range); !i.limitReached(range);
        i.inc(range)) {
      out[i] = this->lambda_k(i) * in[i];
    }
    ///here we dont care about the infinite part, it is handled in the moveParamsFunction
  }

  /**
   * Whole vector field.
   */
  inline void operator()(const ModesContainerType& in, ModesContainerType& out,
      int range = full) {
    L(in, out);
    N(in, tmc, range);

    out += tmc;

    out += yc; ///add y_c
    out += forcing; ///add forcing

    //assign to the output the calculated farTail
    out.farTail = tmc.farTail;
  }

  inline void CalculatePerturbationsUsingFFT(const PolyBdType& in,
      VectorType& out) {
    tpb.copyFinitePartFrom(in);
    tpb.cleanTail();
    tpb2.copyTailPartFrom(in);
    tpb2.cleanFinitePart();
    generalDebug2 << "tpb(finite):\n" << tpb << "\n";
    generalDebug2 << "tpb2(tail):\n" << tpb2 << "\n";
    fft3.fastTransform(tpb, tg3);
    fft3.fastTransform(tpb2, tg3_2);
    addL2Bound(in, tg3_2);
    tg3_3.multiply(tg3_2, tg3_2, tg3_2);
    generalDebug2 << "tg3_3:\n" << tg3_3 << "\n";
    tg3_4.multiply(tg3, tg3_2, tg3_2);
    generalDebug2 << "tg3_4:\n" << tg3_4 << "\n";
    tg3_5.multiply(tg3, tg3, tg3_2);
    generalDebug2 << "tg3_5:\n" << tg3_5 << "\n";
    tg3_4 *= 3.; //in general this the Newton symbol
    tg3_5 *= 3.;
    tg3_3 += tg3_4;
    tg3_3 += tg3_5;

    fft2.fastInverseTransform(tg3_3, tpb3);
    generalDebug2 << "tpb3:\n" << tpb3 << "\n";
    this->multiplyComponents(tpb3);
    ComplexScalarType constant = Ncoeff();
    for (int j = 0; j < this->m_r; j++)
      constant *= ComplexScalarType::i();
    tpb3 *= constant;
    copyFinitePart(tpb3, out);
  }

  void perturbations(const PolyBdType& in, VectorType& out) {
    if (useFFT) {
      CalculatePerturbationsUsingFFT(in, out);
    } else {
      CalculatePerturbationsUsingFFT(in, out);
    }
  }

  using BaseClass::tg1;
  using BaseClass::tg1_2;
  using BaseClass::tg2;
  using BaseClass::tg2_2;
  using BaseClass::tg3;
  using BaseClass::tg3_2;
  using BaseClass::useFFT;
  using BaseClass::tpb;
  using BaseClass::forcing;
  using BaseClass::ir_m;
  using BaseClass::ir_M;
  using BaseClass::firstModeIndex;
  using BaseClass::yc;
  using BaseClass::Ncoeff;
  using BaseClass::sdft;
  using BaseClass::rhs;
  using BaseClass::fft1;
  using BaseClass::fft2;
  using BaseClass::fft3;
  using BaseClass::m;
  using BaseClass::M;
  using BaseClass::ir_overFft3;
  using BaseClass::lambda_k;
  using BaseClass::estimatePolynomialBoundForTheInfinitePart;
  using BaseClass::addBound;
};

template<class EquationT, class FFTT, int COMPONENTS, int ORD,
    class VNormT = capd::jaco::MaximumNorm<typename FFTT::IndexType> >
class DPDE2HighDim: public EquationT {
};

/**A specialization of DPDE2 class for Burgers n dimensional, having nonlinearity (U\cdot\nabla)U,
 * where the solution U is a vector of n components U=(u_1, u_2, ..., u_n).
 *
 * Specialization for n = 2 (2D Burgers).
 */
template<class EquationT, class FFTT, int ORD, class VNormT>
class DPDE2HighDim<EquationT, FFTT, 2, ORD, VNormT> : public EquationT {
public:

  ///przebudowac pod przypadek 2D

  typedef EquationT EquationType;
  typedef FFTT FFTType;
  typedef VNormT VNormType;

  typedef typename FFTType::DFTGridType DFTGridType; //this is not taken from FFT class, as solution U has several components

  typedef typename FFTType::ModesContainerType ModesContainerType; //this is not taken from FFT class, as solution U has several components

  typedef capd::vectalg::Vector<DFTGridType, ORD> GridsContainerType;
  typedef capd::vectalg::Vector<ModesContainerType, ORD> ModesContainerContainerType;

  typedef typename FFTType::IndexType IndexType;
  typedef typename FFTType::IndexRangeType IndexRangeType;
  typedef typename FFTType::ScalarType ScalarT; ///< this is actual scalar, ScalarType typename is reserved for enclosure functions
  typedef typename EquationType::ComplexType ComplexScalarType;
  typedef typename EquationType::RealType RealType;
  typedef RealType ScalarType; ///< ScalarType has to be this, because enclosure functions are using this
  typedef capd::vectalg::Vector<RealType, 0> VectorType;
  ///przebudowac podprz 2D, koniec

  int m, M;

  int dftPts1; ///<number of discrete points for FFT, such that aliasing is omitted
  int& dftPts; ///<an alias for dftPts1
  FFTType fft1; ///<used for calculating only the finite dimensional part (takes m elements and transforms them)

  IndexRangeType ir_m, ir_finiteTail;

  DFTGridType tg1, tg1_2, sdft; ///<auxiliary variables, remember to initialize them!
  ModesContainerType tmc; ///<auxiliary variables, remember to initialize them!


  //DFTGridType r_delta_u, r_abs_mu, r_delta_v, r_abs_mv, r_mu, r_mv, s1, s2, s3, s4;

  DFTGridType empty, ///<auxiliary variable
      rhs; ///<auxiliary variable
  bool useFFT;
  ModesContainerType td, tl;
  VectorType yc, ///<this is y_c ([\delta] midpoint)
      forcing; ///<this is an optional forcing

  RealType pi;
  RealType torusScaling;

  DPDE2HighDim(int m_, int dftPts1_, RealType nu_, RealType pi_, int order) :
      /*EquationType(m_, m_, nu_), m(m_), dftPts1(dftPts1_), dftPts(dftPts1), fft1(
          m, dftPts1, pi_), ir_m(IndexType::zero()), tg1(dftPts1), tg1_2(
          dftPts1), tg1_3(dftPts1), sdft(dftPts1), tmc(m), tmc2(m), tmc3(m), empty(
          dftPts1), rhs(dftPts1), useFFT(true), td(m), tl(m), yc(
          ModesContainerType::modes2realArraySizeStatic(m_)), forcing(
          ModesContainerType::modes2realArraySizeStatic(m_)), pi(pi_), torusScaling(1.) */
      EquationType(m_, m_, nu_), m(m_), dftPts1(dftPts1_), dftPts(dftPts1), fft1(
          m, dftPts1, pi_), ir_m(IndexType::zero()), tg1(dftPts1), tg1_2(
          dftPts1), sdft(dftPts1), tmc(m), rhs(dftPts1), useFFT(true), yc(
          ModesContainerType::modes2realArraySizeStatic(m_)), forcing(
          ModesContainerType::modes2realArraySizeStatic(m_)), pi(pi_), torusScaling(1.){
    ir_m.setRange(0, strong, m, weak);

    /*if( !(typeid(capd::jaco::Index2D)==typeid(typename FFTType::IndexType) and typeid(capd::jaco::Index2D)==typeid(typename ModesContainerType::IndexType)) ){
     std::cerr << "Used classes of Incompatible indexes in DPDE2HighDim class\n";
     throw std::runtime_error("Used classes of Incompatible indexes in DPDE2HighDim class\n");
     }*/
  }

  //calculates the scalar product (U \cdot \nabla)U
  inline void CalculateNonlinearTermUsingFFT(const ModesContainerType& in,
      ModesContainerType& out) {
    if (!in.infiniteDimensional) {
      if (ScalarT::initialConditionIsRealValued())
        ScalarT::switchToComplexValued();
      fft1.fastTransform(in, tg1); //calculate dft of the projected part

      fft1.scalarProduct(tg1, tg1_2);

      //solutions << "tg1_2=\n" << tg1_2 << "\n";
      fft1.inverseExtendedTransform(tg1_2, out, 0);
      //solutions << "out=\n" << out << "\n";
      fft1.inverseExtendedTransform(tg1_2, out, 1);

    //  fft1.project(out, out);

      if (ScalarT::initialConditionIsRealValued())
        ScalarT::switchToRealValued();
    } else {
      std::cerr
          << "DPDE2HighDim.CalculateNonlinearTermUsingFFT is implemented only for finite dimension.";
      throw std::runtime_error(
          "DPDE2HighDim.CalculateNonlinearTermUsingFFT is implemented only for finite dimension.");
    }
  }

  //calculates 'reduced'  (U \cdot \nabla)U (see N1, N2 description)
  inline void CalculateNonlinearTermUsingFFT_component1(const ModesContainerType& in,
      ModesContainerType& out) {
    if (!in.infiniteDimensional) {
      if (ScalarT::initialConditionIsRealValued())
        ScalarT::switchToComplexValued();
      fft1.partialFastTransform4(in, tg1); //calculate dft of the projected part

      fft1.scalarProduct(tg1, tg1_2);

      //solutions << "tg1_2=\n" << tg1_2 << "\n";
      fft1.inverseExtendedTransform(tg1_2, out, 0);
      //solutions << "out=\n" << out << "\n";
      fft1.inverseExtendedTransform(tg1_2, out, 1);
      if (ScalarT::initialConditionIsRealValued())
        ScalarT::switchToRealValued();
    } else {
      std::cerr
          << "DPDE2HighDim.CalculateNonlinearTermUsingFFT is implemented only for finite dimension.";
      throw std::runtime_error(
          "DPDE2HighDim.CalculateNonlinearTermUsingFFT is implemented only for finite dimension.");
    }
  }

  //calculates 'reduced'  (U \cdot \nabla)U (see N1, N2 description)
  inline void CalculateNonlinearTermUsingFFT_component2(const ModesContainerType& in,
      ModesContainerType& out) {
    if (!in.infiniteDimensional) {
      if (ScalarT::initialConditionIsRealValued())
        ScalarT::switchToComplexValued();
      fft1.partialFastTransform1(in, tg1); //calculate dft of the projected part

      fft1.partialScalarProduct2(tg1, tg1_2);

      //solutions << "tg1_2=\n" << tg1_2 << "\n";
      fft1.inverseExtendedTransform(tg1_2, out, 0);
      //solutions << "out=\n" << out << "\n";
      fft1.inverseExtendedTransform(tg1_2, out, 1);
      if (ScalarT::initialConditionIsRealValued())
        ScalarT::switchToRealValued();
    } else {
      std::cerr
          << "DPDE2HighDim.CalculateNonlinearTermUsingFFT is implemented only for finite dimension.";
      throw std::runtime_error(
          "DPDE2HighDim.CalculateNonlinearTermUsingFFT is implemented only for finite dimension.");
    }
  }


  inline void rightHandSide(int i, const GridsContainerType& grids,
      const ModesContainerContainerType& modes, ModesContainerType& rhsSeries,
      DFTGridType& rhsFunctionSpace, bool calculateRhsFunctionSpace = true) {
    int j;
    //if the solution is real valued and even/odd then there are a lot of zeros in the jets, and the additions/multiplication by
    //zeros should be avoided
    //if(ScalarT::initialConditionIsRealValued()) ///optimizing thing
    //  ScalarT::switchToRealValuedL2(); //if solution is real valued then imaginary part of jets is zero
    for (j = 0; j <= i; ++j) {
      ///remark: there is already multiplication by k inside this procedure
      ///< calculates \sum_{k_1\in a Projection}{(k, a_k)\cdot a_{k-k_1}}
      fft1.scalarProduct(grids[j], grids[i - j], sdft);
      if (j == 0)
        rhs = sdft; //this is needed in order to s have the same DPDEContainer
      else
        rhs += sdft;
    }
    fft1.inverseExtendedTransform(rhs, rhsSeries, 0);
    fft1.inverseExtendedTransform(rhs, rhsSeries, 1);

    generalDebug2 << "rhsSeries:\n" << rhsSeries << "\n";
    if (ScalarT::initialConditionIsRealValued())
      ScalarT::switchToRealValued();
    IndexType index;
    for (index = this->firstModeIndex(ir_m); !index.limitReached(ir_m);
        index.inc(ir_m)) {
      //here we have minus - because rhsSeries = Pu\cdot\nabla u
      rhsSeries.set(index, rhsSeries[index] + this->lambda_k(index) * modes[i][index]);
    }
    if (i == 0) { //vector field is calculated, and thus y_c and forcing have to be added
      rhsSeries += yc;
      rhsSeries += forcing;
    }
    rhsSeries *= RealType(1) / RealType(i + 1);
    //we switch back to complex valued, because the FFT is complex
    ScalarT::switchToComplexValued();
    if (calculateRhsFunctionSpace)
      fft1.fastTransform(rhsSeries, rhsFunctionSpace);
  }

  inline void rightHandSide(int i, const ModesContainerContainerType& modes,
      ModesContainerType& rhsSeries) {
    std::cerr << "DPDE2HighDim.rightHandSide is not implemented.";
    throw std::runtime_error("DPDE2HighDim.rightHandSide is not implemented.");
  }

  /**Nonlinear part of the vector field.
   */

  inline virtual void N(const ModesContainerType& in, ModesContainerType& out,
      int range = full) {
    if (useFFT) {
      CalculateNonlinearTermUsingFFT(in, out);
    } else {
      std::cerr << "DPDE2HighDim.N is implemented only for FFT.";
      throw std::runtime_error("DPDE2HighDim.N is implemented only for FFT.");
    }
  }

  /** Returns 'reduced' Nonlinearity vectors, first is only with x-derivative,
   *  second is only with y-derivative
   */
  inline virtual void N1(const ModesContainerType& in, ModesContainerType& out,
        int range = full) {
      if (useFFT) {
        CalculateNonlinearTermUsingFFT_component1(in, out);
      } else {
        std::cerr << "DPDE2HighDim.N is implemented only for FFT.";
        throw std::runtime_error("DPDE2HighDim.N is implemented only for FFT.");
      }
    }
  inline virtual void N2(const ModesContainerType& in, ModesContainerType& out,
        int range = full) {
      if (useFFT) {
        CalculateNonlinearTermUsingFFT_component2(in, out);
      } else {
        std::cerr << "DPDE2HighDim.N is implemented only for FFT.";
        throw std::runtime_error("DPDE2HighDim.N is implemented only for FFT.");
      }
    }


  inline virtual void L(const ModesContainerType& in,
      ModesContainerType& out) const {
    IndexType i;
    const IndexRangeType& range(ir_m);

    for (i = this->firstModeIndex(range); !i.limitReached(range);
        i.inc(range)) {
      ///TEMP: torus scaling
      out[i] = this->torusScaling * this->lambda_k(i) * in[i];
    }
    ///here we dont care about the infinite part, it is handled in the moveParamsFunction
  }

  inline void operator()(const ModesContainerType& in, ModesContainerType& out,
      int range = full) {
    L(in, out);
    N(in, tmc, range);

    out += tmc;
    out += yc; ///add y_c

    out += forcing; ///add forcing

  }

  inline void calculateGrids(int i, const ModesContainerContainerType& modes,
      GridsContainerType& grids) {
    std::cerr
        << "DPDE2HighDim.calculateGrids is  not implemented for DPDE2HighDim.";
    throw std::runtime_error(
        "DPDE2HighDim.calculateGrids is not implemented for DPDE2HighDim.");
  }

};

/**
 * This will be class representing 2D NS equations - stream function formulation
 */
template<class EquationT, class FFTT, int ORD>
class NSStream: public DPDE2<EquationT, FFTT, ORD> {

};

}
}

#endif /* EQUATIONS_H_ */
