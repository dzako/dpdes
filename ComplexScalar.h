/*
 * ComplexInterval.h
 *
 *  Created on: Aug 8, 2011
 *      Author: cyranka
 */

#ifndef COMPLEXINTERVAL_H_
#define COMPLEXINTERVAL_H_

#include "DPDEContainer.h"

namespace capd{
namespace jaco{

///Box representation of complex intervals
//TODO: dokładniejszy opis!!!
template<typename ScalarT>
class ComplexScalar{
public:
  typedef ScalarT ScalarType;

  ScalarType re, im;

  inline ComplexScalar() : re(0), im(0){}
  inline ComplexScalar(const ScalarType& r) : re(r), im(0){}
//  ComplexScalar(double d) : re(d), im(0){}
  inline ComplexScalar(const ScalarType& r, const ScalarType& i) : re(r), im(i){}
  inline ComplexScalar(const ComplexScalar& ci) : re(ci.re), im(ci.im){}
  ///Dummy function
  static void switchToLocalOptimization(){}
  ///Dummy function
  static void switchToGlobalOptimization(){}
  ///Dummy function
  static void switchToRealValued(){}
  ///Dummy function
  static void switchToComplexValued(){}
  ///Dummy function
  static void switchToRealValuedL2(){}
  ///Dummy function
  static bool initialConditionIsRealValued(){return 0;}
  ///Dummy function
  static void setContainer(capd::jaco::DPDEContainer& c){}
  
  static const ComplexScalar i(){
    return ComplexScalar(ScalarType(0), ScalarType(1));
  }

  inline bool isImUnit() const{
    if(re == 0 && im != 0)
      return true;
    return false;
  }

  inline bool isFullComplex() const{
    if(re != 0 && im != 0)
      return true;
    return false;
  }

  inline ScalarType squareNorm() const{
    return power(re, 2)+power(im, 2);
  }

  inline ComplexScalar mid() const{
    ComplexScalar r;
    r.re = re.mid();
    r.im = im.mid();
    return r;
  }
  
  inline void split(ComplexScalar& c, ComplexScalar& r) const{
    (*this).re.split(c.re, r.re);
    (*this).im.split(c.im, r.im);
  }

  inline ComplexScalar abs() const{
    ComplexScalar r;
    r.re = capd::abs(re);
    r.im = capd::abs(im);
    return r;
  }

  inline ComplexScalar supremum() const{
    ComplexScalar r;
    r.re = rightBound(re);
    r.im = rightBound(im);
    return r;
  }

  inline ComplexScalar abs_supremum() const{
    ComplexScalar r;
    r.re = rightBound(capd::abs(re));
    r.im = rightBound(capd::abs(im));
    return r;
  }

  inline ComplexScalar inverse() const{
    ComplexScalar inv;
    inv.re = re / squareNorm();
    inv.im = -im / squareNorm();
    return inv;
  }

  inline ComplexScalar& conjugate(){
    im = -im;
    return *this;
  }

  inline void setImaginaryPartToZero(){
    im = 0;
  }

  inline void setRealPartToZero(){
    re = 0;
  }

  inline void projectOntoImaginarySpace(){
    re = 0;
  }

  inline void projectOntoRealSpace(){
    im = 0;
  }

  inline ComplexScalar& operator=(const ComplexScalar& ci){ re=ci.re; im=ci.im; return (*this); }

  inline ComplexScalar& operator=(const ScalarType& i){ re = i; im = 0; return (*this); }

  inline ComplexScalar& operator+=(const ComplexScalar& ci){ re+=ci.re; im+=ci.im; if(__COUNT_OPERATIONS__){CCadditionsSum++;} return (*this); }

  inline ComplexScalar& operator-=(const ComplexScalar& ci){ re-=ci.re; im-=ci.im; if(__COUNT_OPERATIONS__){CCadditionsSum++;} return (*this); }

  inline ComplexScalar& operator*=(const ComplexScalar& c){
    ScalarType tRe = re,
               tIm = im;
    re = tRe * c.re - tIm * c.im;
    im = tRe * c.im + c.re * tIm;
    if(__COUNT_OPERATIONS__){CCmultiplicationsSum++;}
    
    return (*this);
  }

  inline ComplexScalar& operator*=(const ScalarType& s){
    re *= s;
    im *= s;    
    if(__COUNT_OPERATIONS__){CRmultiplicationsSum++;}
    return *this;
  }

  inline bool operator==(const ScalarType& s) const{
    if(re == s && im == 0)
      return true;
    return false;
  }

//  inline ComplexScalar& operator=(const double& d){ re=d; im=0; return (*this); }

  inline void setConjugate(const ComplexScalar& c){
    re = c.re;
    im = -c.im;
  }

  ///returns a bound for the complex number norm
  inline ScalarType norm() const{
    return sqrt( ScalarType(2.) ) * this->normMax();
  }

  inline ScalarType normMax() const{
    ScalarType max = rightBound(capd::abs(re)), t;
    if((t = rightBound(capd::abs(im))) > max) max = t;
    return max;
  }

  /**A wrapper function needed by the integrators. A dummy function, is doing nothing.
   */
  inline void setVariationalPartToId(int position){
  }

  /**A wrapper function needed by the integrators. A dummy function, is doing nothing.
   */
  inline void setVariationalPartToIdConjugate(int position){
  }

  /**A wrapper function needed by the integrators. A dummy function, returns 0.
   */
  inline int variationalPart(int position) const{
    return 0;
  }

  /**A wrapper function needed by the integrators.
    */
  inline const ComplexScalar& value() const{
    return *this;
  }

  ComplexScalar& value(){
    return *this;
  }


  /**A wrapper function needed by the integrators.
      */
  inline const ComplexScalar& secondFreeCoeff() const{
    return *this;
  }
  /**A wrapper function needed by the integrators.
      */
//  inline const ComplexScalar& thirdFreeCoeff() const{
//    return *this;
//  }

  inline void setFreeCoeff(const ComplexScalar& val){
    *this = val;
  }

  /**A wrapper function needed by the integrators.
      */
  inline void setSecondFreeCoeff(const ComplexScalar& val){
  }
  /**A wrapper function needed by the integrators.
      */
//  inline void setThirdFreeCoeff(const ComplexScalar& val){
//  }

//  inline Interval diam() const{
//    return ( diam(re) > diam(im) ? diam(re) : diam(im) );
//  }

  friend std::ostream& operator<<(std::ostream& out, const ComplexScalar& c){
    //out << "(" << fadbad::val(c.re.val()) << "," << fadbad::val(c.im.val()) << ")"; //TODO: for debug purpose
    out << "(" << c.re << "," << c.im << ") diam (" << diam(c.re) << ", " << diam(c.im) << ")";
    return out;
  }

};

template<class ScalarT>
inline ComplexScalar<ScalarT> conjugate(const ComplexScalar<ScalarT>& cs){
  ComplexScalar<ScalarT> r(cs);
  r.im = -r.im;
  return r;
}

template<typename ScalarT>
inline const ComplexScalar<ScalarT> operator*(const ComplexScalar<ScalarT>& ci1, const ComplexScalar<ScalarT>& ci2){
  ComplexScalar<ScalarT> r;
  r.re = ci1.re*ci2.re - ci1.im*ci2.im;
  r.im = ci1.re*ci2.im + ci1.im*ci2.re;
  if(__COUNT_OPERATIONS__){CCmultiplicationsSum++;}
  return r;
}

template<typename ScalarT>
inline const ComplexScalar<ScalarT> operator/(const ComplexScalar<ScalarT>& ci1, const ComplexScalar<ScalarT>& ci2){
  ComplexScalar<ScalarT> r, t;
  t = ci2.inverse();
  r.re = ci1.re*t.re-ci1.im*t.im;
  r.im = ci1.re*t.im+ci1.im*t.re;
  if(__COUNT_OPERATIONS__){CCmultiplicationsSum++;}
  return r;
}

template<typename ScalarT>
inline const ComplexScalar<ScalarT> operator/(const ComplexScalar<ScalarT>& ci1, const ScalarT& ci2){
  ComplexScalar<ScalarT> r = ci1;
  r.re /= ci2;
  r.im /= ci2;
  if(__COUNT_OPERATIONS__){CRmultiplicationsSum++;}
  return r;
}

template<typename ScalarT>
inline const ComplexScalar<ScalarT> operator+(const ComplexScalar<ScalarT>& ci1, const ComplexScalar<ScalarT>& ci2){
  ComplexScalar<ScalarT> r;
  r.re=ci1.re+ci2.re;
  r.im=ci1.im+ci2.im;
  if(__COUNT_OPERATIONS__){CCadditionsSum++;}
  return r;
}

template<typename ScalarT>
inline const ComplexScalar<ScalarT> operator-(const ComplexScalar<ScalarT>& ci1, const ComplexScalar<ScalarT>& ci2){
  ComplexScalar<ScalarT> r;
  r.re=ci1.re-ci2.re;
  r.im=ci1.im-ci2.im;
  if(__COUNT_OPERATIONS__){CCadditionsSum++;}
  return r;
}


template<typename ScalarT>
inline const ComplexScalar<ScalarT> operator*(const ScalarT& d, const ComplexScalar<ScalarT>& ci){
  ComplexScalar<ScalarT> r;
  r.re=d*ci.re;
  r.im=d*ci.im;
  if(__COUNT_OPERATIONS__){CRmultiplicationsSum++;}
  return r;
}

//updated 08.02.14 now the program compiles without problems for both of the non-rigorous and rigorous version
template<typename ScalarT>
inline const ComplexScalar<ScalarT> operator*(typename ScalarT::BoundType d, const ComplexScalar<ScalarT>& ci){
  ComplexScalar<ScalarT> r;
  r.re = ScalarT(d) * ci.re;
  r.im = ScalarT(d) * ci.im;
  if(__COUNT_OPERATIONS__){CRmultiplicationsSum++;}
  return r;
}

}
}

#endif /* COMPLEXINTERVAL_H_ */
