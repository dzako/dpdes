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
 * Pair.h
 *
 *  Created on: Sep 26, 2011
 *      Author: cyranka
 */

#ifndef PAIR_H_
#define PAIR_H_

#include "capd/vectalg/Matrix.hpp"
#include "DPDEContainer.h"

namespace capd{
namespace jaco{

/* Class representing partial derivatives matrix (derivative of a_k with respect to a_i)
 *  $$
 *  \left[
 *    \begin{array}{cc}
 *      \partial Re{a_k} / \partial Re{a_i} & \partial Re{a_k} / \partial Im{a_i} //
 *      \partial Im{a_k} / \partial Re{a_i} & \partial Im{a_k} / \partial Im{a_i}
 *    \end{array}
 *  \right]
 *  $$
 */
template<typename ScalarT>
class ComplexDerivativePair{
public:
  typedef ScalarT ScalarType;
  typedef typename ScalarType::ScalarType InternalType;
  typedef typename capd::vectalg::Matrix<InternalType, 2, 2> MatrixType;

  ScalarType partialWRespToRe; ///<a complex number - derivative of a complex with respect to a real part $(\partial Re{a_k} / \partial Re{a_i}, \partial Im{a_k} / \partial Re{a_i})$
  ScalarType partialWRespToIm; ///<a complex number - derivative of a complex with respect to an imaginary part $(\partial Re{a_k} / \partial Im{a_i}, \partial Im{a_k} / \partial Im{a_i})$

  inline ComplexDerivativePair() : partialWRespToRe(0), partialWRespToIm(0){
  }

  inline ComplexDerivativePair(double d) : partialWRespToRe(d), partialWRespToIm(d){}

  inline ComplexDerivativePair(const ScalarType& first, const ScalarType& second) : partialWRespToRe(first), partialWRespToIm(second){
  }

  inline void setToMatrix(const MatrixType& m){
    partialWRespToRe.re = m[0][0];
    partialWRespToRe.im = m[0][1];
    partialWRespToIm.re = m[1][0];
    partialWRespToIm.im = m[1][1];
  }

  inline ComplexDerivativePair& operator=(const MatrixType& m){
    setToMatrix(m);
    return *this;
  }

  inline ComplexDerivativePair& operator=(const ComplexDerivativePair& cdp){
      partialWRespToRe = cdp.partialWRespToRe;
      partialWRespToIm = cdp.partialWRespToIm;

    return *this;
  }

  inline const ComplexDerivativePair& operator+=(const ComplexDerivativePair& cdp){
      partialWRespToRe += cdp.partialWRespToRe;
      partialWRespToIm += cdp.partialWRespToIm;
    return *this;
  }

  inline const ComplexDerivativePair& operator-=(const ComplexDerivativePair& cdp){
      partialWRespToRe -= cdp.partialWRespToRe;
      partialWRespToIm -= cdp.partialWRespToIm;
    return *this;
  }

  inline const ComplexDerivativePair& operator*=(const ScalarType& s){
      partialWRespToRe *= s;
      partialWRespToIm *= s;
    return *this;
  }

  /* Returns partial derivative matrix, it is important to perform assingments in the order. The correct square matrix is in this way,
   * see the class description.
   *
   */
  inline operator MatrixType() const{
    MatrixType r;
    r[0][0] = partialWRespToRe.re;
    r[0][1] = partialWRespToIm.re;
    r[1][0] = partialWRespToRe.im;
    r[1][1] = partialWRespToIm.im;
    return r;
  }

   inline void setImaginaryPartToZero(){
     partialWRespToRe.setImaginaryPartToZero();
     partialWRespToIm.setImaginaryPartToZero();
   }

   inline void setRealPartToZero(){
      partialWRespToRe.setRealPartToZero();
      partialWRespToIm.setRealPartToZero();
    }

   inline void projectOntoImaginarySpace(){
     partialWRespToRe.setRealPartToZero();
     partialWRespToIm.setRealPartToZero();
   }

   inline void projectOntoRealSpace(){
      partialWRespToIm.setImaginaryPartToZero();
      partialWRespToRe.setImaginaryPartToZero();
    }

   inline void setToId(const MatrixType& m = MatrixType::Identity(2)){
     setToMatrix(m);
   }

   inline void setToIdConjugate(const MatrixType& m = MatrixType::Identity(2)){
     setToMatrix(m);
     (*this).conjugate();
   }

   inline const ComplexDerivativePair& conjugate(){
     partialWRespToIm.conjugate();
     partialWRespToRe.conjugate();
     return *this;
   }

  friend std::ostream& operator<<(std::ostream& out, const ComplexDerivativePair& d){
    out << d.partialWRespToRe << "\n" << d.partialWRespToIm << "\n";
    return out;
  }

};

template<typename ScalarT>
inline ComplexDerivativePair<ScalarT> conjugate(const ComplexDerivativePair<ScalarT>& p){
  ComplexDerivativePair<ScalarT> r(p);
    r.partialWRespToRe.conjugate();
    r.partialWRespToIm.conjugate();
  return r;
}

template<typename ScalarT>
inline const ComplexDerivativePair<ScalarT> operator+(const ComplexDerivativePair<ScalarT>& p1, const ComplexDerivativePair<ScalarT>& p2){
  ComplexDerivativePair<ScalarT> r;
    r.partialWRespToRe = p1.partialWRespToRe + p2.partialWRespToRe;
    r.partialWRespToIm = p1.partialWRespToIm + p2.partialWRespToIm;
  return r;
}

template<typename ScalarT>
inline const ComplexDerivativePair<ScalarT> operator-(const ComplexDerivativePair<ScalarT>& p1, const ComplexDerivativePair<ScalarT>& p2){
  ComplexDerivativePair<ScalarT> r;
    r.partialWRespToRe = p1.partialWRespToRe - p2.partialWRespToRe;
    r.partialWRespToIm = p1.partialWRespToIm - p2.partialWRespToIm;
  return r;
}

template<typename ScalarT>
inline const ComplexDerivativePair<ScalarT> operator*(const ScalarT& s, const ComplexDerivativePair<ScalarT>& p){
  ComplexDerivativePair<ScalarT> r;
    r.partialWRespToRe = s * p.partialWRespToRe;
    r.partialWRespToIm = s * p.partialWRespToIm;
  return r;
}

//template<typename ScalarT>
//inline const ComplexDerivativePair<ScalarT> operator*(const typename ComplexDerivativePair<ScalarT>::InternalType& s, const ComplexDerivativePair<ScalarT>& p){
//  ComplexDerivativePair<ScalarT> r;
//  r.partialWRespToRe = s * p.partialWRespToRe;
//  r.partialWRespToIm = s * p.partialWRespToIm;
//  return r;
//}

template<typename ScalarT>
inline ComplexDerivativePair<ScalarT> operator/(const ComplexDerivativePair<ScalarT>& p, const ScalarT& s){
  ComplexDerivativePair<ScalarT> r;
    r.partialWRespToRe = p.partialWRespToRe / s;
    r.partialWRespToIm = p.partialWRespToIm / s;
  return r;
}

template<typename ScalarT>
inline const ComplexDerivativePair<ScalarT> operator*(double d, const ComplexDerivativePair<ScalarT>& p){
  ComplexDerivativePair<ScalarT> r;
    r.partialWRespToRe = d * p.partialWRespToRe;
    r.partialWRespToIm = d * p.partialWRespToIm;
  return r;
}

}
}

#endif /* PAIR_H_ */
