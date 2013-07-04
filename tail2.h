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
 * tail2.h
 *
 *  Created on: Sep 9, 2011
 *      Author: Jacek Cyranka
 */

#ifndef _CAPD_JACO_TAIL2_H
#define _CAPD_JACO_TAIL2_H

#include "ComplexScalar.h"
#include "dissipative_enclosure.h"

namespace capd {
namespace jaco {

/////All tails are dedicated to store complex valued tails and that are considered for symmetric Galerkin
/////projection.
/////We use operator()(int index, bool re) when we want to obtain value of i-th mode (index), real or imaginary part boolean re,
/////We use operator[](int index) when explicit location in a table is provided by index.

///infinite part of a tail \prod_{|k|>M}
template <typename ComplexScalarT, typename DimensionT>
class FarTail2 : public DimensionT {
public:
  typedef typename ComplexScalarT::ScalarType ScalarType;
  typedef typename DimensionT::IndexType IndexType;
  typedef ComplexScalarT ComplexScalarType;

  ScalarType m_c;
  int m_s;

  FarTail2(const FarTail2& ft2) :
    DimensionT(ft2.m, ft2.M) {
    this->m_c = ft2.m_c;
    this->m_s = ft2.m_s;
  }

  FarTail2(const ScalarType& c, const ScalarType& s, int m, int M) :
    DimensionT(m, M) {
    this->m_c = c;
    this->m_s = s;
  }

  FarTail2(int m, int M) :
    DimensionT(m, M) {
  }

  FarTail2() {
  }

  inline ScalarType value(int i) const {
    double j = (i > 0 ? i : -i);
    return (m_c / power(ScalarType(j), m_s)) * ScalarType(-1, 1);
  }

  ///Returns a value that is at index i in the internal data storage, used for iterating.
  inline ScalarType operator[](const int& i) const {
    return value(this->array2modeTail(i));
  }

  ///if re==true returns Re{a_i},
  ///else returns Im{a_i}.
  inline ScalarType operator()(const int i, const bool re) const {
    return value(i);
  }

//  inline ComplexScalarType operator()(const IndexType& index) const{
//    return (m_c / power(ScalarType(index.squareNorm()), m_s/2.)) * ScalarType(-1, 1);
//  }

  inline ComplexScalarType operator[](const IndexType& index) const{
    ScalarType s;
    if(index.d() == 1)
      s = (m_c / power(ScalarType((index[0] > 0 ? index[0] : -index[0])), m_s)) * ScalarType(-1, 1);
    else
      s = (m_c / power(ScalarType(index.squareEuclNorm()), m_s/2.)) * ScalarType(-1, 1);
    return farTailMode(s);
  }

  inline void setC(const ScalarType& C) {
    m_c = C;
  }

  inline void setS(int s) {
    this->m_s = s;
  }

  inline ScalarType getC() const {
    return m_c;
  }

  inline int getS() const {
    return m_s;
  }

  inline int getM() const {
    return this->M;
  }

  inline void setM(int M) {
    this->M = M;
  }

  inline bool inFarTail(int index) const;

  inline bool inFarTail(const IndexType& index) const;

  ///Checks if the far tail is subset of a second far tail.
  inline bool subset(const FarTail2& ft2) const;

  ///Returns \prod_{|k|>M}{0}.
  inline FarTail2 mid() const {
    return FarTail2(0., this->m_s, this->m, this->M);
  }

  ///Calculates sum of all elements absolute supremum in the tail.
  inline ScalarType sum() const {
    return sum(this->M + 1);
  }

  ///Calculates sum of elements absolute supremum in the tail, starting at an index i.
  inline ScalarType sum(int i) const;

  ///Euclidean norm of the far tail.
  inline ScalarType euclNorm() const {
    return this->getC() / ((2 * this->getS() - 1) * power(ScalarType(this->M + 1), 2 * this->getS() - 1));
  }

  ///Takes intersection of the tail with another tail ft2, assumes homogeneity of the powers s. We do not allow to store a tail
  ///with varying exponents.
  inline void intersect(const FarTail2& ft2);

  friend void operator+=(FarTail2& ft1, const FarTail2& ft2){
    if(ft1.m_s != ft2.m_s){
      std::cerr << "operator+= of the FarTail class works only with tails of equal exponents.\n";
      throw std::runtime_error("operator+= of the FarTail class works only with tails of equal exponents.\n");
    }
    if(ft1.m_c < ft2.m_c)
      ft1.m_c = ft2.m_c;
  }

  friend std::ostream& operator<<(std::ostream& out, const FarTail2& t) {
    out << "far tail (k>" << t.M << "): \n|a_k| <= " << t.getC() << " / |k|^" << t.getS() << ", a_" << t.M+1 << " =" << t[IndexType(t.M+1)] << "\n";
    return out;
  }
  
  using DimensionT::farTailMode;

};


template <typename ScalarT, typename DimensionT>
inline bool FarTail2<ScalarT, DimensionT>::inFarTail(int index) const {
  if(index > this->M || index < -this->M)
    return true;
  return false;
}

template <typename ScalarT, typename DimensionT>
inline bool FarTail2<ScalarT, DimensionT>::inFarTail(const IndexType& index) const {
  if(index.squareEuclNorm() > this->M*this->M)
    return true;
  return false;
}

template <typename ScalarT, typename DimensionT>
inline bool FarTail2<ScalarT, DimensionT>::subset(const FarTail2& ft2) const {
  if(ft2.getS() > this->getS())
    return false;
  int maxM = (this->M > ft2.getM() ? this->M : ft2.getM());
  if(this->operator()(maxM + 1, 1).subset(ft2(maxM + 1, 1)))
    return true;
  return false;
}

template <typename ScalarT, typename DimensionT>
inline typename FarTail2<ScalarT, DimensionT>::ScalarType FarTail2<ScalarT, DimensionT>::sum(int i) const {
  if(!inFarTail(i)) {
    std::cerr << "Requested to sum modes from the far tail, starting at index i=" << i << " that is out of the far tail.\n";
    throw std::runtime_error("Requested to sum modes from the far tail, starting at index that is out of the far tail.\n");
  }
  return m_c / (POW(ScalarType(i > 0 ? i : -i), m_s - 1) * (m_s - 1));
}

template <typename ScalarT, typename DimensionT>
inline void FarTail2<ScalarT, DimensionT>::intersect(const FarTail2<ScalarT, DimensionT>& ft2){
  if(m_s!=ft2.m_s){
    std::cerr<< "We do not allow to intersect two tails with different exponents.\n";
    throw std::runtime_error("We do not allow to intersect two tails with different exponents.\n");
  }//otherwise we intersect the far tails
  m_c=min(m_c, ft2.m_c);
}


}
}
#endif

