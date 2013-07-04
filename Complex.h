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
 * Complex.h
 *
 *  Created on: Dec 7, 2011
 *      Author: cyranka
 */

#ifndef COMPLEX_H_
#define COMPLEX_H_

#include "Index.h"
#include "ComplexScalar.h"
#include "norms.h"

namespace capd{
namespace jaco{

/**For Complex vectors in which real and imaginary part are kept separately.
 */

///this class is base for all classes that operate on the Fourier modes of one dimensional
///complex functions with periodic bd. conditions.
///Modes (Galerkin projection, near tail) are always stored in an array, and this class
///provides methods to retrieve from this array indicated modes and vice versa.

/**TODO: this class is important and not tested
 *
 */
template<class ScalarT, class IndexT, class NormT,
  class IndexRangeT = capd::jaco::FourierConvolutionIndexRange<IndexT, NormT>,
  class VectorT = capd::vectalg::Vector<ScalarT, _D> >
class Complex
{
public:
  typedef typename capd::jaco::Complex<ScalarT, IndexT, NormT, IndexRangeT, VectorT> Class;
  typedef VectorT VectorType;
  typedef ScalarT ScalarType;
  typedef IndexT IndexType;
  typedef IndexRangeT IndexRangeType;
  typedef NormT NormType;

  int m, M;

  Complex(int m_, int M_) : m(m_), M(M_){}

  ///constructor that uses default dimensions
  Complex() : m(_m), M(_M){}

  inline int mode2array(const IndexType& k, bool re) const{
    IndexType i(k);
    i[t.d()-1] += m;
    return i.mode2array(2*m, re);
  }

  inline int mode2arrayTail(const IndexType& k, bool re) const
  {
    ///TODO: this has to be thinked through
//    IndexType t(k);
//    t[t.d()-1] -= m;
//    if(k.upperHalfspace)
//      return mode2array(t, re);
//    else{
//      t[t.d()-1] += m;
//      return mode2array(t, re);
//    }

  }

  inline IndexType array2modeIndex(int i) const{
    return IndexType::array2modeIndex(m, i);
  }

  inline const IndexType firstModeIndex() const{
    return array2modeIndex(1);
  }

  ///Returns Index of the first mode which is stored and is within IndexRange ir
  inline const IndexType firstModeIndex(const IndexRangeType& ir, int l = 0) const{
    return firstWithinRange(ir);
  }

  ///Returns index of a first mode within range ir
  inline const IndexType firstWithinRange(const IndexRangeType& ir) const{
    IndexType idx;
    int i;
    for(i=0; i<idx.d(); i++)
      idx[i]=ir.returnIndex();
    while(!ir.withinRange(idx)){
      idx.inc(ir);
    }
    return idx;
  }

  ///TODO: finish, gives the smallest element on first and second position for this condition. For example for 2D indices
  ///the used ordering is (-m, 0), ... , (m, 0), (-m, 1), ... . Therefore if the second component is increased by 1, the first
  ///is reseted to -m.

  ///TODO: should be more sharp
  inline const int modes2arraySize(int m) const{
    IndexType idx;
    int i;
    for(i = 0; i<idx.d(); i++)
      idx[i] = m;
    idx.l = idx.d()-1;
    return idx.mode2array(m, 0)+1;
  }

  ///TODO: should be more sharp
  inline const int modes2arraySize() const{
    return modes2arraySize(m);
  }

  ///Modes with the most important index > 0 are stored . e.g. for 2D index k=(k_0, k_1) those with k_1 > 0 are stored in an array.
  ///The rest (a_{-k}) is obtained by conjugating a_k.
  template <typename AVector>
  inline typename capd::jaco::ComplexScalar<typename AVector::ScalarType> mode(const IndexType& k, const AVector& vec) const{
    typedef capd::jaco::ComplexScalar<typename AVector::ScalarType> ComplexScalar;
    if(!k.integer()){
      std::cerr << "k="<<k<<"\nOperation of obtaining a mode with index k with non integer k is not defined.\n";
      throw std::runtime_error("Operation of obtaining a mode with index k with non integer k is not defined.\n");
    }
    if(k.isZero()){
      return ComplexScalar();
    }else{
      if(k.upperHalfspace()){
        return ComplexScalar(vec[this->mode2array(k, 1)], vec[this->mode2array(k, 0)]);
      }
      return ComplexScalar(vec[this->mode2array(-k, 1)], -vec[this->mode2array(-k, 0)]);
    }
  }

  template <typename AVector, typename ComplexScalar>
  inline ComplexScalar mode(const IndexType& k, const AVector& vec) const{
    if(!k.integer()){
      std::cerr << "k="<<k<<"\nOperation of obtaining a mode with index k with non integer k is not defined.\n";
      throw std::runtime_error("Operation of obtaining a mode with index k with non integer k is not defined.\n");
    }
    if(k.isZero()){
      return ComplexScalar();
    }else{
      if(k.upperHalfspace()){
        return ComplexScalar(vec[this->mode2array(k, 1)], vec[this->mode2array(k, 0)]);
      }
      return ComplexScalar(vec[this->mode2array(-k, 1)], -vec[this->mode2array(-k, 0)]);
    }
  }

  ///Stores calculated a_k's mode value in internal representation
  template<typename AVector>
  inline void setMode(const IndexType& k, const typename capd::jaco::ComplexScalar<typename AVector::ScalarType>& mode, AVector& vec) const{
    vec[Class::mode2array(k, 1)] = mode.re;
    vec[Class::mode2array(k, 0)] = mode.im;
  }

  ///returns value of k-th mode in the near tail,
  ///if re==true return Re{a_k},
  ///else return Im{a_k}.
  ///Value of a_k for k=m+1,\dots,M is stored explicitly,
  ///whereas value of -a_k is obtained by taking conjugate of a_k.
  template<typename AVector>
  inline typename capd::jaco::ComplexScalar<typename AVector::ScalarType> closeTailMode( const IndexType& k, const AVector& vec) const
  {
    typedef capd::jaco::ComplexScalar<typename AVector::ScalarType> ComplexScalar;
    if(!k.integer()){
      std::cerr << "Operation of obtaining a mode with index k with non integer k is not defined.\n";
      throw std::runtime_error("Operation of obtaining a mode with index k with non integer k is not defined.\n");
    }
    if(k.storedExplicitly()){
      return ComplexScalar(vec[this->mode2arrayTail(k, 1)], vec[this->mode2arrayTail(k, 0)]);
    }
    return ComplexScalar(vec[this->mode2arrayTail(-k, 1)], -vec[this->mode2arrayTail(-k, 0)]);
  }


  inline bool inTail(const IndexType& index) const{
    if(index.squareNorm() > this->m*this->m)
      return true;
    return false;
  }

  template <typename VectorType, typename TailType>
  inline typename capd::jaco::ComplexScalar<typename VectorType::ScalarType> mode(const IndexType& index, const VectorType& u, const TailType& tail) const{
    typename capd::jaco::ComplexScalar<typename VectorType::ScalarType> r;
    if(inTail(index)){
      r.re = tail(index).re;
      r.im = tail(index).im;
      return r;
    }else
      return mode(index, u);
  }


  inline void printModesIndex(const VectorType& v, capd::auxil::OutputStream& out) const
  {
    IndexRangeType ir(IndexType::zero());
    ir.setRange(0, strong, this->m, weak);
    IndexType index;
    for(index = firstModeIndex(ir); !index.limitReached(ir); index.inc(ir)){
      out << "Re(a_"<<index<<"): ["<<v[mode2array(index, 1)].leftBound()<<","<<v[mode2array(index, 1)].rightBound()<<"] "<<
          " Im(a_"<<index<<"): ["<<v[mode2array(index, 0)].leftBound()<<","<<v[mode2array(index, 0)].rightBound()<<"]\n";
    }
  }
};

}
}
#endif /* COMPLEX_H_ */
