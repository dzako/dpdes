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
 * Odd.h
 *
 *  Created on: Sep 14, 2011
 *      Author: cyranka
 */

#ifndef ODD_H_
#define ODD_H_

#include "Index.h"
#include "ComplexScalar.h"
#include "Real.h"

namespace capd{
namespace jaco{

///Real, odd functions subspace indexed by IndexT.
///this class is base for all classes that operate on the Fourier modes of
///real functions with periodic bd. conditions (invariant subspace a_k=\overline{a_{-k}})
///satisfying oddness condition f(-x)=-f(x).
///Modes (Galerkin projection, near tail) are always stored in an array, and this class
///provides methods to retrieve from this array indicated modes and vice versa.
template<class ScalarT, class IndexT, class NormT,  
  class MatrixT = capd::vectalg::Matrix<ScalarT, 0, 0>, 
  class IndexRangeT = capd::jaco::FourierConvolutionIndexRange<IndexT, NormT> >
class Odd : public capd::jaco::Real<ScalarT, IndexT, NormT, MatrixT, IndexRangeT>
{
public:
  typedef Odd<ScalarT, IndexT, NormT, MatrixT, IndexRangeT> Class;
  typedef capd::jaco::Real<ScalarT, IndexT, NormT, MatrixT, IndexRangeT> BaseClass;
  typedef typename BaseClass::MatrixType MatrixType;
  typedef typename BaseClass::VectorType VectorType;
  typedef typename BaseClass::ScalarType ScalarType;
  typedef typename BaseClass::IndexType IndexType;
  typedef typename BaseClass::IndexRangeType IndexRangeType;

  using BaseClass::m;
  using BaseClass::M;

  Odd(int m_, int M_) : BaseClass(m_, M_){BaseClass::a0IsConstant = true;}

  Odd() : BaseClass(){BaseClass::a0IsConstant = true;}

  ///One dimensional bound.
  inline typename capd::jaco::ComplexScalar<ScalarType> bound(const ScalarType& C, const ScalarType& s, const IndexType& k) const{
    if(k.d == 1){
      return typename capd::jaco::ComplexScalar<ScalarType>(0., 2. * C * C * (1. / (2. * s - 1.)) * power(1. / ScalarType(((M + k[0]) * M)), s - 0.5) * ScalarType(-1, 1));
    }else{
      std::cerr << "Bound can be obtained only for one dimension.\n";
      throw std::runtime_error("Bound can be obtained only for one dimension.\n");
    }
  }

  inline int mode2array(const IndexType& k, bool re) const{
    return (k.mode2array(m, re)) / 2;
  }

  inline int mode2arrayTail(const IndexType& k, bool re) const
  {
    IndexType t(k);
    t[t.d()-1]-=m;
    return (t.mode2array(m, re)) / 2;
  }

  inline IndexType array2modeIndex(int i) const{
    //here we multiply i by 2 in order to use array2modeIndex procedure, which is used by real containers (twice as many modes as here)
    return IndexType::array2modeIndex(m, 2 * i);
  }

  inline virtual IndexType array2modeIndex(int i, bool& re) const{
    //here we multiply i by 2 in order to use array2modeIndex procedure, which is used by real containers (twice as many modes as here)
    return IndexType::array2modeIndex(m, 2 * i + 1, re);
  }

  inline bool array2realPart(int i) const{
    return false;
  }

  inline const int modes2arraySize(int m) const{
    IndexType idx;
    int i;
    for(i = 0; i<idx.d(); i++)
      idx[i] = m;
    idx.l = idx.d()-1;
    return (idx.mode2array(m, 0)+1) / 2;
  }

  inline const int modes2arraySize() const{
    return modes2arraySize(m);
  }

  inline static int modes2arraySizeStatic(int m){
    IndexType idx;
    int i;
    for(i = 0; i<idx.d(); i++)
      idx[i] = m;
    idx.l = idx.d()-1;
    return (idx.mode2array(m, 0)+1)/2;
  }

  inline static int modes2realArraySizeStatic(int m){
    IndexType idx;
    int i;
    for(i = 0; i<idx.d(); i++)
      idx[i] = m;
    idx.l = idx.d()-1;
    return (idx.mode2array(m, 0) + 1) / 2;
  }

  template<class AVector>
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
        return ComplexScalar(0, vec[this->mode2array(k, 0)]);
      }
      return ComplexScalar(0, -vec[this->mode2array(-k, 0)]);
    }
  }

  template<class AVector, class ComplexScalarT>
   inline ComplexScalarT mode(const IndexType& k, const AVector& vec) const{
     if(!k.integer()){
       std::cerr << "k="<<k<<"\nOperation of obtaining a mode with index k with non integer k is not defined.\n";
       throw std::runtime_error("Operation of obtaining a mode with index k with non integer k is not defined.\n");
     }
     if(k.isZero()){
       return ComplexScalarT();
     }else{
       if(k.upperHalfspace()){
         return ComplexScalarT(0, vec[this->mode2array(k, 0)]);
       }
       return ComplexScalarT(0, -vec[this->mode2array(-k, 0)]);
     }
   }

  template<typename AVector>
  inline typename capd::jaco::ComplexScalar<typename AVector::ScalarType> closeTailMode( const IndexType& k, const AVector& vec) const
  {
    typedef capd::jaco::ComplexScalar<typename AVector::ScalarType> ComplexScalar;
    if(!k.integer()){
      std::cerr << "Operation of obtaining a mode with index k with non integer k is not defined.\n";
      throw std::runtime_error("Operation of obtaining a mode with index k with non integer k is not defined.\n");
    }
    if(k.storedExplicitly()){
      return ComplexScalar(0, vec[this->mode2arrayTail(k, 0)]);
    }
    return ComplexScalar(0, -vec[this->mode2arrayTail(-k, 0)]);
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

  typename capd::jaco::ComplexScalar<typename VectorType::ScalarType> farTailMode(const ScalarType& s) const{
    return typename capd::jaco::ComplexScalar<ScalarType>(ScalarType(0.), s);
  }

  template<typename AVector>
  inline void setMode(const IndexType& k, const typename capd::jaco::ComplexScalar<typename AVector::ScalarType>& mode, AVector& vec) const{
    vec[mode2array(k, 0)] = mode.im;
  }

  template<class DerivativePairT, class AMatrixT>
  inline void setDerivative(int i, int j, const DerivativePairT& p, AMatrixT& m){
    m[i][j] = p[1][1];
  }

  ///returns two by two identity matrix.
  template<class AMatrixT>
  inline MatrixT id() const{
    AMatrixT m = AMatrixT::Identity(2);
    m[0][0] = 0;
    return m;
  }


  inline void printModesIndex(const VectorType& v, int n, capd::auxil::OutputStream& out) const
  {
    IndexRangeType ir(IndexType::zero());
    ir.setRange(0, strong, n, weak);
    IndexType index;
    for(index = firstModeIndex(ir); !index.limitReached(ir); index.inc(ir)){
      out << "Im(a_"<<index<<"): "<<v[mode2array(index, 0)]<<"\n";
    }
  }

  using BaseClass::firstModeIndex;

};

}
}
#endif /* ODD_H_ */
