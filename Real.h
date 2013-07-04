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
 * DReal.h
 *
 *  Created on: Aug 8, 2011
 *      Author: cyranka
 */

#ifndef DREAL_H_
#define DREAL_H_

#include "Index.h"
#include "ComplexScalar.h"
#include "norms.h"
#include "capd/vectalg/Matrix.hpp"

namespace capd{
namespace jaco{

/**For Real vectors in which real and imaginary part are kept separately.
 */

///this class is base for all classes that operate on the Fourier modes of one dimensional
///real functions with periodic bd. conditions (invariant subspace a_k=\overline{a_{-k}}).
///Modes (Galerkin projection, near tail) are always stored in an array, and this class
///provides methods to retrieve from this array indicated modes and vice versa. Moreover it
///assumes that arrays stores modes of positive index only, whereas those of neagative index
///are retrieved by taking the conjugate of corresponding positive indexed mode.
template<class ScalarT, class IndexT, class NormT,
  class MatrixT = capd::vectalg::Matrix<ScalarT, 0, 0> ,
  class IndexRangeT = capd::jaco::FourierConvolutionIndexRange<IndexT, NormT> >
class Real
{
public:
  typedef capd::jaco::Real<ScalarT, IndexT, NormT, MatrixT, IndexRangeT> Class;
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef ScalarT ScalarType;
  typedef IndexT IndexType;
  typedef IndexRangeT IndexRangeType;
  typedef NormT NormType;
  //typedef typename ScalarType::BoundType BoundType;

  int m, M;
  bool a0IsConstant;

  Real(int m_, int M_) : m(m_), M(M_), a0IsConstant(false){}

  Real(const Real& r) : a0IsConstant(false){
    m = r.m;
    M = r.M;
  }
  ///constructor that uses default dimensions
  Real() : m(0), M(0), a0IsConstant(false){}

  inline int mode2array(const IndexType& k, bool re) const{
    return k.mode2array(m, re);
  }

  inline int mode2arrayTail(const IndexType& k, bool re) const
  {
    IndexType t(k);
    t[t.d()-1]-=m;
    return mode2array(t, re);
  }

  inline IndexType array2modeIndex(int i) const{
    return IndexType::array2modeIndex(m, i);
  }

  virtual inline IndexType array2modeIndex(int i, bool& re) const{
    return IndexType::array2modeIndex(m, i, re);
  }

  inline bool array2realPart(int i) const{
    if(i % 2 == 0)
      return true;
    return false;
  }

  inline const IndexType firstModeIndex() const{
    return array2modeIndex(0);
  }

  ///Returns Index of the first mode which is stored and is within IndexRange ir
  inline const IndexType firstModeIndex(const IndexRangeType& ir, int l = 0) const{
    IndexType idx=IndexType::zero();
    idx.l = l;
//    int i;
//    for(i=0; i<idx.d()-1; i++)
//      idx[i]=ir.returnIndex();
    while(!ir.withinRange(idx)){
      idx.inc(ir);
    }
    return idx;
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

  inline const IndexType lastModeIndex(const IndexRangeType& ir, int l = 0) const{
    IndexType idx = firstModeIndex(ir, l), t;
    while(ir.withinRange(idx)){
      t = idx;
      idx.inc(ir);
    }
    return t;
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
    return (idx.mode2array(m, 0)+1) / 2;
  }

  ///TODO: should be more sharp
  inline const int modes2arraySize() const{
    return modes2arraySize(m);
  }

  inline static int modes2realArraySizeStatic(int m){
    IndexType idx;
    int i;
    for(i = 0; i<idx.d(); i++)
      idx[i] = m;
    idx.l = idx.d()-1;
    return idx.mode2array(m, 0) + 1;
  }

  inline static int modes2arraySizeStatic(int m){
    IndexType idx;
    int i;
    for(i = 0; i<idx.d(); i++)
      idx[i] = m;
    idx.l = idx.d()-1;
    return (idx.mode2array(m, 0) + 1) / 2;
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
//    if(k.isZero()){
//      return ComplexScalar();
//    }else{
      if(k.upperHalfspace()){
        return ComplexScalar(vec[this->mode2array(k, 1)], vec[this->mode2array(k, 0)]);
      }
      return ComplexScalar(vec[this->mode2array(-k, 1)], -vec[this->mode2array(-k, 0)]);
//    }
  }

  template <typename AVector, typename ComplexScalar>
  inline ComplexScalar mode(const IndexType& k, const AVector& vec) const{
    if(!k.integer()){
      std::cerr << "k="<<k<<"\nOperation of obtaining a mode with index k with non integer k is not defined.\n";
      throw std::runtime_error("Operation of obtaining a mode with index k with non integer k is not defined.\n");
    }
//    if(k.isZero()){
//      return ComplexScalar();
//    }else{
      if(k.upperHalfspace()){
        return ComplexScalar(vec[this->mode2array(k, 1)], vec[this->mode2array(k, 0)]);
      }
      return ComplexScalar(vec[this->mode2array(-k, 1)], -vec[this->mode2array(-k, 0)]);
//    }
  }

  typename capd::jaco::ComplexScalar<typename VectorType::ScalarType> farTailMode(const ScalarType& s) const{
    return typename capd::jaco::ComplexScalar<ScalarType>(s, s);
  }

  ///Stores calculated a_k's mode value in internal representation
  template<typename AVector>
  inline void setMode(const IndexType& k, const typename capd::jaco::ComplexScalar<typename AVector::ScalarType>& mode, AVector& vec) const{
    vec[mode2array(k, 1)] = mode.re;
    vec[mode2array(k, 0)] = mode.im;
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

  template<class ComplexScalarT>
  inline void projectOntoImaginarySpace(ComplexScalarT& c){
    c.projectOntoImaginarySpace();
  }

  template<class ComplexScalarT>
  inline void projectOntoRealSpace(ComplexScalarT& c){
    c.projectOntoRealSpace();
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

  template<class DerivativePairT, class AMatrixT>
  inline void setDerivative(int i, int j, const DerivativePairT& p, AMatrixT& m){
    m[2*i][2*j] = p[0][0];
    m[2*i][2*j+1] = p[0][1];
    m[2*i+1][2*j] = p[1][0];
    m[2*i+1][2*j+1] = p[1][1];
  }

  ///returns two by two identity matrix.
  template<class AMatrixT>
  inline AMatrixT id() const{
    AMatrixT m = AMatrixT::Identity(2);
    return m;
  }

  ///calculates sum of all modes supremum with given range,
  ///returns \sum_{k=-m}^{k=m}{\sup{|a_k|}} without a_0 mode.
  template<typename AVector, typename TailType>
  inline typename AVector::ScalarType sumOfSup(const IndexRangeType& ir, const AVector& u, const TailType& tail) const
  {
    ScalarType r=0;
    IndexType index;
    for(index=firstModeIndex(ir); !index.limitReached(ir); index.inc(ir)){
      r += 2*(sqrt(rightBound(mode(index, u, tail).squareNorm())));
    }
    return rightBound(r);
  }

  virtual inline void printModesIndex(const VectorType& v, int n, capd::auxil::OutputStream& out) const
  {
    IndexRangeType ir(IndexType::zero());
    ir.setRange(0, strong, n, weak);
    IndexType index;
    for(index = firstModeIndex(ir); !index.limitReached(ir); index.inc(ir)){
      out << "Re(a_"<<index<<"): ["<<leftBound(v[mode2array(index, 1)])<<","<<rightBound(v[mode2array(index, 1)])<<"] "<<
          " Im(a_"<<index<<"): ["<<leftBound(v[mode2array(index, 0)])<<","<<rightBound(v[mode2array(index, 0)])<<"]\n";
    }
  }

  inline void printModes(const VectorType& v, capd::auxil::OutputStream& out) const{
    int i;
    bool re;
    IndexType index;
    for(i=0; i < v.size(); ++i){
      index = array2modeIndex(i, re);
      if(re){
        out << "Re(a_" << index << ")=";
      }else{
        out << "Im(a_" << index << ")=";
      }
      out << v[i] << ", diam=" << diam(v[i]) << "\n";
    }
  }

  virtual inline void printModesIndex(const VectorType& v, capd::auxil::OutputStream& out) const
  {
    printModesIndex(v, m, out);
  }


  ///translates index in an array storing modes from a Galerkin projection into index of a corresponding mode,
  ///i.e. returns k: Re{a_k} or Im{a_k} is stored in table[i].
  inline int array2mode(int k) const{
    if(k % 2 == 0) return (k + 2) / 2;
      return (k + 1) / 2;
  }

  ///translates index in an array storing modes from a Galerkin projection into index of a corresponding mode,
  ///i.e. returns k: a_k is stored in table[i]. Regarding that we store real and imaginary parts separately,
  ///information whether table[i] stores Re{a_k} is passed into the variable re.
  inline int array2mode(int k, bool& real) const{
    if(k % 2 == 0){
      real = true;
      return (k + 2) / 2;
    }
    real = false;
    return (k + 1) / 2;
  }
};

}}

#endif /* DREAL_H_ */
