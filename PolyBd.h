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
 * PolyBd.h
 *
 *  Created on: Dec 6, 2011
 *      Author: cyranka
 */

#ifndef POLYBD_H_
#define POLYBD_H_

#include "capd/vectalg/Container.h"
#include "capd/vectalg/Matrix.h"
#include "Real.h"
#include "Odd.h"
#include "norms.h"
#include "DPDEContainer.h"
#include "tail2.h"

namespace capd{
namespace jaco{


/**Container for l_2 coefficients of real-valued solutions (only modes from the upper halfspace are stored in the container, modes
 * from the lower halfspace are obtained by conjugating modes from the upper halfspace - i.e. the condition a_k=\overline{a_{-k}} ).
 * The Vector stores  modes from the projection plus the modes from the finite part of the tail range (0,...,N), additionally
 * the vector stores modes from the redundant range (N+1,...,2*N), which are used when calculating polynomial bound of the convolution.
 *
 * Using convention that when calling operator[](int index) real parts are on even positions, whereas imaginary parts are on odd
 * positions (when index is even the real part of a mode is returned, when index is odd the imaginary part of a mode is
 * returned). Calling operator[](IndexType) returns the mode as the ComplexScalar.
 */
template<class ComplexT, class IndexT, int DIM,
         class SubspaceT = capd::jaco::Real<typename ComplexT::ScalarType, IndexT, capd::jaco::MaximumNorm<IndexT> > >
class RealPolynomialBound : public capd::vectalg::Vector<ComplexT, DIM>, public SubspaceT,
                            public capd::jaco::DPDEContainer{
public:
    typedef ComplexT ComplexScalarType;
    typedef typename ComplexScalarType::ScalarType RealType;
    typedef RealType ScalarType;
    typedef IndexT IndexType;
    typedef capd::vectalg::Vector<ComplexT, DIM> ContainerType;
    typedef SubspaceT SubspaceType;
    typedef typename SubspaceType::NormType NormType;
    typedef typename SubspaceType::IndexRangeType IndexRangeType;
    typedef typename SubspaceType::MatrixType MatrixType;
    typedef typename MatrixType::RowVectorType RealContainerType;    
    typedef capd::jaco::DPDEContainer DPDEContainerType;
    typedef capd::jaco::FarTail2<ComplexScalarType, SubspaceType> FarTailType;

    ///this has to be changed to myiterator
    typedef typename ContainerType::iterator iterator;
    typedef typename ContainerType::const_iterator const_iterator;
    typedef typename ContainerType::iterator ContainerIterator;
    typedef typename ContainerType::const_iterator const_ContainerIterator;

    static const bool storesLowerHalfspaceIndependently = false;
    int n; ///<this is LARGEST POSSIBLE MAXIMUM NORM OF A INDEX
    int d; ///<this is dimension - number of places needed to store all modes with indexed by set of indices having maximum norm bounded by n
    int N;
    int D;
    FarTailType farTail;
    IndexRangeType irProjection; ///< range of indices from the projection
    IndexRangeType irFiniteTail; ///< range of indices from the finite tail
    IndexRangeType irProjectionPlusFiniteTail; ///< range of indices from the projection plus finite tail
    IndexRangeType irRedundantRange; ///< range of indices from the doubled range of the projection plus finite tail (0,...,2*N)
    IndexRangeType irFull; ///<projection plus finiteTail plus redundantRange

    RealContainerType realContainer;

    bool infiniteDimensional; ///< if container is infinite dimensional

    RealPolynomialBound(){}

//    RealPolynomialBound(int n_, int d_) : ContainerType(d_), SubspaceType(d_, d_), n(n_), d(d_), farTail(n, n),
//                                          realContainer(modes2realArraySizeStatic(n_)){
//      infiniteDimensional = false;
//      irProjection.setRange(0, capd::jaco::weak, n, capd::jaco::weak);
//      irFiniteTail.setRange(n, capd::jaco::strong, n, capd::jaco::weak);
//      irProjectionPlusFiniteTail.setRange(0, capd::jaco::weak, n, capd::jaco::weak);
//      irFull.setRange(0, capd::jaco::weak, n, capd::jaco::weak);
//    }

    RealPolynomialBound(int n_) : ContainerType(modes2arraySizeStatic(n_)), SubspaceType(modes2arraySizeStatic(n_),
                                  modes2arraySizeStatic(n_)), n(n_), d(modes2arraySizeStatic(n_)), farTail(n, n),
                                  realContainer(modes2realArraySizeStatic(n_)){
      infiniteDimensional = false;
      irProjection.setRange(0, capd::jaco::weak, n, capd::jaco::weak);
      irFiniteTail.setRange(n, capd::jaco::strong, n, capd::jaco::weak);
      irProjectionPlusFiniteTail.setRange(0, capd::jaco::weak, n, capd::jaco::weak);
      irFull.setRange(0, capd::jaco::weak, n, capd::jaco::weak);
    }

    RealPolynomialBound(int n_, int N_) : ContainerType(modes2arraySizeStatic(2 * N_)),
                                          SubspaceType(modes2arraySizeStatic(n_), modes2arraySizeStatic(N_)), n(n_),
                                          d(modes2arraySizeStatic(n_)), N(N_), D(modes2arraySizeStatic(N_)),
                                          farTail(n_, N_), realContainer(modes2realArraySizeStatic(n_)){
      infiniteDimensional = true;
      irProjection.setRange(0, capd::jaco::weak, n, capd::jaco::weak);
      irFiniteTail.setRange(n, capd::jaco::strong, N, capd::jaco::weak);
      irProjectionPlusFiniteTail.setRange(0, capd::jaco::weak, N, capd::jaco::weak);
      irRedundantRange.setRange(N, capd::jaco::strong, 2 * N, capd::jaco::weak);
      irFull.setRange(0, capd::jaco::weak, 2 * N, capd::jaco::weak);
    }
    
    /**
     * @param container this DPDEContainer defines the type of function that is stored in the bound (odd, even or none) 
     * @return
     */
    RealPolynomialBound(int n_, int N_, const DPDEContainerType& container) : 
                                                       ContainerType(modes2arraySizeStatic(2 * N_)), SubspaceType(modes2arraySizeStatic(n_), 
                                                       modes2arraySizeStatic(N_)), DPDEContainerType(container),
                                                       n(n_), d(modes2arraySizeStatic(n_)), N(N_), D(modes2arraySizeStatic(N_)),
                                                       farTail(n_, N_), realContainer(modes2realArraySizeStatic(n_)){
      infiniteDimensional = true;
      irProjection.setRange(0, capd::jaco::weak, n, capd::jaco::weak);
      irFiniteTail.setRange(n, capd::jaco::strong, N, capd::jaco::weak);
      irProjectionPlusFiniteTail.setRange(0, capd::jaco::weak, N, capd::jaco::weak);
      irRedundantRange.setRange(N, capd::jaco::strong, 2 * N, capd::jaco::weak);
      irFull.setRange(0, capd::jaco::weak, 2 * N, capd::jaco::weak);
    }

    RealPolynomialBound(const RealPolynomialBound& pb){
      *this = pb;
    }


    /**Calculates sum of maximum norms of all modes within the container (but without the far tail)
     */
    inline RealType sumOfNorms() const{
      RealType r(0);
      IndexType index;
      for(index = firstModeIndex(irProjectionPlusFiniteTail); !index.limitReached(irProjectionPlusFiniteTail); index.inc(irProjectionPlusFiniteTail)){
        r += 2 * RealPolynomialBound::operator[](index).normMax(); //multiply by two, because modes are real
      }
      return r;
    }

    ///wrapper function returns index in the array of modes from the lower subspace
    inline int mode2array(const IndexType& k) const{
      bool re=false;
      if(k.upperHalfspace())
        return k.mode2array(n, re) / 2;
      else
        return (-k).mode2array(n, re) / 2;
    }

    ///This is used instead of const operator[] to obtain value of a mode in the redundant range
    inline const ComplexScalarType redundantMode(const IndexType& index) const{
      if(withinRedundantRange(index)){
        if(index.upperHalfspace())
          return ContainerType::operator[](mode2array(index));
        return conjugate(ContainerType::operator[](mode2array(index)));
      }else{
        std::cerr << "RealPolynomialBound.redundantMode called with index ("<<index<<") which is not in the redundant range.\n";
        throw std::runtime_error("RealPolynomialBound.redundantMode called with index which is not in the redundant range.\n");
      }
    }

    /**If the given index is within redundant range (modes that are in far tail), but may be calculated explicitely, as when calculating
     * convolution of two polynomial bounds.
     */
    inline bool withinRedundantRange(const IndexType& index) const{
      if(infiniteDimensional){
        if(index.squareEuclNorm() > 4. * N * N)
          return false;
        return true;
      }else{
        if(index.squareEuclNorm() > n * n)
          return false;
        return true;
      }
    }

    inline bool withinFiniteTail(const IndexType& index) const{
      if(irFiniteTail.withinRange(index))
        return true;
      return false;
    }

    inline const ComplexScalarType operator[](const IndexType& index) const{
      if(!farTail.inFarTail(index)){
        if(index.upperHalfspace()){
          return ContainerType::operator[](mode2array(index));
        }else
            return conjugate(ContainerType::operator[](mode2array(index)));
      }else{
        if(!infiniteDimensional)
          return ComplexScalarType(0.);
        return farTail[index];
      }
    }

    /**Not const version, can assign objects, but only indexed by an index in upper halfspace (which is stored explicitely, the
     * rest is obtained by the conjugacy condition a_k=\overline{a_{-k}}.
     */
    inline ComplexScalarType& operator[](const IndexType& index){
      if(withinRedundantRange(index)){
        if(index.upperHalfspace()){
          return ContainerType::operator[](mode2array(index));
        }
        else{
//          if(index.isZero()){
//            return zero[index.l];
//          }else{
            std::cerr << "RealPolynomialBound.operator[]\nForbidden call with index from the lower halfspace.\n";
            throw std::runtime_error("RealPolynomialBound.operator[]\nForbidden call with index from the lower halfspace.\n");
//          }
        }
      }else{
        std::cerr << "RealPolynomialBound.operator[]\nForbidden call with index out of the range (Index=" << index << ".\n";
        throw std::runtime_error("RealPolynomialBound.operator[]\nForbidden call with index out of the range.\n");
      }
    }

    inline const RealType operator[](int i) const{
      IndexType index = array2modeIndex(i);
      bool re = array2realPart(i);
      if(re){ ///REAL PART is returned
        return (ContainerType::operator[](mode2array(index))).re;
      }else{ ///IMAGINARY PART is returned
        return (ContainerType::operator[](mode2array(index))).im;
      }
    }

    inline RealType& operator[](int i){
      IndexType index = array2modeIndex(i);
      bool re = array2realPart(i);
      if(re){ ///REAL PART is returned
        return (ContainerType::operator[](mode2array(index))).re;
      }else{ ///IMAGINARY PART is returned
        return (ContainerType::operator[](mode2array(index))).im;
      }
    }

    inline RealPolynomialBound& operator*=(const ComplexScalarType& c){
      if(c.isImUnit())
        multiplyByImUnit();
      (ContainerType&)*this *= c;
      if(infiniteDimensional)
        farTail.m_c *= c.normMax();
      return *this;
    }

    inline RealPolynomialBound& operator*=(const RealType& c){
      (ContainerType&)*this *= c;
      if(infiniteDimensional)
        farTail.m_c *= rightBound(abs(c));
      return *this;
    }

    inline RealPolynomialBound& operator+=(const RealPolynomialBound& pb){
      (DPDEContainerType&)*this += (DPDEContainerType&)pb;
      if(infiniteDimensional == pb.infiniteDimensional)
        (ContainerType&)*this += (ContainerType&)pb;
      else{
        IndexType index;
        for(index = firstModeIndex(irProjection); !index.limitReached(irProjection); index.inc(irProjection)){
          (*this)[index] += pb[index];
        }
      }
      return *this;
    }


    inline RealPolynomialBound& operator+=(const RealContainerType& rct){
      IndexType index;
      for(index = firstModeIndex(irProjection); !index.limitReached(irProjection); index.inc(irProjection)){
        (*this)[index] += mode(index, rct);
      }
      return *this;
    }

    //
    inline RealPolynomialBound& operator=(const RealPolynomialBound& rct){
      (DPDEContainer&)*this = (DPDEContainer&)rct;
      (ContainerType&)*this = (ContainerType&)rct;
      n = rct.n;
      N = rct.N;
      d = rct.d;
      D = rct.D;
      irProjection = rct.irProjection;
      irProjectionPlusFiniteTail = rct.irProjectionPlusFiniteTail;
      irFiniteTail = rct.irFiniteTail;
      irRedundantRange = rct.irRedundantRange;
      irFull = rct.irFull;
      realContainer = rct.realContainer;
      infiniteDimensional = rct.infiniteDimensional;
      farTail = rct.farTail;
      return *this;
    }

    inline RealPolynomialBound& operator=(const RealContainerType& rct){
      IndexType index;
      for(index = firstModeIndex(irProjection); !index.limitReached(irProjection); index.inc(irProjection)){
//        if(!index.isZero())
          (*this)[index] = mode(index, rct);
//        else
//          (*this)[index] = 0.;
      }
      return *this;
    }

    inline RealPolynomialBound& copyFinitePartFrom(const RealPolynomialBound& pb){
      (DPDEContainer&)*this = (DPDEContainer&)pb;
      IndexType index;      
      for(index = firstModeIndex(irProjection); !index.limitReached(irProjection); index.inc(irProjection)){
        (*this)[index] = pb[index];
      }
      return *this;
    }

    inline RealPolynomialBound& copyTailPartFrom(const RealPolynomialBound& pb){
      (DPDEContainer&)*this = (DPDEContainer&)pb;
      IndexType index;
      for(index = firstModeIndex(irFiniteTail); !index.limitReached(irFiniteTail); index.inc(irFiniteTail)){
        (*this)[index] = pb[index];
      }
      setC(*this, C(pb));
      setS(*this, s(pb));
      return *this;
    }

    inline RealPolynomialBound& operator=(const ComplexScalarType& c){
      (ContainerType&)*this = c;
      return *this;
    }

    inline void set(const IndexType& index, const ComplexScalarType& val){
      if(withinRedundantRange(index)){
        if(index.upperHalfspace()){
          ContainerType::operator[](mode2array(index)) = val;
        }else{
//          if(index.isZero()){
//            zero[index.l] = val;
//          }else{
            ContainerType::operator[](mode2array(index)) = conjugate(val);
//          }
        }
      }else{
        std::cerr << "PolynomialBound.set - Forbidden call\nIndex is out of the range" << index << ", N=" << N << "\n";
        throw std::runtime_error("PolynomialBound.set - Forbidden call\nIndex is out of the range.\n");
      }
    }

    virtual inline void assignRealContainer(){
      IndexType index;
      for(index = firstModeIndex(irProjection); !index.limitReached(irProjection); index.inc(irProjection)){
        setMode(index, this->operator[](index), realContainer);
      }
    }

    inline operator const RealContainerType&(){
      assignRealContainer();
      return realContainer;
    }

    inline void projectOntoSubspace(){
      IndexType index;

      if(baseReZero){
        for(index = firstModeIndex(irFull); !index.limitReached(irFull); index.inc(irFull)){
          projectOntoImaginarySpace((*this)[index]);
        }
      }
      if(baseImZero){
        for(index = firstModeIndex(irFull); !index.limitReached(irFull); index.inc(irFull)){
          projectOntoRealSpace((*this)[index]);
        }
      }
    }

    ///TODO: check what to return in case of infinite dimensional container
    inline const int dimension() const{
      return d;
    }

    void cleanFinitePart(){
      IndexType index;
      for(index = firstModeIndex(irProjection); !index.limitReached(irProjection); index.inc(irProjection)){
        this->operator[](index) = 0;
      }
    }

    void cleanTail(){
      IndexType index;
      for(index = firstModeIndex(irFiniteTail); !index.limitReached(irFiniteTail); index.inc(irFiniteTail)){
        this->operator[](index) = 0;
      }
      for(index = firstModeIndex(irRedundantRange); !index.limitReached(irRedundantRange); index.inc(irRedundantRange)){
        this->operator[](index) = 0;
      }
      if(infiniteDimensional) setC(*this, 0.);
    }

    void clean(){
      (ContainerType&)*this = ComplexScalarType(0.);
      if(infiniteDimensional) setC(*this, 0.);
    }

    ///this returns index of the position in internal array on which the finite tail starts
    inline int finiteTailBegin() const{
      IndexType i = firstModeIndex(irFiniteTail);
      return mode2array(i, 1);
    }

    ///this returns index of the position in internal array on which the finite tail ends
    inline int finiteTailEnd() const{
      IndexType i = lastModeIndex(irFiniteTail);
      return mode2array(i, 1);
    }

    inline bool subsetFar(const RealPolynomialBound& pb) const{
      if(farTail.subset(pb.farTail)) return true;
      return false;
    }

    inline bool subset(const RealPolynomialBound& pb) const{
      IndexType index;
      for(index = firstModeIndex(irProjectionPlusFiniteTail); !index.limitReached(irProjectionPlusFiniteTail); index.inc(irProjectionPlusFiniteTail)){
        if(! ((*this)[index]).re.subset(pb[index].re))
          return false;
        if(! ((*this)[index]).im.subset(pb[index].im))
          return false;
      }
      return subsetFar(pb);
    }

    friend std::ostream& operator<<(std::ostream& out, const RealPolynomialBound& pb){
      IndexType index;
      out << "RealPolynomialBound\n" << (DPDEContainer&)pb << "\n";
      if(!pb.infiniteDimensional){
        for(index = pb.firstModeIndex(pb.irProjection); !index.limitReached(pb.irProjection); index.inc(pb.irProjection)){
          out << index << ": " << pb[index] << "\n";
        }
      }else{
        for(index = pb.firstModeIndex(pb.irProjectionPlusFiniteTail); !index.limitReached(pb.irProjectionPlusFiniteTail); index.inc(pb.irProjectionPlusFiniteTail)){
          out << index << ": " << pb[index] << "\n";
        }
        out << "\nFarTail(" << pb.n << " - " << pb.N << "):\n" << pb.farTail << "\n";
//      for debug purposes the redundant range is printed out (the range of modes from the far tail that were calculated explicitely)
//      out << "Redundant range(" << pb.N << " - 2*" << pb.N << "):\n";
//      for(index = pb.firstModeIndex(pb.irRedundantRange); !index.limitReached(pb.irRedundantRange); index.inc(pb.irRedundantRange)){
//        out << index << ": " << pb.redundantMode(index) << "\n";
//      }
      }
      return out;
    }

    inline friend const RealType C(const RealPolynomialBound& pb){
      return pb.farTail.m_c;
    }

    inline friend int s(const RealPolynomialBound& pb){
      return pb.farTail.m_s;
    }

    inline friend void setC(RealPolynomialBound& pb, const RealType& c){
      if(!pb.infiniteDimensional){
        std::cerr << "RealPolynomialBound.setC() - Forbidden call\nWorks only with infinite dimensional containers.\n";
        throw std::runtime_error("RealPolynomialBound.setC() - Forbidden call\nWorks only with infinite dimensional containers.\n");
      }
      pb.farTail.m_c = c;
    }

    inline friend void setCLarger(RealPolynomialBound& pb, const RealType& c){
      RealType newC = RealType(c.rightBound()) + RealType(0., diam(c).rightBound());
      pb.farTail.m_c = newC;
    }

    inline friend void setS(RealPolynomialBound& pb, int s){
      if(!pb.infiniteDimensional){
        std::cerr << "RealPolynomialBound.setS() - Forbidden call\nWorks only with infinite dimensional containers.\n";
        throw std::runtime_error("RealPolynomialBound.setS() - Forbidden call\nWorks only with infinite dimensional containers.\n");
      }
      pb.farTail.m_s = s;
    }

    template<typename DoubleT>
    inline void setLeftBound( int i, DoubleT d){
      (*this)[i].setLeftBound(d);
    }

    inline friend void copyFinitePart(const RealPolynomialBound& pb, RealContainerType& r){
      r = RealContainerType(RealPolynomialBound::modes2realArraySizeStatic(pb.n));
      IndexType index;
      for(index = pb.firstModeIndex(pb.irProjection); !index.limitReached(pb.irProjection); index.inc(pb.irProjection)){
        pb.setMode(index, pb[index], r);
      }
    }

    inline friend void copyFinitePart(const RealContainerType& r, RealPolynomialBound& pb){
      IndexType index;
      for(index = pb.firstModeIndex(pb.irProjection); !index.limitReached(pb.irProjection); index.inc(pb.irProjection)){
        pb[index] = pb.mode(index, r);
      }
    }

    template<typename DoubleT>
    inline void setRightBound( int i, DoubleT d){
      (*this)[i].setRightBound(d);
    }

    using SubspaceType::array2realPart;
    using SubspaceType::modes2realArraySizeStatic;
    using SubspaceType::modes2arraySizeStatic;
    using SubspaceType::array2modeIndex;
    using SubspaceType::mode2array;
    using ContainerType::begin;
    using ContainerType::end;
    using ContainerType::size;
    using ContainerType::resize;
    using ContainerType::operator();
    using SubspaceType::firstModeIndex;
    using SubspaceType::setMode;
    using SubspaceType::mode;
    using SubspaceType::projectOntoImaginarySpace;
    using SubspaceType::projectOntoRealSpace;
    using SubspaceType::lastModeIndex;    

};


}
}

#endif /* POLYBD_H_ */
