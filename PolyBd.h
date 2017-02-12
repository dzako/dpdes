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
    int w; ///<this container is capable of storing G(wM) modes
    FarTailType farTail;
    IndexRangeType irProjection; ///< range of indices from the projection
    IndexRangeType irFiniteTail; ///< range of indices from the finite tail
    IndexRangeType irProjectionPlusFiniteTail; ///< range of indices from the projection plus finite tail
    IndexRangeType irRedundantRange; ///< range of indices from the doubled range of the projection plus finite tail (0,...,2*N)
    IndexRangeType irFull; ///<projection plus finiteTail plus redundantRange

    RealContainerType realContainer;

    bool infiniteDimensional; ///< if container is infinite dimensional

    bool useAbsValues; ///<tells if the absolute values should be returned by this container, !used only by FFT!

    RealPolynomialBound(){}


    RealPolynomialBound(int n_) : ContainerType(modes2arraySizeStatic(n_)), SubspaceType(n_, n_), n(n_),
                                  d(modes2arraySizeStatic(n_)), N(n_), D(modes2arraySizeStatic(n_)), w(0),
                                  //w should be 0 here , because the size of the container, and indices should be calculated for N ( equation is (w+1)*N )
                                  farTail(n, n), realContainer(modes2realArraySizeStatic(n_)),
                                  useAbsValues(false){
      infiniteDimensional = false;
      irProjection.setRange(0, capd::jaco::weak, n, capd::jaco::weak);
      irFiniteTail.setRange(n, capd::jaco::strong, n, capd::jaco::weak);
      irProjectionPlusFiniteTail.setRange(0, capd::jaco::weak, n, capd::jaco::weak);
      irFull.setRange(0, capd::jaco::weak, n, capd::jaco::weak);
    }

    RealPolynomialBound(int n_, int N_, int w_ = 2) : ContainerType(modes2arraySizeStatic((w_+1) * N_)),
                                          SubspaceType(n_, N_), n(n_),
                                          d(modes2arraySizeStatic(n_)), N(N_), D(modes2arraySizeStatic(N_)), w(w_),
                                          farTail(n_, N_), realContainer(modes2realArraySizeStatic(n_)),
                                          useAbsValues(false){
      infiniteDimensional = true;
      irProjection.setRange(0, capd::jaco::weak, n, capd::jaco::weak);
      resetRanges();
    }
    
    /**
     * @param container this DPDEContainer defines the type of function that is stored in the bound (odd, even or none) 
     * @return
     */
    RealPolynomialBound(int n_, int N_, const DPDEContainerType& container, int w_ = 1) :
                                                       ContainerType(modes2arraySizeStatic((w_+1) * N_)), SubspaceType(modes2arraySizeStatic(n_),
                                                       modes2arraySizeStatic(N_)), DPDEContainerType(container),
                                                       n(n_), d(modes2arraySizeStatic(n_)), N(N_), D(modes2arraySizeStatic(N_)), w(w_),
                                                       farTail(n_, N_), realContainer(modes2realArraySizeStatic(n_)),
                                                       useAbsValues(false){
      infiniteDimensional = true;
      irProjection.setRange(0, capd::jaco::weak, n, capd::jaco::weak);
      resetRanges();
    }

    RealPolynomialBound(const RealPolynomialBound& pb){
      *this = pb;
    }

    void resetRanges(){
      irFiniteTail.setRange(n, capd::jaco::strong, N, capd::jaco::weak);
      irProjectionPlusFiniteTail.setRange(0, capd::jaco::weak, N, capd::jaco::weak);
      irRedundantRange.setRange(N, capd::jaco::strong, (w+1) * N, capd::jaco::weak);
      irFull.setRange(0, capd::jaco::weak, (w+1) * N, capd::jaco::weak);
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

    /**Calculates the maximum over norms of the modes within this container.
     */
    inline RealType maxOfNorms() const{
      RealType r(0);
      IndexType index;
      for(index = firstModeIndex(irProjectionPlusFiniteTail); !index.limitReached(irProjectionPlusFiniteTail); index.inc(irProjectionPlusFiniteTail)){
        if(r < RealPolynomialBound::operator[](index).normMax())
          r = RealPolynomialBound::operator[](index).normMax();
      }
      return r;
    }

    ///wrapper function returns index in the array of modes from the lower subspace
    inline int mode2array(const IndexType& k, int dim_ = 0) const{
      int dim = (dim_ == 0 ? (w+1) * N : dim_); //the factor here must be (w+1), because the redundant range dimension is (w+1)*N      m
      bool re=false;
      if(k.upperHalfspace())
        return k.mode2array(dim, re) / 2;
      else
        return (-k).mode2array(dim, re) / 2;
    }

    ///This is used instead of const operator[] to obtain value of a mode in the redundant range
    inline const ComplexScalarType redundantMode(const IndexType& index) const{
      if(irFull.withinRange(index)){
        if(index.upperHalfspace())
          return ContainerType::operator[](mode2array(index));
        return conjugate(ContainerType::operator[](mode2array(index)));
      }else{
        std::cerr << "RealPolynomialBound.redundantMode called with index ("<<index<<") which is not in the redundant range, squareEuclNorm="<<index.squareEuclNorm()<<".\n";
        throw std::runtime_error("RealPolynomialBound.redundantMode called with index which is not in the redundant range.\n");
      }
    }

    inline const ComplexScalarType operator[](const IndexType& index) const{
      ComplexScalarType r;
      if(!farTail.inFarTail(index)){
        if(index.upperHalfspace()){
          r = ContainerType::operator[](mode2array(index));
        }else
            r = conjugate(ContainerType::operator[](mode2array(index)));
      }else{
        r = farTail[index];
      }

      if((*this).useAbsValues){
        return r.abs_supremum();
      }else
        return r;
    }

    /**Not const version, can assign objects, but only indexed by an index in upper halfspace (which is stored explicitly, the
     * rest is obtained by the conjugacy condition a_k=\overline{a_{-k}}.
     */
    inline ComplexScalarType& operator[](const IndexType& index){
      if(irFull.withinRange(index)){
        if(index.upperHalfspace()){
          return ContainerType::operator[](mode2array(index));
        }
        else{
//          if(index.isZero()){
//            return zero[index.l];
//          }else{
            std::cerr << "RealPolynomialBound.operator[]\nForbidden call with index from the lower halfspace.\nIndex=" << index << "\n";
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

    //07.08.2013 CHECKED THIS IS THE SAME AS DEFAULT OPERATOR= CONSTRUCTOR
    //
    //IMPORTANT - this is exactly what the C++ predefined operator is doing
    //check if the C++ predefined operator is calling (DPDEContainer&)*this = (DPDEContainer&)hmc; 
    //and (SubspaceType&)*this = (SubspaceType&)hmc; 
    /*inline RealPolynomialBound& operator=(const RealPolynomialBound& rct){
      (DPDEContainer&)*this = (DPDEContainer&)rct;
      (ContainerType&)*this = (ContainerType&)rct;
      (SubspaceType&)*this = (SubspaceType&)rct;
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
      useAbsValues = rct.useAbsValues;
      return *this;
    }*/

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
      if(irFull.withinRange(index)){
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
        std::cerr << "PolynomialBound.set - Forbidden call\nIndex is out of the range" << index << ", squareEuclNorm="<<index.squareEuclNorm()<<" N=" << N << "\n";
        throw std::runtime_error("PolynomialBound.set - Forbidden call\nIndex is out of the range.\n");
      }
    }


    /**auxilliary function used when casting PolyBd object onto RowVector
     * 
     */
    virtual inline void assignRealContainer(){      
      IndexType index;
      for(index = firstModeIndex(irProjection); !index.limitReached(irProjection); index.inc(irProjection)){        
        setMode(index, this->operator[](index), realContainer);        
      }
    }

    /**
     * The cast operator onto RowVector      
     * 
     */
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
        if(! ((*this)[index]).re.subset(pb[index].re)){
          return false;
        }
        if(! ((*this)[index]).im.subset(pb[index].im)){
          return false;
        }
      }
      return subsetFar(pb);
    }

    int calculateExponent(ScalarType s)const{
      int i=0;
      DOUBLE  t = rightBound(s);
      if(t!=0){
        while(abs(t)<1){
          t*=10;
          i--;
        }
      }
      return i;
    }

    friend std::ostream& operator<<(std::ostream& out, const RealPolynomialBound& pb){
      IndexType index;
      //out << "RealPolynomialBound\n" << (DPDEContainer&)pb << "\n";
      //out << pb.dimension() << "\n";

#if __LATEX_OUT__
        int eRe, eC, eCRe;
        ScalarType re, cRe, rRe, Ct;
        //170205 the latex output implemented here only for purely real modes

        out<< "\\begin{array}{|c|c|}\\hline\\mathbf{k} & \\mathbf{a_k} \\\\ \\hline\\hline\n";
        for(index = pb.firstModeIndex(pb.irProjectionPlusFiniteTail); !index.limitReached(pb.irProjectionPlusFiniteTail); index.inc(pb.irProjectionPlusFiniteTail)){
          re = pb[index].re;
          cRe = mid(re); rRe = re - cRe;
          eRe = pb.calculateExponent(rRe);
          eCRe = pb.calculateExponent(cRe);
          if(eCRe < -2)
            out << index << " & " << rightBound( cRe * pow(10, -eCRe) ) << "\\cdot 10^{" << eCRe << "}";
          else
            out << index << " & " << rightBound( cRe );
          if( eRe < -1)
            out << "+" << rRe * pow(10, -eRe) << "10^{" << eRe << "}\\\\ \n";
          else
            out << "+" << rRe << "\\\\ \n";
        }
        Ct = pb.farTail.getC();
        eC = pb.calculateExponent( Ct );
        if(eC < -1)
          out << "\\geq " << index << " & <" << rightBound(Ct) * pow(10, -eC) << "\\cdot 10^{" << eC << "}/k^{" << rightBound( pb.farTail.getS() ) << "}\\\\\\hline\\end{array}\n\n";
        else
          out << "\\geq " << index << " & \\leq" << rightBound(Ct) << "/k^{" << rightBound( pb.farTail.getS() ) << "}\\\\\\hline\\end{array}\n\n";

#else
      out << "mode index: interval of values\n";
      if(!pb.infiniteDimensional){
        for(index = pb.firstModeIndex(pb.irProjection); !index.limitReached(pb.irProjection); index.inc(pb.irProjection)){
          if(__OUTPUT_MODES__ == 2){
            out << index << "(" << pb.mode2array(index) << ": " << pb[index] << "\n";
          }else{
            if(__OUTPUT_MODES__ == 1){
              out << index << "(" << pb.mode2array(index) << ": " << pb[index].im << "\n";
            }else{
              out << index << "(" << pb.mode2array(index) << ": " << pb[index].re << "\n";
            }
          }
        }
      }else{
        //printing out the whole redundant range
        for(index = pb.firstModeIndex(pb.irProjectionPlusFiniteTail); !index.limitReached(pb.irProjectionPlusFiniteTail); index.inc(pb.irProjectionPlusFiniteTail)){
          if(__OUTPUT_MODES__ == 2){
            out << index << ": " << pb[index] << "\n";
          }else{
            if(__OUTPUT_MODES__ == 1){
              out << index << ": " << pb[index].im << "\n";
            }else{
              out << index << ": " << pb[index].re << "\n";
            }
          }
        }
        /*out << "Redundant range:\n";
        for(index = pb.firstModeIndex(pb.irRedundantRange); !index.limitReached(pb.irRedundantRange); index.inc(pb.irRedundantRange)){
          out << index << ": " << pb.redundantMode(index) << "\n";
        }*/
        out  << pb.farTail << "\n";
//      for debug purposes the redundant range is printed out (the range of modes from the far tail that were calculated explicitely)
//      out << "Redundant range(" << pb.N << " - 2*" << pb.N << "):\n";
//      for(index = pb.firstModeIndex(pb.irRedundantRange); !index.limitReached(pb.irRedundantRange); index.inc(pb.irRedundantRange)){
//        out << index << ": " << pb.redundantMode(index) << "\n";
//      }
      }
#endif
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


    inline const void split(RealPolynomialBound& c, RealPolynomialBound& r) const{
      IndexType index;
      for(index = firstModeIndex(irProjectionPlusFiniteTail); !index.limitReached(irProjectionPlusFiniteTail); index.inc(irProjectionPlusFiniteTail)){
        (*this)[index].split(c[index], r[index]);
        r[index] = r[index].supremum();
      }
      (DPDEContainer&)c = (DPDEContainer&)*this;
      (SubspaceType&)c = (SubspaceType&)*this;
      (DPDEContainer&)r = (DPDEContainer&)*this;
      (SubspaceType&)r = (SubspaceType&)*this;
      c.useAbsValues = false;
      r.useAbsValues = true;
      r.farTail = (*this).farTail;
    }

    inline void abs_supremum(RealPolynomialBound& abs){
      IndexType index;
      for(index = firstModeIndex(irProjectionPlusFiniteTail); !index.limitReached(irProjectionPlusFiniteTail); index.inc(irProjectionPlusFiniteTail)){
        abs[index] = (*this)[index].abs_supremum();
      }
    }


    inline int changeM(int newM, bool guard = false){
      IndexType index;
      if(index.d() == 1){
        if(this->N > newM) { /// new dimension is smaller than current dimension
          ScalarType norm, ntp, max = 0;
          int i;
          for(i = this->N; i >= newM + 1; --i) {
            norm = (*this)[IndexType(i)].norm();
            ntp = norm * power(ScalarType(i), farTail.m_s);
            if(guard) {
              if(!(ntp <= max)){
                if(max < 2 * farTail.m_c)
                  max = ntp;
                else {
                  i = i + 1;
                  break;
                }
              }
            } else {
              if(!(ntp <= max))
                max = ntp;
            }
          }

          newM = i;
          if(newM <= this->N) {
            //we have to estimate rest by newC/k^s
            ContainerType tmp( *this );
            this->resize( modes2arraySizeStatic( (w+1) * newM ) );
            for(int i = 0; i < ((ContainerType)(*this)).dimension(); i++)
              ((ContainerType&)(*this))[i] = tmp[i];
            farTail.setC(max);
            farTail.setM(newM);
            N = newM;
            this->M = modes2arraySizeStatic(newM);
            resetRanges();
          }
          return newM;
        }
        if(this->N < newM) { /// new dimension is larger than current dimension
          ContainerType tmp( *this );
          this->resize( modes2arraySizeStatic( (w+1) * newM) );
          for(int i = 0; i < tmp.dimension(); i++)
            ((ContainerType&)(*this))[i] = tmp[i];
          int oldM = N;
          N = newM;
          this->M = modes2arraySizeStatic(newM);
          resetRanges();
          int i;
          for(i = oldM + 1; i <= newM; ++i) {
            (*this)[  IndexType(i) ] = farTail[ IndexType(i) ];
          }
          farTail.setM(newM);
        }
        return newM;
      }else{
        std::cerr << "PolyBd.changeM not implemented for dimensions higher than 1.\n";
        throw std::runtime_error("PolyBd.changeM not implemented for dimensions higher than 1.\n");
      }
    }


    using SubspaceType::array2realPart;
    using SubspaceType::modes2realArraySizeStatic;
    using SubspaceType::modes2arraySizeStatic;
    using SubspaceType::array2modeIndex;
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
