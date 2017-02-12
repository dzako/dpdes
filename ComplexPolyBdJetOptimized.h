/*
 * HalfspaceModesContainer.h
 *
 *  Created on: Oct 28, 2011
 *      Author: cyranka
 */

#ifndef HALFSPACEMODESCONTAINER_H_
#define HALFSPACEMODESCONTAINER_H_

#include "capd/vectalg/Vector.h"
#include "capd/vectalg/Matrix.h"
#include "Real.h"
#include "DPDEContainer.h"

namespace capd{
namespace jaco{

/**
 * Finite dimensional vector (optimized for jets)
 * D is an overestimate of the vector size.
 */
template< class IntervalT, class ScalarT, class IndexT, /*int N,*/ int D,          
          class SubspaceT = capd::jaco::Real<IntervalT, IndexT, capd::jaco::MaximumNorm<IndexT> > >
class ComplexPolyBdJetOptimized : public SubspaceT, public capd::jaco::DPDEContainer{

public:
  typedef IntervalT IntervalType;
  typedef ScalarT ScalarType;
  typedef IndexT IndexType;
  typedef SubspaceT SubspaceType;
  typedef typename SubspaceType::MatrixType RealMatrixType;
  typedef capd::jaco::DPDEContainer DPDEContainerType;
  typedef typename SubspaceType::IndexRangeType IndexRangeType;
  typedef capd::vectalg::Vector<ScalarType, D> VectorType;
  typedef typename RealMatrixType::RowVectorType RealVectorType;
  typedef capd::vectalg::Matrix<IntervalType, 2, 2> QuadMatrixType;
  typedef typename VectorType::iterator Iterator;

  static const bool storesLowerHalfspaceIndependently = true;

  int n; ///<this is LARGEST POSSIBLE MAXIMUM NORM OF A INDEX
  int d; ///<this is dimension - number of places needed to store all modes having maximum norm bounded by n
  VectorType m_lowerHalfspace; ///<lower subspace of modes
  VectorType m_upperHalfspace; ///<upper subspace of modes
  RealVectorType realContainer;
  Iterator iu;
  Iterator il;

  ScalarType empty;///< this is an empty scalar (with zeros on all elements), this is useful, when ScalarType is a Jet scalar
  IndexRangeType irProjection;
  IndexRangeType irFiniteTail;
  IndexRangeType irProjectionPlusFiniteTail;
  IndexRangeType irRedundantRange;
  IndexRangeType irFull;
  IndexRangeType ir;
  bool infiniteDimensional;

  ComplexPolyBdJetOptimized() : infiniteDimensional(false){
  }

  ComplexPolyBdJetOptimized(int n_, int d_) : SubspaceType(n_, n_), n(n_), d(d_), m_lowerHalfspace(d_), m_upperHalfspace(d_),
      realContainer(modes2realArraySizeStatic(n_)), infiniteDimensional(false){
    irProjection.setRange(0, capd::jaco::weak, n, capd::jaco::weak);
    irFiniteTail.setRange(n, capd::jaco::strong, n, capd::jaco::weak);
    irProjectionPlusFiniteTail.setRange(0, capd::jaco::weak, n, capd::jaco::weak);
    irRedundantRange.setRange(n, capd::jaco::strong, n, capd::jaco::weak);
    irFull.setRange(0, capd::jaco::weak, n, capd::jaco::weak);
    ir = irProjection;
  }

  ComplexPolyBdJetOptimized(int n_) : SubspaceType(n_, n_), n(n_), d(modes2arraySizeStatic(n_)),
      m_lowerHalfspace(modes2arraySizeStatic(n_)), m_upperHalfspace(modes2arraySizeStatic(n_)), realContainer(modes2realArraySizeStatic(n_)), infiniteDimensional(false){
    irProjection.setRange(0, capd::jaco::weak, n, capd::jaco::weak);
    irFiniteTail.setRange(n, capd::jaco::strong, n, capd::jaco::weak);
    irProjectionPlusFiniteTail.setRange(0, capd::jaco::weak, n, capd::jaco::weak);
    irRedundantRange.setRange(n, capd::jaco::strong, n, capd::jaco::weak);
    irFull.setRange(0, capd::jaco::weak, n, capd::jaco::weak);
    ir = irProjection;
  }

  ///dummy function
  inline IntervalType sumOfNorms() const{
    return 0;
  }

  ///dummy function
  inline IntervalType maxOfNorms() const{
    return 0;
  }

  ///wrapper function returns index in the array of modes from the lower subspace
  inline int mode2array(const IndexType& k) const{
    bool re=false;
    if(k.upperHalfspace())
      return k.mode2array(n, re) / 2;
    else
      return (-k).mode2array(n, re) / 2;
  }

  inline const VectorType& lowerHalfspace() const{
    return m_lowerHalfspace;
  }

  inline const VectorType& upperHalfspace() const{
    return m_upperHalfspace;
  }

  inline const ScalarType& operator[](const IndexType& index) const{
    if(index.upperHalfspace())
      return m_upperHalfspace[mode2array(index)];
    else
//      if(index.isZero())
//        return m_zero[index.l]; ///temporary solution with zero indexed mode
//      else
        return m_lowerHalfspace[mode2array(index)];
  }

  inline ScalarType& operator[](const IndexType& index){
    if(index.upperHalfspace())
      return m_upperHalfspace[mode2array(index)];
    else
//      if(index.isZero())
//        return m_zero[index.l]; ///temporary solution with zero indexed mode
//      else
        return m_lowerHalfspace[mode2array(index)];
  }

  ScalarType redundantMode(const IndexType& index){
    return ScalarType(0);
  }

  inline ComplexPolyBdJetOptimized& operator*=(const ScalarType& s){    
    m_upperHalfspace *= s;
    m_lowerHalfspace *= conjugate(s);
    if(s.isImUnit())
      multiplyByImUnit();
    return *this;
  }

  inline ComplexPolyBdJetOptimized& operator*=(const IntervalType& i){
    iu = m_upperHalfspace.begin(),
    il = m_lowerHalfspace.begin();
    while(iu != m_upperHalfspace.end()){
      *iu = i * (*iu);
      *il = i * (*il);
      iu++; il++;
    }
    return *this;
  }

  /*inline ComplexPolyBdJetOptimized& operator*=(double d){
    iu = m_upperHalfspace.begin(),
    il = m_lowerHalfspace.begin();
    while(iu != m_upperHalfspace.end()){
      *iu = d * (*iu);
      *il = d * (*il);
      iu++; il++;
    }
    return *this;
  }*/

  inline ComplexPolyBdJetOptimized& operator+=(const RealVectorType& v){    
    IndexType index;    
    for(index = firstModeIndex(ir); !index.limitReached(ir); index.inc(ir)){
      (*this)[index] += mode(index, v);
    }
    return *this;
  }

  inline ComplexPolyBdJetOptimized& operator+=(const ComplexPolyBdJetOptimized& hmc){    
    (DPDEContainerType&)*this += (DPDEContainerType&)hmc;
    m_upperHalfspace += hmc.m_upperHalfspace;
    m_lowerHalfspace += hmc.m_lowerHalfspace;
    return *this;
  }

  inline ComplexPolyBdJetOptimized& operator=(const ScalarType& s){
    m_upperHalfspace = s;
    m_lowerHalfspace = conjugate(s);
    return *this;
  }

  template<class ComplexScalarT>
  inline ComplexPolyBdJetOptimized& operator*=(const ComplexScalarT& c){
    iu = m_upperHalfspace.begin(),
    il = m_lowerHalfspace.begin();
    while(iu != m_upperHalfspace.end()){
      *iu = c * (*iu);
      *il = c * (*il);
      iu++; il++;
    }
    if(c.isImUnit())
      multiplyByImUnit();
    return *this;
  }

//  inline HalfspaceModesContainer& operator=(const IntervalType& i){
//    iu = m_upperHalfspace.begin(),
//    il = m_lowerHalfspace.begin();
//    while(iu != m_upperHalfspace.end()){
//      *iu = i;
//      *il = i;
//      iu++; il++;
//    }
//    return *this;
//  }

  /**This function casts the current instance of ComplexPolyBdJetOptimized to a RealVectorType 
   * (CAPD compatible format - containing intervals) stored in the 'realContainer' variable
   * 
   * IMPORTANT: this function from each jet takes out .secondFreeCoeff() value, which is 
   * \Phi(m([x])), i.e. the value of the numerical method calculated on the middle-point,
   * this value does not have any influence on the variational part \frac{\Phi}{\partial x}.
   *    
   * Then this operator is used in FFTDynSys.current() in order to assign the current \Phi(m([x]))
   * stored in jets to a vector o_phi <- the internal CAPD representation, which is then input 
   * to the Lohner algorithm.  
   * 
   * FFTDynSys is overloading the encloseMap procedure which is providing all input
   * (\Phi(m([x])) and \frac{\Phi}{\partial x}) to the Lohner algorithm.
   *    
   */
  inline void assignRealContainer(){
    IndexType index;
    for(index = firstModeIndex(ir); !index.limitReached(ir); index.inc(ir)){      
      setMode(index, (*this)[index].secondFreeCoeff(), realContainer);
    }
  }

  inline operator const RealVectorType&(){
    assignRealContainer();    
    return realContainer;
  }
  
  inline void monodromyMatrix(RealMatrixType& m){
    int i, j;
    QuadMatrixType p;
    for(i=0; i < d; ++i)
      for(j=0; j < d; ++j){
        if(a0IsConstant && array2modeIndex(i).isZero() && array2modeIndex(j).isZero() && array2modeIndex(i).l == array2modeIndex(j).l){
        //if(i == j && j == 0 && a0IsConstant){
          p[0][0] = 1; p[0][1] = 0; p[1][0] = 0; p[1][1] = 1;
        }else
          p = m_upperHalfspace[i].variationalPart(j);
        setDerivative(i, j, p, m); 
        //method from the SubspaceType must be called 
        //\frac{\partial a_i}{\partial a_j} is set here for matrix m
        //(for only real modes \frac{\partial a_i}{\partial a_j} is a 2x2 matrix, for odd,even 
        //\frac{\partial a_i}{\partial a_j} is only a 1x1 matrix)
      }
  }

  inline void monodromyMatrix2(const RealMatrixType& m){
    int i, j, indexr, indexc;
    QuadMatrixType p;
    for(i = -d + 1; i < d; ++i){
      for(j=0; j < d; ++j){
        p = operator[](IndexType(i)).variationalPart(j);
        indexr = (i >= 0 ? i : m.numberOfRows()/2 + i);
        m[2*indexr][2*j] = p[0][0];
        m[2*indexr][2*j+1] = p[0][1];
        m[2*indexr+1][2*j] = p[1][0];
        m[2*indexr+1][2*j+1] = p[1][1];
      }
      for(j = - d + 1; j < 0; ++j){
        p = operator[](IndexType(i)).variationalPart( - j);
        indexr = (i >= 0 ? i : m.numberOfRows()/2 + i);
        indexc =m.numberOfColumns()/2 + j;
        m[2*indexr][2*indexc] = p[0][0];
        m[2*indexr][2*indexc+1] = p[0][1];
        m[2*indexr+1][2*indexc] = p[1][0];
        m[2*indexr+1][2*indexc+1] = p[1][1];
      }
    }
  }

  inline ComplexPolyBdJetOptimized& operator=(const RealVectorType& rvt){
    IndexType index;
    for(index = firstModeIndex(ir); !index.limitReached(ir); index.inc(ir)){        
        empty.setFreeCoeff(mode(index, rvt));
        set(index, empty); //here we need to set a ScalarType with zero variational part (empty is ok)
        if(index.isZero()){
          if(!a0IsConstant) setVariationalPartToId(index, index);
        }else
          setVariationalPartToId(index, index);
    }
    return *this;
  }

  //check if the C++ predefined operator is calling (DPDEContainer&)*this = (DPDEContainer&)hmc; 
  //and (SubspaceType&)*this = (SubspaceType&)hmc; 
  /*inline ComplexPolyBdJetOptimized& operator=(const ComplexPolyBdJetOptimized& hmc){
    (DPDEContainer&)*this = (DPDEContainer&)hmc;
    (SubspaceType&)*this = (SubspaceType&)hmc;
    n = hmc.n;
    d = hmc.d;
    irProjection = hmc.irProjection;
    irFiniteTail = hmc.irFiniteTail;
    irProjectionPlusFiniteTail = hmc.irProjectionPlusFiniteTail;
    irRedundantRange = hmc.irRedundantRange;
    irFull = hmc.irFull;
    ir = hmc.ir;
    realContainer = hmc.realContainer;
    m_upperHalfspace = hmc.m_upperHalfspace;
    m_lowerHalfspace = hmc.m_lowerHalfspace;
    return *this;
  }*/

  /**Sets a mode, determined by index, and its conjugate to the new value.
   *
   * @param index index of a mode
   * @param val new value of modes
   */
  inline void set(const IndexType& index, const ScalarType& val){
    if(index.upperHalfspace()){
      m_upperHalfspace[mode2array(index)] = val;
      m_lowerHalfspace[mode2array(index)] = conjugate(val);
    }else{
//      if(index.isZero()){
//        m_zero[index.l] = val;
//      }else{
        m_upperHalfspace[mode2array(index)] = conjugate(val);
        m_lowerHalfspace[mode2array(index)] = val;
//      }
    }
  }

  template<class ComplexScalarT>
  inline void setFreeCoeff(const IndexType& index, const ComplexScalarT& val){
    if(index.upperHalfspace()){
      m_upperHalfspace[mode2array(index)].setFreeCoeff(val);
      m_lowerHalfspace[mode2array(index)].setFreeCoeff(conjugate(val));
    }else{
//      if(index.isZero()){
//        m_zero[index.l].setFreeCoeff(val);
//      }else{
        m_upperHalfspace[mode2array(index)].setFreeCoeff(conjugate(val));
        m_lowerHalfspace[mode2array(index)].setFreeCoeff(val);
//      }
    }
  }

  template<class ComplexScalarT>
  inline void setSecondFreeCoeff(const IndexType& index, const ComplexScalarT& val){
    if(index.upperHalfspace()){
      m_upperHalfspace[mode2array(index)].setSecondFreeCoeff(val);
      m_lowerHalfspace[mode2array(index)].setSecondFreeCoeff(conjugate(val));
    }else{
//      if(index.isZero()){
//        m_zero[index.l].setSecondFreeCoeff(val);
//      }else{
        m_upperHalfspace[mode2array(index)].setSecondFreeCoeff(conjugate(val));
        m_lowerHalfspace[mode2array(index)].setSecondFreeCoeff(val);
//      }
    }
  }

  inline void setSecondFreeCoeffs(const RealVectorType& rvt){
    IndexType index;
    for(index = firstModeIndex(ir); !index.limitReached(ir); index.inc(ir)){
      setSecondFreeCoeff(index, mode(index, rvt));
    }
  }


  inline void setConjugate(const IndexType& index, const ScalarType& val){
    if(index.upperHalfspace()){
      m_upperHalfspace[mode2array(index)] = conjugate(val);
      m_lowerHalfspace[mode2array(index)] = val;
    }else{
//      if(index.isZero()){
//        m_zero[index.l] = val;
//      }else{
        m_upperHalfspace[mode2array(index)] = val;
        m_lowerHalfspace[mode2array(index)] = conjugate(val);
//     }
    }
  }

  /**Sets the variational part to Identity (assigns two-dimensional Identity matrix at the position (index1, index2)).
   * ASSUMES that the variational part is zero when calling this function.
   */
  inline void setVariationalPartToId(const IndexType& index1, const IndexType& index2){
    if(index1.upperHalfspace() && index2.upperHalfspace()){
      m_upperHalfspace[mode2array(index1)].setVariationalPartToId(mode2array(index2)/*, SubspaceType::template id<QuadMatrixType>()*/);
      m_lowerHalfspace[mode2array(index1)].setVariationalPartToIdConjugate(mode2array(index2)/*, SubspaceType::template id<QuadMatrixType>()*/);
      return;
    }
    if(!index1.upperHalfspace() && !index2.upperHalfspace()){
      m_upperHalfspace[mode2array(index1)].setVariationalPartToIdConjugate(mode2array(index2)/*, SubspaceType::template id<QuadMatrixType>()*/);
      m_lowerHalfspace[mode2array(index1)].setVariationalPartToId(mode2array(index2)/*, SubspaceType::template id<QuadMatrixType>()*/);
      return;
    }
    std::cerr << "the setVariationalPart function is not defined for the given position.\n";
    throw std::runtime_error("the setVariationalPart function is not defined for the given position.\n");
  }

  inline void projectOntoSubspace(){
    iu = m_upperHalfspace.begin(),
    il = m_lowerHalfspace.begin();
    if(solutionType == realValued){
      if(subspaceType == odd){
        while(iu != m_upperHalfspace.end()){
          projectOntoImaginarySpace(*iu);
          projectOntoImaginarySpace(*il);
          iu++; il++;
        }
      }
      if(subspaceType == even){
        while(iu != m_upperHalfspace.end()){
          projectOntoRealSpace(*iu);
          projectOntoRealSpace(*il);
          iu++; il++;
        }
      }
    }
  }

  inline void setImaginaryPartToZero(){
    iu = m_upperHalfspace.begin(),
    il = m_lowerHalfspace.begin();
    while(iu != m_upperHalfspace.end()){
      (*iu).setImaginaryPartToZero();
      (*il).setImaginaryPartToZero();
      iu++; il++;
    }
  }

  inline const int dimension() const{
    return d;
  }

  inline void split(ComplexPolyBdJetOptimized& c, ComplexPolyBdJetOptimized& r){
    r.m_upperHalfspace = vectalg::midVector(m_upperHalfspace) - m_upperHalfspace;
    c.m_upperHalfspace = vectalg::midVector(m_upperHalfspace);

    r.m_lowerHalfspace = vectalg::midVector(m_lowerHalfspace) - m_lowerHalfspace;
    c.m_lowerHalfspace = vectalg::midVector(m_lowerHalfspace);
  }

  inline void split(RealVectorType& c, RealVectorType& r){
    int i;
    for(i=0; i < d; i++){
      r[2 * i] = mid(m_upperHalfspace[i].value().re) - m_upperHalfspace[i].value().re;
      c[2 * i] = mid(m_upperHalfspace[i].value().re);

      r[2 * i + 1] = mid(m_upperHalfspace[i].value().im) - m_upperHalfspace[i].value().im;
      c[2 * i + 1] = mid(m_upperHalfspace[i].value().im);
    }
    for(i = -d + 1; i < 0; i++){
      r[r.dimension() + 2 * i] = mid(m_lowerHalfspace[-i].value().re) - m_lowerHalfspace[-i].value().re;
      c[c.dimension() + 2 * i] = mid(m_lowerHalfspace[-i].value().re);

      r[r.dimension() + 2 * i + 1] = mid(m_lowerHalfspace[-i].value().im) - m_lowerHalfspace[-i].value().im;
      c[c.dimension() + 2 * i + 1] = mid(m_lowerHalfspace[-i].value().im);
    }
  }

  inline bool subset(ComplexPolyBdJetOptimized& pb2){
    bool r = true;
    int i;
    std::cout << "d=" << d << "\n";
    for(i=0; i < d; i++){
      if(!m_upperHalfspace[i].re.subset(pb2.m_upperHalfspace[i].re)){
        r = false;
        std::cout << m_upperHalfspace[i].re << " " << pb2.m_upperHalfspace[i].re << "\n";
      }
      if(!m_upperHalfspace[i].im.subset(pb2.m_upperHalfspace[i].im)){
        r = false;
        std::cout << m_upperHalfspace[i].im << " " << pb2.m_upperHalfspace[i].im << "\n";
      }
      if(!m_lowerHalfspace[i].re.subset(m_lowerHalfspace[i].re)){
        r = false;
        std::cout << m_upperHalfspace[i].re << " " << pb2.m_upperHalfspace[i].re << "\n";
      }
      if(!m_lowerHalfspace[i].im.subset(m_lowerHalfspace[i].im)){
        r = false;
        std::cout << m_upperHalfspace[i].im << " " << pb2.m_upperHalfspace[i].im << "\n";
      }
    }
    return r;
  }
  
  inline void abs_supremum(ComplexPolyBdJetOptimized& abs){
    int i;
    for(i=0; i < d; i++){
      abs.m_upperHalfspace[i] = m_upperHalfspace[i].abs_supremum();
      abs.m_lowerHalfspace[i] = m_lowerHalfspace[i].abs_supremum();
    }
  }
  
  void clean(){}

  /**Returns maximal diameter of the modes stored within the container.
   */
  inline IntervalType maxDiam(int& coord){
    int i;
    coord = 0;
    IntervalType maxDiam = 0;
    for(i=0; i < dimension(); i++){
      if(diam(m_upperHalfspace[i].re) >= maxDiam){
        coord = i;
        maxDiam = diam(m_upperHalfspace[i].re);
      }
      if(diam(m_upperHalfspace[i].im) >= maxDiam){
        coord = i;
        maxDiam = diam(m_upperHalfspace[i].im);
      }
    }
    return maxDiam;
  }

  /**Returns the average diameter of the modes stored within the container.
   */
  inline IntervalType avgDiam(){
    int i;
    IntervalType avgDiam = 0;
    for(i=0; i < dimension(); i++){
      if(baseReZero || baseImZero){
        if(baseImZero) avgDiam += diam(m_upperHalfspace[i].re);        
        if(baseReZero) avgDiam += diam(m_upperHalfspace[i].im);
      }else{
        avgDiam += diam(m_upperHalfspace[i].re);
        avgDiam += diam(m_upperHalfspace[i].im);
      }
    }
    if(baseReZero || baseImZero)
      return avgDiam / dimension();
    else
      return avgDiam / (2 * dimension());
  }  
  
  friend std::ostream& operator<<(std::ostream& out, const ComplexPolyBdJetOptimized& c){
    out << "ComplexPolyBdJetOptimized\n" << (DPDEContainer&)c << "\n";
    IndexType index;
    for(index = c.firstModeIndex(c.irProjection); !index.limitReached(c.irProjection); index.inc(c.irProjection)){      
      out << index << "," << c.mode2array(index) << " : " << c[index] << "\n";
    }
    out << "conjugate part:\n";
    for(index = c.firstModeIndex(c.irProjection); !index.limitReached(c.irProjection); index.inc(c.irProjection)){      
      out << -index << ": " << c[-index] << "\n";
    }
    return out;
  }

  ///dummy function
  friend const IntervalType C(const ComplexPolyBdJetOptimized& pb){
    return 0;
  }

  ///dummy function
  friend const int s(const ComplexPolyBdJetOptimized& pb){
    return 0;
  }

  ///dummy function
  friend void setS(const ComplexPolyBdJetOptimized& pb, int newS){
  }

  ///dummy function
  friend void setC(const ComplexPolyBdJetOptimized& pb, IntervalType newS){
  }

  using SubspaceType::a0IsConstant;
  using SubspaceType::mode2array;
  using SubspaceType::array2modeIndex;
  using SubspaceType::firstModeIndex;
  using SubspaceType::modes2realArraySizeStatic;
  using SubspaceType::modes2arraySizeStatic;
  using SubspaceType::projectOntoRealSpace;
  using SubspaceType::projectOntoImaginarySpace;
  using SubspaceType::mode;
  using SubspaceType::setMode;
  using SubspaceType::setDerivative;

};


}
}

#endif /* HALFSPACEMODESCONTAINER_H_ */
