/// @addtogroup diffIncl2
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file InclRect2Set.h
///
/// @author Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2007 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.wsb-nlu.edu.pl/ for details. 

#ifndef _CAPD_DIFFINCL2_INCLRECT2SET_H_
#define _CAPD_DIFFINCL2_INCLRECT2SET_H_

//#include "dynset/inclc0set.h"
#include "capd/dynset/C0Rect2Set.h"
#include "capd/dynset/ReorganizedSet.h"
#include "capd/vectalg/Norm.h"

namespace capd{
namespace diffIncl2{

///////////////////////////////////////////////////////////////////////////
// InclRect2Set
/// 
/// Set representation for differential inclusions based on capd::dynset::Rect2Set class
///
///   set is represented as: x + C*r0 + B*r   where
///       C*r0 - basic 'Lipschitz part'
///       B*r - QR-decomposition of the remaining errors
///   
///   it also contains estimates for perturbations (tail modes) Y

template<typename MatrixT, typename PerturbationT>
class InclRect2Set : public capd::dynset::ReorganizedSet< capd::dynset::C0Rect2Set<MatrixT> >{

public:
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef capd::vectalg::Norm<VectorType,MatrixType> NormType;
  typedef capd::dynset::ReorganizedSet< capd::dynset::C0Rect2Set<MatrixT> > BaseSet;

  typedef PerturbationT PerturbationType;

  // constr
  explicit InclRect2Set(int dimension);
  //explicit InclRect2Set(const VectorType& X);
  InclRect2Set(const VectorType& X, const PerturbationType& Y);
  InclRect2Set(PerturbationType& Y);
  InclRect2Set(const BaseSet & X, const PerturbationType& Y);

  capd::dynset::C0Set<MatrixType>* clone(void) const;
  capd::dynset::C0Set<MatrixType>* fresh(void) const;
  capd::dynset::C0Set<MatrixType>* cover(const VectorType& v) const;
  
  template<typename DiffIncl>
  void move( DiffIncl& dynsys);
  template<typename DiffIncl>
  void move( DiffIncl & dynsys, InclRect2Set& result);
  
  std::vector<VectorType> getCorners() const;

  using BaseSet::get_x;
  using BaseSet::get_r;
  using BaseSet::get_r0;
  using BaseSet::get_B;
  using BaseSet::get_C;
  using BaseSet::operator VectorType;
  using BaseSet::show;
  using BaseSet::affineTransformation;

  BaseSet & getBaseSet(){
    return static_cast<BaseSet &>(*this);
  }
  const BaseSet & getBaseSet() const {
      return static_cast<const BaseSet &>(*this);
  }
  void setBaseSet(const BaseSet & newBaseSet){
    BaseSet::operator = (newBaseSet);
  }

  PerturbationType & getPerturbationParams(){
      return m_Y;
  }
  const PerturbationType & getPerturbationParams() const{
      return m_Y;
  }
  void setPerturbationParams(const PerturbationType & newY){
    m_Y = newY;
  }
//
//  PerturbationType & get_Y(){
//    return m_Y;
//  }
//  const PerturbationType & get_Y() const{
//    return m_Y;
//  }
//  void set_Y(const PerturbationType & newY){
//    m_Y = newY;
//  }
protected:
  using BaseSet::m_x;
  using BaseSet::m_r;
  using BaseSet::m_r0;
  using BaseSet::m_B;
  using BaseSet::m_C;
  PerturbationType m_Y;

};

template<typename MatrixType, typename PerturbationType>
std::vector<typename MatrixType::VectorType> getCorners(const InclRect2Set<MatrixType, PerturbationType> & set) ;

// inline definitions
////////////////////////////////////////////////////////////////////////////////
/// Constructors
template<typename MatrixType, typename PerturbationT>
inline InclRect2Set<MatrixType, PerturbationT>::InclRect2Set(int dim)
  :  BaseSet(dim) {
}


template<typename MatrixType, typename PerturbationT>
inline InclRect2Set<MatrixType, PerturbationT>::InclRect2Set(const VectorType& X, const PerturbationType& Y)
  :  BaseSet(X), m_Y(Y) {
  if(!Y.infiniteDimensional){
    std::cerr << "An infinite dimensional PolyBd is expected as the PerturbationType (second parameter)\n";
    throw std::runtime_error("An infinite dimensional PolyBd is expected as the PerturbationType (second parameter)\n");
  }
}

template<typename MatrixType, typename PerturbationT>
inline InclRect2Set<MatrixType, PerturbationT>::InclRect2Set(PerturbationType& Y)
  : BaseSet(Y), m_Y(Y) {
  if(!Y.infiniteDimensional){
    std::cerr << "An infinite dimensional PolyBd is expected as the PerturbationType (first parameter)\n";
    throw std::runtime_error("An infinite dimensional PolyBd is expected as the PerturbationType (first parameter)\n");  
  }
}
template<typename MatrixType, typename PerturbationT>
inline InclRect2Set<MatrixType, PerturbationT>::InclRect2Set(const BaseSet & X, const PerturbationType& Y)
  :  BaseSet(X), m_Y(Y){
  if(!Y.infiniteDimensional){
    std::cerr << "An infinite dimensional PolyBd is expected as the PerturbationType (second parameter)\n";
    throw std::runtime_error("An infinite dimensional PolyBd is expected as the PerturbationType (second parameter)\n");
  }
}

//
//template<typename MatrixType, typename PerturbationT>
//inline InclRect2Set<MatrixType, PerturbationT>::InclRect2Set(const VectorType& the_x,const VectorType& the_r0, const PerturbationType& Y)
//  :  BaseSet(the_x, the_r0), m_Y(Y){
//}
//
//template<typename MatrixType, typename PerturbationT>
//inline InclRect2Set<MatrixType, PerturbationT>::InclRect2Set(
//      const VectorType& the_x,
//      const MatrixType& the_C,
//      const VectorType& the_r0,
//      const PerturbationType& Y
//   )
//  : BaseSet(the_x, the_C, the_r0), m_Y(Y){
//}
//
//template<typename MatrixType, typename PerturbationT>
//inline InclRect2Set<MatrixType, PerturbationT>::InclRect2Set(
//      const VectorType& the_x,
//      const MatrixType &the_C,
//      const VectorType& the_r0,
//      const VectorType& the_r,
//      const PerturbationType& Y
//   ): BaseSet(the_x, the_C, the_r0, the_r), m_Y(Y){
//}
//
//template<typename MatrixType, typename PerturbationT>
//inline InclRect2Set<MatrixType, PerturbationT>::InclRect2Set(
//      const VectorType& the_x,
//      const MatrixType& the_C,
//      const VectorType& the_r0,
//      const MatrixType& the_B,
//      const VectorType& the_r,
//      const PerturbationType& Y
//   ): BaseSet(the_x, the_C, the_r0, the_B, the_r), m_Y(Y){
//}

////////////////////////////////////////////////////////////////////////////////
/// C0set interface overriding
   
template<typename MatrixType, typename PerturbationT>
inline typename capd::dynset::C0Set<MatrixType>* InclRect2Set<MatrixType, PerturbationT>::clone() const{
   return new InclRect2Set(*this);
}

template<typename MatrixType, typename PerturbationT>
inline typename capd::dynset::C0Set<MatrixType>* InclRect2Set<MatrixType, PerturbationT>::fresh() const{
   return new InclRect2Set<MatrixType, PerturbationT>(m_x.dimension());
}

template<typename MatrixType, typename PerturbationT>
inline typename capd::dynset::C0Set<MatrixType>* InclRect2Set<MatrixType, PerturbationT>::cover(const VectorType& v) const{
   return new InclRect2Set(midVector(v),v-midVector(v));
}
////////////////////////////////////////////////////////////////////////////////
template<typename MatrixType, typename PerturbationT>
std::vector<typename MatrixType::VectorType> getCorners(const capd::diffIncl2::InclRect2Set<MatrixType, PerturbationT> & set) ;

}} // namespace capd::diffIncl2

#endif // _CAPD_DIFFINCL2_INCLRECT2SET_H_

/// @}
