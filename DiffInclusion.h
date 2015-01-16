/// @addtogroup diffIncl2
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file DiffInclusion.h
///
/// @author Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.wsb-nlu.edu.pl/ for details.

/* Author: Tomasz Kapela, 2007 */

#ifndef _CAPD_DIFFINCL2_DIFFINCLUSION_H_
#define _CAPD_DIFFINCL2_DIFFINCLUSION_H_

#include <string>
#include <map>
#include "capd/vectalg/Norm.h"
#include "capd/dynsys/Taylor.h"
#include "capd/dynsys/enclosure.h"
#include "capd/dynsys/StepControl.h"
namespace capd{
/// A rigorous integration of the differential inclusions
namespace diffIncl2{

/**
 * Base class for rigorous integration of differential inclusions.
 *
 * Template arguments:
 * - MapT     - MultiMap that stores RHS of the differential inclusion in the form : selection + 'error bounds'
 *              (we assume that it implements all methods that class capd::diffIncl2::MultiMap has).
 * - DynSysT  - numerical method for ODE integration
 *
 * \see capd::diffIncl2::DiffInclusionLN, \see capd::diffIncl2::DiffInclusionCW
 */
template<typename MapT, typename DynSysT = capd::dynsys::Taylor< typename MapT::MapType> >
class DiffInclusion  : public capd::dynsys::StepControlInterface<capd::dynsys::NoStepControl> {

public:
  typedef MapT MultiMapType;
  typedef MultiMapType MapType;
//  typedef typename MapT::MapType VectorFieldType;
//  typedef typename MapT::FunctionType FunctionType;
  typedef typename MapType::MatrixType MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MapType::ParamType ParamType;

  typedef capd::vectalg::Norm<VectorType, MatrixType> NormType;
  typedef DynSysT DynSysType;

  DiffInclusion(   MultiMapType& diffIncl,               // map of the form f(x)+ g(x,e)
      int order,                            // order for integration
      ScalarType const & step,              // time step for integration
      NormType const & norm                 // norm
  );
  
  DiffInclusion(   
      int m,
      int dftPts,
      int M, 
      int dftPts2,
      const ScalarType & pi,
      const ScalarType & nu,
      int order,                            // order for integration
      ScalarType const & step,              // time step for integration
      NormType const & norm                 // norm
  );
  
  ~DiffInclusion();

  /// eclosure of solution of diff. inclusion during one time step starting at x
  virtual ParamType diffInclusionEnclosure(const ParamType & y) const=0;

  /// eclosure of solution of selected ODE during one time step starting at x
  virtual VectorType dynamicalSystemEnclosure(const VectorType& x) ;

  /// eclosure of solution of diff. inclusion during one time step starting at x
  template <typename SetType>
  VectorType enclosure(const SetType & x);

  /// Bounds for perturbation of solution of selected ODE after one time step
  virtual VectorType perturbations(const VectorType & x, ParamType & W_y) = 0;

  virtual void moveParams(const VectorType& x, ParamType & params){}

  template <typename SetType>
  void computeImageAndPerturbation(
      SetType & initialSet,
      SetType & image,
      VectorType & perturbation
  ){
//     clock_t start = clock();
     perturbation = this->perturbations(initialSet.getBaseSet(), initialSet.getPerturbationParams());
//     clock_t end = clock();
//     std::cout << "perturbation time: " << end-start << "\n";
//     start = clock();
     initialSet.getBaseSet().move(getDynamicalSystem(), image.getBaseSet());
//     end = clock();
//    std::cout << "move time: " << end-start << "\n";
     moveParams(image.getBaseSet(), image.getPerturbationParams());
  }

  /// returns RHS of a diff. inclusion
  const MapType& getField() const;
  /// returns RHS of a diff. inclusion
  MapType& getField();

  /// returns order of numerical method
  int getOrder() const;
  /// sets order of numerical method
  void setOrder(int newOrder);

  /// returns current time step
  const ScalarType& getStep() const;
  /// sets currect time step
  void setStep(const ScalarType& newStep);

  /// returns dynamical system (numerical ODE integrator)
  DynSysType & getDynamicalSystem();
  /// returns dynamical system (numerical ODE integrator)
  const DynSysType & getDynamicalSystem() const;

  void clearCoefficients(){
     m_dynamicalSystem.clearCoefficients();
   }

   template <class SetType>
   ScalarType computeNextTimeStep(const SetType& x, const ScalarType& maxStep) {
     return this->m_stepControl.computeNextTimeStep(*this,x,maxStep);
   }

   template <class SetType>
   ScalarType getFirstTimeStep(const SetType& x, const ScalarType& maxStep) {
     return this->m_stepControl.getFirstTimeStep(*this, x, maxStep);
   }


protected:
  //  void operator=(const DiffInclusion & a_t){}
  // DiffInclusion(const DiffInclusion & diffIncl2) : m_norm(m_norm.clone()), m_diffIncl(diffIncl2.m_diffIncl){}

  NormType * m_norm;                  ///<  norm used in perturbation estimations
  DynSysType m_dynamicalSystem;       ///<  dynamical system of selected ODE (numerical integrator)
  MultiMapType & m_diffIncl;          ///<  RHS of differential inclusion

};


// --------------- inline definitions -----------------


template <typename MapT, typename DynSysT>
inline const MapT& DiffInclusion<MapT, DynSysT>::getField() const {
  return m_diffIncl;
}

template <typename MapT, typename DynSysT>
inline MapT & DiffInclusion<MapT, DynSysT>::getField() {
  return m_diffIncl;
}


template <typename MapT, typename DynSysT>
inline int DiffInclusion<MapT, DynSysT>::getOrder() const {
  return m_dynamicalSystem.getOrder();
}

template <typename MapT, typename DynSysT>
inline void DiffInclusion<MapT, DynSysT>::setOrder(int newOrder) {
  m_dynamicalSystem.setOrder(newOrder);
}


template <typename MapT, typename DynSysT>
inline const typename DiffInclusion<MapT, DynSysT>::ScalarType& DiffInclusion<MapT, DynSysT>::getStep() const {
  return m_dynamicalSystem.getStep();
}

template <typename MapT, typename DynSysT>
inline void DiffInclusion<MapT, DynSysT>::setStep(const ScalarType& newStep) {
  m_dynamicalSystem.setStep(newStep);
}

///**
// * Computes enclosure of image of given set for differential inclusion during whole time step
// */
//template <typename MapT, typename DynSysT>
//inline typename DiffInclusion<MapT, DynSysT>::VectorType DiffInclusion<MapT, DynSysT>::diffInclusionEnclosure(const VectorType& x, const ParamType & y) const{
//
//  m_diffIncl.setParameters(y);
//  return capd::dynsys::enclosure(m_diffIncl, x, getStep());
//}

template <typename MapT, typename DynSysT>
inline typename DiffInclusion<MapT, DynSysT>::VectorType DiffInclusion<MapT, DynSysT>::dynamicalSystemEnclosure(const VectorType & x)  {
  return m_dynamicalSystem.enclosure(x);
}

template <typename MapT, typename DynSysT>
template <typename SetType>
inline typename DiffInclusion<MapT, DynSysT>::VectorType DiffInclusion<MapT, DynSysT>::enclosure(const SetType & set) {

  return diffInclusionEnclosure(VectorType(set.getBaseSet(), set.getPerturbationsParams()));
}

template <typename MapT, typename DynSysT>
inline typename DiffInclusion<MapT, DynSysT>::DynSysType & DiffInclusion<MapT, DynSysT>::getDynamicalSystem(){
  return m_dynamicalSystem;
}

template <typename MapT, typename DynSysT>
inline const typename DiffInclusion<MapT, DynSysT>::DynSysType&  DiffInclusion<MapT, DynSysT>::getDynamicalSystem() const{
  return m_dynamicalSystem;
}


//###########################################################//


}} // namespace capd::diffIncl2

#endif // _CAPD_DIFFINCL2_DIFFINCLUSION_H_

/// @}
