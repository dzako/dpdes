/// @addtogroup diffIncl2
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file DiffInclusion.hpp
///
/// @author Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.wsb-nlu.edu.pl/ for details. 

/* Author: Tomasz Kapela, 2007 */

#ifndef _CAPD_DIFFINCL2_DIFFINCLUSION_HPP_
#define _CAPD_DIFFINCL2_DIFFINCLUSION_HPP_

#include <sstream>
#include <string>
#include <stdexcept>

#include "capd/dynsys/Taylor.hpp"
#include "DiffInclusion.h"

namespace capd{
namespace diffIncl2{


template <typename MapT, typename DynSysT>
DiffInclusion<MapT, DynSysT>::DiffInclusion(
           MultiMapType& diffIncl2,
           int order, 
           const ScalarType& step, 
           const NormType & norm
) 
  : m_norm ( norm.clone()),
    m_dynamicalSystem(diffIncl2.getVectorField(), order, step),
    m_diffIncl(diffIncl2){
             std::cout << m_dynamicalSystem.m << ", " << m_dynamicalSystem.dftPts1 << "\n";
}
           
 template <typename MapT, typename DynSysT>
 DiffInclusion<MapT, DynSysT>::DiffInclusion(
            int m,
            int dftPts,
            int M,
            int dftPts2,
            const ScalarType & pi,
            const ScalarType & nu,
            int order, 
            const ScalarType& step, 
            const NormType & norm
 ) 
   : m_norm ( norm.clone()),
     m_dynamicalSystem(m, dftPts, M, dftPts2, step, order, pi, nu),
     m_diffIncl((MapT&)m_dynamicalSystem.getVectorField()){
 }           
                     

template <typename MapT, typename DynSysT>
DiffInclusion<MapT, DynSysT>::~DiffInclusion(){
  delete m_norm;
}
  

  
}} //namespace capd::diffIncl2`

#endif // _CAPD_DIFFINCL2_DIFFINCLUSION_HPP_

/// @}
