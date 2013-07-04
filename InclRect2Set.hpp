/// @addtogroup diffIncl2
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file InclRect2Set.hpp
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.wsb-nlu.edu.pl/ for details. 

#ifndef _CAPD_DIFFINCL2_INCLRECT2SET_HPP_
#define _CAPD_DIFFINCL2_INCLRECT2SET_HPP_

#include "InclRect2Set.h"
#include "capd/dynset/C0Rect2Set.hpp"
#include "capd/vectalg/Norm.hpp"
#include <iostream>

namespace capd{
  namespace diffIncl2{
    
    
    template<typename MatrixType, typename PerturbationT>
    template<typename DiffIncl>
    void InclRect2Set<MatrixType, PerturbationT>::move(DiffIncl & diffIncl2) {
      this->move(diffIncl2, *this);
    }

    template<typename MatrixType, typename PerturbationT>
    template<typename DiffIncl>
    void InclRect2Set<MatrixType, PerturbationT>::move(DiffIncl & diffIncl2, InclRect2Set<MatrixType, PerturbationT>& result) {

      VectorType Delta;
      diffIncl2.computeImageAndPerturbation(*this, result, Delta);
      
      // Rearrangements
      VectorType x = midVector( result.m_x + Delta);
      VectorType dr = result.m_x + Delta - x;
      MatrixType BT = Transpose(result.m_B);
      
      result.m_r = result.m_r + BT * dr;
      result.m_x = x;
    }
    
    template<typename T>
    void corners(T& head, const T & tail, int i, int dim, std::vector<T> & cor) {
      if(i < dim) {
        head[i] = left(tail[i]);
        corners(head, tail, i+1, dim, cor);
        head[i] = right(tail[i]);
        corners(head, tail, i+1, dim, cor);
        head[i] = tail[i];
      }
      else {
        cor.push_back(head);
      }
    }
    
    
    template<typename SetType>
    std::vector<typename SetType::VectorType> getCorners(const SetType & set) {
      typedef typename SetType::VectorType VectorType;
      std::vector<VectorType> cor;
      VectorType v = set.get_r();
      corners(v, set.get_r(), 0, v.dimension(), cor);
      for(typename std::vector<VectorType>::iterator it = cor.begin(); it != cor.end(); ++it){
        *it = set.get_x() + set.get_C() * set.get_r0() + set.get_B() * *it;
      }
      return cor;
    }
    
    template<typename MatrixType, typename PerturbationType>
    std::vector<typename InclRect2Set<MatrixType, PerturbationType>::VectorType> InclRect2Set<MatrixType, PerturbationType>::getCorners() const {
      return ::capd::diffIncl2::getCorners(*this);
    }
    
  }} // namespace capd::diffIncl2

#endif // _CAPD_DIFFINCL2_INCLRECT2SET_HPP_

/// @}
