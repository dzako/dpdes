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

#ifndef _CAPD_JACO_JACOINCLRECT2SET_H
#define _CAPD_JACO_JACOINCLRECT2SET_H

#include "capd/diffIncl/InclRect2Set.hpp"

namespace capd{
namespace jaco{


///Class describing a set in the context of dissipative PDEs. It is composed of a BaseSet that
///represents the finite part, and the tail, instance of TailT that represents the tail.
template<typename MatrixT, typename TailT>
class DPDEInclRect2Set : public capd::diffIncl::InclRect2Set<MatrixT>{
  
public: 
  typedef MatrixT MatrixType;
  typedef TailT TailType;
  typedef capd::diffIncl::InclRect2Set<MatrixType> BaseSet;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef capd::vectalg::Norm<VectorType,MatrixType> NormType;

  explicit DPDEInclRect2Set(int dimension);
  explicit DPDEInclRect2Set(const VectorType& the_x);
  DPDEInclRect2Set(const VectorType& x, const TailType& t);///<version not used, left for compability reasons
//  DPDEInclRect2Set(const VectorType& the_x, const VectorType& the_r0);
  DPDEInclRect2Set(const VectorType& the_x, const MatrixType& the_C, const VectorType& the_r0);
//  DPDEInclRect2Set(const VectorType& the_x, const MatrixType& the_C,
//               const VectorType& the_r0,
//               const VectorType& the_r
//  );

 
  template<typename DiffIncl>
  void move( DiffIncl& dynsys, int stepNumber, std::ostream& f);

  template<typename DiffIncl>
  void move( DiffIncl & dynsys, DPDEInclRect2Set& result) const;

  using BaseSet::get_x;
  using BaseSet::get_r;
  using BaseSet::get_r0;
  using BaseSet::get_B;
  using BaseSet::get_C;
  using BaseSet::operator VectorType;
  using BaseSet::show;
  using BaseSet::affineTransformation;

  TailT m_tail;

protected:
  using BaseSet::m_x;
  using BaseSet::m_r;
  using BaseSet::m_r0;
  using BaseSet::m_B;
  using BaseSet::m_C;
  
};

///==========================================function definitions====================================================

template<typename MatrixType, typename TailT>
inline DPDEInclRect2Set<MatrixType, TailT>::DPDEInclRect2Set(int dim)
  :  BaseSet(dim) {
}

template<typename MatrixType, typename TailT>
inline DPDEInclRect2Set<MatrixType, TailT>::DPDEInclRect2Set(const VectorType& the_x)
  :  BaseSet(the_x) {
}

template<typename MatrixType, typename TailT>
inline DPDEInclRect2Set<MatrixType, TailT>::DPDEInclRect2Set(const VectorType& the_x, const TailType& t)
  :  BaseSet(the_x) {
}

//template<typename MatrixType, typename TailT>
//inline DPDEInclRect2Set<MatrixType, TailT>::DPDEInclRect2Set(const VectorType& the_x,const VectorType& the_r0)
//  :  BaseSet(the_x, the_r0) {
//}

template<typename MatrixType, typename TailType>
inline DPDEInclRect2Set<MatrixType, TailType>::DPDEInclRect2Set(
      const VectorType& the_x,
      const MatrixType& the_C,
      const VectorType& the_r0
   )
  : BaseSet(the_x, the_C, the_r0){
}

//template<typename MatrixType>
//inline DPDEInclRect2Set<MatrixType>::DPDEInclRect2Set(
//      const VectorType& the_x,
//      const MatrixType &the_C,
//      const VectorType& the_r0,
//      const VectorType& the_r
//   ): BaseSet(the_x, the_C, the_r0, the_r){
//}

template<typename MatrixType, typename TailT>
template<typename DiffIncl>
void DPDEInclRect2Set<MatrixType, TailT>::move(DiffIncl & diffIncl, int stepNumber, std::ostream& f) {
  ///opening a file for debugging purposes

  typedef typename DiffIncl::TailType TailType;
  typedef typename DiffIncl::MultiMapType MultiMapType;

  VectorType x = VectorType(*this);
  VectorType x0=x;

  ///variables used for storing
  TailType T=TailType(); /// T: T([0, h]) \subset T
  TailType N=TailType(); /// N: N(x+T) \subset N
  VectorType W_2; ///[W_2]
  //f<<"Step number "<<stepNumber<<" in progress.\n";
  VectorType y_c, Deltha;

  {///steps 1-3, 5-8 of Algorithm 1
  Deltha = diffIncl.perturbations(x, T, N, W_2, y_c, f);
  }///end steps 1-3, 5-8 of Algorithm 1

  ///in T we have T([0,h])
  ///in N we have N([W_2], T[0, h])


  {///step  4 of Algorithm 1
  diffIncl.getDynamicalSystem().setYc(y_c); ///setting y_c of a unperturbed projection
  BaseSet::BaseSet::move(diffIncl .getDynamicalSystem()); ///computation of an unperturbed trajectory
  diffIncl.getDynamicalSystem().eraseYc(); ///y_c=0, in order to make the algorithm produce good enclosures in the
                                           ///next step
  }///end step 4 of Algorithm 1

  {///step 10 of Algorithm 1. Rearrangements
  x = midVector( m_x + Deltha );
  VectorType dr = m_x + Deltha - x;
  MatrixType BT = Transpose(m_B);
  m_r = m_r + BT * dr;
  m_x = x;
  }///end step 10 of Algorithm 1

  TailType Th;
  {///step 11 of Algorithm 1
  Th=diffIncl.calculateTh(x, T, N, W_2, f); ///calculates T(h)
  }///end step 11 of Algorithm 1

  #if __GENERAL_DEBUG__
    f<<"current x: \n"<<m_x<<"\n";
    f<<"current r: \n"<<m_r<<"\n";
    f<<"current r0: \n"<<m_r0<<"\n";
    f<<"calculated T(h):\n";
    Th.print(f);
  #endif
  ///we have T(h) is T(0) for the next step
  diffIncl.setT0(Th);
}

template<typename MatrixType, typename TailT>
template<typename DiffIncl>
void DPDEInclRect2Set<MatrixType, TailT>::move(DiffIncl & diffIncl, DPDEInclRect2Set<MatrixType, TailT>& result) const {
  std::ofstream f;
#if __GENERAL_DEBUG__
  f.open(__FILE_NAME_DEBUG__, std::ofstream::app);
#else
   f=std::cout;
#endif
   f<<"Wrong move function call in DPDEInclRect2Set class. Use move(diffIncl, stepNumber) instead.\n";
   std::cerr<<"Wrong move function call in DPDEInclRect2Set class. Use move(diffIncl, stepNumber) instead.\n";
   throw std::runtime_error("Wrong move function call in DPDEInclRect2Set class. Use move(diffIncl, stepNumber) instead.\n");
#if __GENERAL_DEBUG__
  f.close();
#endif
}

}}


#endif
