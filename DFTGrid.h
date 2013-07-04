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
 * DFTGrid.h
 *
 *  Created on: Oct 26, 2011
 *      Author: cyranka
 */

#ifndef DFTGRID_H_
#define DFTGRID_H_

#include "capd/vectalg/Vector.h"
#include "capd/vectalg/Matrix.h"
#include "DPDEContainer.h"

namespace capd{
namespace jaco{

/**This class represents a one dimensional grid of L_2 coefficients (function values calculated at uniform discrete points).
 */
template<class IntervalT, class ScalarT, int M>
class DFT1DGrid : public capd::vectalg::Vector<ScalarT, M>, public capd::jaco::DPDEContainer{
public:
  typedef ScalarT ScalarType;
  typedef IntervalT IntervalType;
  typedef capd::vectalg::Vector<ScalarType, M> GridType;
  typedef GridType VectorType;
  typedef capd::vectalg::Matrix<IntervalType, 2, 2> QuadMatrixType;
  typedef typename VectorType::iterator iterator;
  typedef typename VectorType::const_iterator const_iterator;
  int m;

  DFT1DGrid() : m(0){ }

  DFT1DGrid(int m_) : VectorType(m_), m(m_){ }

  DFT1DGrid(int m_, bool f) : VectorType(m_, f), m(m_){ }

  DFT1DGrid(int m_, double d) : VectorType(m_, false), m(m_){ }

  DFT1DGrid(const VectorType& v) : VectorType(v), m(v.dimension()){}


  inline DFT1DGrid& operator=(const DFT1DGrid& g){
    (DPDEContainer&)*this = (DPDEContainer&)g;
    m = g.m;
    (VectorType&)*this = g;
    return *this;
  }

  inline DFT1DGrid& operator+=(const DFT1DGrid& g){
    (DPDEContainer&)*this += (DPDEContainer&)g;
    (VectorType&)*this += g;
    return *this;
  }
  
  inline DFT1DGrid& operator+=(const ScalarType& s){
    int j;
    for(j=0; j < m; ++j){
      VectorType::operator[](j) += s;
    }
    return *this;
  }

  inline DFT1DGrid& operator*=(const ScalarType& s){    
    if(s.isImUnit())
      multiplyByImUnit();
    (VectorType&)*this *= s;
    return *this;
  }

  inline DFT1DGrid& operator*=(const IntervalType& i){    
    int j;
    for(j=0; j < m; ++j){
      VectorType::operator[](j) = i * VectorType::operator[](j);
    }
    return *this;
  }

  inline DFT1DGrid& multiplyAndDivide(const DFT1DGrid& dft1, const DFT1DGrid& dft2){
    DPDEContainer::multiply(dft1, dft2);
    int j;
    IntervalType d = IntervalType(1) / IntervalType(m);
    for(j=0; j < m; ++j){
      VectorType::operator[](j) = d * dft1[j];
      VectorType::operator[](j) = VectorType::operator[](j) * dft2[j];
    }
    return *this;
  }

  inline DFT1DGrid& multiply(const DFT1DGrid& dft1, const DFT1DGrid& dft2){
    multiplyAndDivide(dft1, dft2);
    return *this;
  }

  inline DFT1DGrid& multiplyOnly(const DFT1DGrid& dft1, const DFT1DGrid& dft2){
    DPDEContainer::multiply(dft1, dft2);
    int j;
    for(j=0; j < m; ++j){
      VectorType::operator[](j) = dft1[j] * dft2[j];
    }
    return *this;
  }

  inline DFT1DGrid& multiply(const DFT1DGrid& dft1, const DFT1DGrid& dft2, const DFT1DGrid& dft3){
    (*this).multiply(dft1, dft2);
    (*this).multiplyOnly(*this, dft3);
    return *this;
  }

  inline void projectOntoSubspace(){
    int j;
    if(solutionType == capd::jaco::realValued){
      for(j=0; j < m; ++j){
        VectorType::operator[](j).setImaginaryPartToZero();
      }
    }
  }

  inline void projectOntoRealSpace(){
    int j;
    for(j=0; j < m; ++j){
      VectorType::operator[](j).setImaginaryPartToZero();
    }
  }

  inline void projectOntoImaginarySpace(){
    int j;
    for(j=0; j < m; ++j){
      VectorType::operator[](j).setRealPartToZero();
    }
  }

  inline DFT1DGrid& experimentalMultiply(const DFT1DGrid& dft1, const DFT1DGrid& dft2){
    int j;
    IntervalType d = IntervalType(1) / IntervalType(m);
    for(j=0; j < m; ++j){
      *this[j].re = d * (dft1[j].re * dft2[j].re - dft1[j].im * dft2[j].im);
      *this[j].im = d * (dft1[j].re * dft2[j].im + dft1[j].im * dft2[j].re);
    }
    return *this;
  }

  inline int dimension() const{
    return m;
  }

  inline void setImaginaryPartToZero(){
    int i;
    for(i=0; i < m; ++i){
      VectorType::operator[](i).setImaginaryPartToZero();
    }
  }

  friend void print(std::ostream& out, const DFT1DGrid& c, int top){
    int i;
    for(i=0; i < c.m; ++i){
      out << i*2*3.14159265358979323846264338/double(c.m) << " " << c[i].secondFreeCoeff().re.leftBound() << "\n";
    }
    for(i=0; i < c.m; ++i){
      out << i*2*3.14159265358979323846264338/double(c.m) << " " << c[i].secondFreeCoeff().re.rightBound() << "\n";
    }
    out << "\n";
    return out;
  }
  
  friend std::ostream& operator<<(std::ostream& out, const DFT1DGrid& c){
    int i;
//    for(i=0; i < c.m; ++i){
//      out << c[i];
//      if(i < c.m-1)
//        out << "\n";
//    }
    ///Gnuplot verseion, draws a ``wide strip'' representing the function values therefore fist the lower bound is printed, and then 
    ///the higher bound.
    for(i=0; i < c.m; ++i){
      out << i*2*3.14159265358979323846264338/double(c.m) << " " << c[i].secondFreeCoeff().re.leftBound() << "\n";
    }
    for(i=c.m-1; i >= 0; --i){
      out << i*2*3.14159265358979323846264338/double(c.m) << " " << c[i].secondFreeCoeff().re.rightBound() << "\n";
    }
    out << "\n";
    return out;
  }

};

template<class IntervalT, class ScalarT, int M>
inline const DFT1DGrid<IntervalT, ScalarT, M> operator*(const DFT1DGrid<IntervalT, ScalarT, M>& dft1, const DFT1DGrid<IntervalT, ScalarT, M>& dft2){
  DFT1DGrid<IntervalT, ScalarT, M> r(dft1.m);
  ScalarT d = ScalarT(1) / ScalarT(dft1.m);
  int j;
  for(j=0; j < dft1.m; ++j){
    r[j] = d * dft1[j];
    r[j] = r[j]* dft2[j];
  }

  return r;
}

template<class IntervalT, class ScalarT, int M>
inline const DFT1DGrid<IntervalT, ScalarT, M> operator+(const DFT1DGrid<IntervalT, ScalarT, M>& dft1, const DFT1DGrid<IntervalT, ScalarT, M>& dft2){
  DFT1DGrid<IntervalT, ScalarT, M> r(dft1.m);
  int j;
  for(j=0; j < dft1.m; ++j){
    r[j] = dft1[j] + dft2[j];
  }

  return r;
}

/** 05.12.2011 this is used during the initialization of the FFT2D (in ModesContainterOfGrids this is called)
 *
 */
template<class IntervalT, class ScalarT, int M>
inline const DFT1DGrid<IntervalT, ScalarT, M> conjugate(const DFT1DGrid<IntervalT, ScalarT, M>& dft){
  DFT1DGrid<IntervalT, ScalarT, M> r(dft.m);
  int j;
  for(j=0; j < M; ++j){
    r[j] = conjugate(dft[j]);
  }
  return r;
}


template< class IntervalT, class ScalarT, int M>
class DFT2DGrid{
public:
  typedef ScalarT ScalarType;
  typedef IntervalT IntervalType;
  typedef capd::jaco::DFT1DGrid<IntervalType, ScalarType, M> DFT1DGridType;
  typedef capd::vectalg::Vector<DFT1DGridType, M> GridType;

  int m;
  GridType grid;

  DFT2DGrid() : m(0){ }

  DFT2DGrid(int m_) : m(m_), grid(m, false){
    int i;
    for(i=0; i < m; ++i)
      grid[i] = DFT1DGridType(m);
  }

  inline const DFT1DGridType& operator[](int j) const{
    return grid[j];
  }

  inline DFT1DGridType& operator[](int j){
    return grid[j];
  }

  inline DFT2DGrid& operator*=(ScalarType s){
    int i;
    for(i=0; i < m; ++i){
      grid[i] *= s;
    }
    return *this;
  }

  inline DFT2DGrid& multiply(const DFT2DGrid& dft1, const DFT2DGrid& dft2){
    IntervalType d = IntervalType(1) / IntervalType(m*m);
    int i,j;
    for(i=0; i < m; ++i)
      for(j=0; j < m; ++j){
        grid[i][j] = d * dft1[i][j] * dft2[i][j];
      }
    return *this;
  }

  friend std::ostream& operator<<(std::ostream& out, const DFT2DGrid& c){
    int i, j;
    for(i=0; i < c.m; ++i){
      for(j=0; j < c.m; ++j)
        out << i*2*3.14159265358979323846264338/double(c.m) << " " << j*2*3.14159265358979323846264338/double(c.m)
        << " " << c[i][j] << "\n";//leftBound(c[i][j].re) << " " << rightBound(c[i][j].re) << "\n";
    }
    return out;
  }

};

/**A vector of grids.
 */
template<class DFTGridT, int D>
class ComponentGrid{
public:
  typedef DFTGridT DFTGridType;
  typedef capd::vectalg::Vector<DFTGridType, 0> GridVectorType;

  GridVectorType v;

  inline ComponentGrid() { }

  inline ComponentGrid(int m) : v(D){
    int i;
    for(i=0; i < D; ++i)
      v[i] = DFTGridType(m);
  }

  inline const DFTGridType& operator[](int i) const{
    return v[i];
  }

  inline DFTGridType& operator[](int i){
    return v[i];
  }

  friend std::ostream& operator<<(std::ostream& out, const ComponentGrid& c){
    int i;
    for(i=0; i < D; ++i)
      out << "component #" << i << ":\n" << c[i] << "\n";
    return out;
  }
};

template< class IntervalT, class ScalarT, int M>
inline const DFT2DGrid<IntervalT, ScalarT, M> operator*(const DFT2DGrid<IntervalT, ScalarT, M>& dft1, const DFT2DGrid<IntervalT, ScalarT, M>& dft2){
  DFT2DGrid<IntervalT, ScalarT, M> r(dft1.m);
  ScalarT d = ScalarT(1)/ScalarT(dft1.m*dft1.m);
  int i, j;
  for(i=0; i < dft1.m; ++i){
    for(j=0; j < dft1.m; ++j){
      r[i][j] = d * dft1[i][j] * dft2[i][j];
    }
  }

  return r;
}

template< class IntervalT, class ScalarT, int M>
inline const DFT2DGrid<IntervalT, ScalarT, M> operator+(const DFT2DGrid<IntervalT, ScalarT, M>& dft1, const DFT2DGrid<IntervalT, ScalarT, M>& dft2){
  DFT2DGrid<IntervalT, ScalarT, M> r(dft1.m);
  int i, j;
  for(i=0; i < dft1.m; ++i){
    for(j=0; j < dft1.m; ++j){
      r[i][j] = dft1[i][j] + dft2[i][j];
    }
  }

  return r;
}


}
}

#endif /* DFTGRID_H_ */
