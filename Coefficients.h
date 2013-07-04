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
 * Derivatives.h
 *
 *  Created on: Sep 15, 2011
 *      Author: cyranka
 * 
 * 
 */

#ifndef COEFFICIENTS_H_
#define COEFFICIENTS_H_

#include "capd/vectalg/Vector.h"

namespace capd{
namespace jaco{

/**Class used for calculating coefficients of the polynomial upto the given order.
  *Important information: dedicated for use in ODEs originating from a PDE discretization, thus
  *automatic diffrentation formulas are implemented for multiplication and addition only.
  *Coefficients of order "order" are the last calculated.
**/
template< class ScalarT, int order>
class Coefficients{
public:
  typedef capd::vectalg::Vector<ScalarT, 0> VectorType;
  typedef ScalarT ScalarType;
  typedef typename ScalarType::ScalarType InternalType;
  VectorType d;
  const static int dim=order+1;
  static Coefficients* imaginaryUnit;

  static Coefficients* buffer; ///<pointer to the buffer storing temporary results

  static Coefficients* zero;

  static int currentDegree;///<stores the index of the current degree which is being computed. Degrees are computed recursively

  Coefficients() : d(dim) {
  }

  Coefficients(ScalarType s) : d(dim){
    d[0] = s;
  }

//  Coefficients(InternalType s) : d(dim){
//    d[0] = s;
//  }

  Coefficients(const Coefficients& c2) : d(dim){
    int i;
    for(i=0; i < dim; i++)
      d[i] = c2[i];
  }

  inline static Coefficients& i(){
    return *imaginaryUnit;
  }

  ///Constructs a new Object representing the imaginary unit i.
  inline static Coefficients constructI(){
    return Coefficients(ScalarType::i());
  }

  inline ScalarType x() const{
    return d[0];
  }

///Automatic diffrentation, assigning right-hand side to the left-hand side, which is diffrentiated, therefore the values are
///assigned with switching.
  inline Coefficients& setCoefficients(const Coefficients& c2, int i){
    d[i+1] = (InternalType(1.) / InternalType(i+1.)) * c2.d[i];
    return *this;
  }

  inline Coefficients& operator=(const Coefficients& c2){
    int i;
    for(i=0; i < dim; ++i){
      d[i] = c2.d[i];
    }
    return *this;
  }

  inline Coefficients& operator=(ScalarType s){
    d[0] = s;
    int i;
    for(i=1; i < dim; ++i){
      d[i] = 0;
    }
    return *this;
  }

//  inline Coefficients& operator=(InternalType v){
//    d[0] = v;
//    int i;
//    for(i=1; i <= order; ++i){
//      d[i] = 0;
//    }
//    return *this;
//  }

  inline Coefficients& operator*=(const Coefficients& c){
    int i, j;
    ScalarType r;
    for(i=currentDegree; i >= 0; --i){///the loop has to be 'downto' type, because we cannot modify d[i] during the process
      r = 0;
      for(j=0; j <= i; ++j){
        r += d[j] * c.d[i-j];
      }
      d[i] = r;
    }
    return *this;
  }

  inline const ScalarType& operator[](int i) const{
    return d[i];
  }

  inline ScalarType& operator[](int i){
    return d[i];
  }

  inline Coefficients& operator+=(const Coefficients& c2){
    int i;
    for(i=0; i <= Coefficients::currentDegree; ++i)
      d[i] += c2[i];
    return *this;
  }

  friend std::ostream& operator<<(std::ostream& out, const Coefficients& c){
    int i;
    for(i=0; i < dim; ++i){
      out << c.d[i];
      if(i < dim-1)
        out << " ";
    }
    return out;
  }

};

template< class ScalarT, int order>
int Coefficients<ScalarT, order>::currentDegree = 0;

template< class ScalarT, int order>
Coefficients<ScalarT, order>* Coefficients<ScalarT, order>::imaginaryUnit = 0;

template< class ScalarT, int order>
Coefficients<ScalarT, order>* Coefficients<ScalarT, order>::buffer = 0;

template< class ScalarT, int order>
Coefficients<ScalarT, order>* Coefficients<ScalarT, order>::zero = 0;

///Automatic diffrentation, the multiplication implementation
template< typename ScalarT, int order>
inline const Coefficients<ScalarT, order> operator*(const Coefficients<ScalarT, order>& d1, const Coefficients<ScalarT, order>& d2){
  int j;
  Coefficients<ScalarT, order> r;
  int i;
  for(i=0; i <= Coefficients<ScalarT, order>::currentDegree; ++i){
    r.d[i] = 0;
    for(j=0; j <= i; ++j){
      r.d[i] += d1.d[j] * d2.d[i-j];
    }
  }
  return r;
}

///Automatic diffrentation, the addition implementation
template< typename ScalarT, int order>
inline const Coefficients<ScalarT, order> operator+(const Coefficients<ScalarT, order>& d1, const Coefficients<ScalarT, order>& d2){
  int i;
  Coefficients<ScalarT, order> r;
  for(i=0; i <= Coefficients<ScalarT, order>::currentDegree; ++i){
    r[i] = d1.d[i] + d2.d[i];
  }
  return r;
}

///Automatic diffrentation, the subtraction implementation
template< typename ScalarT, int order>
inline const Coefficients<ScalarT, order> operator-(const Coefficients<ScalarT, order>& d1, const Coefficients<ScalarT, order>& d2){
  int i;
  Coefficients<ScalarT, order> r;
  for(i=0; i <= Coefficients<ScalarT, order>::currentDegree; ++i){
    r[i] = d1.d[i] - d2.d[i];
  }
  return r;
}

template< typename ScalarT, int order>
inline const Coefficients<ScalarT, order> operator*(double d, const Coefficients<ScalarT, order>& c){
  int i;
  Coefficients<ScalarT, order> r;
  for(i=0; i <= Coefficients<ScalarT, order>::currentDegree; ++i){
    r[i] = d * c[i];
  }
  return r;
}

template< typename ScalarT, int order>
inline const Coefficients<ScalarT, order> operator*(typename Coefficients<ScalarT, order>::InternalType s, const Coefficients<ScalarT, order>& c){
  int i;
  Coefficients<ScalarT, order> r;
  for(i=0; i <= Coefficients<ScalarT, order>::currentDegree; ++i){
    r[i] = s * c[i];
  }
  return r;
}

template< typename ScalarT, int order>
inline const Coefficients<ScalarT, order> conjugate(const Coefficients<ScalarT, order>& coeff){
  Coefficients<ScalarT, order> r;
  int i;
  for(i=0; i < Coefficients<ScalarT, order>::dim; ++i){
    r.d[i] = conjugate(coeff[i]);
  }
  return r;
}

}
}

#endif /* DERIVATIVES_H_ */
