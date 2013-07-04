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
 * DPDEContainer.h
 *
 *  Created on: Jan 15, 2012
 *      Author: cyranka
 */

#ifndef DPDECONTAINER_H_
#define DPDECONTAINER_H_

namespace capd{
namespace jaco{

enum Subspace{non, odd, even};

///this realValuedL2 value is needed, when container is storing L_2 coefficients of a real valued function
enum Solution{complexValued, realValued, realValuedL2};

enum Solution2{other, realValuedOdd, realValuedEven};

/**Used for optimizing multiplication of Jets. For example for real valued and even / odd solutions some of the elements
 * are zero, therefore some of the multiplications can be avoided.
 */
class DPDEContainer{
public:
  int subspaceType;
  int solutionType;
  int solutionType2;
  bool baseReZero; ///<if the real part of the base is zero
  bool baseImZero; ///<if the imag part of the base is zero
  bool partialReReZero; ///<if the real part of partial derivative with resp to Re is zero
  bool partialReImZero; ///<if the imag part of partial derivative with resp to Re is zero
  bool partialImReZero; ///<if the real part of partial derivative with resp to Im is zero
  bool partialImImZero; ///<if the imag part of partial derivative with resp to Im is zero

  DPDEContainer() : subspaceType(non), solutionType(complexValued), solutionType2(other), baseReZero(false), baseImZero(false), partialReReZero(false), 
                    partialReImZero(false), partialImReZero(false), partialImImZero(false){}

  inline int productSubspace(int subspace1, int subspace2) const{
    if(subspace1 == odd && subspace2 == odd)
      return even;
    if(subspace1 == even && subspace2 == even)
      return even;
    if((subspace1 == even && subspace2 == odd) || (subspace1 == odd && subspace2 == even))
      return odd;
    return non;
  }
  inline int productSolutionType(int type1, int type2) const{
    if(type1 == complexValued || type2 == complexValued)
      return complexValued;
    return realValued; //we have realValuedOdd and realValuedEven
  }
  inline int productSolutionType2(int type1, int type2) const{
    if(type1 == realValuedEven && type2 == realValuedEven)
      return realValuedEven;
    if(type1 == realValuedOdd && type2 == realValuedOdd)
      return realValuedOdd;
    return other;
  }
  inline int sumSubspace(int subspace1, int subspace2) const{
    if(subspace1 == odd && subspace2 == odd)
      return odd;
    if(subspace1 == even && subspace2 == even)
      return even;
    return non;
  }
  inline int sumSolutionType(int type1, int type2) const{
    if(type1 == realValued && type2 == realValued)
      return realValued;
    return complexValued;
  }
  inline int sumSolutionType2(int type1, int type2) const{
    if(type1 == realValuedOdd && type2 == realValuedOdd)
      return realValuedOdd;
    if(type1 == realValuedEven && type2 == realValuedEven)
      return realValuedEven;
    return other;
  }
  
  inline DPDEContainer& multiply(const DPDEContainer& c1, const DPDEContainer& c2){
    DPDEContainer r;
    r.baseReZero = (c1.baseReZero || c2.baseReZero) && (c1.baseImZero || c2.baseImZero);
    r.baseImZero = (c1.baseReZero || c2.baseImZero) && (c1.baseImZero || c2.baseReZero);
    r.partialReReZero = (c1.partialReReZero || c2.baseReZero) && (c1.baseReZero || c2.partialReReZero) &&
                      (c1.baseImZero || c2.partialReImZero) && (c1.partialReImZero || c2.baseImZero);
    r.partialImReZero = (c1.partialImReZero || c2.baseReZero) && (c1.baseReZero || c2.partialImReZero) &&
                      (c1.baseImZero || c2.partialImImZero) && (c1.partialImImZero || c2.baseImZero);
    r.partialReImZero = (c1.partialReImZero || c2.baseReZero) && (c1.baseReZero || c2.partialReImZero) &&
                      (c1.baseImZero || c2.partialReReZero) && (c1.partialReReZero || c2.baseImZero);
    r.partialImImZero = (c1.partialImImZero || c2.baseReZero) && (c1.baseReZero || c2.partialImImZero) &&
                      (c1.baseImZero || c2.partialImReZero) && (c1.partialImReZero || c2.baseImZero);
    
    r.subspaceType = productSubspace(c1.subspaceType, c2.subspaceType);
    r.solutionType = productSolutionType(c1.solutionType, c2.solutionType);
    r.solutionType2 = productSolutionType2(c1.solutionType2, c2.solutionType2);
    setSubspaceType(r);
    return *this;   
  }

  inline DPDEContainer& operator*=(const DPDEContainer& c){
    *this = multiply(*this, c);
    return *this;
  }

  inline void setSubspaceType(const DPDEContainer& c){
    baseReZero = c.baseReZero;
    baseImZero = c.baseImZero;
    partialReReZero = c.partialReReZero;
    partialReImZero = c.partialReImZero;
    partialImReZero = c.partialImReZero;
    partialImImZero = c.partialImImZero;
    subspaceType = c.subspaceType;
    solutionType = c.solutionType;
    solutionType2 = c.solutionType2;
  }

  inline DPDEContainer& operator=(const DPDEContainer& c){
    setSubspaceType(c);
    return *this;
  }

  inline void multiplyByImUnit(){
    bool t = baseReZero;
    baseReZero = baseImZero;
    baseImZero = t;
    t = partialReReZero;
    partialReReZero = partialReImZero;
    partialReImZero = t;
    t = partialImReZero;
    partialImReZero = partialImImZero;
    partialImImZero = t;

    if(subspaceType == odd)
      subspaceType = even;
    else
      if(subspaceType == even)
        subspaceType = odd;
  }

  ///TODO: when a multiplication by complex is performed then the container type is updated
  inline void multiplyByFullComplex(){
    bool t = baseReZero && baseImZero;
    baseReZero = t;
    baseImZero = t;
    t = partialReReZero && partialReImZero;
    partialReReZero = t;
    partialReImZero = t;
    t = partialImReZero && partialImImZero;
    partialImReZero = t;
    partialImImZero = t;
    subspaceType = non;
    solutionType = complexValued;
  }

  inline void setToZero(){
    baseReZero = true;
    baseImZero = true;
    partialReReZero = true;
    partialReImZero = true;
    partialImReZero = true;
    partialImImZero = true;
  }

  inline void setBaseToZero(){
    baseReZero = true;
    baseImZero = true;
  }

  inline void setVariationalPartToId(){
    partialReReZero = false;
    partialImImZero = false;
    
    partialReImZero = true;
    partialImReZero = true;
  }

  inline void setRealPartToZero(){
    baseReZero = true;
    partialReReZero = true;
    partialImReZero = true;
    
    baseImZero = false;
    partialReImZero = false;
    partialImImZero = false;
  }

  inline void setImaginaryPartToZero(){
    baseImZero = true;
    partialReImZero = true;
    partialImImZero = true;
    
    baseReZero = false;
    partialImReZero = false;
    partialReReZero = false;
  }

  inline void projectContainerOntoRealSpace(){
    baseImZero = true;
    partialReImZero = true;
    partialImImZero = true;

    baseReZero = false;
    partialImReZero = false;
    partialReReZero = false;
  }

  inline void projectContainerOntoImaginarySpace(){
    baseReZero = true;
    partialReReZero = true;
    partialImReZero = true;

    baseImZero = false;
    partialReImZero = false;
    partialImImZero = false;
  }

  inline DPDEContainer& operator+=(const DPDEContainer& c){
    baseReZero = baseReZero && c.baseReZero;
    baseImZero = baseImZero && c.baseImZero;
    partialReReZero = partialReReZero && c.partialReReZero;
    partialImReZero = partialImReZero && c.partialImReZero;
    partialReImZero = partialReImZero && c.partialReImZero;
    partialImImZero = partialImImZero && c.partialImImZero;

    subspaceType = sumSubspace(subspaceType, c.subspaceType);
    solutionType = sumSolutionType(solutionType, c.solutionType);
    solutionType2 = sumSolutionType2(solutionType2, c.solutionType2);
    return *this;
  }

  inline int isRealValued() const{
    if(solutionType2 == other)
      return false;
    return solutionType;
  }

  inline bool isOdd() const{
    if(subspaceType == odd || solutionType2 == realValuedOdd)
      return true;
    return false;
  }

  inline bool isEven() const{
    if(subspaceType == even || solutionType2 == realValuedEven)
      return true;
    return false;
  }

  inline bool isRealValuedOdd() const{
    if(solutionType2 == realValuedOdd)
      return true;
    return false;
  }

  inline bool isRealValuedEven() const{
    if(solutionType2 == realValuedEven)
      return true;
    return false;
  }

  inline bool isOther() const{
    if(solutionType2 == other)
      return true;
    return false;
  }

  void setToRealValuedOdd(){
    solutionType2 = realValuedOdd;
    solutionType = realValued;
    subspaceType = odd;
    baseImZero = false;
    baseReZero = true;
    partialReReZero = false;
    partialReImZero = true;
    partialImReZero = true;
    partialImImZero = false;
  }
  
  void setToRealValuedEven(){
    solutionType2 = realValuedEven;
    solutionType = realValued;
    subspaceType = even;
    baseImZero = true;
    baseReZero = false;
    partialReReZero = false;
    partialReImZero = true;
    partialImReZero = true;
    partialImImZero = false;
  }
  
  void setToRealValued(){
    solutionType2 = other;
    solutionType = realValued;
    subspaceType = non;
    baseImZero = false;
    baseReZero = false;
    partialReReZero = false;
    partialReImZero = false;
    partialImReZero = false;
    partialImImZero = false;
  }
  
  virtual void projectOntoSubspace(){};

  friend std::ostream& operator<<(std::ostream& out, const DPDEContainer& c){
    out << "[" << c.subspaceType <<", "<<c.solutionType<<", "<<c.solutionType2<<"]\n";
    out <<"("<<c.baseReZero<<", "<<c.baseImZero<<")\n";
    out <<c.partialReReZero<<" "<<c.partialReImZero<<"\n";
    out <<c.partialImReZero<<" "<<c.partialImImZero<<"\n";
    return out;
  }

};

inline DPDEContainer operator+(const DPDEContainer& c1, const DPDEContainer& c2){
  DPDEContainer r = c1;
  r += c2;
  return r;
}

inline DPDEContainer operator*(const DPDEContainer& c1, const DPDEContainer& c2){
  DPDEContainer r = c1;
  r *= c2;
  return r;
}

}
}

#endif /* DPDECONTAINER_H_ */
