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
 * FirstOrderJet.h
 *
 *  Created on: Sep 22, 2011
 *      Author: cyranka
 */

#ifndef FIRSTORDERJET_H_
#define FIRSTORDERJET_H_

#include "capd/vectalg/Vector.hpp"
#include "DPDEContainer.h"

namespace capd{
namespace jaco{

/**This value determines an optimization control of operations performed on Jets (additions and multiplications). We want to avoid
 * additions and multiplications by zero.
 */
enum OptimizationControl{global, local};

///First Order Jet class. Optimized for use in the dpdes integrator.
template< class DerivativeT, int D>
class FirstOrderJet : public capd::jaco::DPDEContainer{
public:
  typedef DerivativeT DerivativeType;
  typedef typename DerivativeType::ScalarType ScalarType;
  typedef typename ScalarType::ScalarType InternalType;
  typedef capd::vectalg::Vector<DerivativeT, D> VectorType;
  typedef typename VectorType::iterator IteratorType;
  typedef typename VectorType::const_iterator ConstIteratorType;
  typedef typename DerivativeType::MatrixType MatrixType;
  typedef capd::jaco::DPDEContainer DPDEContainer;

//  static SubspaceType* subspace; /**<keep the reference to a subspace in order to get the index iterator (needed, cause indices may be
//                             *stored as a map of index -> internal index (index in the array of modes)
//                             */
  static int dim;

  static FirstOrderJet* buffer;

  static capd::jaco::DPDEContainer* initialCondition;
  static capd::jaco::DPDEContainer* currentDPDEContainer;
  
  static int optimizationControl;

  static bool initialized;

  ScalarType val; ///< base part (enclosure on which the variational part is calculated)
  ScalarType val2; ///< second free coefficient (not taken into account when calculating the variational part)
  VectorType ksi; ///< first order coefficients (an array)
  IteratorType iter;
  ConstIteratorType iter2;

  inline FirstOrderJet() : DPDEContainer(*initialCondition), val(), val2(), ksi(FirstOrderJet::hasBeenInitialized(dim)){
//    std::cout << "|";
  }
  inline FirstOrderJet(const ScalarType& s) : DPDEContainer(*initialCondition), val(s), val2(s), ksi(FirstOrderJet::hasBeenInitialized(dim)){
//    std::cout << ".";
    if(s == 0)
      setToZero();
  }
  inline FirstOrderJet(const InternalType& i) : DPDEContainer(*initialCondition), val(i), val2(i), ksi(FirstOrderJet::hasBeenInitialized(dim)){
//    std::cout << ",";
    if(i == 0)
      setToZero();
  }
  
  /**Initializes the first order jets, 
   * @param dim the dimension of the vector
   * @param container the DPDEContainer instance defining the subpspace of the initial condition, e.g. RealValuedEven ,
   *    this is important for the optimizations that are performed while algebraic operations on jets.
   */    
  static void initialize(int dim, DPDEContainer& container){
    FirstOrderJet::dim = dim;
    FirstOrderJet::initialCondition = new DPDEContainer(container);    
    FirstOrderJet::currentDPDEContainer = new DPDEContainer(container);
    FirstOrderJet::initialized = 1;
    FirstOrderJet::buffer = new FirstOrderJet();
  }
  
  static void destroy(){
    delete FirstOrderJet::initialCondition;
    delete FirstOrderJet::currentDPDEContainer;
    delete FirstOrderJet::buffer;
  }
  
  /**Filter through this function in order to check if the static variable dim has been properly initialized.
   */
  static int hasBeenInitialized(int dim){
    if(dim <= 0){
      std::cerr << "FirstOrderJet class (static variable dim) has to be initialized!!!\nCall FisrtOrderJet::initialize()\n";
      throw std::runtime_error("FirstOrderJet class (static variable dim) has to be initialized!!!\nCall FisrtOrderJet::initialize()\n");
    }
    return dim;
  }
  
  static bool initialConditionIsRealValued(){
    if(initialCondition->isRealValued())
      return true;
    return false;
  }
  
  static void switchToRealValued(){
    currentDPDEContainer->solutionType = capd::jaco::realValued;
    if(currentDPDEContainer->isRealValuedOdd()){
      currentDPDEContainer->baseImZero = false;
      currentDPDEContainer->baseReZero = true;
    }
    if(currentDPDEContainer->isRealValuedEven()){      
      currentDPDEContainer->baseImZero = true;
      currentDPDEContainer->baseReZero = false;
    }
    currentDPDEContainer->partialReReZero = false;
    currentDPDEContainer->partialReImZero = true;
    currentDPDEContainer->partialImReZero = true;
    currentDPDEContainer->partialImImZero = false;
  }

  static void switchToComplexValued(){
    //do not touch solutionType2, this stores information if the solution is real-valued and odd/even (only complex component of the 
    //derivative is taken into account)
    currentDPDEContainer->solutionType = capd::jaco::complexValued;
    currentDPDEContainer->baseImZero = false;
    currentDPDEContainer->baseReZero = false;
    currentDPDEContainer->partialReReZero = false;
    currentDPDEContainer->partialReImZero = false;
    currentDPDEContainer->partialImReZero = false;
    currentDPDEContainer->partialImImZero = false;
  }

  static void switchToRealValuedL2(){
    //do not touch solutionType2, this stores information if the solution is real-valued and odd/even (only complex component of the 
    //derivative is taken into account)
    currentDPDEContainer->solutionType = capd::jaco::realValuedL2;
    currentDPDEContainer->baseImZero = true;
    currentDPDEContainer->baseReZero = false;
    currentDPDEContainer->partialReReZero = false;
    currentDPDEContainer->partialReImZero = true;
    currentDPDEContainer->partialImReZero = false;
    currentDPDEContainer->partialImImZero = true;
  }

  static void switchToLocalOptimization(){
    optimizationControl = local;
  }

  static void switchToGlobalOptimization(){
    optimizationControl = global;
  }

  /**Sets the DPDEContainer defining the subspace to which belongs the stored modes to the initial condition's subspace.    
   */
  static void setContainer(DPDEContainer& c){
    c = *FirstOrderJet::initialCondition;
  }
  
  const ScalarType& value() const{
    return val;
  }

  ScalarType& value(){
    return val;
  }

  inline InternalType normMax() const{;
    InternalType max = abs(val.re), t;
    if((t = rightBound(abs(val.im))) > max) max = t;
    return max;
  }

  inline static FirstOrderJet i(){
    FirstOrderJet r;
    r.val = ScalarType::i();
    r.val2 = ScalarType::i();
    return r;
  }

  inline bool isImUnit(){
    if(val.isImUnit())
      return true;
    return false;
  }

  inline FirstOrderJet& conjugate(){
    val.conjugate();
    val2.conjugate();
//    val3.conjugate();
    int i;
    for(i=0; i < FirstOrderJet<DerivativeT, D>::dim; ++i){
      if(optimizationControl == global){
        if(currentDPDEContainer->isRealValuedEven() || currentDPDEContainer->isOther()){
          if(!currentDPDEContainer->partialReImZero)
            ksi[i].partialWRespToRe.conjugate();
        }
        if(currentDPDEContainer->isRealValuedOdd() || currentDPDEContainer->isOther()){
          if(!currentDPDEContainer->partialImImZero)
            ksi[i].partialWRespToIm.conjugate();
        }

//        if(currentDPDEContainer->isEven() || !currentDPDEContainer->subspaceType){
//          if(!currentDPDEContainer->isRealValued()){
//            ksi[i].partialWRespToRe.conjugate();
//          }
//        }
//        if(currentDPDEContainer->isOdd() == capd::jaco::odd || !currentDPDEContainer->subspaceType){
//          if(!currentDPDEContainer->isRealValued())
//            ksi[i].partialWRespToIm.conjugate();
//          else
//            ksi[i].partialWRespToIm.im = -ksi[i].partialWRespToIm.im;
//        }
      }
      if(optimizationControl == local){
        if(isRealValuedEven() || isOther()){
          if(!partialReImZero)
            ksi[i].partialWRespToRe.conjugate();
        }
        if(isRealValuedOdd() || isOther()){
          if(!partialImImZero)
            ksi[i].partialWRespToIm.conjugate();
        }
      }
    }
    return *this;
  }

  inline FirstOrderJet& operator=(const FirstOrderJet& foj){
    (DPDEContainer&)*this = (DPDEContainer&)foj;
    val = foj.val;
    val2 = foj.val2;
//    val3 = foj.val3;
    iter = ksi.begin();
    iter2 = foj.ksi.begin();
    while(iter != ksi.end()){
      if(optimizationControl == global){
        if(currentDPDEContainer->isRealValuedEven() || currentDPDEContainer->isOther()){
          if(!currentDPDEContainer->partialReReZero)
            (*iter).partialWRespToRe.re = (*iter2).partialWRespToRe.re;
          if(!currentDPDEContainer->partialReImZero)
            (*iter).partialWRespToRe.im = (*iter2).partialWRespToRe.im;
        }
        if(currentDPDEContainer->isRealValuedOdd() || currentDPDEContainer->isOther()){
          if(!currentDPDEContainer->partialImReZero)
            (*iter).partialWRespToIm.re = (*iter2).partialWRespToIm.re;
          if(!currentDPDEContainer->partialImImZero)
            (*iter).partialWRespToIm.im = (*iter2).partialWRespToIm.im;
        }
      }
      if(optimizationControl == local){
        if(foj.isRealValuedEven() || foj.isOther()){
          if(!foj.partialReReZero)
            (*iter).partialWRespToRe.re = (*iter2).partialWRespToRe.re;
          else
            (*iter).partialWRespToRe.re = 0;
          if(!foj.partialReImZero)
            (*iter).partialWRespToRe.im = (*iter2).partialWRespToRe.im;
          else
            (*iter).partialWRespToRe.im = 0;
        }
        if(foj.isRealValuedOdd() || foj.isOther()){
          if(!foj.partialImReZero)
            (*iter).partialWRespToIm.re = (*iter2).partialWRespToIm.re;
          else
            (*iter).partialWRespToIm.re = 0;
          if(!foj.partialImImZero)
            (*iter).partialWRespToIm.im = (*iter2).partialWRespToIm.im;
          else
            (*iter).partialWRespToIm.im = 0;
        }
      }
    	
//      if(currentDPDEContainer->isEven() || !currentDPDEContainer->subspaceType)
//        (*iter).partialWRespToRe = (*iter2).partialWRespToRe;
//      if(currentDPDEContainer->isOdd() || !currentDPDEContainer->subspaceType)
//        (*iter).partialWRespToIm = (*iter2).partialWRespToIm;
      iter++; iter2++;
    }
    return *this;
  }

  inline FirstOrderJet& operator=(ScalarType s){
    setFreeCoeff(s);
    val2 = s;
    ksi = 0;
    return *this;
  }

  inline void setFreeCoeff(const ScalarType& s){
    val = s;
    if(s == 0){
      setToZero();
    }else{
      if(s.isImUnit()){
        baseReZero = true;
        baseImZero = false;
      }else{
        if(s.isFullComplex()){
          baseReZero = false;
          baseImZero = false;
        }else{
          baseReZero = false;
          baseImZero = true;
        }
      }
    }
  }

  inline void setSecondFreeCoeff(const ScalarType& s){
    val2 = s;
  }

//  inline void setThirdFreeCoeff(const ScalarType& s){
//    val3 = s;
//  }

  inline const ScalarType& secondFreeCoeff() const{
    return val2;
  }

//  inline const ScalarType& thirdFreeCoeff() const{
//    return val3;
//  }

  inline FirstOrderJet& operator+=(const FirstOrderJet& foj){
//    val += foj.val;
//    val2 += foj.val2;
////    val3 += foj.val3;
//    int i;
//    for(i=0; i < FirstOrderJet::dim; ++i){
//      if(optimizationControl == global){
//        if(currentDPDEContainer->isEven() || !currentDPDEContainer->subspaceType){
//          if(!currentDPDEContainer->isRealValued())
//            ksi[i].partialWRespToRe += foj.ksi[i].partialWRespToRe;
//          else
//            ksi[i].partialWRespToRe.re += foj.ksi[i].partialWRespToRe.re;
//        }
//        if(currentDPDEContainer->isOdd() || !currentDPDEContainer->subspaceType){
//          if(!currentDPDEContainer->isRealValued())
//            ksi[i].partialWRespToIm += foj.ksi[i].partialWRespToIm;
//          else{
//            if(currentDPDEContainer->isRealValued() == capd::jaco::realValuedL2)
//              ksi[i].partialWRespToIm.re += foj.ksi[i].partialWRespToIm.re;
//            else
//              ksi[i].partialWRespToIm.im += foj.ksi[i].partialWRespToIm.im;
//          }
//        }
//      }
//    }
    *this = *this + foj;
    return *this;
  }
  
  inline FirstOrderJet& operator+=(const ScalarType& s){
    val += s;
    val2 += s;
    return *this;
  }

  inline FirstOrderJet& operator-=(const FirstOrderJet& foj){
//    (DPDEContainer&)*this += (DPDEContainer&)foj;
//    val -= foj.val;
//    val2 -= foj.val2;
////    val3 += foj.val3;
//    int i;
//    for(i=0; i < FirstOrderJet::dim; ++i){
//      if(optimizationControl == global){
//        if(currentDPDEContainer->isEven() || !currentDPDEContainer->subspaceType){
//          if(!currentDPDEContainer->isRealValued())
//            ksi[i].partialWRespToRe -= foj.ksi[i].partialWRespToRe;
//          else
//            ksi[i].partialWRespToRe.re -= foj.ksi[i].partialWRespToRe.re;
//        }
//        if(currentDPDEContainer->isOdd() || !currentDPDEContainer->subspaceType){
//          if(!currentDPDEContainer->isRealValued())
//            ksi[i].partialWRespToIm -= foj.ksi[i].partialWRespToIm;
//          else{
//            if(currentDPDEContainer->isRealValued() == capd::jaco::realValuedL2)
//              ksi[i].partialWRespToIm.re -= foj.ksi[i].partialWRespToIm.re;
//            else
//              ksi[i].partialWRespToIm.im -= foj.ksi[i].partialWRespToIm.im;
//          }
//        }
//      }
//    }
    *this = *this - foj;
    return *this;
  }


  inline FirstOrderJet& operator*=(const FirstOrderJet& foj){
//    int i;
//    for(i=0; i < FirstOrderJet::dim; ++i){
//      if(optimizationControl == global){
//        if(currentDPDEContainer->isEven() || !currentDPDEContainer->subspaceType){
//          if(!currentDPDEContainer->isRealValued())
//            ksi[i].partialWRespToRe = (foj.val * ksi[i].partialWRespToRe) + (val * foj.ksi[i].partialWRespToRe);
//          else
//            ksi[i].partialWRespToRe.re = (foj.val.re * ksi[i].partialWRespToRe.re) + (val.re * foj.ksi[i].partialWRespToRe.re);
//        }
//        if(currentDPDEContainer->isOdd() || !currentDPDEContainer->subspaceType){
//          if(!currentDPDEContainer->isRealValued())
//            ksi[i].partialWRespToIm = (foj.val * ksi[i].partialWRespToIm) + (val * foj.ksi[i].partialWRespToIm);
//          else{
//            if(currentDPDEContainer->isRealValued() == capd::jaco::realValuedL2)
//              ksi[i].partialWRespToIm.re = (foj.val.re * ksi[i].partialWRespToIm.re) + (val.re * foj.ksi[i].partialWRespToIm.re);
//            else
//              ksi[i].partialWRespToIm.re = -((foj.val.im * ksi[i].partialWRespToIm.im) + (val.im * foj.ksi[i].partialWRespToIm.im));
//          }
//        }
//      }
//    }
//    val *= foj.val;
//    val2 *= foj.val2;
////    val3 *= foj.val3;
    *this = *this * foj;

    return *this;
  }


  inline FirstOrderJet& operator*=(const ScalarType& s){

    int i;
    ScalarType t;
    for(i=0; i < FirstOrderJet::dim; ++i){
      if(optimizationControl == global){
        if(currentDPDEContainer->isRealValuedEven() || currentDPDEContainer->isOther()){
          t = ksi[i].partialWRespToRe;
          ksi[i].partialWRespToRe = 0;
          if(!currentDPDEContainer->partialReReZero){
            ksi[i].partialWRespToRe.re += t.re * s.re;
            ksi[i].partialWRespToRe.im += t.re * s.im;
          }
          if(!currentDPDEContainer->partialReImZero){
            ksi[i].partialWRespToRe.re += - (t.im * s.im);
            ksi[i].partialWRespToRe.im += t.im * s.re;
          }
        }
        if(currentDPDEContainer->isRealValuedOdd() || currentDPDEContainer->isOther()){
          t = ksi[i].partialWRespToIm;
          ksi[i].partialWRespToIm = 0;
          if(!currentDPDEContainer->partialImReZero){
            ksi[i].partialWRespToIm.re += t.re * s.re;
            ksi[i].partialWRespToIm.im += t.re * s.im;
          }
          if(!currentDPDEContainer->partialImImZero){
            ksi[i].partialWRespToIm.re += - (t.im * s.im);
            ksi[i].partialWRespToIm.im += t.im * s.re;
          }
        }
//        if(currentDPDEContainer->isEven() || !currentDPDEContainer->subspaceType){
//          if(!currentDPDEContainer->isRealValued())
//            ksi[i].partialWRespToRe *= s;
//          else{
//            ksi[i].partialWRespToRe.im = ksi[i].partialWRespToRe.re * s.im;
//            ksi[i].partialWRespToRe.re *= s.re;
//          }
//        }
//        if(currentDPDEContainer->isOdd() || !currentDPDEContainer->subspaceType){
//          if(!currentDPDEContainer->isRealValued())
//            ksi[i].partialWRespToIm *= s;
//          else{
//            if(currentDPDEContainer->isRealValued() == capd::jaco::realValuedL2){
//              ksi[i].partialWRespToIm.im = ksi[i].partialWRespToIm.re * s.im;
//              ksi[i].partialWRespToIm.re *= s.re;
//            }else{
//              ksi[i].partialWRespToIm.re = -ksi[i].partialWRespToIm.im * s.im;
//              ksi[i].partialWRespToIm.im *= s.re;
//            }
//          }
//        }
      }
      if(optimizationControl == local){
        if(isRealValuedEven() || isOther()){
          t = ksi[i].partialWRespToRe;
          ksi[i].partialWRespToRe = 0;
          if(!partialReReZero){
            ksi[i].partialWRespToRe.re += t.re * s.re;
            ksi[i].partialWRespToRe.im += t.re * s.im;
          }
          if(!partialReImZero){
            ksi[i].partialWRespToRe.re += - (t.im * s.im);
            ksi[i].partialWRespToRe.im += t.im * s.re;
          }
        }
        if(isRealValuedOdd() || isOther()){
          t = ksi[i].partialWRespToIm;
          ksi[i].partialWRespToIm = 0;
          if(!partialImReZero){
            ksi[i].partialWRespToIm.re += t.re * s.re;
            ksi[i].partialWRespToIm.im += t.re * s.im;
          }
          if(!partialImImZero){
            ksi[i].partialWRespToIm.re += - (t.im * s.im);
            ksi[i].partialWRespToIm.im += t.im * s.re;
          }
        }
      }
    }
    val *= s;
    val2 *= s;
    if(s == 0)
      setToZero();
    if(s.isImUnit())
      multiplyByImUnit();
    if(s.isFullComplex())
      multiplyByFullComplex();
    return *this;
  }

  inline FirstOrderJet& operator*=(const InternalType& s){
    if(s == 0)
      setToZero();
    int i;
    for(i=0; i < FirstOrderJet::dim; ++i){
      if(optimizationControl == global){
        if(currentDPDEContainer->isRealValuedEven() || currentDPDEContainer->isOther()){
          if(!currentDPDEContainer->partialReReZero){
            ksi[i].partialWRespToRe.re *= s;
          }
          if(!currentDPDEContainer->partialReImZero){
            ksi[i].partialWRespToRe.im *= s;
          }
        }
        if(currentDPDEContainer->isRealValuedOdd() || currentDPDEContainer->isOther()){
          if(!currentDPDEContainer->partialImReZero){
            ksi[i].partialWRespToIm.re *= s;
          }
          if(!currentDPDEContainer->partialImImZero){
            ksi[i].partialWRespToIm.im *= s;
          }
        }
//        if(currentDPDEContainer->isEven() || !currentDPDEContainer->subspaceType){
//          if(!currentDPDEContainer->isRealValued())
//            ksi[i].partialWRespToRe *= s;
//          else
//            ksi[i].partialWRespToRe.re *= s;
//        }
//        if(currentDPDEContainer->isOdd() || !currentDPDEContainer->subspaceType){
//          if(!currentDPDEContainer->isRealValued())
//            ksi[i].partialWRespToIm *= s;
//          else{
//            if(currentDPDEContainer->isRealValued() == capd::jaco::realValuedL2){
//              ksi[i].partialWRespToIm.re *= s;
//            }else{
//              ksi[i].partialWRespToIm.im *= s;
//            }
//          }
//        }
      }
      if(optimizationControl == local){
        if(isRealValuedEven() || isOther()){
          if(!partialReReZero){
            ksi[i].partialWRespToRe.re *= s;
          }
          if(!partialReImZero){
            ksi[i].partialWRespToRe.im *= s;
          }
        }
        if(isRealValuedOdd() || isOther()){
          if(!partialImReZero){
            ksi[i].partialWRespToIm.re *= s;
          }
          if(!partialImImZero){
            ksi[i].partialWRespToIm.im *= s;
          }
        }
      }
    }
    val *= s;
    val2 *= s;
    return *this;
  }


  /**
   * @param position this is the position in the vector where the Identity matrix should be put.
   */
  inline void setVariationalPartToId(int position/*, const MatrixType& m = MatrixType::Identity(2)*/){
    ((DPDEContainer&)*this).setVariationalPartToId();
    ksi[position].setToId();
  }

  /**
   * @param position this is the position in the vector where the conjugate of the Identity matrix should be put.
   */
  inline void setVariationalPartToIdConjugate(int position/*, const MatrixType& m = MatrixType::Identity(2)*/){
    ((DPDEContainer&)*this).setVariationalPartToId();
    ksi[position].setToIdConjugate();
  }

  inline DerivativeType variationalPart(int position) const{
    return ksi[position];
  }

  inline void setConjugate(const FirstOrderJet& foj){
    (DPDEContainer&)*this = (DPDEContainer&)foj;
    val = foj.val;
    val.conjugate();
    val2 = foj.val2;
    val2.conjugate();
//    val3 = foj.val3;
//    val3.conjugate();
    int i;
    for(i=0; i < FirstOrderJet::dim; ++i){
      if(optimizationControl == global){
        if(currentDPDEContainer->isRealValuedEven() || currentDPDEContainer->isOther()){
          ksi[i].partialWRespToRe = foj.ksi[i].partialWRespToRe;
          if(!currentDPDEContainer->partialReImZero)
            ksi[i].partialWRespToRe.conjugate();
        }
        if(currentDPDEContainer->isRealValuedOdd() || currentDPDEContainer->isOther()){
          ksi[i].partialWRespToIm = foj.ksi[i].partialWRespToIm;
          if(!currentDPDEContainer->partialImImZero)
            ksi[i].partialWRespToIm.conjugate();
        }
//        if(currentDPDEContainer->isEven() || !currentDPDEContainer->subspaceType){
//          ksi[i].partialWRespToRe = foj.ksi[i].partialWRespToRe;
//          ksi[i].partialWRespToRe.conjugate();
//        }
//        if(currentDPDEContainer->isOdd() || !currentDPDEContainer->subspaceType){
//          ksi[i].partialWRespToIm = foj.ksi[i].partialWRespToIm;
//          ksi[i].partialWRespToIm.conjugate();
//        }
      }
      if(optimizationControl == local){
        ksi[i].partialWRespToRe = 0;
        if(foj.isRealValuedEven() || foj.isOther()){
          ksi[i].partialWRespToRe = foj.ksi[i].partialWRespToRe;
          if(!foj.partialReImZero)
            ksi[i].partialWRespToRe.conjugate();
        }
        ksi[i].partialWRespToIm = 0;
        if(foj.isRealValuedOdd() || foj.isOther()){
          ksi[i].partialWRespToIm = foj.ksi[i].partialWRespToIm;
          if(!foj.partialImImZero)
            ksi[i].partialWRespToIm.conjugate();
        }
      }
    }

  }

  inline void setImaginaryPartToZero(){
    ((DPDEContainer&)*this).setImaginaryPartToZero();
    int i;
    for(i=0; i < FirstOrderJet::dim; ++i)
      ksi[i].setImaginaryPartToZero();
    val.setImaginaryPartToZero();
    val2.setImaginaryPartToZero();
//    val3.setImaginaryPartToZero();
  }

  /**this is not acually used*/
  inline void setRealPartToZero(){
    ((DPDEContainer&)*this).setRealPartToZero();
    int i;
    for(i=0; i < FirstOrderJet::dim; ++i)
      ksi[i].setRealPartToZero();
    val.setRealPartToZero();
    val2.setRealPartToZero();
  }

  /**This is different than setting real part to zero. In case of complex derivatives this leaves only \partial{Im}/\partial{Im},
   * while setRealPartToZero leaves \partial{Im}/\partial{Im} and \partial{Re}/\partial{Im}
   */
  inline void projectOntoImaginarySpace(){
    ((DPDEContainer&)*this).projectContainerOntoImaginarySpace();
    int i;
    for(i=0; i < FirstOrderJet::dim; ++i)
      ksi[i].projectOntoImaginarySpace();
    val.projectOntoImaginarySpace();
    val2.projectOntoImaginarySpace();
  }

  inline void projectOntoRealSpace(){
    ((DPDEContainer&)*this).projectContainerOntoRealSpace();
    int i;
    for(i=0; i < FirstOrderJet::dim; ++i)
      ksi[i].projectOntoRealSpace();
    val.projectOntoRealSpace();
    val2.projectOntoRealSpace();
  }

  friend std::ostream& operator<<(std::ostream& out, const FirstOrderJet& foj){
    out << (DPDEContainer&)foj << "\n";
    out << "val:" << foj.val << ", val2:" << foj.val2 << "\nksi:\n";
    int i;
    for(i=0; i < FirstOrderJet::dim; ++i){
      out << foj.ksi[i] << "\n";
    }
    out << "\n";
    return out;
  }

};

//template<class DerivativeT, int D>
//SubspaceT* FirstOrderJet<DerivativeT, D>::subspace = 0;

template<class DerivativeT, int D>
int FirstOrderJet<DerivativeT, D>::optimizationControl = 0;

template<class DerivativeT, int D>
int FirstOrderJet<DerivativeT, D>::dim = -1;

template<class DerivativeT, int D>
capd::jaco::DPDEContainer* FirstOrderJet<DerivativeT, D>::initialCondition = 0;

template<class DerivativeT, int D>
capd::jaco::DPDEContainer* FirstOrderJet<DerivativeT, D>::currentDPDEContainer = 0;

template<class DerivativeT, int D>
FirstOrderJet<DerivativeT, D>* FirstOrderJet<DerivativeT, D>::buffer = 0;

template<class DerivativeT, int D>
bool FirstOrderJet<DerivativeT, D>::initialized = 0;

///TODO: try iterators
///19.11.11 important change (The order here is very important, otherwise if one of the foj1, foj2 is the buffer then the val has to be calculated at the end,
///because otherwise the new value (wrong value) of val will be used in ksi calculations.
template< class DerivativeT, int D>
inline const FirstOrderJet<DerivativeT, D>& operator*(const FirstOrderJet<DerivativeT, D>& foj1, const FirstOrderJet<DerivativeT, D>& foj2){
  //FirstOrderJet<DerivativeT, SubspaceT, D> r;
  ///The order here is very important, otherwise if one of the foj1, foj2 is the buffer then the val has to be calculated at the end,
  ///because otherwise the new value (wrong value) of val will be used in ksi calculations.
  (typename FirstOrderJet<DerivativeT, D>::DPDEContainer&)*FirstOrderJet<DerivativeT, D>::buffer = 
    (typename FirstOrderJet<DerivativeT, D>::DPDEContainer&)foj1 * (typename FirstOrderJet<DerivativeT, D>::DPDEContainer&)foj2;
  DerivativeT t;
  int i;
  for(i=0; i < FirstOrderJet<DerivativeT, D>::dim; ++i){
    if(FirstOrderJet<DerivativeT, D>::optimizationControl == global){
      if(FirstOrderJet<DerivativeT, D>::currentDPDEContainer->isRealValuedEven() || FirstOrderJet<DerivativeT, D>::currentDPDEContainer->isOther()){
        FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToRe = 0;
        if(!FirstOrderJet<DerivativeT, D>::currentDPDEContainer->partialReReZero){
          if(!FirstOrderJet<DerivativeT, D>::currentDPDEContainer->baseReZero)
            FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToRe.re += foj1.ksi[i].partialWRespToRe.re * foj2.val.re + foj2.ksi[i].partialWRespToRe.re * foj1.val.re;
          if(!FirstOrderJet<DerivativeT, D>::currentDPDEContainer->baseImZero)  
            FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToRe.im += foj1.ksi[i].partialWRespToRe.re * foj2.val.im + foj2.ksi[i].partialWRespToRe.re * foj1.val.im;
        }
        if(!FirstOrderJet<DerivativeT, D>::currentDPDEContainer->partialReImZero){
          if(!FirstOrderJet<DerivativeT, D>::currentDPDEContainer->baseReZero)
            FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToRe.im += foj1.ksi[i].partialWRespToRe.im * foj2.val.re + foj2.ksi[i].partialWRespToRe.im * foj1.val.re; 
          if(!FirstOrderJet<DerivativeT, D>::currentDPDEContainer->baseImZero)
            FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToRe.re += -(foj1.ksi[i].partialWRespToRe.im * foj2.val.im + foj2.ksi[i].partialWRespToRe.im * foj1.val.im);
        }
      }
      if(FirstOrderJet<DerivativeT, D>::currentDPDEContainer->isRealValuedOdd() || FirstOrderJet<DerivativeT, D>::currentDPDEContainer->isOther()){
        FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToIm = 0;
        if(!FirstOrderJet<DerivativeT, D>::currentDPDEContainer->partialImReZero){
          if(!FirstOrderJet<DerivativeT, D>::currentDPDEContainer->baseReZero)
            FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToIm.re += foj1.ksi[i].partialWRespToIm.re * foj2.val.re + foj2.ksi[i].partialWRespToIm.re * foj1.val.re;
          if(!FirstOrderJet<DerivativeT, D>::currentDPDEContainer->baseImZero)  
            FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToIm.im += foj1.ksi[i].partialWRespToIm.re * foj2.val.im + foj2.ksi[i].partialWRespToIm.re * foj1.val.im;
        }
        if(!FirstOrderJet<DerivativeT, D>::currentDPDEContainer->partialImImZero){
          if(!FirstOrderJet<DerivativeT, D>::currentDPDEContainer->baseReZero)
            FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToIm.im += foj1.ksi[i].partialWRespToIm.im * foj2.val.re + foj2.ksi[i].partialWRespToIm.im * foj1.val.re; 
          if(!FirstOrderJet<DerivativeT, D>::currentDPDEContainer->baseImZero)
            FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToIm.re += -(foj1.ksi[i].partialWRespToIm.im * foj2.val.im + foj2.ksi[i].partialWRespToIm.im * foj1.val.im);
        }
      }
//      if(FirstOrderJet<DerivativeT, D>::currentDPDEContainer->isEven() || !FirstOrderJet<DerivativeT, D>::currentDPDEContainer->subspaceType){
//        if(!FirstOrderJet<DerivativeT, D>::currentDPDEContainer->isRealValued())
//          FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToRe = (foj2.val * foj1.ksi[i].partialWRespToRe) + (foj1.val * foj2.ksi[i].partialWRespToRe);
//        else
//          FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToRe.re = (foj2.val.re * foj1.ksi[i].partialWRespToRe.re) + (foj1.val.re * foj2.ksi[i].partialWRespToRe.re);
//      }
//      if(FirstOrderJet<DerivativeT, D>::currentDPDEContainer->isOdd() || !FirstOrderJet<DerivativeT, D>::currentDPDEContainer->subspaceType){
//        if(!FirstOrderJet<DerivativeT, D>::currentDPDEContainer->isRealValued())
//          FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToIm = (foj2.val * foj1.ksi[i].partialWRespToIm) + (foj1.val * foj2.ksi[i].partialWRespToIm);
//        else{
//          if(FirstOrderJet<DerivativeT, D>::currentDPDEContainer->isRealValued() == capd::jaco::realValuedL2)
//            FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToIm.re = (foj2.val.re * foj1.ksi[i].partialWRespToIm.re) + (foj1.val.re * foj2.ksi[i].partialWRespToIm.re);
//          else
//            FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToIm.re = -((foj2.val.im * foj1.ksi[i].partialWRespToIm.im) + (foj1.val.im * foj2.ksi[i].partialWRespToIm.im));
//        }
//      }
    }
    if(FirstOrderJet<DerivativeT, D>::optimizationControl == local){
      if(FirstOrderJet<DerivativeT, D>::buffer->isRealValuedEven() || FirstOrderJet<DerivativeT, D>::buffer->isOther()){
        FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToRe = 0;
        if(!foj1.partialReReZero || !foj2.partialReReZero){
          if(!foj1.baseReZero || !foj2.baseReZero){
            FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToRe.re += foj1.ksi[i].partialWRespToRe.re * foj2.val.re + foj2.ksi[i].partialWRespToRe.re * foj1.val.re;
          }if(!foj1.baseImZero || !foj2.baseImZero){
            FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToRe.im += foj1.ksi[i].partialWRespToRe.re * foj2.val.im + foj2.ksi[i].partialWRespToRe.re * foj1.val.im;
          }
        }
        if(!foj1.partialReImZero || !foj2.partialReImZero){
          if(!foj1.baseReZero || !foj2.baseReZero){
            FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToRe.im += foj1.ksi[i].partialWRespToRe.im * foj2.val.re + foj2.ksi[i].partialWRespToRe.im * foj1.val.re; 
          }if(!foj1.baseImZero || !foj2.baseImZero){
            FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToRe.re += -(foj1.ksi[i].partialWRespToRe.im * foj2.val.im + foj2.ksi[i].partialWRespToRe.im * foj1.val.im);
          }
        }
      }
      if(FirstOrderJet<DerivativeT, D>::buffer->isRealValuedOdd() || FirstOrderJet<DerivativeT, D>::buffer->isOther()){
        FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToIm = 0;
        if(!foj1.partialImReZero || !foj2.partialImReZero){
          if(!foj1.baseReZero || !foj2.baseReZero){
            FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToIm.re += foj1.ksi[i].partialWRespToIm.re * foj2.val.re + foj2.ksi[i].partialWRespToIm.re * foj1.val.re;
          }
          if(!foj1.baseImZero || !foj2.baseImZero)  {
            FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToIm.im += foj1.ksi[i].partialWRespToIm.re * foj2.val.im + foj2.ksi[i].partialWRespToIm.re * foj1.val.im;
          }
        }
        if(!foj1.partialImImZero || !foj2.partialImImZero){
          if(!foj1.baseReZero || !foj2.baseReZero){
            FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToIm.im += foj1.ksi[i].partialWRespToIm.im * foj2.val.re + foj2.ksi[i].partialWRespToIm.im * foj1.val.re;
          }
          if(!foj1.baseImZero || !foj2.baseImZero){            
            FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToIm.re += -(foj1.ksi[i].partialWRespToIm.im * foj2.val.im + foj2.ksi[i].partialWRespToIm.im * foj1.val.im);            
          }
        }
      }
    }

    //r.ksi[i] = (foj2.val * foj1.ksi[i]) + (foj1.val * foj2.ksi[i]);
  }
  FirstOrderJet<DerivativeT, D>::buffer->val = foj1.val * foj2.val;
  FirstOrderJet<DerivativeT, D>::buffer->val2 = foj1.val2 * foj2.val2;
//  FirstOrderJet<DerivativeT, D>::buffer->val3 = foj1.val3 * foj2.val3;
  //r.val = foj1.val * foj2.val;
//  *FirstOrderJet<DerivativeT, D>::buffer = foj1;
//  *FirstOrderJet<DerivativeT, D>::buffer *= foj2;
  return *FirstOrderJet<DerivativeT, D>::buffer;
  //return r;
}

template< class DerivativeT, int D>
inline const FirstOrderJet<DerivativeT, D>& operator+(const FirstOrderJet<DerivativeT, D>& foj1, const FirstOrderJet<DerivativeT, D>& foj2){
  //FirstOrderJet<DerivativeT, SubspaceT, D> r;
  (typename FirstOrderJet<DerivativeT, D>::DPDEContainer&)*FirstOrderJet<DerivativeT, D>::buffer = 
    (typename FirstOrderJet<DerivativeT, D>::DPDEContainer&)foj1 + (typename FirstOrderJet<DerivativeT, D>::DPDEContainer&)foj2;
  FirstOrderJet<DerivativeT, D>::buffer->val = foj1.val + foj2.val;
  FirstOrderJet<DerivativeT, D>::buffer->val2 = foj1.val2 + foj2.val2;
//  FirstOrderJet<DerivativeT, D>::buffer->val3 = foj1.val3 + foj2.val3;
  int i;
  for(i=0; i < FirstOrderJet<DerivativeT, D>::dim; ++i){
    if(FirstOrderJet<DerivativeT, D>::optimizationControl == global){
      if(FirstOrderJet<DerivativeT, D>::currentDPDEContainer->isRealValuedEven() || FirstOrderJet<DerivativeT, D>::currentDPDEContainer->isOther()){
        if(!FirstOrderJet<DerivativeT, D>::currentDPDEContainer->partialReReZero){          
          FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToRe.re = foj1.ksi[i].partialWRespToRe.re + foj2.ksi[i].partialWRespToRe.re;
        }
        if(!FirstOrderJet<DerivativeT, D>::currentDPDEContainer->partialReImZero){          
          FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToRe.im = foj1.ksi[i].partialWRespToRe.im + foj2.ksi[i].partialWRespToRe.im; 
        }
      }
      if(FirstOrderJet<DerivativeT, D>::currentDPDEContainer->isRealValuedOdd() || FirstOrderJet<DerivativeT, D>::currentDPDEContainer->isOther()){
        if(!FirstOrderJet<DerivativeT, D>::currentDPDEContainer->partialImReZero){
            FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToIm.re = foj1.ksi[i].partialWRespToIm.re + foj2.ksi[i].partialWRespToIm.re;
        }
        if(!FirstOrderJet<DerivativeT, D>::currentDPDEContainer->partialImImZero){
          FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToIm.im = foj1.ksi[i].partialWRespToIm.im + foj2.ksi[i].partialWRespToIm.im; 
        }
      }
    }
//      if(FirstOrderJet<DerivativeT, D>::currentDPDEContainer->isEven() || !FirstOrderJet<DerivativeT, D>::currentDPDEContainer->subspaceType){
//        if(!FirstOrderJet<DerivativeT, D>::currentDPDEContainer->isRealValued())
//          FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToRe = foj1.ksi[i].partialWRespToRe + foj2.ksi[i].partialWRespToRe;
//        else
//          FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToRe.re = foj1.ksi[i].partialWRespToRe.re + foj2.ksi[i].partialWRespToRe.re;
//      }
//      if(FirstOrderJet<DerivativeT, D>::currentDPDEContainer->isOdd() || !FirstOrderJet<DerivativeT, D>::currentDPDEContainer->subspaceType){
//        if(!FirstOrderJet<DerivativeT, D>::currentDPDEContainer->isRealValued())
//          FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToIm = foj1.ksi[i].partialWRespToIm + foj2.ksi[i].partialWRespToIm;
//        else{
//          if(FirstOrderJet<DerivativeT, D>::currentDPDEContainer->isRealValued() == capd::jaco::realValuedL2)
//            FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToIm.re = foj1.ksi[i].partialWRespToIm.re + foj2.ksi[i].partialWRespToIm.re;
//          else
//            FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToIm.im = foj1.ksi[i].partialWRespToIm.im + foj2.ksi[i].partialWRespToIm.im;
//        }
//      }      
    if(FirstOrderJet<DerivativeT, D>::optimizationControl == local){
      if(FirstOrderJet<DerivativeT, D>::buffer->isRealValuedEven() || FirstOrderJet<DerivativeT, D>::buffer->isOther()){
        if(!FirstOrderJet<DerivativeT, D>::buffer->partialReReZero){          
          FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToRe.re = foj1.ksi[i].partialWRespToRe.re + foj2.ksi[i].partialWRespToRe.re;
        }
        if(!FirstOrderJet<DerivativeT, D>::buffer->partialReImZero){          
          FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToRe.im = foj1.ksi[i].partialWRespToRe.im + foj2.ksi[i].partialWRespToRe.im; 
        }
      }
      if(FirstOrderJet<DerivativeT, D>::buffer->isRealValuedOdd() || FirstOrderJet<DerivativeT, D>::buffer->isOther()){
        if(!FirstOrderJet<DerivativeT, D>::buffer->partialImReZero){
            FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToIm.re = foj1.ksi[i].partialWRespToIm.re + foj2.ksi[i].partialWRespToIm.re;
        }
        if(!FirstOrderJet<DerivativeT, D>::buffer->partialImImZero){
          FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToIm.im = foj1.ksi[i].partialWRespToIm.im + foj2.ksi[i].partialWRespToIm.im; 
        }
      }
    }
  }
  return *FirstOrderJet<DerivativeT, D>::buffer;
//  FirstOrderJet<DerivativeT, D> r;
//  r = foj1;
//  r += foj2;
//
//  return r;
}

template< class DerivativeT, int D>
inline const FirstOrderJet<DerivativeT, D>& operator-(const FirstOrderJet<DerivativeT, D>& foj1, const FirstOrderJet<DerivativeT, D>& foj2){
  //FirstOrderJet<DerivativeT, SubspaceT, D> r;
  (typename FirstOrderJet<DerivativeT, D>::DPDEContainer&)*FirstOrderJet<DerivativeT, D>::buffer = 
    (typename FirstOrderJet<DerivativeT, D>::DPDEContainer&)foj1 + (typename FirstOrderJet<DerivativeT, D>::DPDEContainer&)foj2;
  FirstOrderJet<DerivativeT, D>::buffer->val = foj1.val - foj2.val;
  FirstOrderJet<DerivativeT, D>::buffer->val2 = foj1.val2 - foj2.val2;
//  FirstOrderJet<DerivativeT, D>::buffer->val3 = foj1.val3 - foj2.val3;
  int i;
  for(i=0; i < FirstOrderJet<DerivativeT, D>::dim; ++i){
    if(FirstOrderJet<DerivativeT, D>::optimizationControl == global){
      if(FirstOrderJet<DerivativeT, D>::currentDPDEContainer->isRealValuedEven() || FirstOrderJet<DerivativeT, D>::currentDPDEContainer->isOther()){
        if(!FirstOrderJet<DerivativeT, D>::currentDPDEContainer->partialReReZero){          
          FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToRe.re = foj1.ksi[i].partialWRespToRe.re - foj2.ksi[i].partialWRespToRe.re;
        }
        if(!FirstOrderJet<DerivativeT, D>::currentDPDEContainer->partialReImZero){          
          FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToRe.im = foj1.ksi[i].partialWRespToRe.im - foj2.ksi[i].partialWRespToRe.im; 
        }
      }
      if(FirstOrderJet<DerivativeT, D>::currentDPDEContainer->isRealValuedOdd() || FirstOrderJet<DerivativeT, D>::currentDPDEContainer->isOther()){
        if(!FirstOrderJet<DerivativeT, D>::currentDPDEContainer->partialImReZero){
            FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToIm.re = foj1.ksi[i].partialWRespToIm.re - foj2.ksi[i].partialWRespToIm.re;
        }
        if(!FirstOrderJet<DerivativeT, D>::currentDPDEContainer->partialImImZero){
          FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToIm.im = foj1.ksi[i].partialWRespToIm.im - foj2.ksi[i].partialWRespToIm.im; 
        }
      }
//      if(FirstOrderJet<DerivativeT, D>::currentDPDEContainer->isEven() || !FirstOrderJet<DerivativeT, D>::currentDPDEContainer->subspaceType){
//        if(!FirstOrderJet<DerivativeT, D>::currentDPDEContainer->isRealValued())
//          FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToRe = foj1.ksi[i].partialWRespToRe - foj2.ksi[i].partialWRespToRe;
//        else
//          FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToRe.re = foj1.ksi[i].partialWRespToRe.re - foj2.ksi[i].partialWRespToRe.re;
//      }
//      if(FirstOrderJet<DerivativeT, D>::currentDPDEContainer->isOdd() || !FirstOrderJet<DerivativeT, D>::currentDPDEContainer->subspaceType){
//        if(!FirstOrderJet<DerivativeT, D>::currentDPDEContainer->isRealValued())
//          FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToIm = foj1.ksi[i].partialWRespToIm - foj2.ksi[i].partialWRespToIm;
//        else{
//          if(FirstOrderJet<DerivativeT, D>::currentDPDEContainer->isRealValued() == capd::jaco::realValuedL2)
//            FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToIm.re = foj1.ksi[i].partialWRespToIm.re - foj2.ksi[i].partialWRespToIm.re;
//          else
//            FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToIm.im = foj1.ksi[i].partialWRespToIm.im - foj2.ksi[i].partialWRespToIm.im;
//        }
//      }
    }
    if(FirstOrderJet<DerivativeT, D>::optimizationControl == local){
      if(foj1.isRealValuedEven() || foj1.isOther() || foj2.isRealValuedEven() || foj2.isOther()){
        if(!foj1.partialReReZero || !foj2.partialReReZero){          
          FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToRe.re = foj1.ksi[i].partialWRespToRe.re - foj2.ksi[i].partialWRespToRe.re;
        }
        if(!foj1.partialReImZero || !foj2.partialReImZero){          
          FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToRe.im = foj1.ksi[i].partialWRespToRe.im - foj2.ksi[i].partialWRespToRe.im; 
        }
      }
      if(foj1.isRealValuedOdd() || foj1.isOther() || foj2.isRealValuedOdd() || foj2.isOther()){
        if(!foj1.partialImReZero || !foj2.partialImReZero){
            FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToIm.re = foj1.ksi[i].partialWRespToIm.re - foj2.ksi[i].partialWRespToIm.re;
        }
        if(!foj1.partialImImZero || !foj2.partialImImZero){
          FirstOrderJet<DerivativeT, D>::buffer->ksi[i].partialWRespToIm.im = foj1.ksi[i].partialWRespToIm.im - foj2.ksi[i].partialWRespToIm.im; 
        }
      }
    }
  }
  return *FirstOrderJet<DerivativeT, D>::buffer;
//  FirstOrderJet<DerivativeT, D> r;
//  r = foj1;
//  r -= foj2;
//
//  return r;
}

template< class DerivativeT, int D>
inline const FirstOrderJet<DerivativeT, D>& conjugate(const FirstOrderJet<DerivativeT, D>& foj){
////  FirstOrderJet<DerivativeT, SubspaceT, D> r;
////  r.val = conjugate(foj.val);
//  FirstOrderJet<DerivativeT, D>::buffer->val = conjugate(foj.val);
//  FirstOrderJet<DerivativeT, D>::buffer->val2 = conjugate(foj.val2);
////  FirstOrderJet<DerivativeT, D>::buffer->val3 = conjugate(foj.val3);
//  int i;
//  for(i=0; i < FirstOrderJet<DerivativeT, D>::dim; ++i){
//    FirstOrderJet<DerivativeT, D>::buffer->ksi[i]  = conjugate(foj.ksi[i]);
//  }

  *FirstOrderJet<DerivativeT, D>::buffer = foj;
  FirstOrderJet<DerivativeT, D>::buffer->conjugate();

  return *FirstOrderJet<DerivativeT, D>::buffer;
  //return r;
}

template< class DerivativeT, int D>
inline const FirstOrderJet<DerivativeT, D>& operator*(const typename FirstOrderJet<DerivativeT, D>::ScalarType& s, const FirstOrderJet<DerivativeT, D>& foj){
//  //FirstOrderJet<DerivativeT, SubspaceT, D> r;
//  //r.val = s * foj.val;
//  FirstOrderJet<DerivativeT, D>::buffer->val = s * foj.val;
//  FirstOrderJet<DerivativeT, D>::buffer->val2 = s * foj.val2;
////  FirstOrderJet<DerivativeT, D>::buffer->val3 = s * foj.val3;
//  int i;
//  for(i=0; i < FirstOrderJet<DerivativeT, D>::dim; ++i){
//    FirstOrderJet<DerivativeT, D>::buffer->ksi[i] = s * foj.ksi[i];
//    //r.ksi[i] = s * foj.ksi[i];
//  }

  *FirstOrderJet<DerivativeT, D>::buffer = foj;
  *FirstOrderJet<DerivativeT, D>::buffer *= s;

  return *FirstOrderJet<DerivativeT, D>::buffer;
  //return r;
}

//template< class DerivativeT, int D>
//inline const FirstOrderJet<DerivativeT, D>& operator/(const FirstOrderJet<DerivativeT, D>& foj, const typename FirstOrderJet<DerivativeT, D>::ScalarType& s){
//  //FirstOrderJet<DerivativeT, SubspaceT, D> r;
//  FirstOrderJet<DerivativeT, D>::buffer->val = foj.val / s;
//  FirstOrderJet<DerivativeT, D>::buffer->val2 = foj.val2 / s;
////  FirstOrderJet<DerivativeT, D>::buffer->val3 = foj.val3 / s;
//  int i;
//  for(i=0; i < FirstOrderJet<DerivativeT, D>::dim; ++i)
//    FirstOrderJet<DerivativeT, D>::buffer->ksi[i] = foj.ksi[i] / s;
//  return *FirstOrderJet<DerivativeT, D>::buffer;
//}

template< class DerivativeT, int D>
inline const FirstOrderJet<DerivativeT, D>& operator*(double d, const FirstOrderJet<DerivativeT, D>& foj){
//  //FirstOrderJet<DerivativeT, SubspaceT, D> r;
//  //r.val = d * foj.val;
//  FirstOrderJet<DerivativeT, D>::buffer->val = d * foj.val;
//  FirstOrderJet<DerivativeT, D>::buffer->val2 = d * foj.val2;
////  FirstOrderJet<DerivativeT, D>::buffer->val3 = d * foj.val3;
//  int i;
//  for(i=0; i < FirstOrderJet<DerivativeT, D>::dim; ++i)
//    FirstOrderJet<DerivativeT, D>::buffer->ksi[i] = d * foj.ksi[i];
//    //r.ksi[i] = d * foj.ksi[i];

  *FirstOrderJet<DerivativeT, D>::buffer = foj;
  *FirstOrderJet<DerivativeT, D>::buffer *= d;

  return *FirstOrderJet<DerivativeT, D>::buffer;
  //return r;
}

}
}

#endif /* FIRSTORDERJET_H_ */
