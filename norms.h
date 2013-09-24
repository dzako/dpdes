/*
 * norms.h
 *
 *  Created on: Oct 25, 2011
 *      Author: cyranka
 */

#ifndef NORMS_H_
#define NORMS_H_

/**This file contains some classes defining norms, in a sense providing a method calculating norm of an index that is used by
 * the Subspace classes (Real and so on).
 */

namespace capd{
namespace jaco{

template<class IndexT>
class EuclideanNorm{
public:
  typedef IndexT IndexType;

  template <class ScalarType>
  static ScalarType norm(const IndexType& k){
    return sqrt(k.squareEuclNorm());
  }

  static int squareNorm(const IndexType& k){
    return k.squareEuclNorm();
  }

  template <class ScalarType>
  static ScalarType harmonic1DBound(const ScalarType& start, int s){
    return (1. / (s - 1.)) * power(1. / start, s - 1.);
  }

  template <class ScalarType>
  static ScalarType harmonic2DBound(const ScalarType& start, int s){
    ScalarType sqrt2 = ScalarType(1.414213562373093, 1.414213562373095),
               cd = 2 * ScalarType::pi(); //a constant depending on the dimension
    return power((1 + sqrt2 /(2 * start)), s) * (cd / (s-2)) * power((start - sqrt2 / 2.), -(s-2)); //this line is hard-coded for 2D
  }

};


template<class IndexT>
class MaximumNorm{
public:
  typedef IndexT IndexType;

  static int norm(const IndexType& k){
    int max = 0;
    int t;
    if((t = abs(k[0])) > max)
      max = t;
    if((t = abs(k[1])) > max)
      max = t;
    if((t = abs(k[2])) > max)
      max = t;
    return max;
  }

  static int squareNorm(const IndexType& k){
    int k0sq = k[0] * k[0],
           k1sq = k[1] * k[1],
           k2sq = k[2] * k[2],
           max=0;
    if(k0sq > max)
      max = k0sq;
    if(k1sq > max)
      max = k1sq;
    if(k2sq > max)
      max = k2sq;
    return max;
  }



  template <class ScalarType>
  static ScalarType harmonic1DBound(const ScalarType& start, int s){
    return (1. / (s - 1.)) * power(1. / start, s - 1.);
  }

  template <class ScalarType>
  static ScalarType harmonic2DBound(const ScalarType& start, int s){
    if(s <= 2){
      std::cerr << "harmonic2DBound in MaximumNorm class can not be called with s <= 2\n";
      throw std::runtime_error("harmonic2DBound in MaximumNorm class can not be called with s <= 2\n");
    }
    return (4 * 2) / ((s - 2) * power((start - 1), s - 2)) ;
  }

};

}
}

#endif /* NORMS_H_ */
