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

  static double squareNorm(const IndexType& k){
    return k.squareEuclNorm();
  }

};


template<class IndexT>
class MaximumNorm{
public:
  typedef IndexT IndexType;

  static double squareNorm(const IndexType& k){
    double k0sq = k[0] * k[0],
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
};

}
}

#endif /* NORMS_H_ */
