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

  const double squareNorm(const IndexType& k) const{
    return k.squareEuclNorm();
  }

};


template<class IndexT>
class MaximumNorm{
public:
  typedef IndexT IndexType;

  const double squareNorm(const IndexType& k) const{
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
