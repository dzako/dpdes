/*
 * Index.h
 *
 *  Created on: Aug 8, 2011
 *      Author: cyranka
 */

#ifndef INDEX_H_
#define INDEX_H_

namespace capd{
namespace jaco{

enum Inequality{
  weak, strong
};

template<typename IndexT>
class IndexRange{
public:
  IndexRange(){}


  bool isWithin(double what, int leftBound, Inequality left, int rightBound, Inequality right) const{
    if( what<leftBound || what>rightBound )
      return false;
    if( (left==strong && what==leftBound) || (right==strong && what==rightBound) )
      return false;
    return true;
  }

  virtual int moduloM1Index() const=0;///<each coordinate of the index is taken modulo (this value + 1)
  virtual int returnIndex() const=0;///<when modulo value is reached at a coordinate of the index, this coordinate is reseted to this value
  virtual bool withinRange(const IndexT& i) const=0;
};


///15.08.2011
///Some kind of ITERATOR.
///Class used for generating indices required to calculate the convolution sum \sum{a_k a_{k-k_1}}.
///Responsible for the part of the sum. That is, returns {k_1} such that the norms |k_1| and |k-k_1| are within required range.

template<class IndexT, class NormT>
class FourierConvolutionIndexRange : public IndexRange<IndexT>, NormT{
public:
  typedef IndexT IndexType;
  typedef NormT NormType;
  IndexType k;
  int k_1NormLeft; ///<assumes that the bounds of the norms are given with integers, NEGATIVE NORM MEANS INFINITY
  int k_1NormRight;
  int kmk_1NormLeft;
  int kmk_1NormRight;
  Inequality k_1IneqLeft;
  Inequality k_1IneqRight;
  Inequality kmk_1IneqLeft;
  Inequality kmk_1IneqRight;
  bool k_1orKmk_1;//withinRange returns true if k_1 or k-k_1 is withinRange

  using IndexRange<IndexT>::isWithin;
  using NormType::squareNorm;

  FourierConvolutionIndexRange(const IndexType& k_ = IndexType::zero()) : k(k_), k_1orKmk_1(false){}
  void setRange(int NL, Inequality L, int NR, Inequality R){
    if(NL > NR && NR > 0){
      std::cerr << "Range is wrong NL("<<NL<<") > NR("<<NR<<")\n";
      throw std::runtime_error("Range is wrong NL > NR\n");
    }
    k_1NormLeft=NL;
    k_1NormRight=NR;
    k_1IneqLeft=L;
    k_1IneqRight=R;
    kmk_1NormLeft=NL;
    kmk_1NormRight=NR;
    kmk_1IneqLeft=L;
    kmk_1IneqRight=R;
    k_1orKmk_1=false;
  }
  void setK_1Range(int k_1NL, Inequality k_1L, int k_1NR, Inequality k_1R){
    if(k_1NL > k_1NR && k_1NR > 0){
      std::cerr << "Range is wrong k_1NL("<<k_1NL<<") > k_1NR("<<k_1NR<<")\n";
      throw std::runtime_error("Range is wrong k_1NL > k_1NR\n");
    }
    k_1NormLeft=k_1NL;
    k_1NormRight=k_1NR;
    k_1IneqLeft=k_1L;
    k_1IneqRight=k_1R;
    k_1orKmk_1=false;
  }
  void setKmk_1Range(int kmk_1NL, Inequality kmk_1L, int kmk_1NR, Inequality kmk_1R){
    if(kmk_1NL > kmk_1NR && kmk_1NR > 0){
      std::cerr << "Range is wrong kmk_1NL("<<kmk_1NL<<") > kmk_1NR("<<kmk_1NR<<")\n";
      throw std::runtime_error("Range is wrong kmk_1NL > kmk_1NR\n");
    }
    kmk_1NormLeft=kmk_1NL;
    kmk_1NormRight=kmk_1NR;
    kmk_1IneqLeft=kmk_1L;
    kmk_1IneqRight=kmk_1R;
    k_1orKmk_1=false;
  }

  void setK_1orKmK_1Range(int NL, Inequality L, int NR, Inequality R){
    if(NR < 0){
      std::cerr << "FourierConvolutionIndexRange.setK_1orKmK_1Range() Wrong range\nRange of k_1 or k-k_1 within (...infty) results" <<
          "in infinite loops.\n";
      throw std::runtime_error("FourierConvolutionIndexRange.setK_1orKmK_1Range() Wrong range\nRange of k_1 or k-k_1 within (...infty) results in infinite loops.\n");
    }
    if(NL > NR){
      std::cerr << "Range is wrong NL("<<NL<<") > NR("<<NR<<")\n";
      throw std::runtime_error("Range is wrong NL > NR\n");
    }
    k_1NormLeft=NL;
    k_1NormRight=NR;
    k_1IneqLeft=L;
    k_1IneqRight=R;
    kmk_1NormLeft=NL;
    kmk_1NormRight=NR;
    kmk_1IneqLeft=L;
    kmk_1IneqRight=R;
    k_1orKmk_1=true;
  }

  int moduloM1Index() const{
    return k_1NormRight;
  }
  int returnIndex() const{
    return -k_1NormRight;
  }

  bool k_1UptoInfinity() const{
    if(k_1NormRight<0)
      return true;
    return false;
  }

  bool kmk_1UptoInfinity() const{
    if(kmk_1NormRight<0)
      return true;
    return false;
  }

  bool withinRange(const IndexType& k_1) const{
    bool r=true;
    double k_1SquareNorm = squareNorm(k_1),
               kmk_1SquareNorm = squareNorm(k-k_1);

    if(k_1orKmk_1){
      if( isWithin(kmk_1SquareNorm, kmk_1NormLeft*kmk_1NormLeft, kmk_1IneqLeft, kmk_1NormRight*kmk_1NormRight, kmk_1IneqRight) ||
        isWithin(k_1SquareNorm, k_1NormLeft*k_1NormLeft, k_1IneqLeft, k_1NormRight*k_1NormRight, k_1IneqRight) )
        return true;
      return false;
    }else{
      if(k_1UptoInfinity()){
        if(k_1SquareNorm < k_1NormLeft*k_1NormLeft)
          r = false;
        if(k_1IneqLeft == strong && k_1SquareNorm == k_1NormLeft*k_1NormLeft)
          r = false;
      }else{
        if(!isWithin(k_1SquareNorm, k_1NormLeft*k_1NormLeft, k_1IneqLeft, k_1NormRight*k_1NormRight, k_1IneqRight))
          r = false;
      }
      if(kmk_1UptoInfinity()){
        if(kmk_1SquareNorm < kmk_1NormLeft*kmk_1NormLeft)
          r = false;
        if(kmk_1IneqLeft == strong && kmk_1SquareNorm == kmk_1NormLeft*kmk_1NormLeft)
          r = false;
      }else{
        if(!isWithin(kmk_1SquareNorm, kmk_1NormLeft*kmk_1NormLeft, kmk_1IneqLeft, kmk_1NormRight*kmk_1NormRight, kmk_1IneqRight))
          r = false;
      }
    }
    return r;
  }
};

enum IncludeZero{withoutZero, withZero};

class Index{
public:
  double k[3];///<coordinates are doubles, because we want to compare integer indices with double indices, a result of division
              ///of a index which is not an integer itself.
  double l;    
  
  Index() : l(0){
    k[0]=0; k[1]=0; k[2]=0; }
  inline Index(const Index& i2) : l(i2.l){
    k[0]=i2.k[0]; k[1]=i2.k[1]; k[2]=i2.k[2]; }
  inline bool isZero() const{
      if(k[0]==0 && k[1]==0 && k[2]==0)
        return true;
      else
        return false;
    }

  void cast(){
    //02.11.2011 conceptional mistake, casting negative should return first integer smaller or equal to the value that is being casted
    if(k[0] != static_cast<int>(k[0]))
      k[0] = (k[0] >= 0 ? static_cast<int>(k[0]) : static_cast<int>(k[0]-1));
    if(k[1] != static_cast<int>(k[1]))
      k[1] = (k[1] >= 0 ? static_cast<int>(k[1]) : static_cast<int>(k[1]-1));
    if(k[2] != static_cast<int>(k[2]))
      k[2] = (k[2] >= 0 ? static_cast<int>(k[2]) : static_cast<int>(k[2]-1));
  }

  bool integer() const{
    if(k[0]!=static_cast<int>(k[0]))
      return false;
    if(k[1]!=static_cast<int>(k[1]))
      return false;
    if(k[2]!=static_cast<int>(k[2]))
      return false;
    return true;
  }

  virtual int d() const=0;

  ///Number of components stored in the index
  virtual int components() const=0; 
  
  inline double squareEuclNorm() const{
    double r=k[0]*k[0];
    if(d()>1) r+=k[1]*k[1];
    if(d()>2) r+=k[2]*k[2];
    return r;
  }

  Index& operator=(const Index& i2){
    l=i2.l;
    k[0]=i2.k[0];
    k[1]=i2.k[1];
    k[2]=i2.k[2];
    return *this;
  }

  virtual double& operator[](int l){
    return (k[l]);
  }

  virtual const double& operator[](int l) const{
    return k[l];
  }

  bool divisibleBy2() const{
    bool f=true;
    if(static_cast<int>((*this)[0])%2 != 0)
      f=false;
    if(d() > 1 && static_cast<int>((*this)[1])%2 != 0)
      f=false;
    if(d() > 2 && static_cast<int>((*this)[2])%2 != 0)
      f=false;
    return f;
  }

  bool isDivisibleBy2() const{
    return divisibleBy2();
  }

  ///returns a norm of the index which sets the order
  virtual int orderNorm() const=0;

};


class Index1D : public Index{

public:
  inline Index1D() : Index(){}
  inline Index1D(const Index1D& i2) : Index(i2){}

  explicit inline Index1D(int i) : Index(){
    k[0] = i;
  }

  inline virtual int d() const{return 1;}

  inline virtual int components() const{ return 1;}
  
  using Index::operator[];

  ///increases the current index by one, and returns true if current index is within range ir.
  template<class IndexRangeT>
  inline bool inc(const IndexRangeT& ir, bool unconditional = false){
    cast();
    do{
      k[0]++;
    }while(!ir.withinRange(*this) && k[0] <= ir.k_1NormRight);
    if(ir.withinRange(*this)){
      return true;
    }else{
      return false;
    }
  }

  ///if index last component is nonnegative
  inline bool upperHalfspace() const{
    if(k[0]>=0)
      return true;
    return false;
  }

  inline static Index1D firstInUpperHalfspace(){
    Index1D index;
    index[0] = 1;
    return index;
  }

  template<class IndexRangeT>
  inline bool limitReached(const IndexRangeT& ir){
    if(ir.k_1orKmk_1)
      return abs(k[0]) > ir.k_1NormRight && abs(ir.k[0]-k[0]) > ir.kmk_1NormRight; //both have to be out of range because k_1orKmk_1 is true
    return (k[0] > ir.k_1NormRight && !ir.k_1UptoInfinity()) || (ir.k[0]-k[0]) < -ir.kmk_1NormRight;
  }

  inline static Index1D zero(){
    Index1D i1D;
    i1D[0]=0;
    return i1D;
  }

  inline Index1D operator/(double divisor) const{
    Index1D r(*this);
    r.k[0]=k[0]/divisor;
    return r;
  }

  inline Index1D operator-() const{
    Index1D r(*this);
    r.k[0]=-r.k[0];
    return r;
  }

  inline int orderNorm() const{
    return k[0];
  }

  /**returns
    *\f[
    *\sum_{k\in{range}}{\frac{1}{|k_1|^s}}
    *\f]
  */
  template <class ScalarType, class IndexRangeT>
  static ScalarType harmonicSumK_1(const IndexRangeT& range, int s){
    if(!range.k_1UptoInfinity()){
      std::cerr << "forbidden call of harmonicSumK_1 with range which is not upto infinity.\n";
      throw std::runtime_error("forbidden call of harmonicSumK_1 with range which is not upto infinity.\n");
    }
    ScalarType start = (range.k_1IneqLeft == strong ? range.k_1NormLeft : range.k_1NormLeft - 1);
    return (1. / (s - 1.)) * power(1. / ScalarType(start), s - 1.);
  }

  /**returns
    *\f[
    *\sum_{k\in{range}}{\frac{1}{|k-k_1|^s}}
    *\f]
  */
  template <class ScalarType, class IndexRangeT>
  static ScalarType harmonicSumKmk_1(const IndexRangeT& range, int s){
    if(!range.kmk_1UptoInfinity()){
      std::cerr << "forbidden call of harmonicSumKmk_1 with range which is not upto infinity.\n";
      throw std::runtime_error("forbidden call of harmonicSumKmk_1 with range which is not upto infinity.\n");
    }
    ScalarType start = (range.kmk_1IneqLeft == strong ? range.kmk_1NormLeft : range.kmk_1NormLeft - 1);
    return (1. / (s - 1.)) * power(1. / ScalarType(start), s - 1.);
  }

  inline int mode2array(int m, bool re, int includeZero = withZero) const{
    int t = (includeZero == withoutZero ? -1 : 0);
    if(re)
      return 2 * (k[0] + t);
    return 2 * (k[0] + t) + 1;
  }

  inline static Index1D array2modeIndex(int m, int i, int includeZero = withZero){
    Index1D r;
    i /= 2;
    int t = (includeZero == withoutZero ? 1 : 0);
    r[0] = i + t;
    return r;
  }

  inline static Index1D array2modeIndex(int m, int i, bool& re, int includeZero = withZero){
    Index1D r;
    if(i % 2 == 0)
      re = true;
    else
      re = false;
    i /= 2;
    int t = (includeZero == withoutZero ? 1 : 0);
    r[0] = i + t;
    return r;
  }

  friend std::ostream& operator<<(std::ostream& out, const Index1D& i) // output
  {
    out<<i.k[0];
    return out;
  }

};

inline Index1D operator-(const Index1D& i1, const Index1D& i2){
  Index1D r(i1);
  r.k[0]=i1.k[0]-i2.k[0];
  return r;
}

inline Index1D operator+(const Index1D& i1, int i){
  Index1D r(i1);
  r.k[0]=i1.k[0]+i;
  return r;
}

//inline bool operator<=(const Index1D& i, int k){
//  return i[0]<=k;
//}

inline bool operator<=(const Index1D& i1, const Index1D& i2){
  return i1[0]<=i2[0];
}

inline bool operator<(const Index1D& i1, const Index1D& i2){
  return i1[0]<i2[0];
}

inline bool operator>=(const Index1D& i1, const Index1D& i2){
  return i1[0]>=i2[0];
}

inline bool operator>(const Index1D& i1, const Index1D& i2){
  return i1[0]>i2[0];
}


/**2D Index for one-component functions.
 * 
 */
class Index2D : public Index{
public:
  Index2D() : Index(){}
  Index2D(int k0, int k1){
    k[0] = k0;
    k[1] = k1;
  }
  inline Index2D(const Index2D& i2) : Index(i2){}
  inline virtual int d() const{return 2;}

  inline virtual int components() const{ return 1;}
  
  template<class IndexRangeT>
  inline bool inc(const IndexRangeT& ir, bool unconditional = false){
    Index2D t(*this);
    cast();
    do{
      if(k[0]<ir.moduloM1Index())
        k[0]++;
      else{
        k[1]++;
        k[0]=ir.returnIndex();
      }
      //the last inequality used when noninteger index is increased
    }while((!ir.withinRange(*this) && (k[1] <= ir.k_1NormRight)) || (*this <= t));    
    
    if(ir.withinRange(*this)){
      return true;
    }else{
      return false;
    }
  }

  ///TODO: verify this, 31.10.11 probably was error here not (ir.k[1]-k[1])>ir.kmk_1NormRight but (ir.k[1]-k[1])<-ir.kmk_1NormRight
  template<class IndexRangeT>
  inline bool limitReached(const IndexRangeT& ir){
    if(ir.k_1orKmk_1)
      return k[1] > ir.k_1NormRight && (ir.k[1] - k[1]) < -ir.kmk_1NormRight; //both have to be out of range because k_1orKmk_1 is true
    return k[1] > ir.k_1NormRight || (ir.k[1] - k[1]) < -ir.kmk_1NormRight;
  }

  inline bool operator<=(Index2D& i2) const{
    if((*this)[1] == i2[1])
      return (*this)[0] <= i2[0];
    return (*this)[1] <= i2[1];
  }

  inline Index2D operator/(double divisor) const{
    Index2D r(*this);
    r.k[0]=k[0]/divisor;
    r.k[1]=k[1]/divisor;
    return r;
  }

  inline Index2D operator-() const{
    Index2D r(*this);
    r.k[0]=-r.k[0];
    r.k[1]=-r.k[1];
    return r;
  }

  static Index2D zero(){
    Index2D i2D;
    i2D[0]=0;
    i2D[1]=0;
    return i2D;
  }

  inline int orderNorm() const{
    return 0;//k[0]+max+(2*max+1)*k[1];
  }

  inline int mode2array(int m, bool re, int includeZero = withZero) const{
    int t = (includeZero == withoutZero ? -1 : 0);
    
    int r;
    if(k[1] == 0)
      r = k[0] + t;
    else
      r = k[1] * (2 * m + 1) + k[0] + t; //this is mysterious, but works OK
      
    r *= 2;
    
    if(!re)
      r++;
    return r;    
  }

  /**Function reverse to mode2array , i.e. array index i is provided, and Index2D is calculated
   * 
   */
  inline static Index2D array2modeIndex(int m, int i, int includeZero = withZero){
    int t = (includeZero == withoutZero ? -1 : 0);
        
    Index2D r;
    if(i <= 2 * (m + t) + 1){
      r[1] = 0;
      r[0] = i / 2 + t;
    }else{
      int s = i;
      s /= 2;
      s -= t;
      int u = s / (2 * m + 1);
      int p =  s - u * (2 * m + 1);
      if(p <= m)
        r[0] = p;
      else
        r[0] = s - (u + 1) * (2 * m + 1);
      
      r[1] = (s - r[0]) / (2 * m + 1); 
    }
    return r;
  }

  inline static Index2D array2modeIndex(int m, int i, bool& re, int includeZero = withZero){
    if(i % 2 == 0)
      re = true;
    else
      re = false;
    return array2modeIndex(m, i, includeZero);
  }

  friend std::ostream& operator<<(std::ostream& out, const Index2D& i) // output
  {
    out<<"("<<i.k[0]<<", "<<i.k[1]<<"),"<<i.l;
    return out;
  }

  /**That is how we define the upper halfspace, i.e. { (k_1, k_2) | k_2 > 0 or k_2 = 0 and k_1 >= 0 }
   *
   * The upper-halfspace has to start at the index (0,0), because this mode has to be stored explicitly
   */
  inline bool upperHalfspace() const{
    if(k[1] > 0)
      return true;
    if(k[1] == 0 && k[0] >= 0)
      return true;
    return false;
  }

  inline Index2D firstInUpperHalfspace(int m){
    Index2D index;
    index[1] = 0;
    index[0] = -m;
    return index;
  }

};

inline Index2D operator-(const Index2D& i1, const Index2D& i2){
  Index2D r(i1);
  r.l = i2.l; //this is important for the fadbad algorithms to work
  r[0]=i1[0]-i2[0];
  r[1]=i1[1]-i2[1];
  return r;
}


//inline bool operator<=(const Index2D& i, int k){
//  return i[1]<=k;
//}


inline bool operator<=(const Index2D& i1, const Index2D& i2){
  if(i1[1]==i2[1])
    return i1[0]<=i2[0];
  return i1[1]<=i2[1];
}

inline bool operator<(const Index2D& i1, const Index2D& i2){
  if(i1[1]==i2[1])
    return i1[0]<i2[0];
  return i1[1]<i2[1];
}

inline bool operator>=(const Index2D& i1, const Index2D& i2){
  if(i1[1]==i2[1])
    return i1[0]>=i2[0];
  return i1[1]>=i2[1];
}

inline bool operator>(const Index2D& i1, const Index2D& i2){
  if(i1[1]==i2[1])
    return i1[0]>i2[0];
  return i1[1]>i2[1];
}









class Index2DTwoComponents : public Index{
public:
  Index2DTwoComponents() : Index(){}
  inline Index2DTwoComponents(const Index2DTwoComponents& i2) : Index(i2){}
  inline virtual int d() const{return 2;}
  
  inline virtual int components() const{ return 2;}

  /**Switches l (component of a mode) during increasing
   *
   * @param unconditional if true then the l of the index is not being changed, otherwise it changes cyclically
   */
  template<class IndexRangeT>
  inline bool inc(const IndexRangeT& ir, bool unconditional = false){
    Index2DTwoComponents t(*this);
    cast();
    if(l == 0 && !unconditional)
      l++;
    else{
      if(!unconditional) l=0;
      do{
        if(k[0]<ir.moduloM1Index())
          k[0]++;
        else{
          k[1]++;
          k[0]=ir.returnIndex();
        }
        //the last inequality used when noninteger index is increased
      }while((!ir.withinRange(*this) && (k[1] <= ir.k_1NormRight)) || (*this <= t));
    }
    if(ir.withinRange(*this)){
      return true;
    }else{
      return false;
    }
  }

  ///TODO: verify this, 31.10.11 probably was error here not (ir.k[1]-k[1])>ir.kmk_1NormRight but (ir.k[1]-k[1])<-ir.kmk_1NormRight
  template<class IndexRangeT>
  inline bool limitReached(const IndexRangeT& ir){
    if(ir.k_1orKmk_1)
      return k[1]>ir.k_1NormRight && (ir.k[1]-k[1])<-ir.kmk_1NormRight; //both have to be out of range because k_1orKmk_1 is true
    return k[1]>ir.k_1NormRight || (ir.k[1]-k[1])<-ir.kmk_1NormRight;
  }

  inline bool operator<=(Index2DTwoComponents& i2) const{
    if((*this)[1]==i2[1])
      return (*this)[0]<=i2[0];
    return (*this)[1]<=i2[1];
  }

  inline Index2DTwoComponents operator/(double divisor) const{
    Index2DTwoComponents r(*this);
    r.k[0]=k[0]/divisor;
    r.k[1]=k[1]/divisor;
    return r;
  }

  inline Index2DTwoComponents operator-() const{
    Index2DTwoComponents r(*this);
    r.k[0]=-r.k[0];
    r.k[1]=-r.k[1];
    return r;
  }

  static Index2DTwoComponents zero(){
    Index2DTwoComponents i2D;
    i2D[0]=0;
    i2D[1]=0;
    return i2D;
  }

  inline int orderNorm() const{
    return 0;//k[0]+max+(2*max+1)*k[1];
  }

  inline int mode2array(int m, bool re, int includeZero = withZero) const{
    if(re)
      return 2*(k[0]+m+(2*m+1)*k[1]+(2*m*m+3*m+1)*l); //index of the first mode stored is (i_0,i_1)=(-m,0) (upper halfspace is for i_1>=0)
    return 2*(k[0]+m+(2*m+1)*k[1]+(2*m*m+3*m+1)*l)+1;
  }

  inline static Index2DTwoComponents array2modeIndex(int m, int i, int includeZero = withZero){
    Index2DTwoComponents r;
    i /= 2;
    r.l = i / (2*m*m+3*m+1);
    i = i % (2*m*m+3*m+1);
    r[1] = i / (2*m+1);
    i = i % (2*m+1);
    r[0] = i - m;
    return r;
  }

  inline static Index2DTwoComponents array2modeIndex(int m, int i, bool& re, int includeZero = withZero){
    if(i % 2 == 0)
      re = true;
    else
      re = false;
    Index2DTwoComponents r;
    i /= 2;
    r.l = i / (2*m*m+3*m+1);
    i = i % (2*m*m+3*m+1);
    r[1] = i / (2*m+1);
    i = i % (2*m+1);
    r[0] = i - m;
    return r;
  }

  friend std::ostream& operator<<(std::ostream& out, const Index2DTwoComponents& i) // output
  {
    out<<"("<<i.k[0]<<", "<<i.k[1]<<"),"<<i.l;
    return out;
  }

  /**That is how we define the upper halfspace, i.e. { (k_1, k_2) | k_2 > 0 or k_2 = 0 and k_1 >= 0 }
   *
   * The upper-halfspace has to start at the index (0,0), because this mode has to be stored explicitly
   */
  inline bool upperHalfspace() const{
    if(k[1] > 0)
      return true;
    if(k[1] == 0 && k[0] >= 0)
      return true;
    return false;
  }

  inline Index2DTwoComponents firstInUpperHalfspace(int m){
    Index2DTwoComponents index;
    index[1] = 0;
    index[0] = -m;
    return index;
  }

};

inline Index2DTwoComponents operator-(const Index2DTwoComponents& i1, const Index2DTwoComponents& i2){
  Index2DTwoComponents r(i1);
  r.l = i2.l; //this is important for the fadbad algorithms to work
  r[0]=i1[0]-i2[0];
  r[1]=i1[1]-i2[1];
  return r;
}

inline bool operator<=(const Index2DTwoComponents& i1, const Index2DTwoComponents& i2){
  if(i1[1]==i2[1])
    return i1[0]<=i2[0];
  return i1[1]<=i2[1];
}

inline bool operator<(const Index2DTwoComponents& i1, const Index2DTwoComponents& i2){
  if(i1[1]==i2[1])
    return i1[0]<i2[0];
  return i1[1]<i2[1];
}

inline bool operator>=(const Index2DTwoComponents& i1, const Index2DTwoComponents& i2){
  if(i1[1]==i2[1])
    return i1[0]>=i2[0];
  return i1[1]>=i2[1];
}

inline bool operator>(const Index2DTwoComponents& i1, const Index2DTwoComponents& i2){
  if(i1[1]==i2[1])
    return i1[0]>i2[0];
  return i1[1]>i2[1];
}























//========================================================================================================================


///TODO: in construction
class Index3D : public Index{
public:
  Index3D() : Index(){}
  inline Index3D(const Index3D& i2) : Index(i2){}
  inline virtual int d() const{return 3;}

  inline virtual int components() const{ return 3;}
  
  inline Index3D operator/(double divisor) const{
    Index3D r(*this);
    r.k[0]=k[0]/divisor;
    r.k[1]=k[1]/divisor;
    r.k[2]=k[2]/divisor;
    return r;
  }
  inline Index3D operator-() const{
    Index3D r(*this);
    r.k[0]=-r.k[0];
    r.k[1]=-r.k[1];
    r.k[2]=-r.k[2];
    return r;
  }

  inline static Index3D zero(){
    Index3D i3D;
    i3D[0]=0;
    i3D[1]=0;
    i3D[2]=0;
    return i3D;
  }

  inline int orderNorm() const{
    return 0;//k[0]+max+(2*max+1)*k[1];
  }

  //TODO: hasn't been tested
  inline int mode2array(int m, bool re) const{
    if(re)
      return 2*(k[0]+m+2*(m+1)*(k[1]+m)+(4*m*m+6*m+2)*k[2]+(4*m*m*m+10*m*m+7*m+2)*l);
    return 2*(k[0]+m+2*(m+1)*(k[1]+m)+(4*m*m+6*m+2)*k[2]+(4*m*m*m+10*m*m+7*m+2)*l)-1;
  }

  //TODO: not done
  inline static Index3D firstInUpperHalfspace(){
    Index3D index;
    return index;
  }

  friend std::ostream& operator<<(std::ostream& out, const Index3D& i) // output
  {
    out<<"("<<i.k[0]<<", "<<i.k[1]<<", "<<i.k[2]<<"), "<<i.l;
    return out;
  }

  //TODO: not implemented
  inline bool storedExplicitly() const{
    return true;
  }
};

inline Index3D operator-(const Index3D& i1, const Index3D& i2){
  Index3D r(i1);
  r.k[0]=i1.k[0]-i2.k[0];
  r.k[1]=i1.k[1]-i2.k[1];
  r.k[2]=i1.k[2]-i2.k[2];
  return r;
}

inline Index3D operator+(const Index3D& i1, int i){
  Index3D r(i1);
  r.k[0]=i1.k[0]+i;
  return r;
}

inline bool operator<=(const Index3D& i, int k){
  return i.k[2]<=k;
}


}
}

#endif /* INDEX_H_ */
