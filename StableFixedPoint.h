/**
 * 24.06.2015 this is a class dedicated for automatic proving of a stable fixed point, and
 * rigorously estimating the basin of attraction.
 *
 * The code is mainly taken from the code used in the paper (context of the Burgers equation)
 *
 * Existence of globally attracting fixed points of viscous Burgers equation with constant forcing. A computer assisted proof
 *  Topological Methods in Nonlinear Analysis, Vol. 45, 2 (2015), 655–697
 *
 */



#include "dissipative_enclosure.h"

#include "capd/alglib/hsschur.h"
#include "capd/alglib/hessenberg.h"


#define __BOX_FIND_INFLATE_C__ 1.01

namespace capd {
namespace jaco {


///The class represents a trapping region, like in the Definition 13 in the paper.
///A trapping region is enclosing a fixed point m_x0.
///In suitable coordinate frame Q vector field is pointing inwards.
///Coordinates from the set \mathcal{I} are kept in the m_discs table.
///We assume that in block coordinates for each (i)\in\mathcal{I} exists index j such that  (i)=(i_1, i_2)
///and i_1=j, i_2=j+1.
template<class VectorT, class MatrixT, class TailT>
class Box{
public:
  typedef VectorT VectorType;
  typedef MatrixT MatrixType;
  typedef typename VectorT::ScalarType ScalarType;
  typedef TailT TailType;
  int m_n;
  VectorType m_v; ///the trapping region in BLOCK COORDINATES
  VectorType m_x0; ///the approximate fixed point location, the center point of the trapping region
  VectorType m_r; ///radiuses of discs
  capd::vectalg::Vector<bool, 0> m_disc; ///on which coordinates are discs, we assume that in coordinates Q , discs are attached to coordinates next to each other
  MatrixType m_Q; ///coordinate frame
  MatrixType m_Qinv;
  TailType m_tail; ///tail that satisfies all assumptions, and together with finite dimensional box defines trapping region around a fixed point

  Box(int dim) : m_n(dim), m_v(dim), m_x0(dim), m_r(dim), m_disc(dim){}

  Box(int dim, const MatrixType& Q, const MatrixType& Qinv) : m_n(dim), m_v(dim), m_x0(dim), m_r(dim), m_disc(dim), m_Q(Q), m_Qinv(Qinv){}

  Box(int dim, const TailType& tail, const MatrixType& Q, const MatrixType& Qinv) : m_n(dim), m_v(dim), m_x0(dim), m_r(dim), m_disc(dim),
      m_Q(Q), m_Qinv(Qinv), m_tail(tail){}

  Box(const VectorType& v, const VectorType& x0, const MatrixType& Q, const MatrixType& Qinv) : m_n(v.size()), m_x0(x0), m_v(m_n),
      m_r(m_n), m_disc(m_n), m_Q(Q), m_Qinv(Qinv){ setVector(v);}

  Box(const VectorType& v, const VectorType& x0, const TailType& tail, const MatrixType& Q, const MatrixType& Qinv) : m_n(v.size()), m_v(m_n),
      m_x0(x0), m_r(m_n), m_disc(m_n), m_Q(Q), m_Qinv(Qinv), m_tail(tail){ setVector(v);}

  Box(ScalarType s, const VectorType& x0, int dim, const MatrixType& Q, const MatrixType& Qinv) : m_n(dim), m_v(dim, s), m_x0(x0),
      m_r(dim), m_disc(dim), m_Q(Q), m_Qinv(Qinv){}

  Box(const Box& box2) : m_n(box2.m_n), m_v(box2.m_v), m_x0(box2.m_x0), m_r(box2.m_r), m_disc(box2.m_disc),
      m_Q(box2.m_Q), m_Qinv(box2.m_Qinv), m_tail(box2.m_tail){}

  ///After finding good tail with isolation we save it here.
  inline void setTail(const TailType& tail){ m_tail=tail;}

  inline TailType& getTail() { return m_tail;}

  inline bool empty() const;

  ///Wraps canonical vector into the box coordinate frame.
  inline VectorType expressInCoordinates(const VectorType& v) const{ return m_Q*v; }


  ///Checks if the box contains in its interior a box in CANONICAL COORDINATES.
  inline bool contains(const VectorType& v) const;

  ///Checks if the box contains in its interior a box in another coordinates, along with corresponding tail.
  inline bool contains( const Box& box2) const;

  ///Checks if the box contains in its interior a box in canonical coordinates, along with corresponding tail.
  inline bool contains( const Box& box2, const VectorType& fp) const;

  ///Takes intersection with another box.
  inline void intersect( const Box& box2);

  inline bool subset( const Box& box2) const;

  void setRadius(ScalarType& r);

  ///Sets radius of the 2-dim disc at coordinates i and ip1.
  void setRadius(int i, int ip1, typename ScalarType::BoundType r);

  inline ScalarType& radius(int i){return m_r[i];}

  inline void setCoordinateFrame(const MatrixType& Q, const MatrixType& Qinv);

  inline MatrixType& getQ(){ return m_Q;}

  inline MatrixType& getQinv(){ return m_Qinv;}

  ///Checks if i=i_1 or i=i_2 where (i_1, i_2)\in\mathcal{I}.
  inline bool isDisc(int i) const{ return m_disc[i];}

  inline void setDisc(int i, int ip1);

  inline ScalarType& operator[](int i);

  //Set current box in coordinates Q.
  void setVector(const VectorType& v);

  ///Function wraps vector m_x expressed in coordinates Q into coordinates canonical coordinates minus m_x0.
  inline VectorType wrap() const{ return m_Qinv*m_v;}

  ///Function wraps vector m_x expressed in coordinates Q into coordinates canonical coordinates, using the provided x0.
  inline VectorType wrap(const VectorType& x0) const{ return m_Qinv*m_v+x0;}

  ///Function wraps vector m_x + provided x_0 in coordinates Q into canonical coordinates.
  inline VectorType wrapAffine() const{ return m_Qinv*(m_v)+m_x0;}

  inline VectorType r() const{ return m_v-midVector(m_v);}

  inline ScalarType last() const{ return m_v[m_n-1];}

  inline int dimension() const{ return m_n;}

  ///Calculates the Euclidean norm of the box.
  ScalarType euclNorm() const{ return diam(m_v).euclNorm()+m_tail.euclNorm();}

  void inflate(double c);

  ///Inflates the smallest element in the box.
  void inflateMinElement(double c);

  ///Inflates the element at given coordinate in the box.
  void inflateElement(int index, double c);

  ///Returns a box right neighbourhood of the box that has radius c.
  Box getNeighbourhood(double radius);

  ///Sums modes absolute supremum, starting at index i up to infinity.
  ScalarType sum(int i, ScalarType a0, bool re) const;

  VectorType getX0() const{ return m_x0;}

  void setX0(const VectorType& x0){ m_x0=x0;}

  ///Returns block from diagonal of matrix T either one or two dimensional (depending if eigenvalue is real of complex)
  ///at coordinate i, complex=true indicates if this block is complex.
  MatrixType getBlock(const MatrixType& T, int i, bool& complex) const;

  ///Calculates norm of the block.
  ScalarType blockNorm(const MatrixType& block) const;

  ///Calculates \sum_{j=1}^{m}{\sup{\frac{\partial N_{(i)}}{\partial x_{(j)}}(x)}}, where i=m+1
  ScalarType partialSum(ScalarType c, ScalarType a0, int i) const;

  ///Calculates the Lipshitz constant in the trapping region in block infinity norm. Needed to prove uniqueness and stability of
  ///a fixed point.
  template<class FadMapT>
  VectorType calculateL(const FadMapT& map, std::ostream& out) const;

  ///Returns the doubleton representation of the box, see (4) in the paper. It is an interval set
  ///that is larger than original box, due to wrapping effect.
  void doubleton(VectorType& x, MatrixType& C, VectorType& r0) const;

  void printRaw(std::ostream& stream) const;

  void print(std::ostream& out) const;

  ///Used by print2Latex function.
  int calculateExponent(const ScalarType& s)const;

  ///Prints the box in latex format.
  void print2Latex(std::ostream& out);

  ///Returns maximal diameter of an element in the box.
  ScalarType maxDiam() const;

  //OPERATORS
  //returns box in its coordinates Q
  operator VectorType& () { return m_v; }

  ///Sums current box with another box (takes interval hull).
  friend void operator+=(Box& box1, const Box& box2){
    if(!box1.empty()){
      box1.m_r=intervalHull(box1.m_r, box2.m_r);
      box1.m_v=intervalHull(box1.m_v, box2.m_v);
      box1.m_tail += box2.m_tail;
    }else{
      box1.m_r=box2.m_r;
      box1.m_v=box2.m_v;
      box1.m_tail=box2.m_tail;
    }
  }

  friend std::ostream& operator<<(std::ostream& out, const Box& box) // output
  {
    int i;
    out<<"\n{";
    for(i=0; i<box.m_n; ++i){
      out<<(i>0 ? ", \n" : "");
      if(box.m_disc[i]){
        out << "ball_{"<<i<<", "<<(i+1)<<"}, r="<<box.m_r[i];
        i++; continue;
      }else{
        out << box.m_v[i];
      }
    }
    out<<"}\n";
    return out;
  }

};

///==========================================function definitions====================================================

template<class VectorT, class MatrixT, class TailT>
inline bool Box<VectorT, MatrixT, TailT>::empty() const{
  VectorType t(m_n);
  if(m_v==t)
    return true;
  else
    return false;
}

template<class VectorT, class MatrixT, class TailT>
inline bool Box<VectorT, MatrixT, TailT>::contains(const VectorType& v) const{
  VectorType vInCoordinatesQ=expressInCoordinates(v);

  if(capd::vectalg::subset(vInCoordinatesQ, m_v)){
    //verify if on coordinates, where we have isolating discs is OK
    int i;
    for(i=0; i<m_n; i++){
      if(m_disc[i]){
        if(!(sqrt(power(vInCoordinatesQ[i], 2)+power(vInCoordinatesQ[i+1], 2)).subset(m_r[i]))){
          return false;
        }
        i++;
      }
    }
    return true;
  }else{
    return false;
  }
}

template<class VectorT, class MatrixT, class TailT>
inline bool Box<VectorT, MatrixT, TailT>::contains( const Box& box2) const{
  //first, we check if coordinate frames of boxes is the same
  //wraps box2 into canonical coordinates (=box') and then check if (box') in coordinates Q is inside the box
  if(contains(box2.wrap()) && box2.m_tail.subset(m_tail))
    return true;
  return false;
}

template<class VectorT, class MatrixT, class TailT>
inline bool Box<VectorT, MatrixT, TailT>::contains( const Box& box2, const VectorType& fp) const{
  //first, we check if coordinate frames of boxes is the same
  //wraps box2 into canonical coordinates (=box') and then check if (box') in coordinates Q is inside the box
  if(contains(box2.wrap(-fp)) && box2.m_tail.subset(m_tail))
    return true;
  return false;
}

template<class VectorT, class MatrixT, class TailT>
inline void Box<VectorT, MatrixT, TailT>::intersect( const Box& box2){
  m_r=intersection(m_r, box2.m_r);
  m_v=intersection(m_v, box2.m_v);
  m_tail.intersect(box2.m_tail);
}

template<class VectorT, class MatrixT, class TailT>
inline bool Box<VectorT, MatrixT, TailT>::subset( const Box& box2) const{
  if(capd::vectalg::subset(m_v, box2.m_v) && m_tail.subset(box2.m_tail))
    return true;
  else
    return false;
}

template<class VectorT, class MatrixT, class TailT>
void Box<VectorT, MatrixT, TailT>::setRadius(ScalarType& r){
  int i;
  for(i=0; i<m_n; ++i){
    m_r[i]=ScalarType(0, r);
    m_v[i]=ScalarType(-r, r);
  }
}

template<class VectorT, class MatrixT, class TailT>
void Box<VectorT, MatrixT, TailT>::setRadius(int i, int ip1, typename ScalarType::BoundType r){
  m_r[i]=ScalarType(0, r);
  m_v[i]=ScalarType(-r, r);
  m_r[ip1]=ScalarType(0, r);
  m_v[ip1]=ScalarType(-r, r);
}

template<class VectorT, class MatrixT, class TailT>
inline void Box<VectorT, MatrixT, TailT>::setCoordinateFrame(const MatrixType& Q, const MatrixType& Qinv){
  m_Q=Q;
  m_Qinv=Qinv;
}

template<class VectorT, class MatrixT, class TailT>
inline void Box<VectorT, MatrixT, TailT>::setDisc(int i, int ip1){
  m_disc[i]=true;
  m_disc[ip1]=true;
}

template<class VectorT, class MatrixT, class TailT>
inline typename Box<VectorT, MatrixT, TailT>::ScalarType& Box<VectorT, MatrixT, TailT>::operator[](int i){
  if(i<0 || i>=m_n){
    std::cerr<<"Index out of bounds in operator[] function of the class Box.\n";
    throw std::runtime_error("Index out of bounds in operator[] function of the class Box.\n");
  }
  return m_v[i];
}

template<class VectorT, class MatrixT, class TailT>
void Box<VectorT, MatrixT, TailT>::setVector(const VectorType& v){
  int i;
  for(i=0; i<m_n; ++i){
    m_v[i]=v[i];
    m_r[i]=ScalarType(0, rightBound(abs(v[i])));
  }
}

template<class VectorT, class MatrixT, class TailT>
void Box<VectorT, MatrixT, TailT>::inflate(double c){
  int i;
  for(i=0; i<m_n; i++){
    m_v[i]=capd::jaco::inflate<VectorType>(m_v[i], c);
    if(m_disc[i]){
      m_r[i]=m_v[i];
      m_r[i].setLeftBound(0.);
    }
  }
  //m_tail.inflate(c);
}

template<class VectorT, class MatrixT, class TailT>
void Box<VectorT, MatrixT, TailT>::inflateMinElement(double c){
  int i,
      min=0;
  bool disc=false;
  ScalarType minElem=diam(m_v[0]);
  for(i=0; i<m_n; i++){
    if(diam(m_v[i])<minElem){
      min=i;
      minElem=diam(m_v[i]);
    }
    if(m_disc[i])
      disc=true;
    else
      disc=false;
  }
  if(disc){
    m_v[min]=capd::jaco::inflate<VectorType>(m_v[min], c);
    m_v[min+1]=capd::jaco::inflate<VectorType>(m_v[min+1], c);
    m_r[min]=m_v[min];
    m_r[min].setLeftBound(0.);
    m_r[min+1]=m_v[min+1];
    m_r[min+1].setLeftBound(0.);
  }else{
    m_v[min]=capd::jaco::inflate<VectorType>(m_v[min], c);
  }
}

template<class VectorT, class MatrixT, class TailT>
void Box<VectorT, MatrixT, TailT>::inflateElement(int index, double c){
  bool disc=false;
  if(m_disc[index])
    disc=true;
  else
    disc=false;
  if(disc){
    m_v[index]=capd::jaco::inflate<VectorType>(m_v[index], c);
    m_v[index+1]=capd::jaco::inflate<VectorType>(m_v[index+1], c);
    m_r[index]=m_v[index];
    m_r[index].setLeftBound(0.);
    m_r[index+1]=m_v[index+1];
    m_r[index+1].setLeftBound(0.);
  }else{
    m_v[index]=capd::jaco::inflate<VectorType>(m_v[index], c);
  }
}

template<class VectorT, class MatrixT, class TailT>
Box<VectorT, MatrixT, TailT> Box<VectorT, MatrixT, TailT>::getNeighbourhood(double radius){
  Box newBox((*this));
  VectorType r(m_n);
  r=ScalarType(-radius, radius);
  r+=rightVector(m_v);
  newBox.setVector(r);
  return newBox;
}

template<class VectorT, class MatrixT, class TailT>
typename Box<VectorT, MatrixT, TailT>::ScalarType Box<VectorT, MatrixT, TailT>::sum(int i, ScalarType a0, bool re) const{
  VectorType v=wrapAffine();
  if(abs(i)<=m_tail.getm()){
    int start=(i>0 ? i : -i),
        end=m_tail.getm(),
        j;
    ScalarType sum=0;
    if(start==0){
      sum+=rightBound(abs(a0));
      start=1;
    }
    for(j=start; j<=end; j++){
      sum+=rightBound(abs(v[m_tail.mode2array(j, re)]));
    }
    return sum+m_tail.sum(re);
  }else{
    return m_tail.sum(i, re);
  }
}

template<class VectorT, class MatrixT, class TailT>
typename Box<VectorT, MatrixT, TailT>::MatrixType Box<VectorT, MatrixT, TailT>::getBlock(const MatrixType& T, int i, bool& complex) const{
  MatrixType r;
  if(m_disc[i]){
    complex=true;
    r=MatrixType(2,2);
    r[0][0]=T[i][i];
    r[0][1]=T[i][i+1];
    r[1][0]=T[i+1][i];
    r[1][1]=T[i+1][i+1];
  }else{
    complex=false;
    r=MatrixType(1,1);
    r[0][0]=T[i][i];
  }
  return r;
}

template<class VectorT, class MatrixT, class TailT>
typename Box<VectorT, MatrixT, TailT>::ScalarType Box<VectorT, MatrixT, TailT>::blockNorm(const MatrixType& block) const{
  if(block.numberOfRows()==2){
    return capd::max(block[0][0], block[1][1])+0.5*abs(block[0][1]+block[1][0]);
  }
  if(block.numberOfRows()==1){
    return block[0][0];
  }
  std::cerr<<"Function calculating the block norm is not working with blocks of size different than 1 or 2.\nBlock="<<block<<".\n";
  throw std::runtime_error("Function calculating the block norm is not working with blocks of size different than 1 or 2.\n");
}

template<class VectorT, class MatrixT, class TailT>
typename Box<VectorT, MatrixT, TailT>::ScalarType Box<VectorT, MatrixT, TailT>::partialSum(ScalarType c, ScalarType a0, int i) const{
  VectorType v=wrapAffine();
  int l, lbar, j;
  bool re;
  ScalarType termS,
             term1Re, term1Im,
             term2Re, term2Im,
             sum=0;
  for(j=0; j<m_n; j++){
    for(lbar=0; lbar<m_n; lbar++){
      l=m_tail.array2mode(lbar, re);
      if(abs(i+l)<=m_tail.getm()){
        if(i+l==0){
          term2Im=0;
          term2Re=a0;
        }else{
          term2Im=v[m_tail.mode2array(i+l, false)];
          term2Re=v[m_tail.mode2array(i+l, true)];
        }
      }else{
        term2Im=m_tail(i+l, false);
        term2Im=m_tail(i+l, true);
      }
      if(abs(i-l)<=m_tail.getm()){
        if(i-l==0){
          term1Im=0;
          term1Re=a0;
        }else{
          term1Im=v[m_tail.mode2array(i-l, false)];
          term1Re=v[m_tail.mode2array(i-l, true)];
        }
      }else{
        term1Im=m_tail(i-l, false);
        term1Re=m_tail(i-l, true);
      }

      if(re){
        sum+=abs(2*c*i*(term1Im+term2Im)*m_Qinv[lbar][j])+abs(2*c*i*(term1Re+term2Re)*m_Qinv[lbar][j]);
      }else{ //im
        sum+=abs(2*c*i*(term1Re-term2Re)*m_Qinv[lbar][j])+abs(2*c*i*(term1Im-term2Im)*m_Qinv[lbar][j]);
      }
    }
  }
  return sum;
}

template<class VectorT, class MatrixT, class TailT>
template<class FadMapT>
typename Box<VectorT, MatrixT, TailT>::VectorType Box<VectorT, MatrixT, TailT>::calculateL(const FadMapT& map, std::ostream& out) const{
  int i,j,k, kbar;
  bool complex;
  VectorType r(m_n+1);
  VectorType v=wrapAffine();
  MatrixType dFtilde=m_Q*map[v]*m_Qinv;
#if __CALCULATE_L_DEBUG__
  out<<"d\\tilde{F}(W) matrix (the Jacobian matrix calculated on the finite part of the trapping region): \n"<<dFtilde<<"\n";
#endif
  ScalarType c=map.getNCoeff();
  ScalarType a0=map.getA0();
  ScalarType lambda=map.lambda_k(map.getm() + 1);

#if __CALCULATE_L_DEBUG__
   out<<"wrapped affine="<<v<<"\n";
   out<<"dFtilde="<<dFtilde<<"\n";
#endif

  for(i=0; i<m_n; i++){
#if __CALCULATE_L_DEBUG__
    out<<"i="<<i<<"\n\n";
    out<<"block="<<getBlock(dFtilde, i, complex)<<"\n";
    out<<"block norm="<<blockNorm(getBlock(dFtilde, i, complex))<<"\n";
#endif
    r[i]=blockNorm(getBlock(dFtilde, i, complex));
    ///first, we collect the finite number of terms related to the projection coordinates
    if(complex){
      for(j=0; j<m_n; j++){
        if(i!=j && i+1!=j)
          r[i]+=rightBound(sqrt(power(dFtilde[i][j], 2)+power(dFtilde[i+1][j], 2)));
      }
    }else{
      for(j=0; j<m_n; j++){
        if(i!=j)
          r[i]+=rightBound(abs(dFtilde[i][j]));
      }
    }
    ///second, we sum the infinite number of terms related to the coordinates outside the projection
    for(kbar=m_n-1; kbar>=0; kbar--){
      k=m_tail.array2mode(kbar);
#if __CALCULATE_L_DEBUG__
      out<<"kbar="<<kbar<<"\n";
      out<<"k="<<k<<"\n";
      int m=m_tail.getm();
      out<<"Q[i][kbar]="<<abs(m_Q[i][kbar])<<"\n";
      out<<"Im(sum("<<k<<"-"<<m<<"-1))="<<sum(k-m_tail.getm()-1, a0, false)<<"\n";
      out<<"Im(sum("<<k<<"+"<<m<<"+1))="<<sum(k+m_tail.getm()+1, a0, false)<<"\n";
      out<<"Re(sum("<<k<<"-"<<m<<"-1))="<<sum(k-m_tail.getm()-1, a0, true)<<"\n";
      out<<"Re(sum("<<k<<"+"<<m<<"+1))="<<sum(k+m_tail.getm()+1, a0, true)<<"\n";
#endif

      if(complex){ //we use formulas for two dimensional blocks from [ZAKS]
        r[i]+=(abs(m_Q[i][kbar])+abs(m_Q[i+1][kbar]))*2*c*k*(sum(k-m_tail.getm()-1, a0, false)+sum(k+m_tail.getm()+1, a0, false)+
            sum(k-m_tail.getm()-1, a0, true)+sum(k+m_tail.getm()+1, a0, true));
      }else{ //we use formulas for one dimensional blocks from [ZAKS]
        r[i]+=abs(m_Q[i][kbar])*2*c*k*(sum(k-m_tail.getm()-1, a0, false)+sum(k+m_tail.getm()+1, a0, false)+sum(k-m_tail.getm()-1, a0, true)+
            sum(k+m_tail.getm()+1, a0, true));
      }
    }
#if __CALCULATE_L_DEBUG__
    out<<"r["<<i<<"]="<<r[i]<<"\n";
#endif
    if(complex){
      r[i+1]=r[i];
      i++;
    }
  }
  //looking for maximum in m_Qinv matrix
  ScalarType max=0;
  for(i=0; i<m_n; i++)
    for(j=0; j<m_n; j++){
      if(abs(m_Qinv[i][j])>max)
        max=abs(m_Qinv[i][j]);
    }

#if __CALCULATE_L_DEBUG__
  out<<"ibar="<<i<<"\n";
#endif
  i=m_tail.array2mode(i);
#if __CALCULATE_L_DEBUG__
  out<<"partialSum="<<partialSum(c, a0, i)<<"\n";
  out<<"i="<<i<<"\n";
  int m=m_tail.getm();
  out<<"Re(sum("<<i<<"-"<<m<<"))="<<sum(i-m_tail.getm(), a0, true)<<"\n";
  out<<"Im(sum("<<i<<"-"<<m<<"))="<<sum(i-m_tail.getm(), a0, false)<<"\n";
  out<<"Re(sum("<<i<<"+1))="<<sum(i+1, a0, true)<<"\n";
  out<<"Im(sum("<<i<<"+1))="<<sum(i+1, a0, false)<<"\n";
  out<<"Im(sum("<<i<<"-"<<m<<"-2))="<<sum(i-m_tail.getm()-2, a0, false)<<"\n";
  out<<"Im(sum("<<i<<"+"<<m<<"+2))="<<sum(i+m_tail.getm()+2, a0, false)<<"\n";
  out<<"Re(sum("<<i<<"-"<<m<<"-2))="<<sum(i-m_tail.getm()-2, a0, true)<<"\n";
  out<<"Re(sum("<<i<<"+"<<m<<"+2))="<<sum(i+m_tail.getm()+2, a0, true)<<"\n";
#endif

  int m = m_tail.getm();
  r[m_n] = lambda;
  r[m_n] += 2*c*i*m * max * (sum(i-m+1, a0, true) + sum(i-m+1, a0, false) + sum(i+2, a0, true) + sum(i+2, a0, false));
  r[m_n] += 2*c*i*(sum(i-m_tail.getm()-1, a0, false)+sum(i+m_tail.getm()+1, a0, false)+sum(i-m_tail.getm()-1, a0, true)+sum(i+m_tail.getm()+1, a0, true));
#if __CALCULATE_L_DEBUG__
    out<<"r["<<m_n<<"]="<<r[m_n]<<"\n";
#endif
  return r;
}

template<class VectorT, class MatrixT, class TailT>
void Box<VectorT, MatrixT, TailT>::printRaw(std::ostream& stream) const{
  int i;
  stream<<"box:\n";
  for(i=0; i<m_n; ++i){
    stream << m_v[i]<<"\n";
  }
  m_tail.printRaw(stream);
}

template<class VectorT, class MatrixT, class TailT>
void Box<VectorT, MatrixT, TailT>::print(std::ostream& out) const{
  out << "The finite part \\prod_{|k|<=m} a_k:\n";
  m_tail.printModes(m_v, out);
  out << "The tail T=\\prod_{|k|>m} a_k:\n";
  m_tail.print(out);
  out << "\n";
}

template<class VectorT, class MatrixT, class TailT>
int Box<VectorT, MatrixT, TailT>::calculateExponent(const ScalarType& s)const{
  int i=0;
  ScalarType t=s;
  if(t!=0){
    while(abs(t)<1){
      t*=10;
      i--;
    }
  }
  return i;
}

template<class VectorT, class MatrixT, class TailT>
void Box<VectorT, MatrixT, TailT>::doubleton(VectorType& x, MatrixType& C, VectorType& r0) const{
  x=m_x0;
  C=m_Q;
  r0=wrap(); //here is wrapping effect
}

template<class VectorT, class MatrixT, class TailT>
void Box<VectorT, MatrixT, TailT>::print2Latex(std::ostream& out){
   int first=1,
       last=m_tail.lastModeInTail(),
       i, eRe, eIm, eC, eCRe, eCIm;
   ScalarType re, im, cIm, rIm, cRe, rRe;
   out<< "\\begin{array}{|c|c|c|}\\hline\\mathbf{k} & \\mathbf{\\re{a_k}} & \\mathbf{\\im{a_k}}\\\\\\hline\\hline\n";
   for(i=first; i<last; i++){
     if(i<m_tail.array2mode(m_n))
       re=(*this)[m_tail.mode2array(i, 1)];
     else
       re=m_tail(i, 1);
     if(i<m_tail.array2mode(m_n))
       im=(*this)[m_tail.mode2array(i, 0)];
     else
       im=m_tail(i, 0);
     cRe=mid(re); rRe=re-cRe;
     eRe=calculateExponent(rRe);
     cIm=mid(im); rIm=im-cIm;
     eIm=calculateExponent(rIm);
     eCRe=calculateExponent(cRe);
     eCIm=calculateExponent(cIm);
     if(eCRe<-2)
       out<<i<<" & "<<rightBound(cRe*pow(10, -eCRe))<<"\\cdot 10^{"<<eCRe<<"}";
     else
       out<<i<<" & "<<rightBound(cRe);
     if(eRe<-1)
       out<<"+"<<rRe*pow(10, -eRe)<<"10^{"<<eRe<<"} & ";
     else
       out<<"+"<<rRe<<" & ";
     if(eCIm<-2)
       out<<rightBound(cIm*pow(10, -eCIm))<<"\\cdot 10^{"<<eCIm<<"}";
     else
       out<<rightBound(cIm);
     if(eIm<-1)
       out<<"+"<<rIm*pow(10, -eIm)<<"10^{"<<eIm<<"}\\\\\n";
     else
       out<<"+"<<rIm<<"\\\\\n";
   }
   eC=calculateExponent(C(m_tail));
   if(eC<-1)
     out<<"<"<<i<<" & <\\multicolumn{2}{|c|}{"<<rightBound(C(m_tail))*pow(10, -eC)<<"10^{"<<eC<<"}/k^{"<<rightBound(s(m_tail))<<"}}\\\\\\hline\\end{array}\n\n";
   else
     out<<">"<<i<<" & \\multicolumn{2}{|c|}{\\leq"<<rightBound(C(m_tail))<<"/k^{"<<rightBound(s(m_tail))<<"}}\\\\\\hline\\end{array}\n\n";
}

template<class VectorT, class MatrixT, class TailT>
typename Box<VectorT, MatrixT, TailT>::ScalarType Box<VectorT, MatrixT, TailT>::maxDiam() const{
  int i;
  ScalarType maxDiam = 0;
  for(i = 0; i < m_v.size(); ++i)
    if(diam(m_v[i]) > maxDiam) maxDiam = diam(m_v[i]);
  return maxDiam;
}



///Class dedicated for finding a trapping region enclosing a given fixed point.
template <class MultiMapT, class JetMultiMapT, class DoubleT, int D_>
class BoxFinder {
public:
  typedef MultiMapT    MultiMapType;
  typedef JetMultiMapT JetMultiMapType;
  typedef DoubleT DoubleType;
  typedef typename MultiMapType::ScalarType ScalarType;
  typedef typename MultiMapType::VectorType VectorType;
  typedef typename VectorType::iterator iterator;
  typedef typename MultiMapType::MatrixType MatrixType;
  typedef typename MultiMapType::PolyBdType TailType;
  typedef capd::vectalg::Matrix< DoubleType, _D, _D > DoubleMatrixType;
  typedef capd::vectalg::Vector< DoubleType, _D > DoubleVectorType;
  typedef Box<VectorType, MatrixType, TailType> BoxType;

  MultiMapType& m_multiMap;
  JetMultiMapType& m_jetMultiMap;

  TailType m_tailIn;
  TailType m_tailOut;
  TailType m_tailInfiniteIn;
  TailType m_tailInfiniteOut;
  int m_dim;

  ///parameters are: mm-is the map, whose fixed points we are looking for
  ///fp-is the approximate fixed point location (from the nonrigorous integraion or the Newton iterations),
  ///ib-is the initial box enclosing the fixed point,
  ///t-is the tail of our box, we will verify if this tail is in fact isolating.
  BoxFinder(MultiMapT& mm, JetMultiMapT& jetMm ) :
    m_multiMap(mm), m_dim(mm.dimension()), m_tailIn( mm.m ), m_tailOut( mm.m ) ,
    m_tailInfiniteIn( mm.m, mm.M ), m_tailInfiniteOut( mm.m, mm.M ),
    m_jetMultiMap( jetMm ) {
  }

  ///Function finds an isolating tail for a given box. The result is box, an interval set satisfying C1, C2, C3 and C4a of
  ///Definition 2.
  bool findIsolatingTail(VectorType& box, TailType& tail, std::ostream& out) ;

  ///Function takes calculated trapping region and inflates it until largest possible trapping region, satisfying $l<0$ enclosing
  ///an approximate fixed point, is obtained.
  bool inflateTrappingRegion( BoxType& trappingRegion, const VectorType& eigenvaluesRe, const MatrixType& T, DoubleVectorType& fp,
      const VectorType& FatFP, std::ostream& out );

  ///Checks if a given interval set is a trapping region for each Galerkin projection l>m, it assumes that Qbox satisfies
  ///C1, C2, C3 and C4a of Definition 2 (has isolating tail).
  bool vectorFieldOnBoundaryPointingInwards(BoxType& Qbox, const VectorType& eigenvaluesRe, const MatrixType& T, DoubleVectorType& fp,
      const VectorType& FatFP) ;

  ///Function translates complex change of basis into real change of basis, as a result almost upper triangular matrix is obtained.
  DoubleMatrixType translateComplexChOfCoordIntoRealOne(const DoubleMatrixType& rVec, const DoubleMatrixType& iVec) const ;

  ///Method refining candidate for inverse of matrix A, by calculating the Krawczyk operator. Ainv is approximate inverse of A.
  ///It refines only one column of A at a time, x is the columnIndex-th column of [Ainv]. We assume that A is a square matrix.
  bool krawczykRefineColumn( DoubleMatrixType& A, const DoubleMatrixType& Ainv, VectorType& x, int columnIndex) const;

  ///Method refining candidate for inverse of matrix A, Ainv by calculating the Krawczyk operator.
  ///Returns true if Ainv contains the rigorous inverse A^{-1} of A, and the result i.e. K([x], m, C)\cap Ainv is returned in
  ///Ainv.
  ///Returns false if Ainv does not contain A^{-1}, as a consequence of the intersection of the Krawczyk operator K([x],m, C)
  ///and the candidate Ainv being empty.
  bool krawczykRefineMatrix(DoubleMatrixType& A, const DoubleMatrixType& approxAinv, MatrixType& Ainv) const;

  ///Function validating given initial box, finds a box enclosing approximate fixed point location such that the vector field is
  ///pointing inwards.
  ///Parameters are:
  ///fp - approximate fixed point location,
  ///ib - initial box, surrounding given fixed point,
  ///t - the tail
  bool validateBox(DoubleVectorType& fp, DoubleMatrixType& der,
      VectorType& zeroCenteredBox, ///a zero centred box (box around fixed point - fixed point)
      BoxType& Qbox, TailType& tail, DoubleMatrixType& E, VectorType& eigenvaluesRe, MatrixType& T,
      VectorType& FevaluatedAtFP, std::ostream& out);

  bool schurDecompose(const DoubleMatrixType& A, DoubleMatrixType& T, DoubleMatrixType& S) const;

  ///Calculates the Schur decomposition of the matrix A.
  bool schurMatrix(const DoubleMatrixType& A, DoubleMatrixType& T, DoubleMatrixType& S, std::ostream& out) const;


  void printMatrix(const MatrixType& m, std::ostream& out) const;

  ///Used by printMatrix2Latex function.
  int calculateExponent(const ScalarType& s) const;

  ///Prints a matix in latex format.
  void printMatrix2Latex(const MatrixType& m, std::ostream& out) const;

};

///==========================================function definitions====================================================


template <class MultiMapT, class JetMultiMapT, class DoubleT, int D_>
bool BoxFinder<MultiMapT, JetMultiMapT, DoubleT, D_>::schurDecompose(const DoubleMatrixType& A, DoubleMatrixType& T, DoubleMatrixType& S) const{
   if(A.numberOfRows() != A.numberOfColumns())
      throw std::invalid_argument("schurDecompose works only for square matrices");
    ap::real_2d_array a;
    int n = A.numberOfRows();
    a.setbounds(0,n-1, 0, n-1) ;
    int i,j,k;
    for(i =0; i < n; i++)
      for(j=0; j< n; j++){
        a(i,j) = A[i][j];
      }
    ap::real_2d_array t;
    ap::real_2d_array s;
    ap::real_1d_array tau;
    ap::real_2d_array h;
    ap::real_2d_array h1;
    ap::real_2d_array s1;
    ap::real_2d_array s2;
    ap::real_2d_array s3;
    h1.setbounds(1, n, 1, n);
    s1.setbounds(1, n, 1, n);
    s2.setbounds(1, n, 1, n);
    s3.setbounds(1, n, 1, n);

    rmatrixhessenberg(a, n, tau);

    rmatrixhessenbergunpackq(a, n, tau, s);

    rmatrixhessenbergunpackh(a, n, h);

    //changing indexing of h and s from [0...n-1] to [1...n]
    for(i=0; i<n; i++)
      for(j=0; j<n; j++){
        h1(i+1, j+1)=h(i,j);
        s1(i+1, j+1)=s(i,j);
      }

    if(! upperhessenbergschurdecomposition(h1, n, s2))
        throw std::runtime_error("algorithm for schur decomposition did not converge!");

    ///multiplying orthogonal s1 by orthogonal s2
    for(i=1; i<=n; i++)
      for(j=1; j<=n; j++){
        s3(i, j)=0;
        for(k=1; k<=n; k++){
          s3(i, j)+=s1(i, k)*s2(k, j);
        }
      }

    for(int i =0; i < n; i++){
      for(int j=0; j< n; j++){
        T[j][i]=h1(j+1,i+1);
        S[j][i]=s3(j+1,i+1);
      }
    }
    return true;
}


template <class MultiMapT, class JetMultiMapT, class DoubleT, int D_>
bool BoxFinder<MultiMapT, JetMultiMapT, DoubleT, D_>::schurMatrix(const DoubleMatrixType& A, DoubleMatrixType& T, DoubleMatrixType& S, std::ostream& out) const{
  bool r=schurDecompose(A, T, S);
  return r;
}


template <class MultiMapT, class JetMultiMapT, class DoubleT, int D_>
bool BoxFinder<MultiMapT, JetMultiMapT, DoubleT, D_>::findIsolatingTail(VectorType& box, TailType& tail, std::ostream& out){
  int s_ = s(tail);
  TailType tailBak(tail);
  int i, steps = 0, ss;
  TailType n(tail);
  ScalarType diff;
  bool isolation = false;
  bool tailFound = false;

  while(s_ >= m_multiMap.m_sufficientlyLarge && !tailFound) {
#if __BOX_DEBUG__
    out << "Tail finding iteration, s=" << ss << "\n";
    out << "box=" << box << "\n";
    out << "tail=" << tail << "\n";
    out.flush();
#endif
    tailFound = true;
    steps = 0;
    while(!isolation && tailFound && steps++ < 50) {
      isolation = true;

      m_tailInfiniteIn = tail;
      copyFinitePart( box, m_tailInfiniteIn );
#if __BOX_DEBUG__
      out << "m_tailInfiniteIn=" << m_tailInfiniteIn << "\n";
#endif
      m_multiMap.N(m_tailInfiniteIn, n);

#if __BOX_DEBUG__
      out << "Improving isolation step# "<<steps<<"\n current n: " << n << "\n";
      out.flush();
#endif
      for(i = n.finiteTailBegin() ; i <= n.finiteTailEnd(); i++) { //for close tail

        if(leftBound(n[i]) != leftBound(n[i]) || rightBound(n[i]) != rightBound(n[i])) {
          tailFound = false;
          break;
        }
#if __BOX_DEBUG__
        out << "i=" << i << "\n";
        out << "n[" << i << "]=" << n[i] << "\n";
        out << "T[" << i << "]=" << tail[i] << "\n";
        out << "lambda(" << i << ")=" << m_multiMap.lambda_k(tail.array2modeIndex(i)) << "\n";
#endif

        if(!((diff = -tail[i].leftBound() + leftBound((n[i] / -m_multiMap.lambda_k( tail.array2modeIndex(i) )))) > 0)) {
          isolation = false;
          tail.setLeftBound( i, leftBound(tail[i]) + leftBound(diff) );
          tail.set( tail.array2modeIndex(i), capd::jaco::inflate<MultiMapType>(tail[i], __BOX_FIND_INFLATE_C__) );
#if __BOX_DEBUG__
          out << "! left, diff=" << diff << "\n";
          out << "new leftBound=" << leftBound(tail[i]) << "\n";
#endif
        }
        if(!((diff = -tail[i].rightBound() + rightBound((n[i] / -m_multiMap.lambda_k( tail.array2modeIndex(i) )))) < 0)) {
          isolation = false;
          tail.setRightBound( i, rightBound(tail[i]) + rightBound(diff) );
          tail.set( tail.array2modeIndex(i), capd::jaco::inflate<MultiMapType>(tail[i], __BOX_FIND_INFLATE_C__) );
#if __BOX_DEBUG__
          out << "! right, diff=" << diff << "\n";
          out << "new rightBound=" << rightBound(tail[i]) << "\n";
#endif
        }
      }
#if __BOX_DEBUG__
      out << "i from farTail\n";
      out << "array2modeT=" << tail.array2mode(i) << "\n";
      out << "C(T)=" << C(tail) << "\n";
      out << "C(n)=" << C(n) << "\n";
      out << "lambda(" << i << ")=" << m_multiMap.lambda_k(m_multiMap.M + 1) << "\n";
#endif
      //one more check (regarding far tail), because of symmetry we check only one condition
      if(!((diff = C(tail) * m_multiMap.lambda_k(m_multiMap.M + 1) + C(n) * power(ScalarType(m_multiMap.M + 1), s(tail)
          - s(n))) < 0)) {
#if __BOX_DEBUG__
        out << "!, diff=" << diff << "\n";
        out << "new C=" << 1.01 * (C(tail) + rightBound(diff / -m_multiMap.lambda_k(m_multiMap.M + 1))) << "\n";
#endif
        isolation = false;
        setC(tail, 1.01 * (C(tail) + rightBound(diff / -m_multiMap.lambda_k(m_multiMap.M + 1))));
      }
    }
    if(!tailFound) {
      ///tail has not been found successfully, we increase s and start right again
      s_ -= 1;
      tail = tailBak;
      m_multiMap.changeS(tail, (int)s_);
    }
  }
#if __BOX_DEBUG__
  //check if found tail is isolating in fact
  if(tailFound && isolation){

    copyFinitePart(box, tail);

    m_multiMap.N( tail, n );
    for(i = n.finiteTailBegin() ; i <= n.finiteTailEnd(); i++) {
      if(!((diff=(-tail[i].leftBound() + leftBound((n[i] / -m_multiMap.lambda_k(tail.array2modeIndex(i)))))) > 0)) {
        out << "Found tail at coordinate i=" << i << " is not isolating, -tail[i]^- + b[i]^-=" << diff << "\n";
        std::cerr << "Found tail at coordinate i=" << i << " is not isolating, -tail[i]^- + b[i]^-=" << diff << "\n";
        throw std::runtime_error("Found tail is not isolating!\n");
      }
      if(!((diff=(-tail[i].rightBound() + rightBound((n[i] / -m_multiMap.lambda_k(tail.array2modeIndex(i)))))) < 0)) {
        out << "Found tail at coordinate i=" << i << " is not isolating, -tail[i]^+ + b[i]^+=" << diff << "\n";
        std::cerr << "Found tail at coordinate i=" << i << " is not isolating, -tail[i]^+ + b[i]^+=" << diff << "\n";
        throw std::runtime_error("Found tail is not isolating!\n");
      }
    }
  }
  out << "Found tail:" << tail << "\n";
#endif
  return tailFound && isolation;
}

template <class MultiMapT, class JetMultiMapT, class DoubleT, int D_>
bool BoxFinder<MultiMapT, JetMultiMapT, DoubleT, D_>::inflateTrappingRegion(BoxType& trappingRegion, const VectorType& eigenvaluesRe, const MatrixType& T, DoubleVectorType& fp,
    const VectorType& FatFP, std::ostream& out ) {


  TailType tail(trappingRegion.getTail());
  BoxType potentialTR(trappingRegion);
  VectorType box,
             l;
  int iterations = 0;
  double c=1.5;
  //int indexOfLargestMargin=-1; //index at which there is the largest margin, the vector field value at this coordinate is the smallest
  while(true) {
    //potentialTR.inflateMinElement(c);
    if(potentialTR.maxDiam() > 10e-04) //if diameter of the largest element increases we decrease c
      c = 1.1;
    potentialTR.inflate(c);

    box = potentialTR.wrap(VectorType(fp));
    if(!findIsolatingTail(box, tail, out)){
      break;
    }
    else {
      potentialTR.setTail(tail);

      if(!vectorFieldOnBoundaryPointingInwards(potentialTR, eigenvaluesRe, T, fp, FatFP)){ //there is no isolation on current box

        break;
      }
      trappingRegion = potentialTR;
    }
    iterations++;
  }
  if(iterations > 0)
    return true;
  return false;
}

template <class MultiMapT, class JetMultiMapT, class DoubleT, int D_>
bool BoxFinder<MultiMapT, JetMultiMapT, DoubleT, D_>::vectorFieldOnBoundaryPointingInwards(BoxType& Qbox, const VectorType& eigenvaluesRe, const MatrixType& T, DoubleVectorType& fp,
    const VectorType& FatFP){
  bool r = true;
  int i, j;
  ScalarType diff;
  VectorType a(Qbox.dimension()), rest( Qbox.dimension() );
  VectorType box = Qbox.wrap(VectorType(fp)),
             boxMfp = Qbox.wrap();

             m_tailInfiniteIn = box;
             m_tailInfiniteIn.copyTailPartFrom(Qbox.getTail());
             m_multiMap.perturbations(m_tailInfiniteIn, rest);///<Q_mF(x) (the rest which is not affected by change of coordinates and linearization)

             /// !!!!!! 27.06 TODO DOROBIC HESIAN
  VectorType sum = FatFP /*+ 0.5 * m_multiMap.hessian(boxMfp)*/ + rest,
             Qsum = Qbox.getQ() * sum;

  for(i = box.dimension() - 1; i > 0; i--) {  //TODO 11.08.15 DALEM i > 0 zamiast i >= 0 , bo 1 wspolrzedna to a_0

    //check if on the projection of the box on the i-th m_dimension v.f. is pointing inwards, first end
    if(!Qbox.isDisc(i)) { //purely real eigenvalue, coordinate belongs to the cube coordinates
      a[i] = 0;
      for(j = box.dimension() - 1; j >= 0; j--) {
        if(j != i)
          a[i] += T[i][j] * Qbox[j];
      }
      if(!((diff = eigenvaluesRe[i] * rightBound(Qbox[i]) + a[i] + Qsum[i]) < 0)) {
        r = false;
        std::cout << "FAILED ON i=" << i << "\n" ;
      }
      //check if on the projection of the box on the i-th m_dimension v.f. is pointing inwards, second end
      if(!((diff = eigenvaluesRe[i] * leftBound(Qbox[i]) + a[i] + Qsum[i]) > 0)) {
        r = false;
        std::cout << "FAILED ON i=" << i << "\n" ;
      }

    } else { //complex conjugate pair of eigenvalues
      a[i] = 0;
      for(j = box.dimension() - 1; j >= 0; j--) {
        if(j != i && j != i - 1) //we exclude j=i-1 also, because there is a disc corresponding to (i-1,i) coordinates
          a[i] += T[i][j] * Qbox[j];
      }
      //have to additionally calculate a[i-1]
      a[i - 1] = 0;
      for(j = m_multiMap.dimension() - 1; j >= 0; j--) {
        if(j != i - 1 && j != i) //we exclude j=i also, because there is a disc corresponding to (i-1,i) coordinates
          a[i - 1] += T[i - 1][j] * Qbox[j];
      }

      if(!((diff = (eigenvaluesRe[i]) * rightBound(Qbox.radius(i)) + sqrt(power((a[i] + Qsum[i]), 2) + power((a[i - 1] + Qsum[i - 1]), 2))) < 0)) {
        r = false;
      }
      i--;
    }
  }
  return r;
}

template <class MultiMapT, class JetMultiMapT, class DoubleT, int D_>
typename BoxFinder<MultiMapT, JetMultiMapT, DoubleT, D_>::DoubleMatrixType BoxFinder<MultiMapT, JetMultiMapT, DoubleT, D_>::translateComplexChOfCoordIntoRealOne(const DoubleMatrixType& rVec, const DoubleMatrixType& iVec) const {
  //we scan imaginary matrix iVec looking for blocks 2x2 corresponding to i*w and -i*w part from eigenvectors v+i*w, v-i*w
  int i, j;
  DoubleMatrixType r(m_dim, m_dim);
  for(i = 0; i < m_dim - 1; i += 2)
    for(j = 0; j < m_dim; j++) {
      if(iVec[i][j] != 0 || iVec[i + 1][j] != 0) { //block with eigenvectors in form v+iw, v-iw is detected
        //instead of block [v+iw,v-iw] we take block [v,w]
        r[i][j] = rVec[i][j]; //v_1
        r[i + 1][j] = rVec[i + 1][j]; //v_2
        r[i][j + 1] = iVec[i][j]; //w_1
        r[i + 1][j + 1] = iVec[i + 1][j]; //w_2
        j++;
      } else { //we rewrite real eigenvectors
        r[i][j] = rVec[i][j];
        r[i + 1][j] = rVec[i + 1][j];
      }
    }
  //in case of odd dimension we check last row
  if(m_dim % 2 == 1){
    for(j = 0; j < m_dim; j++) {
      if(iVec[i][j] != 0){
        r[i][j + 1] = iVec[i][j];
        r[i][j] = rVec[i][j];
        j++;
      }else{
        r[i][j]=rVec[i][j];
      }
    }
  }
  return r;
}

template <class MultiMapT, class JetMultiMapT, class DoubleT, int D_>
bool BoxFinder<MultiMapT, JetMultiMapT, DoubleT, D_>::krawczykRefineColumn( DoubleMatrixType& A, const DoubleMatrixType& Ainv, VectorType& x, int columnIndex) const {
  int dim = A.numberOfRows();
  VectorType m(dim), r(dim);
  MatrixType der = MatrixType(A); /// F(x)=A*x, dF(x)=A
  DoubleMatrixType C=Ainv; /// C=dF(x)^{-1}=A^{-1}
  MatrixType Id = MatrixType::Identity(dim);
  capd::vectalg::split(x, m, r);
  VectorType fm = MatrixType(A) * m - capd::vectalg::transpose(Id)[columnIndex];
  VectorType K = m - MatrixType(C) * fm  + ( Id - MatrixType(C) * der ) * r;
  if(capd::vectalg::intersection(x, K, x))
    return true;
  else
    return false;
}

template <class MultiMapT, class JetMultiMapT, class DoubleT, int D_>
bool BoxFinder<MultiMapT, JetMultiMapT, DoubleT, D_>::krawczykRefineMatrix(DoubleMatrixType& A, const DoubleMatrixType& approxAinv, MatrixType& Ainv) const {
  int i;
  VectorType x;
  int dim=Ainv.numberOfColumns();
  MatrixType AinvT=capd::vectalg::transpose(Ainv),
             candidate;
  for(i=0; i<dim; i++){
    x=AinvT[i];
    if(!krawczykRefineColumn( A, approxAinv, x, i ))
      return false;
    AinvT[i]=x;
  }
  candidate=capd::vectalg::transpose(AinvT);
  intersection(Ainv, candidate, Ainv);
  return true;
}

template <class MultiMapT, class JetMultiMapT, class DoubleT, int D_>
bool BoxFinder<MultiMapT, JetMultiMapT, DoubleT, D_>::validateBox(DoubleVectorType& fp, DoubleMatrixType& dDer,
    VectorType& zeroCenteredBox, ///a zero centred box (box around fixed point - fixed point)
    BoxType& Qbox, TailType& tail, DoubleMatrixType& E, VectorType& eigenvaluesRe, MatrixType& T,
    VectorType& FevaluatedAtFP, std::ostream& out) {

  DoubleMatrixType bMatrix(m_dim, m_dim), dQ(m_dim, m_dim), dQinv(m_dim, m_dim), dS(m_dim, m_dim), dT(m_dim, m_dim), R(m_dim, m_dim), rVec(m_dim,
      m_dim), //real part of eigenvectors
      iVec(m_dim, m_dim), rVecInv(m_dim, m_dim), dSinv(m_dim, m_dim) //imaginary part of eigenvectors
      , dq(m_dim, m_dim), dr(m_dim, m_dim);
  DoubleVectorType rV(m_dim), iV(m_dim); //nonrigorous eigenvalues

  MatrixType der(m_dim, m_dim), S(m_dim, m_dim), Q(m_dim, m_dim), Qinv(m_dim, m_dim), eigenvectors(m_dim, m_dim), eigenvectorsInv(m_dim, m_dim),
      Sinv(m_dim, m_dim);
  der = MatrixType(dDer);

  //Q and Qinv is our coordinate change
  schurMatrix(dDer, dT, dS, out);
  dS=capd::vectalg::transpose(dS);
  computeEigenvaluesAndEigenvectors(dT, rV, iV, rVec, iVec);
  //we have a complex change of coordinates in two matrices rVec the real part and iVec the imaginary part
  //we unite them into one matrix which generates block diagonal form of T
  rVec = translateComplexChOfCoordIntoRealOne(rVec, -iVec);
  //rVec = capd::vectalg::transpose(rVec);

  int j;

#if __BOX_DEBUG__
  out << "rVec: " << rVec << "\n";
  out << "iVec: " << iVec << "\n";
  out << "Re eigenvalue: "<<rV<<"\n";
  out << "Im eigenvalue: "<<iV<<"\n";
#endif

  eigenvectors = MatrixType(rVec);
  rVecInv = capd::matrixAlgorithms::inverseMatrix(rVec);
  eigenvectorsInv = MatrixType(rVecInv);

  der = m_jetMultiMap.jacobian( VectorType(fp) );

  S = MatrixType(dS);
  //matrices S together with eigenvectors is our change of coordinates, inverses are calculated rigorously

  dSinv = capd::matrixAlgorithms::inverseMatrix(dS);
  Sinv = MatrixType(dSinv);
  std::cout<<"Refining inverse matrices using the Krawczyk operator.\n";
  //inflating matrix Sinv
  Sinv+=ScalarType(-1, 1);
  //refining it using the Krawczyk operator
  if(!krawczykRefineMatrix(dS, dSinv, Sinv)){
    out << "Procedure of rigorously inverting a matrix using the Krawczyk operator failed.\nTry to rerun program with different"
        << " projection size m.\n";
    std::cerr << "Procedure of rigorously inverting a matrix using the Krawczyk operator failed.\nTry to rerun program with different"
              << " projection size m.\n";
    throw std::runtime_error("Procedure of rigorously inverting a matrix using the Krawczyk operator failed.\nTry to rerun program with different projection size m.\n");
  }

  //the same with rVec matrix
  eigenvectorsInv += ScalarType(-1, 1);

  if(!krawczykRefineMatrix(rVec, rVecInv, eigenvectorsInv)){
    out << "Procedure of rigorously inverting a matrix using the Krawczyk operator failed.\nTry to rerun program with different"
        << " projection size m.\n";
    std::cerr << "Procedure of rigorously inverting a matrix using the Krawczyk operator failed.\nTry to rerun program with different"
              << " projection size m.\n";
    throw std::runtime_error("Procedure of rigorously inverting a matrix using the Krawczyk operator failed.\nTry to rerun program with different projection size m.\n");
  }

  Q = eigenvectorsInv * S;

  Qinv = Sinv * eigenvectors;

  std::cout<<"Matrices refined.\n";
  std::cout<<"Looking for a trapping region enclosing the fixed point.\n";
  //almost upper triangular matrix T,
  //with eigenvalues on diagonal, for real eigenvalues, and blocks of size 2x2 for complex conjugate eigenvalues
  T = Q * der * Qinv; //rigorous evaluation of T
  //dT represents approximate T, we use it, because it has zeros under diagonal blocks

#if __BOX_DEBUG__
  out << "eigenvectorsInv=" << eigenvectorsInv << "\n";
  out << "eigenvectors*eigenvectorsInv=" << eigenvectors * eigenvectorsInv << "\n";
  out << "rVec: " << rVec << "\n";
  out << "rVecInv: " << rVecInv << "\n";
  out << "rVec*T*rVecInv: " << rVecInv * dT * rVec << "\n";
  out << "S: " << dS << "\n";
  out << "Q:" << Q << "\n";
  out << "Qinv:" << Qinv << "\n";
  out << "Q*Qinv: " << Q * Qinv << "\n";
  out << "dT:" << dT << "\n";
#endif
#if __PRINT_DATA_IN_LATEX_FORMAT__
  out.precision(3);
  out << "the Jacobian at \\overline{x} dF(\\overline{x}): \n";
  printMatrix2Latex(der, out);
  out << "calculated change of basis:\n";
  out << "A: \n";
  printMatrix2Latex(Q, out);
  out << "Ainv: \n";
  printMatrix2Latex(Qinv, out);
  out << "calculated matrix T=Q*dF(\\ overline{x})*Qinv (the Jacobian in new coordinates)\n";
  out << "Almost diagonal form with blocks 1x1 and 2x2 on the diagonal, blocks 1x1 contain rigorous values of the real "
         "eigenvalues, whereas blocks 2x2 contain rigorous values of the complex eigenvalues: \n";
  printMatrix2Latex(T, out);
#else
  out << "the Jacobian at \\overline{x} dF(\\overline{x}): \n"<< dDer <<"\n";
  out << "calculated change of basis:\n";
  out << "A: \n"<< Q <<"\n";
  out << "Ainv: \n" << Qinv << "\n";
  out << "calculated matrix [D]=A*dF(\\ overline{x})*Ainv (the Jacobian in new coordinates)\n";
  out << "Almost diagonal form with blocks 1x1 and 2x2 on the diagonal, blocks 1x1 contain rigorous values of the real "
         "eigenvalues, whereas blocks 2x2 contain rigorous values of the complex eigenvalues: \n" << T << "\n";
#endif
  out << "approximate eigenvalues (complex numbers): \n" <<  eigenvaluesToString( rV, iV ) << "\n";

  Qbox.setX0(VectorType(fp));
  Qbox.setCoordinateFrame(Q, Qinv); //Qbox is a candidate for a trapping region around a fixed point

  VectorType box(m_dim), boxMfp(m_dim);

  m_tailIn = VectorType(fp) ;

  m_multiMap( m_tailIn, m_tailOut );

  copyFinitePart(m_tailOut, FevaluatedAtFP);

  //m_multiMap.setT0(tail);

  VectorType rest(m_dim), sum(m_dim), Qsum(m_dim), a(m_dim); //sum of elements off diagonal

  bool validated = false, wasUpdated = false;
  capd::vectalg::Vector<bool, 0> updated(m_dim);
  ScalarType diff;
  int steps, lsteps = 0, refinementSteps=0; //steps for counting how many times we could not find a trapping region at given box

  int i;
  //Initialising box, for each coordinate choosing if it belongs to cube or a ball
  for(i = 0; i < m_dim - 1; i++) {
    if(abs(dT[i + 1][i]) > abs(0.01*dT[i][i])) { //if we have at coordinates i, i-1 pair of complex conjugate eigenvalues we inform
                                       //about it the Qbox. We check if the imaginary part is sufficiently large to create a
                                       //disc at coordinate (i, i+1).
      Qbox.setDisc(i, i + 1);
      i++;
      continue;
    }
  }
  //calculating rigorous values of real part of eigenvalues
  for(i = 0; i < m_dim; i++) {
    if(dT[i + 1][i] != 0) {
      eigenvaluesRe[i] = (T[i][i] + T[i + 1][i + 1]) / 2.;
      eigenvaluesRe[i + 1] = eigenvaluesRe[i];
      i++;
    } else {
      eigenvaluesRe[i] = T[i][i];
    }
  }

  TailType tailBak = tail;
  //local trapping region finding algorithm, starts at some large box enclosing a fixed point,
  //if it does not manage to find trapping region in 100 steps it takes a larger box.

  while(!validated && lsteps++ < 100) {
//      std::cout<<"iteration #"<<lsteps<<"\n";
    steps = 0;
    ///setting box
    Qbox.setVector(zeroCenteredBox);
    boxMfp = Qbox.wrap(); //this is box in canonical coordinates minus fixed point
    box = Qbox.wrapAffine(); //this is box around fixed point in canonical coordinates
#if __BOX_DEBUG__
    out << "Step #" << lsteps << "\n";
    out << "Trying to find trapping region in a box.\n";
    out << "Current box: " << box << "\n";
    out << "Current Qbox: " << Qbox << "\n";
    out.flush();
#endif

    if(findIsolatingTail(box, tail, out)) {
      Qbox.setTail(tail);
      while(!validated && steps < 10) {

        validated = true;
        wasUpdated = false;
#if __BOX_DEBUG__
        out << "step #" << steps << ".\n";
#endif
        boxMfp = Qbox.wrap(); //this is a box in canonical coordinates minus fixed point
        box = Qbox.wrapAffine(); //this is a box around fixed point in canonical coordinates


        m_tailIn = VectorType(fp) ;
          m_multiMap( m_tailIn, m_tailOut );
          copyFinitePart(m_tailOut, FevaluatedAtFP);

        m_tailInfiniteIn = box;
        m_tailInfiniteIn.copyTailPartFrom(tail);
        m_multiMap.perturbations(m_tailInfiniteIn, rest);

        /// !!!!!! 27.06 TODO DOROBIC HESIAN
        sum = FevaluatedAtFP /*+ 0.5 * m_multiMap.hessian(boxMfp)*/ + rest;
        Qsum = Q * sum;
#if __BOX_DEBUG__
        out << "m_tailIn=" << m_tailIn << "\n";
        out << "FatFP=" << FevaluatedAtFP << "\n";
        out << "box-fp=" << boxMfp << "\n";
        out << "box=" << box << "\n";
        out << "sum=" << sum << "\n";
        out << "rest=" << rest << "\n";
        out << "Qsum=" << Qsum << "\n";
        out << "Qbox=" << Qbox << "\n";
#endif
        for(i = m_dim - 1; i > 0; i--) { //TODO 11.08.15 DALEM i > 0 zamiast i >= 0 , bo 1 wspolrzedna to a_0
          //check if on the projection of the box on the i-th m_dimension v.f. is pointing inwards, the first end
          if(eigenvaluesRe[i] > 0) {
            std::cerr << "Error. One of the eigenvalues is positive, the fixed point is probably not stable.\n";
            throw std::runtime_error("Error. One of the eigenvalues is positive, the fixed point is probably not stable.\n");
          }

          if(!Qbox.isDisc(i)) { //purely real eigenvalue, coordinate belongs to the cube coordinates
            a[i] = 0;
            for(j = m_dim - 1; j >= 0; j--) {
              if(j != i)
                a[i] += T[i][j] * Qbox[j];
            }
#if __BOX_DEBUG__
            out << "Qbox[" << i << "]=" << Qbox[i] << "\n";
            out << "Qsum[" << i << "]=" << Qsum[i] << "\n";
            out << "a[" << i << "]=" << a[i] << "\n";
            out << "\neigenvalues[" << i << "]*rightBound(Qbox[" << i << "])=" << eigenvaluesRe[i] * rightBound(Qbox[i]) << "\n";
#endif
            if(!((diff = eigenvaluesRe[i] * rightBound(Qbox[i]) + a[i] + Qsum[i]) < 0)) {
#if __BOX_DEBUG__
              out << "!eigenvalues[" << i << "]*rightBound(Qbox[" << i << "])+a[" << i << "]+Qsum[" << i << "]<0\n";
              out << "diff=" << diff << "\n";
              out << "validated false\n\n";
#endif
              validated = false;
            } else {
#if __BOX_DEBUG__
              out << "decreasing Qbox[" << i << "].rightBound=" << rightBound(Qbox[i]) << "\n";
              out << "new rightBound=" << rightBound((a[i] + Qsum[i]) / -eigenvaluesRe[i]) << "\n\n";
#endif
              Qbox[i].setRightBound(rightBound((a[i] + Qsum[i]) / -eigenvaluesRe[i]));
              updated[i] = true;
              wasUpdated = true;
            }
#if __BOX_DEBUG__
            out << "\neigenvalues[" << i << "]*leftBound(Qbox[" << i << "])=" << eigenvaluesRe[i] * leftBound(Qbox[i]) << "\n";
            out << "Qsum[" << i << "]=" << Qsum[i] << "\n";
            out << "a[" << i << "]=" << a[i] << "\n";
#endif
            //check if on the projection of the box on the i-th m_dimension v.f. is pointing inside, second end
            if(!((diff = eigenvaluesRe[i] * leftBound(Qbox[i]) + a[i] + Qsum[i]) > 0)) {
#if __BOX_DEBUG__
              out << "!eigenvalues[" << i << "]*leftBound(Qbox[" << i << "])+a[" << i << "]+Qsum[" << i << "]>0\n";
              out << "diff=" << diff << "\n\n";
              out << "validated false\n\n";
#endif
              validated = false;
            } else {
#if __BOX_DEBUG__
              out << "decreasing Qbox[" << i << "].leftBound=" << leftBound(Qbox[i]) << "\n";
              out << "new leftBound=" << leftBound((a[i] + Qsum[i]) / -eigenvaluesRe[i]) << "\n\n";
#endif
              Qbox[i].setLeftBound(leftBound((a[i] + Qsum[i]) / -eigenvaluesRe[i]));
              updated[i] = true;
              wasUpdated = true;
            }
          } else { //complex conjugate pair of eigenvalues
            a[i] = 0;

            for(j = m_dim - 1; j >= 0; j--) {
              if(j != i && j != i - 1) //we exclude j=i-1 also, because there is a disc corresponding to (i-1,i) coordinates
                a[i] += T[i][j] * Qbox[j];
            }
            //have to calculate a[i-1]
            a[i - 1] = 0;
            for(j = m_dim - 1; j >= 0; j--) {
              if(j != i - 1 && j != i) //we exclude j=i-1 also, because there is a disc corresponding to (i-1,i) coordinates
                a[i - 1] += T[i - 1][j] * Qbox[j];
            }
#if __BOX_DEBUG__
            out << "r=" << Qbox[i] << "\n";
            out << "a[i]=" << a[i] << ", a[i-1]=" << a[i-1] << "\n";
            out << "alpha("<<eigenvaluesRe[i]<<")*r=" << (eigenvaluesRe[i]) * rightBound(Qbox.radius(i)) << " rest=" << (sqrt(power((a[i] + Qsum[i]), 2) + power((a[i - 1]
                + Qsum[i - 1]), 2))) << "\n";
#endif
            if(!((diff = eigenvaluesRe[i] * rightBound(Qbox.radius(i)) + sqrt(power((a[i] + Qsum[i]), 2) + power((a[i - 1] + Qsum[i - 1]), 2))) < 0)) {
#if __BOX_DEBUG__
              out << "!, diff=" << diff << "\n\n";
#endif
              validated = false;
            } else {
#if __BOX_DEBUG__
              out << "decreasing ball at i=" << i << "\n";
              out << "new r=" << rightBound(sqrt(power((a[i] + Qsum[i]), 2) + power((a[i - 1] + Qsum[i - 1]), 2)) / -(eigenvaluesRe[i])) << "\n\n";
#endif
              Qbox.setRadius(i - 1, i, rightBound(sqrt(power((a[i] + Qsum[i]), 2) + power((a[i - 1] + Qsum[i - 1]), 2)) / -eigenvaluesRe[i]));
              updated[i] = true;
              wasUpdated = true;
              updated[i - 1] = true;
            }
            i--;
          }
        }
        if(wasUpdated) { //calculating new tail and recalculating rest, in case a mode was updated
#if __BOX_DEBUG__
          out << "\n a mode was updated, recalculating tail: \n";
#endif
          tail = tailBak;
          findIsolatingTail(box, tail, out);
          Qbox.setTail(tail);
          boxMfp = Qbox.wrap(); //this is box in canonical coordinates minus fixed point
          box = Qbox.wrap(VectorType(fp)); //this is box around fixed point in canonical coordinates


          m_tailInfiniteIn = box;
          m_tailInfiniteIn.copyTailPartFrom(tail);
          m_multiMap.perturbations(m_tailInfiniteIn, rest);


          /// !!!!!! 27.06 TODO DOROBIC HESIAN
          sum = FevaluatedAtFP /*+ 0.5 * m_multiMap.hessian(boxMfp)*/ + rest;
          Qsum = Q * sum;
        }
        steps++;
      }
    }
    if(!validated) {
      tail = tailBak;
      double factor = (zeroCenteredBox.euclNorm().rightBound() > 0.01 ? 1.1 : 2);
      zeroCenteredBox = factor * zeroCenteredBox; //make a box potentially enclosing fixed point larger
    }
  }


std::cout << "lsteps=" << lsteps << "\n";

//TWO refinement steps
for(refinementSteps=0; refinementSteps<1; refinementSteps++){ // TODO TURNED OFF

  validated = true;
  wasUpdated = false;
#if __BOX_DEBUG__
  out << "Refinement step #"<<refinementSteps<<".\n";
#endif
  boxMfp = Qbox.wrap(); //this is a box in canonical coordinates minus fixed point
  box = Qbox.wrapAffine(); //this is a box around fixed point in canonical coordinates


  m_tailInfiniteIn = box;
  m_tailInfiniteIn.copyTailPartFrom(tail);
  m_multiMap.perturbations(m_tailInfiniteIn, rest);


  /// !!!!!! 27.06 TODO DOROBIC HESIAN
  sum = FevaluatedAtFP /*+ 0.5 * m_multiMap.hessian(boxMfp)*/ + rest;
  Qsum = Q * sum;
#if __BOX_DEBUG__
  out << "T=" << T << "\n";
  out << "FatFP=" << FevaluatedAtFP << "\n";
  out << "box-fp=" << boxMfp << "\n";
  out << "box=" << box << "\n";
  out << "sum=" << sum << "\n";
  out << "rest=" << rest << "\n";
  out << "Qsum=" << Qsum << "\n";
  out << "Qbox=" << Qbox << "\n";
#endif
  for(i = m_dim - 1; i > 0; i--) {   //TODO 11.08.15 DALEM i > 0 zamiast i >= 0 , bo 1 wspolrzedna to a_0
    //check if on the projection of the box on the i-th m_dimension v.f. is pointing inwards, the first end
    if( eigenvaluesRe[i] > 0 ) {
      std::cerr << "Error. One of the eigenvalues is positive, the fixed point is probably not stable.\n";
      throw std::runtime_error("Error. One of the eigenvalues is positive, the fixed point is probably not stable.\n");
    }

    if(!Qbox.isDisc(i)) { //purely real eigenvalue, coordinate belongs to the cube coordinates
      a[i] = 0;

      for(j = m_dim - 1; j >= 0; j--) {
        if(j != i)
          a[i] += T[i][j] * Qbox[j];
      }
#if __BOX_DEBUG__
      out << "Qbox[" << i << "]=" << Qbox[i] << "\n";
      out << "Qsum[" << i << "]=" << Qsum[i] << "\n";
      out << "a[" << i << "]=" << a[i] << "\n";
      out << "\neigenvalues[" << i << "]*rightBound(Qbox[" << i << "])=" << eigenvaluesRe[i] * rightBound(Qbox[i]) << "\n";
#endif
      if(!((diff = eigenvaluesRe[i] * rightBound(Qbox[i]) + a[i] + Qsum[i]) < 0)) {
#if __BOX_DEBUG__
        out << "!eigenvalues[" << i << "]*rightBound(Qbox[" << i << "])+a[" << i << "]+Qsum[" << i << "]<0\n";
        out << "diff=" << diff << "\n\n";
#endif
        validated = false;
      } else {
#if __BOX_DEBUG__
        out << "decreasing Qbox[" << i << "].rightBound=" << rightBound(Qbox[i]) << "\n";
        out << "new rightBound=" << rightBound((a[i] + Qsum[i]) / -eigenvaluesRe[i]) << "\n\n";
#endif
        Qbox[i].setRightBound(rightBound((a[i] + Qsum[i]) / -eigenvaluesRe[i]));
        updated[i] = true;
        wasUpdated = true;
      }
#if __BOX_DEBUG__
      out << "\neigenvalues[" << i << "]*leftBound(Qbox[" << i << "])=" << eigenvaluesRe[i] * leftBound(Qbox[i]) << "\n";
      out << "Qsum[" << i << "]=" << Qsum[i] << "\n";
      out << "a[" << i << "]=" << a[i] << "\n";
#endif
      //check if on the projection of the box on the i-th m_dimension v.f. is pointing inside, second end
      if(!((diff = eigenvaluesRe[i] * leftBound(Qbox[i]) + a[i] + Qsum[i]) > 0)) {
#if __BOX_DEBUG__
        out << "!eigenvalues[" << i << "]*leftBound(Qbox[" << i << "])+a[" << i << "]+Qsum[" << i << "]>0\n";
        out << "diff=" << diff << "\n\n";
#endif
        validated = false;
      } else {
#if __BOX_DEBUG__
        out << "decreasing Qbox[" << i << "].leftBound=" << leftBound(Qbox[i]) << "\n";
        out << "new leftBound=" << leftBound((a[i] + Qsum[i]) / -eigenvaluesRe[i]) << "\n\n";
#endif
        Qbox[i].setLeftBound(leftBound((a[i] + Qsum[i]) / -eigenvaluesRe[i]));
        updated[i] = true;
        wasUpdated = true;
      }
    } else { //complex conjugate pair of eigenvalues
      a[i] = 0;
      for(j = m_dim - 1; j >= 0; j--) {
        if(j != i && j != i-1) //we exclude j=i-1 also, because there is a disc corresponding to (i-1,i) coordinates
          a[i] += T[i][j] * Qbox[j];
      }
      //have to calculate a[i-1]
      a[i - 1] = 0;
      for(j = m_dim - 1; j >= 0; j--) {
        if(j != i - 1 && j != i) //we exclude j=i-1 also, because there is a disc corresponding to (i-1,i) coordinates
          a[i - 1] += T[i - 1][j] * Qbox[j];
      }
#if __BOX_DEBUG__
      out << "r=" << Qbox[i] << "\n";
      out << "a[i]=" << a[i] << ", a[i-1]=" << a[i-1] << "\n";
      out << "alpha("<<eigenvaluesRe[i]<<")*r=" << (eigenvaluesRe[i]) * rightBound(Qbox.radius(i)) << " rest=" << (sqrt(power((a[i] + Qsum[i]), 2) + power((a[i - 1]
          + Qsum[i - 1]), 2))) << "\n";
#endif
      if(!((diff = eigenvaluesRe[i] * rightBound(Qbox.radius(i)) + sqrt(power((a[i] + Qsum[i]), 2) + power((a[i - 1] + Qsum[i - 1]), 2))) < 0)) {
#if __BOX_DEBUG__
        out << "!, diff=" << diff << "\n\n";
#endif
        validated = false;
      } else {
#if __BOX_DEBUG__
        out << "decreasing ball at i=" << i << "\n";
        out << "new r=" << rightBound(sqrt(power((a[i] + Qsum[i]), 2) + power((a[i - 1] + Qsum[i - 1]), 2)) / -(eigenvaluesRe[i])) << "\n\n";
#endif
        Qbox.setRadius(i - 1, i, rightBound(sqrt(power((a[i] + Qsum[i]), 2) + power((a[i - 1] + Qsum[i - 1]), 2)) / -eigenvaluesRe[i]));
        updated[i] = true;
        wasUpdated = true;
        updated[i - 1] = true;
      }
      i--;
    }
  }
  if(wasUpdated) { //calculating new tail and recalculating rest, in case a mode was updated
#if __BOX_DEBUG__
    out << " recalculating tail: \n";
#endif
    tail = tailBak;
    findIsolatingTail(box, tail, out);
    Qbox.setTail(tail);
    boxMfp = Qbox.wrap(); //this is box in canonical coordinates minus fixed point
    box = Qbox.wrap(VectorType(fp)); //this is box around fixed point in canonical coordinates

    m_tailInfiniteIn = box;
    m_tailInfiniteIn.copyTailPartFrom(tail);
    m_multiMap.perturbations(m_tailInfiniteIn, rest);

    /// !!!!!! 27.06 TODO DOROBIC HESIAN
    sum = FevaluatedAtFP /*+ 0.5 * m_multiMap.hessian(boxMfp)*/ + rest;
    Qsum = Q * sum;
  }
}

  if(!validated){
    out << "updated=" << updated << "\n";
    out << "A box enclosing a fixed point was not found. Probably the provided Galerkin projection dimension is too small," <<
           "or there is a complex eigenvalue \\alpha+i \\beta, with \\alpha being too small\n";
    std::cerr<<"A box enclosing a fixed point was not found. Probably the provided Galerkin projection dimension is too small," <<
        "or there is a complex eigenvalue \\alpha+i \\beta, with \\alpha being too small\n";;
    throw std::runtime_error("A box enclosing a fixed point was not found. Probably the provided Galerkin projection dimension is too small or there is a complex eigenvalue \\alpha+i \\beta the difference between \\alpha and \\beta is too small.\n");
  }
  //checking if found box is in fact isolating
  if(!vectorFieldOnBoundaryPointingInwards(Qbox, eigenvaluesRe, T, fp, FevaluatedAtFP)){
    out << "Vector field on the found potential trapping region (Qbox in validateBox function) is not isolating.\n";
    std::cerr << "Vector field on the found potential trapping region (Qbox in validateBox function) is not isolating.\n";
    throw std::runtime_error("Vector field on the found potential trapping region (Qbox in validateBox function) is not isolating.\n");
  }

#if __BOX_DEBUG__
  out << "\nFinished.\n";
  out << "box-fp=" << boxMfp << "\n";
  out << "box=" << box << "\n";
  out << "Qbox=" << Qbox << "\n";
  out << "tail: " << tail << "\n";
#endif
  return validated;
}

template <class MultiMapT, class JetMultiMapT, class DoubleT, int D_>
void BoxFinder<MultiMapT, JetMultiMapT, DoubleT, D_>::printMatrix(const MatrixType& m, std::ostream& out) const{
  MatrixType c=m,
             r=m;
  capd::vectalg::split(c, r);
  int i,j;
  out<<"{";
  for(i=0; i<m.numberOfRows(); i++){
    for(j=0; j<m.numberOfColumns(); j++){
      if(j>0) out<<", ";
      out<<rightBound(m[i][j])<<"+"<<r[i][j];
    }
    out<<";\n";
  }
  out<<"}\n";
}

template <class MultiMapT, class JetMultiMapT, class DoubleT, int D_>
int BoxFinder<MultiMapT, JetMultiMapT, DoubleT, D_>::calculateExponent(const ScalarType& s)const{
  int i=0;
  ScalarType t=s;
  if(t!=0){
    while(abs(t)<1){
      t*=10;
      i--;
    }
  }
  return i;
}

template <class MultiMapT, class JetMultiMapT, class DoubleT, int D_>
void BoxFinder<MultiMapT, JetMultiMapT, DoubleT, D_>::printMatrix2Latex(const MatrixType& m, std::ostream& out) const{
  MatrixType c(m);
  MatrixType r(m);
  DoubleType maxBound;
  int i,j,eM;
  for(i=0; i<m.numberOfRows(); i++){
    for(j=0; j<m.numberOfColumns(); j++){
      c[i][j]=mid(m[i][j]);
      r[i][j]=m[i][j]-c[i][j];
    }
  }
  out<<"m=\n\\left[\\begin{array}{@{\\,}c @{\\,} c @{\\,} c @{\\,} c @{\\,} c @{\\,} c@{\\,}}";
  for(i=0; i<m.numberOfRows(); i++){
    for(j=0; j<m.numberOfColumns(); j++){
      if(j>0) out<<" & ";
      eM=calculateExponent(rightBound(m[i][j]));
      if(eM<-2)
        out<<rightBound(c[i][j]) * power(10, -eM)<<"\\cdot 10^{"<<eM<<"}";
      else
        out<<rightBound(c[i][j]);
    }
    out<<"\\\\\n";
  }
  out<<"\\end{array}\\right].\n";
  int a, b, max;
  out<<"r=\n[-1,1]\\cdot\\left[\\begin{array}{@{\\,}c @{\\,} c @{\\,} c @{\\,} c @{\\,} c @{\\,} c@{\\,}}";
  for(i=0; i<m.numberOfRows(); i++){
    for(j=0; j<m.numberOfColumns(); j++){
      if(j>0) out<<" & ";
      a=calculateExponent(leftBound(r[i][j]));
      b=calculateExponent(rightBound(r[i][j]));

      max=-min(a, b);
      maxBound=abs(leftBound(r[i][j]))>rightBound(r[i][j]) ? abs(leftBound(r[i][j])) : rightBound(r[i][j]);
      if(max!=0)
        out<<maxBound*power(10, max)<<"\\cdot 10^{-"<<max<<"}";
      else
        out<<maxBound;
    }
    out<<"\\\\\n";
  }
  out<<"\\end{array}\\right].\n\n";
}


}
}
