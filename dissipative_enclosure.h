#ifndef _DISSIPATIVEENCLOSURE_H_ 
#define _DISSIPATIVEENCLOSURE_H_ 

#include "config.h"
#include "SetBlowUpException.h"

namespace capd {
namespace jaco {

///finds maximal element in a vector
template <typename VectorT>
typename VectorT::ScalarType max(const VectorT& v);

///inflates given interval z in both directions by c.
template <typename FadMapT>
typename FadMapT::ScalarType inflate(typename FadMapT::ScalarType const & z, typename FadMapT::ScalarType::BoundType const& c);

///inflates given interval z in both directions by c, another template.
template <typename ScalarT>
ScalarT inflate(ScalarT const & z, typename ScalarT::BoundType const& c);

///the function finds an C^0 enclosure for \varphi([0,step],x), using combination of both techniques
///-C^1 rough-enclosure,
///-enclosure algorithm based on isolation designed for dPDEs, see [Z3],
///vField is a vector field, its type is the template parameter, should be derivative of FadMap class.
///Data is logged into out.
template <typename FadMapT>
typename FadMapT::VectorType enclosure(FadMapT & vField, typename FadMapT::MatrixType::RowVectorType const & x,
    typename FadMapT::ScalarType const & step);

///finds an enclosure for the solution \varphi([0,step],x). It is dedicated to the DPDE maps that take into account the tail
///(provided as the input parameter T). It uses OutputStream as data logger.
template <typename FadMapT>
typename FadMapT::VectorType enclosureInputTail(FadMapT & vField, typename FadMapT::MatrixType::RowVectorType const & x,
    typename FadMapT::TailType const & T, typename FadMapT::ScalarType const & step);

///version working with PolyBds
template <typename EquationT>
typename EquationT::PolyBdType enclosure(EquationT & vField, typename EquationT::PolyBdType const & x,
    typename EquationT::RealType const & step);


///=================08.2011 Older versions, left for compability reasons===================


///the function finds an C^0 enclosure for \varphi([0,step],x), using combination of both techniques
///-C^1 rough-enclosure,
///-enclosure algorithm based on isolation designed for dPDEs, see [Z3],
///vField is a vector field, its type is the template parameter, should be derivative of FadMap class.
///Data is logged into out.
template <typename FadMapT>
typename FadMapT::VectorType enclosure(FadMapT & vField, typename FadMapT::MatrixType::RowVectorType const & x,
    typename FadMapT::ScalarType const & step, std::ostream& out);


///like regular enclosure function, but calculates on provided tail, FadMapT should provide functions
///setT0() and getT0().
template <typename FadMapT, typename Tail>
typename FadMapT::VectorType enclosure(FadMapT & vField, const Tail& tail, typename FadMapT::MatrixType::RowVectorType const & x,
    typename FadMapT::ScalarType const & step);

///the function finds enclosure for Jacobian matrix (variational part)
///source- "C^1-Lohner algorithm" by P. Zgliczynski.
///Extract from CAPD library.
template <typename FadMapT, typename NormType>
typename FadMapT::MatrixType jacEnclosure(const FadMapT& vectorField, const typename FadMapT::ScalarType& step,
    const typename FadMapT::VectorType &enc, const NormType &the_norm, typename FadMapT::ScalarType* o_logNormOfDerivative = 0);

///==========================================function definitions====================================================

template <typename VectorT>
typename VectorT::ScalarType max(const VectorT& v){
  int i;
  typename VectorT::ScalarType max=v[0];
  for(i=0; i<v.size(); ++i)
    if(v[i]>max) max=v[i];
  return max;
}

template <typename FadMapT>
typename FadMapT::ScalarType inflate(typename FadMapT::ScalarType const & x, typename FadMapT::ScalarType::BoundType const& c) {
  typedef typename FadMapT::ScalarType::BoundType BoundT;
  BoundT left = x.leftBound();
  BoundT right = x.rightBound();
  BoundT dia = diam(x).leftBound();
  BoundT center = (left + right) / 2;
  BoundT sigma = c * dia / 2;//some small value added to prevent from this function doing nothing because 0
  typename FadMapT::ScalarType r(center - sigma, center + sigma);
  return r;
}

template <typename ScalarT>
ScalarT inflate(ScalarT const & x, typename ScalarT::BoundType const& c) {
  typedef typename ScalarT::BoundType BoundT;
  BoundT left = x.leftBound();
  BoundT right = x.rightBound();
  BoundT dia = diam(x).leftBound();
  BoundT center = (left + right) / 2;
  BoundT sigma = c * dia / 2;//some small value added to prevent from this function doing nothing because 0
  ScalarT r(center - sigma, center + sigma);
  return r;
}

template <typename FadMapT>
typename FadMapT::VectorType enclosure(FadMapT & vField, typename FadMapT::MatrixType::RowVectorType const & x,
    typename FadMapT::ScalarType const & step) {
  typedef typename FadMapT::ScalarType ScalarT;
  typedef typename FadMapT::VectorType VectorT;
  typedef typename ScalarT::BoundType DoubleT;

  int steps = 5, s = 0, dimension = vField.getDimension();    
  
  if(dimension > 5) {
    steps = dimension*2;
  }
  DoubleT c = 1.1;
  DoubleT cd = 0.1;
  VectorT z(dimension), N(dimension), b(dimension), g(dimension), y(dimension), w(dimension), fz(dimension), Nz(dimension), Ftemp(dimension),
      L(dimension);
  std::vector<bool> validated(dimension);
  bool allValidated = false;
  int i, j;
  ScalarT h = ScalarT(0, 1) * step;
  FadMapT& F = vField;
  ScalarT small = ScalarT(-1, 1) * __SMALL__;
  ScalarT sm;
  //initial guess
  VectorT fx = F(x);
  for(i = 0; i < dimension; i++) {
    if(F.lambda(i) < 0) {
      z[i] = x[i];
    } else {
      z[i] = x[i] + h * fx[i];
      z[i] = inflate<FadMapT> (z[i], c);
      //z[i]+=small;//here we increase z by small vector to be sure it has nonempty interior
    }
    z[i] += small;
    validated[i] = false;
  }
  double d_l, d_r;
  ScalarT fzn;
  while(!allValidated && s++ < steps) {
    enclosureDebug<<"step: "<<s<<"\n";
    allValidated = true;
    fz = F(z);
    Nz = F.N(z);
    L = F.L(z);
    enclosureDebug<<"current z: \n"; F.printModesIndex(z, enclosureDebug);
    //Initialization of validated array
    for(i = 0; i < dimension; i++) {
      validated[i] = true;
    }
    for(i = 0; i < dimension; i++) {
      if(F.isDissipative(i)) { //dissipative direction we use the algorithm based on isolation from [Z3]
        enclosureDebug<<"\ndirection "<<F.array2modeIndex(i)<<"is dissipative.\n";
        N[i] = Nz[i];
        b[i] = N[i] / -F.lambda(i);
        enclosureDebug<<"checking if implication from dissipative dir is true\n";
        fzn = L[i] + N[i];
        if(leftBound(b[i]) != leftBound(b[i]) || leftBound(N[i]) != leftBound(N[i]) || leftBound(g[i]) != leftBound(g[i]) || leftBound(fz[i])
            != leftBound(fz[i])) {
          enclosureDebug << "Blow up! Enclosure function.\n";
          enclosureDebug << "b[i]!=b[i] " << (b[i] != b[i]) << ", N[i]!=N[i] " << (N[i] != N[i]) << ", g[i]!=g[i] " << (g[i] != g[i]) << ", fz[i]!=fz[i] "
              << (fz[i] != fz[i]) << "\n";
          std::cerr << "Blow up! Enclosure function.\n";
          std::cerr << "b[i]!=b[i] " << (b[i] != b[i]) << ", N[i]!=N[i] " << (N[i] != N[i]) << ", g[i]!=g[i] " << (g[i] != g[i]) << ", fz[i]!=fz[i] "
              << (fz[i] != fz[i]) << "\n";
          throw capd::jaco::SetBlowUpException<VectorT>(x, "enclosure");
        }
        if(!(fz[i].subset(fzn))) {
          enclosureDebug<<"!(fz[i]("<<fz[i]<<").subset(fzn)("<<fzn<<"))\n";
          d_r = fz[i].rightBound() - fzn.rightBound();
          if(d_r < 0)
            d_r = 0;
          d_l = fzn.leftBound() - fz[i].leftBound();
          if(d_l < 0)
            d_l = 0;
          N[i] += ScalarT(-d_l, d_r);
          fzn += ScalarT(-d_l, d_r);
          enclosureDebug<<"N[i] diam="<<ScalarT(-d_l, d_r)<<"\n";
        }
        if(fz[i] != fz[i]) {
          enclosureDebug << "Blow up! Enclosure function.\n";
          throw capd::jaco::SetBlowUpException<VectorT>(x, "enclosure");
        }
#if __VERIFY_ENCLOSURE__
        if(fz[i].subset(fzn)) {
          //    enclosureDebug<<"OK\n";
        } else {
          enclosureDebug << "SOMETHING WENT WRONG: " << fz[i] << ".subsetInterior(" << fzn << "\n";
          enclosureDebug << "F.lambda(i): " << F.lambda(i) << " z[i]= " << z[i] << " L[i]=" << L[i] << " N[i]=" << N[i] << "\n";
          enclosureDebug << "N(z)[i]=" << F.N(z)[i] << "\n";
          if(fz[i].leftBound() < (F.lambda(i) * z[i] + N[i]).leftBound())
            enclosureDebug << "left bound ! " << fz[i].leftBound() << " < " << fzn.leftBound() << " diam: " << (fzn).leftBound() - fz[i].leftBound() << "\n";
          if(fz[i].rightBound() > (F.lambda(i) * z[i] + N[i]).rightBound())
            enclosureDebug << "right bound ! " << fz[i].rightBound() << " < " << fzn.rightBound() << " diam: " << fz[i].rightBound() - (fzn).rightBound()
                << "\n";
          throw std::runtime_error("WARNING\n Error calculating enclosure, N^-[] and N^+[] aren't proper. \n See temp.txt for details.\n");
        }
#endif
        enclosureDebug<<"lambda(i)="<<F.lambda(i)<<"\n";
        enclosureDebug<<"z[i]="<<z[i]<<"\n"<<"x+h*Fz[i]="<<x[i]+h*fz[i]<<"\n"<<"b[i] for this dissipative direction: "<<b[i]<<"\n";
        enclosureDebug<<"N[i]="<<N[i]<<"\n"<<"F.L[i] +N[i]="<<fzn<<"\n"<<"F(z)[i]="<<fz[i]<<"\n";
        if(!b[i].subsetInterior(x[i])) {
          g[i] = F.g_encl(step, i, x[i], b[i]);
          enclosureDebug<<"g[i]="<<g[i]<<"\n";
          w[i] = diam(intervalHull(g[i], x[i]));
          enclosureDebug<<"w[i]="<<w[i]<<"\n"<<"z[i].right<=b[i].right "<<!(z[i].rightBound()>b[i].rightBound())<<"\n";
          enclosureDebug<<"z[i].right("<<z[i].rightBound()<<")<g[i].right("<<g[i].rightBound()<<") "<<(z[i].rightBound()<g[i].rightBound())<<"\n";
          if(!(z[i].rightBound() > b[i].rightBound()) && z[i].rightBound() < g[i].rightBound()) {
            validated[i] = false;
            allValidated = false;
            z[i].setRightBound(b[i].rightBound());
            if(g[i].rightBound() + cd * w[i].rightBound() < b[i].rightBound()) {
              z[i].setRightBound(g[i].rightBound() + cd * w[i].rightBound());
            }
            enclosureDebug<<"and new z[i]="<<z[i]<<"\n";
          }
          enclosureDebug<<"z[i].left>=b[i].left "<<!(z[i].leftBound()<b[i].leftBound())<<"\n";
          enclosureDebug<<"z[i].left("<<z[i].leftBound()<<")>g[i].left("<<g[i].leftBound()<<") "<<(z[i].leftBound()>g[i].leftBound())<<"\n";
          //symetrically for leftBound
          if(!(z[i].leftBound() < b[i].leftBound()) && z[i].leftBound() > g[i].leftBound()) {
            validated[i] = false;
            allValidated = false;
            z[i].setLeftBound(b[i].leftBound());
            if(g[i].leftBound() - cd * w[i].leftBound() > b[i].leftBound()) {
              z[i].setLeftBound(g[i].leftBound() - cd * w[i].leftBound());
            }
            enclosureDebug<<"and new z[i]="<<z[i]<<"\n";
          }
        }
        y[i] = z[i];
      } else { // direction is not dissipative, we use C^1 rough enclosure algorithm, see [ZLo]
        enclosureDebug<<"\ndirection "<<F.array2modeIndex(i)<<" is NOT dissipative(lambda="<<F.lambda(i)<<", ni="<<F.ni(i)<<").\n";
        y[i] = x[i] + h * fz[i];
        enclosureDebug<<"z[i]="<<z[i]<<"\n"<<"fz[i]="<<fz[i]<<"\n"<<"y[i]="<<y[i]<<"\n";
        if(!(y[i].subsetInterior(z[i]))) {
          enclosureDebug<<"!(y[i].subsetInterior(z[i]))\n";
          validated[i] = false;
          allValidated = false;
          sm = intervalHull(z[i], y[i]);
          z[i] = inflate<FadMapT> (sm, c);
          enclosureDebug<<"new z[i]="<<z[i]<<"\n";
        }
      }
    }
  }
  if(!allValidated) {
    for(i = 0; i < dimension; i++)
      if(!validated[i])
        enclosureDebug << "ith index for i=" << F.array2mode(i) << " not validated\n";
//    enclosureDebug << "DEBUGGING ENCLOSURE:\n";
//    debug_enclosure(vField, x, step, enclosureDebug);
    enclosureDebug << "Problem occurred in dissipative_enclosure function. Possible blow up.\n";
    std::cerr << "Problem occurred in dissipative_enclosure function. Possible blow up.\n";
    throw capd::jaco::SetBlowUpException<VectorT>(x, "enclosure");
  }
  int refinementSteps = __REFINEMENT_STEPS__;
  enclosureDebug<<"calculated candidate for enclosure:\n";
  F.printModesIndex(y, enclosureDebug);
  VectorT Fy;

  //refinement steps
  for(j = 0; j < refinementSteps; j++) {
    Fy = F(y);
    N = F.N(y);
    for(i = 0; i < dimension; i++) {
      if(F.isDissipative(i)) {
        b[i] = N[i] / -F.lambda(i);
        g[i] = F.g_encl(step, i, x[i], b[i]);
        if(x[i].rightBound() >= b[i].rightBound())
          y[i].setRightBound(x[i].rightBound());
        else
          y[i].setRightBound(g[i].rightBound());
        if(x[i].leftBound() <= b[i].leftBound())
          y[i].setLeftBound(x[i].leftBound());
        else
          y[i].setLeftBound(g[i].leftBound());
      } else {
        y[i] = x[i] + h * Fy[i];
      }
    }
    enclosureDebug<<"refinement step: "<<j<<"\n";
    F.printModesIndex(y, enclosureDebug);
  }
#if __VERIFY_ENCLOSURE__
  //checking if found enclosure satisfies the theorem assumption
  //candidate is in y
  N = F.N(z);
  fz = F(z);
  for(i = 0; i < dimension; i++) {
    if(!(F.isDissipative(i))) {
      y[i] = x[i] + h * fz[i];
      enclosureDebug<<"checking enclosure at "<<F.array2modeIndex(i)<<" yi "<<y[i]<<" subset zi "<<z[i]<<" \n";
      if(!(y[i].subsetInterior(z[i]))) {
        enclosureDebug << "Problem occurred in dissipative_enclosure function. Possible blow up."
            "at " << F.array2mode(i) << "th coordinate y[i] (" << y[i] << ") subsetInt (" << z[i] << ")\n";
        std::cerr << "Problem occurred in dissipative_enclosure function. Possible blow up.\n"
            "at " << F.array2mode(i) << "th coordinate y[i] (" << y[i] << ") subsetInt (" << z[i] << ")";
        throw capd::jaco::SetBlowUpException<VectorT>(x, "enclosure");
      }
    } else {
      b[i] = N[i] / -F.lambda(i);
      g[i] = F.g_encl(step, i, x[i], b[i]);
      enclosureDebug<<"checking enclosure at "<<F.array2modeIndex(i)<<" yi="<<y[i]<<" x[i]="<<x[i]<<" b[i]="<<b[i]<<" g[i]="<<g[i]<<" \n";
      if(z[i].rightBound() < b[i].rightBound()) {
        if(!(z[i].rightBound() >= g[i].rightBound())) {
          enclosureDebug << "Problem occurred in dissipative_enclosure function. Possible blow up."
              "at " << F.array2mode(i) << "th coordinate z[i] (" << z[i] << " g[i]=" << g[i] << " b[i]=" << b[i]
              << " w[i]=" << w[i] << " x[i]=" << x[i] << "!(z[i].rightBound()(" <<
              z[i].rightBound() << ")>=g[i].rightBound()(" << g[i].rightBound() << "))\n";
          std::cerr << "Problem occurred in dissipative_enclosure function. Possible blow up.\n"
              "at " << F.array2mode(i) << "th coordinate z[i] (" << z[i] << " g[i]=" << g[i] << " b[i]=" << b[i]
              << " w[i]=" << w[i] << " x[i]=" << x[i] <<
              "!(z[i].rightBound()(" << z[i].rightBound() << ")>=g[i].rightBound()(" << g[i].rightBound() << "))\n";
          //enclosureDebug << "DEBUGGING ENCLOSURE:\n";
          //debug_enclosure(vField, x, step, enclosureDebug);
          throw capd::jaco::SetBlowUpException<VectorT>(x, "enclosure");
        }
      }
      if(z[i].leftBound() > b[i].leftBound()) {
        if(!(z[i].leftBound() <= g[i].leftBound())) {
          enclosureDebug << "Problem occurred in dissipative_enclosure function. Possible blow up."
              "at " << F.array2mode(i) << "th coordinate z[i] (" << z[i] << " g[i]=" << g[i] << " b[i]=" << b[i]
              << " w[i]=" << w[i] << " x[i]=" << x[i] <<
              "!(z[i].leftBound()(" << z[i].leftBound() << ")>=g[i].leftBound()(" << g[i].leftBound() << "))\n";
          std::cerr << "Problem occurred in dissipative_enclosure function. Possible blow up.\n"
              "at " << F.array2mode(i) << "th coordinate z[i] (" << z[i] << " g[i]=" << g[i] << " b[i]=" << b[i]
              << " w[i]=" << w[i] << " x[i]=" << x[i] <<
              "!(z[i].leftBound()(" << z[i].leftBound() << ")>=g[i].leftBound()(" << g[i].leftBound() << "))\n";
          //enclosureDebug << "DEBUGGING ENCLOSURE:\n";
          //debug_enclosure(vField, x, step, enclosureDebug);
          throw capd::jaco::SetBlowUpException<VectorT>(x, "enclosure");
        }
      }
    }
  }
#endif
  return y;
}

template <typename EquationT>
typename EquationT::PolyBdType enclosure(EquationT & vField, typename EquationT::PolyBdType const & x,
    typename EquationT::RealType const & step){
  typedef typename EquationT::ScalarType ScalarT;
  typedef typename EquationT::VectorType VectorT;
  typedef typename ScalarT::BoundType DoubleT;
  typedef typename EquationT::PolyBdType PolyBdT;

  int steps = 5, s = 0, dimension = vField.dimension();  
  
  if(dimension > 5) {
    steps = dimension*2;
  }
  DoubleT c = 1.1;
  DoubleT cd = 0.1;
  int n = x.n;
  PolyBdT z(x), y(x), tpb(x);
  VectorT N(dimension), b(dimension), g(dimension), w(dimension), fz(dimension), Nz(dimension), Ftemp(dimension), L(dimension), fx(dimension);
  std::vector<bool> validated(dimension);
  bool allValidated = false;
  int i, j;
  ScalarT h = ScalarT(0, 1) * step;
  EquationT& F = vField;

  ScalarT small = ScalarT(-1, 1) * __SMALL__;
  ScalarT sm;
  //initial guess  
  F(x, tpb, F.getFinitePart());
  fx = tpb;  
    
  enclosureDebug << "calculated F: " << tpb << "\n";
  
  for(i = 0; i < dimension; i++) {
    if(F.isDissipative(i)) {
      z[i] = x[i];
    } else {
      z[i] = x[i] + h * fx[i];
      z[i] = inflate<EquationT> (z[i], c);
      //z[i]+=small;//here we increase z by small vector to be sure it has nonempty interior
    }
    z[i] += small;
    validated[i] = false;
  }
  double d_l, d_r;
  ScalarT fzn;
  while(!allValidated && s++ < steps) {
    enclosureDebug<<"step: "<<s<<"\n";
    allValidated = true;    
    F.N(z, tpb, F.getFinitePart());
    Nz = tpb;
    F.L(z, tpb);
    L = tpb;
//    F(z, tpb, F.getFinitePart()); //this was redundant
    fz = Nz;
    fz += L;
    enclosureDebug<<"current z: \n"<<z;
    //Initialization of validated array
    for(i = 0; i < dimension; i++) {
      validated[i] = true;
    }
    for(i = 0; i < dimension; i++) {
      if(F.isDissipative(i)) { //dissipative direction we use the algorithm based on isolation from [Z3]
        enclosureDebug<<"\ndirection "<<F.array2modeIndex(i)<<"is dissipative.\n";
        N[i] = Nz[i];
        b[i] = N[i] / -F.lambda(i);
        enclosureDebug<<"checking if implication from dissipative dir is true\n";
        fzn = L[i] + N[i];
        if(leftBound(b[i]) != leftBound(b[i]) || leftBound(N[i]) != leftBound(N[i]) || leftBound(g[i]) != leftBound(g[i]) || leftBound(fz[i])
            != leftBound(fz[i])) {
          enclosureDebug << "Blow up! Enclosure function.\n";
          enclosureDebug << "b[i]!=b[i] " << (b[i] != b[i]) << ", N[i]!=N[i] " << (N[i] != N[i]) << ", g[i]!=g[i] " << (g[i] != g[i]) << ", fz[i]!=fz[i] "
              << (fz[i] != fz[i]) << "\n";
          std::cerr << "Blow up! Enclosure function.\n";
          std::cerr << "b[i]!=b[i] " << (b[i] != b[i]) << ", N[i]!=N[i] " << (N[i] != N[i]) << ", g[i]!=g[i] " << (g[i] != g[i]) << ", fz[i]!=fz[i] "
              << (fz[i] != fz[i]) << "\n";
          throw capd::jaco::SetBlowUpException<PolyBdT>(x, "enclosure");
        }
        if(!(fz[i].subset(fzn))) {
          enclosureDebug<<"!(fz[i]("<<fz[i]<<").subset(fzn)("<<fzn<<"))\n";
          d_r = fz[i].rightBound() - fzn.rightBound();
          if(d_r < 0)
            d_r = 0;
          d_l = fzn.leftBound() - fz[i].leftBound();
          if(d_l < 0)
            d_l = 0;
          N[i] += ScalarT(-d_l, d_r);
          fzn += ScalarT(-d_l, d_r);
          enclosureDebug<<"N[i] diam="<<ScalarT(-d_l, d_r)<<"\n";
        }
        if(fz[i] != fz[i]) {
          enclosureDebug << "Blow up! Enclosure function.\n";
          throw capd::jaco::SetBlowUpException<PolyBdT>(x, "enclosure");
        }
#if __VERIFY_ENCLOSURE__
        if(fz[i].subset(fzn)) {
          //    enclosureDebug<<"OK\n";
        } else {
          enclosureDebug << "SOMETHING WENT WRONG: " << fz[i] << ".subsetInterior(" << fzn << "\n";
          enclosureDebug << "F.lambda(i): " << F.lambda(i) << " z[i]= " << z[i] << " L[i]=" << L[i] << " N[i]=" << N[i] << "\n";
          if(fz[i].leftBound() < (F.lambda(i) * z[i] + N[i]).leftBound())
            enclosureDebug << "left bound ! " << fz[i].leftBound() << " < " << fzn.leftBound() << " diam: " << (fzn).leftBound() - fz[i].leftBound() << "\n";
          if(fz[i].rightBound() > (F.lambda(i) * z[i] + N[i]).rightBound())
            enclosureDebug << "right bound ! " << fz[i].rightBound() << " < " << fzn.rightBound() << " diam: " << fz[i].rightBound() - (fzn).rightBound()
                << "\n";
          throw std::runtime_error("WARNING\n Error calculating enclosure, N^-[] and N^+[] aren't proper. \n See temp.txt for details.\n");
        }
#endif
        enclosureDebug<<"lambda(i)="<<F.lambda(i)<<"\n";
        enclosureDebug<<"z[i]="<<z[i]<<"\n"<<"x+h*Fz[i]="<<x[i]+h*fz[i]<<"\n"<<"b[i] for this dissipative direction: "<<b[i]<<"\n";
        enclosureDebug<<"N[i]="<<N[i]<<"\n"<<"F.L[i] +N[i]="<<fzn<<"\n"<<"F(z)[i]="<<fz[i]<<"\n";
        if(!b[i].subsetInterior(x[i])) {
          g[i] = F.g_encl(step, i, x[i], b[i]);
          enclosureDebug<<"g[i]="<<g[i]<<"\n";
          w[i] = diam(intervalHull(g[i], x[i]));
          enclosureDebug<<"w[i]="<<w[i]<<"\n"<<"z[i].right<=b[i].right "<<!(z[i].rightBound()>b[i].rightBound())<<"\n";
          enclosureDebug<<"z[i].right("<<z[i].rightBound()<<")<g[i].right("<<g[i].rightBound()<<") "<<(z[i].rightBound()<g[i].rightBound())<<"\n";
          if(!(z[i].rightBound() > b[i].rightBound()) && z[i].rightBound() < g[i].rightBound()) {
            validated[i] = false;
            allValidated = false;
            z[i].setRightBound(b[i].rightBound());
            if(g[i].rightBound() + cd * w[i].rightBound() < b[i].rightBound()) {
              z[i].setRightBound(g[i].rightBound() + cd * w[i].rightBound());
            }
            enclosureDebug<<"and new z[i]="<<z[i]<<"\n";
          }
          enclosureDebug<<"z[i].left>=b[i].left "<<!(z[i].leftBound()<b[i].leftBound())<<"\n";
          enclosureDebug<<"z[i].left("<<z[i].leftBound()<<")>g[i].left("<<g[i].leftBound()<<") "<<(z[i].leftBound()>g[i].leftBound())<<"\n";
          //symetrically for leftBound
          if(!(z[i].leftBound() < b[i].leftBound()) && z[i].leftBound() > g[i].leftBound()) {
            validated[i] = false;
            allValidated = false;
            z[i].setLeftBound(b[i].leftBound());
            if(g[i].leftBound() - cd * w[i].leftBound() > b[i].leftBound()) {
              z[i].setLeftBound(g[i].leftBound() - cd * w[i].leftBound());
            }
            enclosureDebug<<"and new z[i]="<<z[i]<<"\n";
          }
        }
        y[i] = z[i];
      } else { // direction is not dissipative, we use C^1 rough enclosure algorithm, see [ZLo]
        enclosureDebug<<"\ndirection "<<F.array2modeIndex(i)<<" is NOT dissipative(lambda="<<F.lambda(i)<<", ni="<<F.ni(i)<<").\n";
        y[i] = x[i] + h * fz[i];
        enclosureDebug<<"z[i]="<<z[i]<<"\n"<<"fz[i]="<<fz[i]<<"\n"<<"y[i]="<<y[i]<<"\n";
        if(!(y[i].subsetInterior(z[i]))) {
          enclosureDebug<<"!(y[i].subsetInterior(z[i]))\n";
          validated[i] = false;
          allValidated = false;
          sm = intervalHull(z[i], y[i]);
          z[i] = inflate<EquationT> (sm, c);
          enclosureDebug<<"new z[i]="<<z[i]<<"\n";
        }
      }
    }
  }
  if(!allValidated) {
    for(i = 0; i < dimension; i++)
      if(!validated[i])
        enclosureDebug << "ith index for i=" << F.array2mode(i) << " not validated\n";
//    enclosureDebug << "DEBUGGING ENCLOSURE:\n";
//    debug_enclosure(vField, x, step, enclosureDebug);
    enclosureDebug << "Problem occurred in dissipative_enclosure function. Possible blow up.\n";
    std::cerr << "Problem occurred in dissipative_enclosure function. Possible blow up.\n";
    throw capd::jaco::SetBlowUpException<PolyBdT>(x, "enclosure");
  }
  int refinementSteps = __REFINEMENT_STEPS__;
  enclosureDebug<<"calculated candidate for enclosure:\n";
  F.printModesIndex(y, n, enclosureDebug);
  PolyBdT Fy(x.n);


  //refinement steps
  for(j = 0; j < refinementSteps; j++) {    
    F.N(y, tpb, F.getFinitePart());
    N = tpb;
    F.L(y, tpb);
    L = tpb;
    //F(y, tpb, F.getFinitePart()); //this was redundant
    Fy = N;
    Fy += L;
    for(i = 0; i < dimension; i++) {
      if(F.isDissipative(i)) {
        b[i] = N[i] / -F.lambda(i);
        g[i] = F.g_encl(step, i, x[i], b[i]);
        if(x[i].rightBound() >= b[i].rightBound())
          y[i].setRightBound(x[i].rightBound());
        else
          y[i].setRightBound(g[i].rightBound());
        if(x[i].leftBound() <= b[i].leftBound())
          y[i].setLeftBound(x[i].leftBound());
        else
          y[i].setLeftBound(g[i].leftBound());
      } else {
        y[i] = x[i] + h * Fy[i];
      }
    }
    enclosureDebug<<"refinement step: "<<j<<"\n";
    F.printModesIndex(y, n, enclosureDebug);
  }
#if __VERIFY_ENCLOSURE__
  //checking if found enclosure satisfies the theorem assumption
  //candidate is in y
  F.N(z, tpb, F.getFinitePart());
  N = tpb;
  F.L(z, tpb);
  L = tpb;
//  F(z, tpb, F.getFinitePart());
  fz = N;
  fz += L;
  for(i = 0; i < dimension; i++) {
    if(!(F.isDissipative(i))) {
      y[i] = x[i] + h * fz[i];
      enclosureDebug<<"checking enclosure at "<<F.array2modeIndex(i)<<" yi "<<y[i]<<" subset zi "<<z[i]<<" \n";
      if(!(y[i].subsetInterior(z[i]))) {
        enclosureDebug << "Problem occurred in dissipative_enclosure function. Possible blow up."
            "at " << F.array2mode(i) << "th coordinate y[i] (" << y[i] << ") subsetInt (" << z[i] << ")\n";
        std::cerr << "Problem occurred in dissipative_enclosure function. Possible blow up.\n"
            "at " << F.array2mode(i) << "th coordinate y[i] (" << y[i] << ") subsetInt (" << z[i] << ")";
        throw capd::jaco::SetBlowUpException<PolyBdT>(x, "enclosure");
      }
    } else {
      b[i] = N[i] / -F.lambda(i);
      g[i] = F.g_encl(step, i, x[i], b[i]);
      enclosureDebug<<"checking enclosure at "<<F.array2modeIndex(i)<<" yi="<<y[i]<<" x[i]="<<x[i]<<" b[i]="<<b[i]<<" g[i]="<<g[i]<<" \n";
      if(z[i].rightBound() < b[i].rightBound()) {
        if(!(z[i].rightBound() >= g[i].rightBound())) {
          enclosureDebug << "Problem occurred in dissipative_enclosure function. Possible blow up."
              "at " << F.array2mode(i) << "th coordinate z[i] (" << z[i] << " g[i]=" << g[i] << " b[i]=" << b[i]
              << " w[i]=" << w[i] << " x[i]=" << x[i] << "!(z[i].rightBound()(" <<
              z[i].rightBound() << ")>=g[i].rightBound()(" << g[i].rightBound() << "))\n";
          std::cerr << "Problem occurred in dissipative_enclosure function. Possible blow up.\n"
              "at " << F.array2mode(i) << "th coordinate z[i] (" << z[i] << " g[i]=" << g[i] << " b[i]=" << b[i]
              << " w[i]=" << w[i] << " x[i]=" << x[i] <<
              "!(z[i].rightBound()(" << z[i].rightBound() << ")>=g[i].rightBound()(" << g[i].rightBound() << "))\n";
          //enclosureDebug << "DEBUGGING ENCLOSURE:\n";
          //debug_enclosure(vField, x, step, enclosureDebug);
          throw capd::jaco::SetBlowUpException<PolyBdT>(x, "enclosure");
        }
      }
      if(z[i].leftBound() > b[i].leftBound()) {
        if(!(z[i].leftBound() <= g[i].leftBound())) {
          enclosureDebug << "Problem occurred in dissipative_enclosure function. Possible blow up."
              "at " << F.array2mode(i) << "th coordinate z[i] (" << z[i] << " g[i]=" << g[i] << " b[i]=" << b[i]
              << " w[i]=" << w[i] << " x[i]=" << x[i] <<
              "!(z[i].leftBound()(" << z[i].leftBound() << ")>=g[i].leftBound()(" << g[i].leftBound() << "))\n";
          std::cerr << "Problem occurred in dissipative_enclosure function. Possible blow up.\n"
              "at " << F.array2mode(i) << "th coordinate z[i] (" << z[i] << " g[i]=" << g[i] << " b[i]=" << b[i]
              << " w[i]=" << w[i] << " x[i]=" << x[i] <<
              "!(z[i].leftBound()(" << z[i].leftBound() << ")>=g[i].leftBound()(" << g[i].leftBound() << "))\n";
          //enclosureDebug << "DEBUGGING ENCLOSURE:\n";
          //debug_enclosure(vField, x, step, enclosureDebug);
          throw capd::jaco::SetBlowUpException<PolyBdT>(x, "enclosure");
        }
      }
    }
  }
#endif
  return y;

}


template <typename FadMapT>
typename FadMapT::VectorType enclosureInputTail(FadMapT & vField, typename FadMapT::MatrixType::RowVectorType const & x,
    typename FadMapT::TailType const & T, typename FadMapT::ScalarType const & step) {
  typedef typename FadMapT::ScalarType ScalarT;
  typedef typename FadMapT::VectorType VectorT;
  typedef typename ScalarT::BoundType DoubleT;

  int steps = 5, s = 0, dimension = vField.getDimension();
  if(dimension > 5) {
    steps = dimension*2;
  }
  DoubleT c = 1.1;
  DoubleT cd = 0.1;
  VectorT z(dimension), N(dimension), b(dimension), g(dimension), y(dimension), w(dimension), fz(dimension), Nz(dimension), Ftemp(dimension), L(
      dimension);
  std::vector<bool> validated(dimension);
  bool allValidated = false;
  int i, j;
  ScalarT h = ScalarT(0, 1) * step;
  FadMapT& F = vField;
  ScalarT small = ScalarT(-1, 1) * __SMALL__;
  ScalarT sm;
  //initial guess
  VectorT fx = F(x, T);
  for(i = 0; i < dimension; i++) {
    if(F.lambda(i) < 0) {
      z[i] = x[i];
    } else {
      z[i] = x[i] + h * fx[i];
      z[i] = inflate<FadMapT> (z[i], c);
      //z[i]+=small;//here we increase z by small vector to be sure it has nonempty interior
    }
    z[i] += small;
    validated[i] = false;
  }
  double d_l, d_r;
  ScalarT fzn;
  while(!allValidated && s++ < steps) {
    enclosureDebug<<"step: "<<s<<"\n";
    allValidated = true;
    fz = F(z, T);
    Nz = F.N(z, T);
    L = F.L(z);
    enclosureDebug<<"current z: \n"; F.printModes(z, enclosureDebug);
    //Initialisation of validated array
    for(i = 0; i < dimension; i++) {
      validated[i] = true;
    }
    for(i = 0; i < dimension; i++) {
      if(F.isDissipative(i)) { //dissipative direction we use the algorithm based on isolation from [Z3]
        enclosureDebug<<"\ndirection "<<F.array2mode(i)<<"is dissipative.\n";
        N[i] = Nz[i];
        b[i] = N[i] / -F.lambda(i);
        enclosureDebug<<"checking if implication from dissipative dir is true\n";
        fzn = L[i] + N[i];
        if(leftBound(b[i]) != leftBound(b[i]) || leftBound(N[i]) != leftBound(N[i]) || leftBound(g[i]) != leftBound(g[i]) || leftBound(fz[i])
            != leftBound(fz[i])) {
          enclosureDebug << "Blow up! Enclosure function.\n";
          enclosureDebug << "b[i]!=b[i] " << (b[i] != b[i]) << ", N[i]!=N[i] " << (N[i] != N[i]) << ", g[i]!=g[i] " << (g[i] != g[i]) << ", fz[i]!=fz[i] "
              << (fz[i] != fz[i]) << "\n";
          std::cerr << "Blow up! Enclosure function.\n";
          std::cerr << "b[i]!=b[i] " << (b[i] != b[i]) << ", N[i]!=N[i] " << (N[i] != N[i]) << ", g[i]!=g[i] " << (g[i] != g[i]) << ", fz[i]!=fz[i] "
              << (fz[i] != fz[i]) << "\n";
          throw capd::jaco::SetBlowUpException<VectorT>(x, "enclosure");
        }
        if(!(fz[i].subset(fzn))) {
          enclosureDebug<<"!(fz[i]("<<fz[i]<<").subset(fzn)("<<fzn<<"))\n";
          d_r = fz[i].rightBound() - fzn.rightBound();
          if(d_r < 0)
            d_r = 0;
          d_l = fzn.leftBound() - fz[i].leftBound();
          if(d_l < 0)
            d_l = 0;
          N[i] += ScalarT(-d_l, d_r);
          fzn += ScalarT(-d_l, d_r);
          enclosureDebug<<"N[i] diam="<<ScalarT(-d_l, d_r)<<"\n";
        }
        if(fz[i] != fz[i]) {
          enclosureDebug << "Blow up! Enclosure function.\n";
          throw capd::jaco::SetBlowUpException<VectorT>(x, "enclosure");
        }
#if __VERIFY_ENCLOSURE__
        if(fz[i].subset(fzn)) {
          //    out<<"OK\n";
        } else {
          enclosureDebug << "SOMETHING WENT WRONG: " << fz[i] << ".subsetInterior(" << fzn << "\n";
          enclosureDebug << "F.lambda(i): " << F.lambda(i) << " z[i]= " << z[i] << " L[i]=" << L[i] << " N[i]=" << N[i] << "\n";
          enclosureDebug << "N(z)[i]=" << F.N(z, T)[i] << "\n";
          if(fz[i].leftBound() < (F.lambda(i) * z[i] + N[i]).leftBound())
            enclosureDebug << "left bound ! " << fz[i].leftBound() << " < " << fzn.leftBound() << " diam: " << (fzn).leftBound() - fz[i].leftBound() << "\n";
          if(fz[i].rightBound() > (F.lambda(i) * z[i] + N[i]).rightBound())
            enclosureDebug << "right bound ! " << fz[i].rightBound() << " < " << fzn.rightBound() << " diam: " << fz[i].rightBound() - (fzn).rightBound()
                << "\n";
          throw std::runtime_error("WARNING\n Error calculating enclosure, N^-[] and N^+[] aren't proper. \n See temp.txt for details.\n");
        }
#endif

        enclosureDebug<<"z[i]="<<z[i]<<"\n"<<"x+h*Fz[i]="<<x[i]+h*fz[i]<<"\n"<<"b[i] for this dissipative direction: "<<b[i]<<"\n";
        enclosureDebug<<"N[i]="<<N[i]<<"\n"<<"F.L[i] +N[i]="<<fzn<<"\n"<<"F(z)[i]="<<fz[i]<<"\n";
        if(!b[i].subsetInterior(x[i])) {
          g[i] = F.g_encl(step, i, x[i], b[i]);
          enclosureDebug<<"g[i]="<<g[i]<<"\n";
          w[i] = diam(intervalHull(g[i], x[i]));
          enclosureDebug<<"w[i]="<<w[i]<<"\n"<<"z[i].right<=b[i].right "<<!(z[i].rightBound()>b[i].rightBound())<<"\n";
          enclosureDebug<<"z[i].right("<<z[i].rightBound()<<")<g[i].right("<<g[i].rightBound()<<") "<<(z[i].rightBound()<g[i].rightBound())<<"\n";
          if(!(z[i].rightBound() > b[i].rightBound()) && z[i].rightBound() < g[i].rightBound()) {
            validated[i] = false;
            allValidated = false;
            z[i].setRightBound(b[i].rightBound());
            if(g[i].rightBound() + cd * w[i].rightBound() < b[i].rightBound()) {
              z[i].setRightBound(g[i].rightBound() + cd * w[i].rightBound());
            }
            enclosureDebug<<"and new z[i]="<<z[i]<<"\n";
          }
          enclosureDebug<<"z[i].left>=b[i].left "<<!(z[i].leftBound()<b[i].leftBound())<<"\n";
          enclosureDebug<<"z[i].left("<<z[i].leftBound()<<")>g[i].left("<<g[i].leftBound()<<") "<<(z[i].leftBound()>g[i].leftBound())<<"\n";
          //symetrically for leftBound
          if(!(z[i].leftBound() < b[i].leftBound()) && z[i].leftBound() > g[i].leftBound()) {
            validated[i] = false;
            allValidated = false;
            z[i].setLeftBound(b[i].leftBound());
            if(g[i].leftBound() - cd * w[i].leftBound() > b[i].leftBound()) {
              z[i].setLeftBound(g[i].leftBound() - cd * w[i].leftBound());
            }
            enclosureDebug<<"and new z[i]="<<z[i]<<"\n";
          }
        }
        y[i] = z[i];
      } else { // direction is not dissipative, we use C^1 rough enclosure algorithm, see [ZLo]
        enclosureDebug<<"\ndirection "<<F.array2mode(i)<<" is NOT dissipative(lambda="<<F.lambda(i)<<", ni="<<F.ni(i)<<").\n";
        y[i] = x[i] + h * fz[i];
        enclosureDebug<<"z[i]="<<z[i]<<"\n"<<"fz[i]="<<fz[i]<<"\n"<<"y[i]="<<y[i]<<"\n";
        if(!(y[i].subsetInterior(z[i]))) {
          enclosureDebug<<"!(y[i].subsetInterior(z[i]))\n";
          validated[i] = false;
          allValidated = false;
          sm = intervalHull(z[i], y[i]);
          z[i] = inflate<FadMapT> (sm, c);
          enclosureDebug<<"new z[i]="<<z[i]<<"\n";
        }
      }
    }
  }
  if(!allValidated) {
    for(i = 0; i < dimension; i++)
      if(!validated[i])
        std::cerr<< "ith index for i=" << F.array2mode(i) << " not validated\n";
//    out << "DEBUGGING ENCLOSURE:\n";
//    debug_enclosure(vField, x, step, out);
    std::cerr << "Problem occurred in dissipative_enclosure function. Possible blow up.\n";
    throw capd::jaco::SetBlowUpException<VectorT>(x, "enclosure");
  }
  int refinementSteps = __REFINEMENT_STEPS__;
  enclosureDebug<<"calculated candidate for enclosure:\n";
  F.printModes(y, enclosureDebug);
  VectorT Fy;

  //refinement steps
  for(j = 0; j < refinementSteps; j++) {
    Fy = F(y, T);
    N = F.N(y, T);
    for(i = 0; i < dimension; i++) {
      if(F.isDissipative(i)) {
        b[i] = N[i] / -F.lambda(i);
        g[i] = F.g_encl(step, i, x[i], b[i]);
        if(x[i].rightBound() >= b[i].rightBound())
          y[i].setRightBound(x[i].rightBound());
        else
          y[i].setRightBound(g[i].rightBound());
        if(x[i].leftBound() <= b[i].leftBound())
          y[i].setLeftBound(x[i].leftBound());
        else
          y[i].setLeftBound(g[i].leftBound());
      } else {
        y[i] = x[i] + h * Fy[i];
      }
    }
    enclosureDebug<<"refinement step: "<<j<<"\n";
    F.printModes(y, enclosureDebug);
  }
#if __VERIFY_ENCLOSURE__
  //checking if found enclosure satisfies the theorem assumption
  //candidate is in y
  N = F.N(z, T);
  fz = F(z, T);
  for(i = 0; i < dimension; i++) {
    if(!(F.isDissipative(i))) {
      y[i] = x[i] + h * fz[i];
      enclosureDebug<<"checking enclosure at "<<F.array2mode(i)<<" yi "<<y[i]<<" subset zi "<<z[i]<<" \n";
      if(!(y[i].subsetInterior(z[i]))) {
        enclosureDebug << "Problem occurred in dissipative_enclosure function. Possible blow up."
            "at " << F.array2mode(i) << "th coordinate y[i] (" << y[i] << ") subsetInt (" << z[i] << ")\n";
        std::cerr << "Problem occurred in dissipative_enclosure function. Possible blow up.\n"
            "at " << F.array2mode(i) << "th coordinate y[i] (" << y[i] << ") subsetInt (" << z[i] << ")";
        throw capd::jaco::SetBlowUpException<VectorT>(x, "enclosure");
      }
    } else {
      b[i] = N[i] / -F.lambda(i);
      g[i] = F.g_encl(step, i, x[i], b[i]);
      enclosureDebug<<"checking enclosure at "<<F.array2mode(i)<<" yi="<<y[i]<<" b[i]="<<b[i]<<" g[i]="<<g[i]<<" \n";
      if(z[i].rightBound() < b[i].rightBound()) {
        if(!(z[i].rightBound() >= g[i].rightBound())) {
          enclosureDebug << "Problem occurred in dissipative_enclosure function. Possible blow up."
              "at " << F.array2mode(i) << "th coordinate z[i] (" << z[i] << " g[i]=" << g[i] << " b[i]=" << b[i]
              << " w[i]=" << w[i] << " x[i]=" << x[i] << "!(z[i].rightBound()(" <<
              z[i].rightBound() << ")>=g[i].rightBound()(" << g[i].rightBound() << "))\n";
          std::cerr << "Problem occurred in dissipative_enclosure function. Possible blow up.\n"
              "at " << F.array2mode(i) << "th coordinate z[i] (" << z[i] << " g[i]=" << g[i] << " b[i]=" << b[i]
              << " w[i]=" << w[i] << " x[i]=" << x[i] <<
              "!(z[i].rightBound()(" << z[i].rightBound() << ")>=g[i].rightBound()(" << g[i].rightBound() << "))\n";
          //out << "DEBUGGING ENCLOSURE:\n";
          //debug_enclosure(vField, x, step, out);
          throw capd::jaco::SetBlowUpException<VectorT>(x, "enclosure");
        }
      }
      if(z[i].leftBound() > b[i].leftBound()) {
        if(!(z[i].leftBound() <= g[i].leftBound())) {
          enclosureDebug << "Problem occurred in dissipative_enclosure function. Possible blow up."
              "at " << F.array2mode(i) << "th coordinate z[i] (" << z[i] << " g[i]=" << g[i] << " b[i]=" << b[i]
              << " w[i]=" << w[i] << " x[i]=" << x[i] <<
              "!(z[i].leftBound()(" << z[i].leftBound() << ")>=g[i].leftBound()(" << g[i].leftBound() << "))\n";
          std::cerr << "Problem occurred in dissipative_enclosure function. Possible blow up.\n"
              "at " << F.array2mode(i) << "th coordinate z[i] (" << z[i] << " g[i]=" << g[i] << " b[i]=" << b[i]
              << " w[i]=" << w[i] << " x[i]=" << x[i] <<
              "!(z[i].leftBound()(" << z[i].leftBound() << ")>=g[i].leftBound()(" << g[i].leftBound() << "))\n";
          //out << "DEBUGGING ENCLOSURE:\n";
          //debug_enclosure(vField, x, step, out);
          throw capd::jaco::SetBlowUpException<VectorT>(x, "enclosure");
        }
      }
    }
  }
#endif
  return y;
}

template <typename FadMapT>
typename FadMapT::VectorType enclosure(FadMapT & vField, typename FadMapT::MatrixType::RowVectorType const & x,
    typename FadMapT::ScalarType const & step, std::ostream& out) {
  typedef typename FadMapT::ScalarType ScalarT;
  typedef typename FadMapT::VectorType VectorT;
  typedef typename ScalarT::BoundType DoubleT;

  int steps = 5, s = 0, dimension = vField.getDimension();
  if(dimension > 5) {
    steps = dimension*2;
  }
  DoubleT c = 1.1;
  DoubleT cd = 0.1;
  VectorT z(dimension), N(dimension), b(dimension), g(dimension), y(dimension), w(dimension), fz(dimension), Nz(dimension), Ftemp(dimension), L(
      dimension);
  std::vector<bool> validated(dimension);
  bool allValidated = false;
  int i, j;
  ScalarT h = ScalarT(0, 1) * step;
  FadMapT& F = vField;
  ScalarT small = ScalarT(-1, 1) * __SMALL__;
  ScalarT sm;
  //initial guess
  VectorT fx = F(x);
  for(i = 0; i < dimension; i++) {
    if(F.lambda(i) < 0) {
      z[i] = x[i];
    } else {
      z[i] = x[i] + h * fx[i];
      z[i] = inflate<FadMapT> (z[i], c);
      //z[i]+=small;//here we increase z by small vector to be sure it has nonempty interior
    }
    z[i] += small;
    validated[i] = false;
  }
#if __DEBUG_ENCLOSURE__
  F.debugMsg(out);
#endif
  double d_l, d_r;
  ScalarT fzn;
  while(!allValidated && s++ < steps) {

#if __DEBUG_ENCLOSURE__
    out<<"step: "<<s<<"\n";
#endif
    allValidated = true;
    fz = F(z);
    Nz = F.N(z);
    L = F.L(z);
#if __DEBUG_ENCLOSURE__
    out<<"current z: \n"; F.printModes(z, out);
#endif
    //Initialisation of validated array
    for(i = 0; i < dimension; i++) {
      validated[i] = true;
    }
    for(i = 0; i < dimension; i++) {
      if(F.isDissipative(i)) { //dissipative direction we use the algorithm based on isolation from [Z3]
#if __DEBUG_ENCLOSURE__
        out<<"\ndirection "<<F.array2mode(i)<<"is dissipative.\n";
#endif
        N[i] = Nz[i];
        b[i] = N[i] / -F.lambda(i);
#if __DEBUG_ENCLOSURE__
        out<<"checking if implication from dissipative dir is true\n";
#endif
        fzn = L[i] + N[i];
        if(leftBound(b[i]) != leftBound(b[i]) || leftBound(N[i]) != leftBound(N[i]) || leftBound(g[i]) != leftBound(g[i]) || leftBound(fz[i])
            != leftBound(fz[i])) {
          out << "Blow up! Enclosure function.\n";
          out << "b[i]!=b[i] " << (b[i] != b[i]) << ", N[i]!=N[i] " << (N[i] != N[i]) << ", g[i]!=g[i] " << (g[i] != g[i]) << ", fz[i]!=fz[i] "
              << (fz[i] != fz[i]) << "\n";
          std::cerr << "Blow up! Enclosure function.\n";
          std::cerr << "b[i]!=b[i] " << (b[i] != b[i]) << ", N[i]!=N[i] " << (N[i] != N[i]) << ", g[i]!=g[i] " << (g[i] != g[i]) << ", fz[i]!=fz[i] "
              << (fz[i] != fz[i]) << "\n";
          throw capd::jaco::SetBlowUpException<VectorT>(x, "enclosure");
        }
        if(!(fz[i].subset(fzn))) {
#if __DEBUG_ENCLOSURE__
          out<<"!(fz[i]("<<fz[i]<<").subset(fzn)("<<fzn<<"))\n";
#endif
          d_r = fz[i].rightBound() - fzn.rightBound();
          if(d_r < 0)
            d_r = 0;
          d_l = fzn.leftBound() - fz[i].leftBound();
          if(d_l < 0)
            d_l = 0;
          N[i] += ScalarT(-d_l, d_r);
          fzn += ScalarT(-d_l, d_r);
#if __DEBUG_ENCLOSURE__
          out<<"N[i] diam="<<ScalarT(-d_l, d_r)<<"\n";
#endif
        }
        if(fz[i] != fz[i]) {
          out << "Blow up! Enclosure function.\n";
          throw capd::jaco::SetBlowUpException<VectorT>(x, "enclosure");
        }
#if __VERIFY_ENCLOSURE__
        if(fz[i].subset(fzn)) {
          //    out<<"OK\n";
        } else {
          out << "SOMETHING WENT WRONG: " << fz[i] << ".subsetInterior(" << fzn << "\n";
          out << "F.lambda(i): " << F.lambda(i) << " z[i]= " << z[i] << " L[i]=" << L[i] << " N[i]=" << N[i] << "\n";
          out << "N(z)[i]=" << F.N(z)[i] << "\n";
          if(fz[i].leftBound() < (F.lambda(i) * z[i] + N[i]).leftBound())
            out << "left bound ! " << fz[i].leftBound() << " < " << fzn.leftBound() << " diam: " << (fzn).leftBound() - fz[i].leftBound() << "\n";
          if(fz[i].rightBound() > (F.lambda(i) * z[i] + N[i]).rightBound())
            out << "right bound ! " << fz[i].rightBound() << " < " << fzn.rightBound() << " diam: " << fz[i].rightBound() - (fzn).rightBound()
                << "\n";
          throw std::runtime_error("WARNING\n Error calculating enclosure, N^-[] and N^+[] aren't proper. \n See temp.txt for details.\n");
        }
#endif

#if __DEBUG_ENCLOSURE__
        out<<"z[i]="<<z[i]<<"\n";
        out<<"x+h*Fz[i]="<<x[i]+h*fz[i]<<"\n";
        out<<"b[i] for this dissipative direction: "<<b[i]<<"\n";
        out<<"N[i]="<<N[i]<<"\n";
        out<<"F.L[i] +N[i]="<<fzn<<"\n";
        out<<"F(z)[i]="<<fz[i]<<"\n";
#endif
        if(!b[i].subsetInterior(x[i])) {
          g[i] = F.g_encl(step, i, x[i], b[i]);
#if __DEBUG_ENCLOSURE__
          out<<"g[i]="<<g[i]<<"\n";
#endif
          w[i] = diam(intervalHull(g[i], x[i]));
#if __DEBUG_ENCLOSURE__
          out<<"w[i]="<<w[i]<<"\n";
          out<<"z[i].right<=b[i].right "<<!(z[i].rightBound()>b[i].rightBound())<<"\n";
          out<<"z[i].right("<<z[i].rightBound()<<")<g[i].right("<<g[i].rightBound()<<") "<<(z[i].rightBound()<g[i].rightBound())<<"\n";
#endif
          if(!(z[i].rightBound() > b[i].rightBound()) && z[i].rightBound() < g[i].rightBound()) {
            validated[i] = false;
            allValidated = false;
            z[i].setRightBound(b[i].rightBound());
            if(g[i].rightBound() + cd * w[i].rightBound() < b[i].rightBound()) {
              z[i].setRightBound(g[i].rightBound() + cd * w[i].rightBound());
            }
#if __DEBUG_ENCLOSURE__
            out<<"and new z[i]="<<z[i]<<"\n";
#endif
          }
#if __DEBUG_ENCLOSURE__
          out<<"z[i].left>=b[i].left "<<!(z[i].leftBound()<b[i].leftBound())<<"\n";
          out<<"z[i].left("<<z[i].leftBound()<<")>g[i].left("<<g[i].leftBound()<<") "<<(z[i].leftBound()>g[i].leftBound())<<"\n";
#endif
          //symetrically for leftBound
          if(!(z[i].leftBound() < b[i].leftBound()) && z[i].leftBound() > g[i].leftBound()) {
            validated[i] = false;
            allValidated = false;
            z[i].setLeftBound(b[i].leftBound());
            if(g[i].leftBound() - cd * w[i].leftBound() > b[i].leftBound()) {
              z[i].setLeftBound(g[i].leftBound() - cd * w[i].leftBound());
            }
#if __DEBUG_ENCLOSURE__
            out<<"and new z[i]="<<z[i]<<"\n";
#endif
          }
        }
        y[i] = z[i];
      } else { // direction is not dissipative, we use C^1 rough enclosure algorithm, see [ZLo]
#if __DEBUG_ENCLOSURE__
        out<<"\ndirection "<<F.array2mode(i)<<" is NOT dissipative(lambda="<<F.lambda(i)<<", ni="<<F.ni(i)<<").\n";
#endif
        y[i] = x[i] + h * fz[i];
#if __DEBUG_ENCLOSURE__
        out<<"z[i]="<<z[i]<<"\n";
        out<<"fz[i]="<<fz[i]<<"\n";
        out<<"y[i]="<<y[i]<<"\n";
#endif
        if(!(y[i].subsetInterior(z[i]))) {
#if __DEBUG_ENCLOSURE__
          out<<"!(y[i].subsetInterior(z[i]))\n";
#endif
          validated[i] = false;
          allValidated = false;
          sm = intervalHull(z[i], y[i]);
          z[i] = inflate<FadMapT> (sm, c);
#if __DEBUG_ENCLOSURE__
          out<<"new z[i]="<<z[i]<<"\n";
#endif
        }
      }
    }
  }
  if(!allValidated) {
    for(i = 0; i < dimension; i++)
      if(!validated[i])
        out << "ith index for i=" << F.array2mode(i) << " not validated\n";
//    out << "DEBUGGING ENCLOSURE:\n";
//    debug_enclosure(vField, x, step, out);
    out << "Problem occurred in dissipative_enclosure function. Possible blow up.\n";
    std::cerr << "Problem occurred in dissipative_enclosure function. Possible blow up.\n";
    throw capd::jaco::SetBlowUpException<VectorT>(x, "enclosure");
  }
  int refinementSteps = __REFINEMENT_STEPS__;
#if __DEBUG_ENCLOSURE__
  out<<"calculated candidate for enclosure:\n";
  F.printModes(y, out);
#endif
  VectorT Fy;

  //refinement steps
  for(j = 0; j < refinementSteps; j++) {
    Fy = F(y);
    N = F.N(y);
    for(i = 0; i < dimension; i++) {
      if(F.isDissipative(i)) {
        b[i] = N[i] / -F.lambda(i);
        g[i] = F.g_encl(step, i, x[i], b[i]);
        if(x[i].rightBound() >= b[i].rightBound())
          y[i].setRightBound(x[i].rightBound());
        else
          y[i].setRightBound(g[i].rightBound());
        if(x[i].leftBound() <= b[i].leftBound())
          y[i].setLeftBound(x[i].leftBound());
        else
          y[i].setLeftBound(g[i].leftBound());
      } else {
        y[i] = x[i] + h * Fy[i];
      }
    }
#if __DEBUG_ENCLOSURE__
    out<<"refinement step: "<<j<<"\n";
    F.printModes(y, out);
#endif
  }
#if __VERIFY_ENCLOSURE__
  //checking if found enclosure satisfies the theorem assumption
  //candidate is in y
  N = F.N(z);
  fz = F(z);
  for(i = 0; i < dimension; i++) {
    if(!(F.isDissipative(i))) {
      y[i] = x[i] + h * fz[i];
#if __DEBUG_ENCLOSURE__
      out<<"checking enclosure at "<<F.array2mode(i)<<" yi "<<y[i]<<" subset zi "<<z[i]<<" \n";
#endif
      if(!(y[i].subsetInterior(z[i]))) {
        out << "Problem occurred in dissipative_enclosure function. Possible blow up."
            "at " << F.array2mode(i) << "th coordinate y[i] (" << y[i] << ") subsetInt (" << z[i] << ")\n";
        std::cerr << "Problem occurred in dissipative_enclosure function. Possible blow up.\n"
            "at " << F.array2mode(i) << "th coordinate y[i] (" << y[i] << ") subsetInt (" << z[i] << ")";
        throw capd::jaco::SetBlowUpException<VectorT>(x, "enclosure");
      }
    } else {
      b[i] = N[i] / -F.lambda(i);
      g[i] = F.g_encl(step, i, x[i], b[i]);
#if __DEBUG_ENCLOSURE__
      out<<"checking enclosure at "<<F.array2mode(i)<<" yi="<<y[i]<<" b[i]="<<b[i]<<" g[i]="<<g[i]<<" \n";
#endif
      if(z[i].rightBound() < b[i].rightBound()) {
        if(!(z[i].rightBound() >= g[i].rightBound())) {
          out << "Problem occurred in dissipative_enclosure function. Possible blow up."
              "at " << F.array2mode(i) << "th coordinate z[i] (" << z[i] << " g[i]=" << g[i] << " b[i]=" << b[i]
              << " w[i]=" << w[i] << " x[i]=" << x[i] << "!(z[i].rightBound()(" <<
              z[i].rightBound() << ")>=g[i].rightBound()(" << g[i].rightBound() << "))\n";
          std::cerr << "Problem occurred in dissipative_enclosure function. Possible blow up.\n"
              "at " << F.array2mode(i) << "th coordinate z[i] (" << z[i] << " g[i]=" << g[i] << " b[i]=" << b[i]
              << " w[i]=" << w[i] << " x[i]=" << x[i] <<
              "!(z[i].rightBound()(" << z[i].rightBound() << ")>=g[i].rightBound()(" << g[i].rightBound() << "))\n";
          //out << "DEBUGGING ENCLOSURE:\n";
          //debug_enclosure(vField, x, step, out);
          throw capd::jaco::SetBlowUpException<VectorT>(x, "enclosure");
        }
      }
      if(z[i].leftBound() > b[i].leftBound()) {
        if(!(z[i].leftBound() <= g[i].leftBound())) {
          out << "Problem occurred in dissipative_enclosure function. Possible blow up."
              "at " << F.array2mode(i) << "th coordinate z[i] (" << z[i] << " g[i]=" << g[i] << " b[i]=" << b[i]
              << " w[i]=" << w[i] << " x[i]=" << x[i] <<
              "!(z[i].leftBound()(" << z[i].leftBound() << ")>=g[i].leftBound()(" << g[i].leftBound() << "))\n";
          std::cerr << "Problem occurred in dissipative_enclosure function. Possible blow up.\n"
              "at " << F.array2mode(i) << "th coordinate z[i] (" << z[i] << " g[i]=" << g[i] << " b[i]=" << b[i]
              << " w[i]=" << w[i] << " x[i]=" << x[i] <<
              "!(z[i].leftBound()(" << z[i].leftBound() << ")>=g[i].leftBound()(" << g[i].leftBound() << "))\n";
          //out << "DEBUGGING ENCLOSURE:\n";
          //debug_enclosure(vField, x, step, out);
          throw capd::jaco::SetBlowUpException<VectorT>(x, "enclosure");
        }
      }
    }
  }
#endif
  return y;
}


template <typename FadMapT, typename Tail>
typename FadMapT::VectorType enclosure(FadMapT & vField, const Tail& T, typename FadMapT::MatrixType::RowVectorType const & x,
    typename FadMapT::ScalarType const & step) {
  FadMapT& nvField = vField;
  Tail T0 = nvField.getT0();
  nvField.setT0(T);
  typename FadMapT::VectorType result = capd::jaco::enclosure(nvField, x, step);
  nvField.setT0(T0);
  return result;
}

template <typename MapType, typename NormType>
typename MapType::MatrixType jacEnclosure(const MapType& vectorField, const typename MapType::ScalarType& step,
    const typename MapType::VectorType &enc, const NormType &the_norm, typename MapType::ScalarType* o_logNormOfDerivative)
{
  typedef typename MapType::MatrixType MatrixType;
  typedef typename MapType::ScalarType ScalarType;

  int dimension = enc.dimension();
  MatrixType der = vectorField[enc], result(dimension, dimension);

  ScalarType l = the_norm(der).rightBound(); // computation of lagarithmic norm
  ScalarType s = ScalarType(0, 1) * step;
  ScalarType w = ScalarType(-1, 1) * exp(s * l);

  MatrixType W(dimension, dimension); // W_3 in paper "C^1 - Lohner algorithm"
  W = w;

  result = MatrixType::Identity(dimension) + s * der * W;

  int i, j;
  for(i = 1; i <= dimension; ++i)
    for(j = 1; j <= dimension; ++j) {
      ScalarType d = result(i, j);
      typename ScalarType::BoundType l = (w.leftBound() > d.leftBound() ? w.leftBound() : d.leftBound()), r =
          (w.rightBound() < d.rightBound() ? w.rightBound() : d.rightBound());
      result(i, j) = ScalarType(l, r);
    }
  if(o_logNormOfDerivative)
    *o_logNormOfDerivative = l;
  return result;
}

}
}

#endif // _DISSIPATIVEENCLOSURE_H_ 
