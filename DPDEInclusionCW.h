/*
 * DPDEInclusionCW3.h
 *
 *  Created on: Mar 10, 2011
 *      Author: cyranka
 */

#ifndef DPDEINCLUSIONCW3_H_
#define DPDEINCLUSIONCW3_H_

#include "DiffInclusion.hpp"


namespace capd{
namespace jaco{

/**Differential Inclusion class. Implements algorithms for solving differential inclusions, see for example
 *
 * T. Kapela, P. Zgliczy\'nski, A Lohner-type algorithm for control systems and ordinary differential inclusions,
 * Discrete Cont. Dyn. Sys. B, vol. 11(2009), 365-385.
 *
 */
template<typename MapT, typename DynSysT = capd::dynsys::Taylor< typename MapT::MapType> >
class DPDEInclusionCW : public capd::diffIncl2::DiffInclusion<MapT, DynSysT>{

public:
  typedef capd::diffIncl2::DiffInclusion<MapT, DynSysT> BaseClass;
  typedef typename BaseClass::MultiMapType   MultiMapType;
//  typedef typename BaseClass::MapType        MapType;
//  typedef typename BaseClass::FunctionType   FunctionType;
  typedef typename BaseClass::MatrixType     MatrixType;
  typedef typename BaseClass::VectorType     VectorType;
  typedef typename BaseClass::ScalarType     ScalarType;
  typedef typename BaseClass::NormType       NormType;
  typedef typename BaseClass::ParamType      ParamType;
  typedef DynSysT                            DynSysType;

  ParamType m_W_2; ///<stores an enclosure for the Tail, in order to use it during different steps
  ParamType m_N; ///<auxillary object, stored, in order to avoid redundant calculations
  ParamType m_g; ///<auxillary object, stored, in order to avoid redundant calculations
  ScalarType m_errorTolerance;   ///< error tolerance in the exp(A) computations

  /// Constructor
  DPDEInclusionCW(
      MultiMapType &     diffInclusion,   ///< map of the form f(x)+ g(x,e)
      int                order,           ///< order of Taylor method for integration
      ScalarType const & step,            ///< time step for integration
      NormType const &   norm,            ///< norm used in computation of perturbations
      ScalarType const & expErrorTolerance =  capd::TypeTraits<ScalarType>::epsilon()   ///< error tolerance in the exp(A) computations
  );

  ///09.05.12 Jaco: I have added this constructor, because it better fits my idea of integrator
  DPDEInclusionCW(
      int m,
      int dftPts,
      int M, 
      int dftPts2,
      const ScalarType & pi,
      const ScalarType & nu,
      int                order,           ///< order of Taylor method for integration
      ScalarType const & step,            ///< time step for integration
      NormType const &   norm,            ///< norm used in computation of perturbations
      ScalarType const & expErrorTolerance =  capd::TypeTraits<ScalarType>::epsilon()   ///< error tolerance in the exp(A) computations
  );
  
  ///calculates the following enclosures, [W_2] and T: T([0,h])\subset T as in the description of Algorithm 1.
  ///Also the nonlinear part N_k(x+T) calculated for the far modes $|k|>m$ is stored in the m_N field. T is stored
  ///in the m_W_2 field.
  void enclosureWithTail(ParamType& x_0){
    int j, i;    

    ScalarType bk, g;
    ParamType bt, gt, T0h; ///TODO: at each step these objects are created, this is redundant
    bool validated = false;
    int stepNr = 0;
    bool validate = true;
    bool farTailValidate = false;
    bool changedLastStep = false;
    ScalarType step = getStep();

    if(!(s(x_0) >= m_diffIncl.m_sufficientlyLarge)) {
      tailDebug << "The assumption regarding dissipative PDE is not satisfied, i.e. s(T0)=" << s(x_0) << "<" << m_diffIncl.m_sufficientlyLarge << "=d+p+1.\n";
      std::cerr << "The assumption regarding dissipative PDE is not satisfied, i.e. s(T0)=" << s(x_0) << "<" << m_diffIncl.m_sufficientlyLarge
          << "=d+r+1.\n";
      throw std::runtime_error("The assumption regarding dissipative PDE is not satisfied, i.e. s(T0)<d+p+1.\n");
    }        
    m_W_2 = x_0;// = this->template startT<VectorType> (x, T0, m_W_2, step, changedLastStep);
    m_N = x_0;
//    clock_t start = clock();    
    m_diffIncl.findS(step, x_0, m_W_2, changedLastStep);
//    clock_t end = clock();
//    std::cout << "findS time: " << end-start << "\n";
    ScalarType previousdL = -1;
    int previousL = -1;
    bool firstStep = true;
    bool recalculateN = true;
    ScalarType guessedC = 0;
    //validation loop
    while(!validated && stepNr++ < __MAX_VALIDATE_STEPS__) {
      //T=this->template guessT<Scalar, VectorType>(x, T0);
      tailDebug<<"current W_2:\n"; tailDebug << m_W_2;
      tailDebug<<"!validated\n";
      m_W_2.copyFinitePartFrom(x_0);
      m_W_2 = diffInclusionEnclosure(m_W_2);
      validated = m_diffIncl.validateT(step, x_0, m_W_2, validate, farTailValidate, guessedC, previousL, previousdL, changedLastStep, firstStep, recalculateN);
    }
    std::cout << "validating steps nr=" << stepNr << "\n";
//    ::TOTAL_ITER += stepNr;

//    std::cout << "validating steps nr=" << stepNr << "\n";

    if(!validated) {
      std::cerr << "Failed to validate a tail, number of attempts exceeded threshold value" << __MAX_VALIDATE_STEPS__ << ".\n" <<
          "Probably the set has grown too large.\n For the usage refer the documentation, for the examples refer Section 7 in the paper."<<"\n";
      throw std::runtime_error("Failed to validate a tail, number of attempts exceeded threshold value.\n Probably the set has grown too large.\n For the usage refer the documentation, for the examples refer Section 7 in the paper.");
    }

    int first = x_0.finiteTailBegin(),
        last = x_0.finiteTailEnd();
    //  ONLY ONE REFINEMENT STEP!!! NEVER CHANGE IT!!!
    // WHEN THIS STEP IS NOT PERFORMED T_0 \SUBSET W_2 MAY NOT BE SATISFIED
    for(j = 0; j < 1; j++) {
      tailDebug<<"refinement step: "<<j<<"\n"<<"calculated W_2=\n"<<m_W_2<<"\n";

      //values are copied from equation (was calculated in last validateT step
      m_N = m_diffIncl.m_Nt;
      //m_diffIncl.N(m_W_2, m_N);
      m_diffIncl.bt(m_N, bt);
      m_diffIncl.gt(step, x_0, m_N, gt);
      for(i = first; i <= last; ++i) {
        tailDebug<<"i= "<<i<<"\n";
        bk = bt[i];
        g = gt[i];
        if(x_0[i].rightBound() >= bk.rightBound()) {
          tailDebug<<" T0.right "<<x_0[i].rightBound()<<">=bk.right "<<bk.rightBound()<<"\n";
          m_W_2.setRightBound(i, x_0[i].rightBound());
          tailDebug<<" T.right=T0.right "<<x_0[i].rightBound()<<"\n";
        } else {
          tailDebug<<" T0.right "<<x_0[i].rightBound()<<"<bk.right "<<bk.rightBound()<<"\n";
          m_W_2.setRightBound( i, g.rightBound());
          tailDebug<<"Re T.right=g.right "<<g.rightBound()<<"\n";
        }
        if(x_0[i].leftBound() <= bk.leftBound()) {
          tailDebug<<"T0.left "<<x_0[i].leftBound()<<">=bk.left "<<bk.leftBound()<<"\n";
          m_W_2.setLeftBound(i, x_0[i].leftBound());
          tailDebug<<"T.left=T0.left "<<x_0[i].leftBound()<<"\n";
        } else {
          tailDebug<<"T0.left "<<x_0[i].leftBound()<<"<bk.left "<<bk.leftBound()<<"\n";
          m_W_2.setLeftBound( i, g.leftBound());
          tailDebug<<"T.left=g.left "<<g.leftBound()<<"\n";
        }
      }
    }

    tailDebug<<"chosen W_2:\n";
    tailDebug << m_W_2;

    ///saving calculated T
    ///saving calculated T
    ///take intersection of the far tail of T with the far tail of T0
    if(!x_0.subsetFar(m_W_2)) {
      setCLarger(m_W_2, C(x_0) * power(m_W_2.M, s(m_W_2) - s(x_0)));
    }

    ///CHECKS IF IN FACT VALIDATED TAIL IS OK, in the sense satisfies T([0,h])\subset T
  #if __VERIFY_TAIL__
    m_W_2.copyFinitePartFrom(x_0);
    m_W_2 = diffInclusionEnclosure(m_W_2);
    m_diffIncl.N(m_W_2, m_N);

    if(!x_0.subset(m_W_2)) {
      tailDebug<<"!x_0.subset(W_2)\n" << x_0 << "\n" << m_W_2 << "\n";
      std::cerr<<"!x_0.subset(W_2)\n";
      throw std::runtime_error("!x_0.subset(W_2)\n");
    }
    T0h = m_W_2;
    m_diffIncl.bt(m_N, bt);
    m_diffIncl.gt(step, x_0, m_N, gt);
    tailDebug<<"Calculated g:\n";
    tailDebug << gt;
    tailDebug<<"Checking if W_2[0, h]subset W_2\n";
    for(i=first; i<=last; ++i) {
      bk=bt[i];
      g=gt[i];
      if(x_0[i].rightBound() >= bk.rightBound()) {
        T0h.setRightBound(i, x_0[i].rightBound());
      } else {
        T0h.setRightBound(i, g.rightBound());
      }
      if(x_0[i].leftBound() <= bk.leftBound()) {
        T0h.setLeftBound(i, x_0[i].leftBound());
      } else {
        T0h.setLeftBound(i, g.leftBound());
      }
      //test
      if(!(T0h[i].subset(m_W_2[i]))) {
        tailDebug<<"tail test failed ! W_2[0,h] "<<T0h[i]<<" subset W_2 "<<m_W_2[i]<<" at coordinate i="<<i<<" re\n";
        if(T0h[i].leftBound() < m_W_2[i].leftBound()) {
          tailDebug<<"T0h(i, 1).leftBound() "<<T0h[i].leftBound()<<" < W_2(i, 1).leftBound() "<<m_W_2[i].leftBound()<<"diff="<<T0h[i].leftBound()-m_W_2[i].leftBound()<<"\n";
        }
        if(T0h[i].rightBound() > m_W_2[i].rightBound()) {
          tailDebug<<"T0h(i, 1).rightBound() "<<T0h[i].rightBound()<<"> W_2(i, 1).rightBound() "<<m_W_2[i].rightBound()<<"diff="<<T0h[i].rightBound()-m_W_2[i].rightBound()<<"\n";
        }
        std::cerr<<"tail test failed ! W_2[0,h] "<<T0h[i]<<" subset W_2 "<<m_W_2[i]<<" at coordinate i="<<i<<" re\n";
        throw std::runtime_error("Tail test failed ! W_2[0,h] subset W_2 \n");
      }
      tailDebug<<"re W_2[0,h]_"<<i<<"="<<T0h[i]<<"im W_2[0,h]_"<<i<<"="<<T0h[i]<<"\n";
    }
//    int L;
    ScalarType L2;
    if(C(x_0)!=0 && C(bt)!=0 && s(bt)-s(x_0) != 0)
      L2 = power(C(bt)/C(x_0), ScalarType(1./(s(bt)-s(x_0)))).leftBound();
    else
      L2=-1;
//    if(C(gt)!=0 && C(m_W_2)!=0 && s(m_W_2)-s(gt) != 0)
//      L = static_cast<int>(power(C(m_W_2)/C(gt), 1./(s(m_W_2)-s(gt))).rightBound())+1;
//    else
//      L = -1;

    //validation of the far tali
    if(s(bt) > s(x_0)) {
      if(L2 > 0) { //L2<\infty
        if(L2 > m_diffIncl.M) {
          if(!(C(m_W_2) >= C(gt) * power(m_diffIncl.M, s(m_W_2)-s(gt)))) {
            tailDebug<<"tail test failed!\nCase s(b)>s(x_0), L2>M+1\nT_{M+1}<g_{M+1}\n";
            std::cerr<<"tail test failed!\nCase s(b)>s(x_0), L2>M+1\nT_{M+1}<g_{M+1}\n";
            throw std::runtime_error("tail test failed! T_{M+1}<g_{M+1}\n");
          }
          if(!(C(m_W_2) >= C(gt) * power(m_diffIncl.M, s(m_W_2)-s(gt)))) {
            tailDebug<<"tail test failed!\nCase s(b)>s(x_0), L2>M+1\nT_L2<g_L2\n";
            std::cerr<<"tail test failed!\nCase s(b)>s(x_0), L2>M+1\nT_L2<g_L2\n";
            throw std::runtime_error("tail test failed! T_L2<g_L2\n");
          }
          if(!(C(m_W_2) >= C(x_0) * power(m_diffIncl.M, s(m_W_2)-s(x_0)))) {
            tailDebug<<"tail test failed!\nCase s(b)>s(x_0), L2>M+1\nC(W_2) "<<C(m_W_2)<<" < C(x_0)*POW "<<C(x_0)*power(L2, s(m_W_2)-s(x_0))<<"\n";
            std::cerr<<"tail test failed!\nCase s(b)>s(x_0), L2>M+1\nC(W_2) "<<C(m_W_2)<<" < C(x_0)*POW "<<C(x_0)*power(L2, s(m_W_2)-s(x_0))<<"\n";
            throw std::runtime_error("tail test failed! C(W_2) < C(x_0)*POW \n");
          }
        } else { //L2<=M+1
          if(s(m_W_2) > s(x_0)) {
            tailDebug<<"Tail test failed!\nCase s(b)>s(x_0)\nUnexpected values of L2=infty and s(x_0)<s(W_2).\n";
            std::cerr<<"Tail test failed!\nCase s(b)>s(x_0)\nUnexpected values of L2=infty and s(x_0)<s(W_2).\n";
            throw std::runtime_error("tail test failed!\nEquation is probably not dissipative.\n");
          }
          if(!(C(m_W_2) >= C(x_0) * power(m_diffIncl.M, s(m_W_2)-s(x_0)))) {
            tailDebug<<"tail test failed!\nCase s(b)>s(x_0), L2<=M+1\nT_{M+1}<T0_{M+1}\n";
            std::cerr<<"tail test failed!\nCase s(b)>s(x_0), L2<=M+1\nT_{M+1}<T0_{M+1}\n";
            throw std::runtime_error("tail test failed!\nCase s(b)>s(x_0), L2<=M+1\nT_{M+1}<T0_{M+1}\n");
          }
        }

      } else { //L2=\infty
        if(s(gt) >= s(m_W_2)) {
          if(!(C(m_W_2) >= C(gt) * power(m_diffIncl.M, s(m_W_2)-s(gt)))) {
            tailDebug<<"tail test failed!\nCase s(b)>s(x_0)\nT_{M+1} < g_{M+1}\n";
            std::cerr<<"tail test failed!\nCase s(b)>s(x_0)\nT_{M+1} < g_{M+1}\n";
            throw std::runtime_error("tail test failed! T_{M+1} < g_{M+1}\n");
          }
        } else {//s(gt)<s(W_2) validation is not possible
          tailDebug<<"Tail test failed!\nCase s(b)>s(x_0)\nUnexpected values of L2=infty and s(g)<s(W_2).\nEquation is probably not dissipative, reason: p<=q.\n";
          std::cerr<<"Tail test failed!\nCase s(b)>s(x_0)\nUnexpected values of L2=infty and s(g)<s(W_2).\nEquation is probably not dissipative, reason: p<=q.\n";
          throw std::runtime_error("tail test failed!\nEquation is probably not dissipative.\n");
        }
      }

    }
    if(s(bt) == s(x_0)) {
      if(C(x_0) >= C(bt)) {
        if(!(C(m_W_2) >= C(x_0) * power(m_diffIncl.M, s(m_W_2)-s(x_0)))) {
          tailDebug<<"tail test failed !\nCase s(b)==s(x_0)\nT_{M+1} < T0_{M+1}\n";
          std::cerr<<"tail test failed !\nCase s(b)==s(x_0)\nT_{M+1} < T0_{M+1}\n";
          throw std::runtime_error("Tail test failed ! T_{M+1} < T0_{M+1}\n");
        }
      } else { //C(b)>C(x_0)
        if(s(gt) < s(m_W_2)) {
          tailDebug<<"Tail test failed!\nCase s(b)==s(x_0)\nUnexpected s(g).\nWith s(g) "<<s(gt)<<" > s(W_2) "<<s(m_W_2)<<" validation is not possible.\n";
          std::cerr<<"Tail test failed!\nCase s(b)==s(x_0)\nUnexpected s(g).\nWith s(g) "<<s(gt)<<" > s(W_2) "<<s(m_W_2)<<" validation is not possible.\n";
          throw std::runtime_error("Tail test failed ! W_2[0,h] subset W_2\nUnexpected s(W_2)!\nWith s(g)>s(W_2) validation is not possible.\n");
        }
        if(!(C(m_W_2) >= C(gt) * power(m_diffIncl.M, s(m_W_2)-s(gt)))) {
          tailDebug<<"tail test failed !\nCase s(b)==s(x_0)\nT_{M+1} "<<C(m_W_2)<<" <g_{M+1} "<<C(gt)*power(m_diffIncl.M, s(m_W_2)-s(gt))<<" \n";
          std::cerr<<"tail test failed !\nCase s(b)==s(x_0)\nT_{M+1} "<<C(m_W_2)<<" <g_{M+1} "<<C(gt)*power(m_diffIncl.M, s(m_W_2)-s(gt))<<" \n";
          //std::cout<<"M(W_2) "<<W_2.getM()<<", M(gt) "<<gt.getM()<<"\n";
          throw std::runtime_error("Tail test failed ! T_{M+1}<g_{M+1}\n");
        }
      }
    }

    if(s(bt) < s(x_0)) {
      if(L2 > 0) { //L2<\infty
        if(L2 > m_diffIncl.M) {
          if(!(C(m_W_2) >= C(x_0) * power(m_diffIncl.M, s(m_W_2)-s(x_0)))) {
            tailDebug<<"tail test failed !\nCase s(b)<s(x_0), L2>M+1\nT_{M+1}<T0_{M+1}\n";
            std::cerr<<"tail test failed !\nCase s(b)<s(x_0), L2>M+1\nT_{M+1}<T0_{M+1}\n";
            throw std::runtime_error("Tail test failed ! T_{M+1}<T0_{M+1}\n");
          }
          if(!(C(m_W_2) >= C(x_0) * power(L2, s(m_W_2)-s(x_0)))) {
            tailDebug<<"tail test failed !\nCase s(b)<s(x_0), L2>M+1\nT_{L2}<T0_{L2}\n";
            std::cerr<<"tail test failed !\nCase s(b)<s(x_0), L2>M+1\nT_{L2}<T0_{L2}\n";
            throw std::runtime_error("Tail test failed ! T_{L2}<T0_{L2}\n");
          }
          if(!(C(m_W_2) >= C(gt) * power(L2, s(m_W_2)-s(gt)))) {
            tailDebug<<"tail test failed!\nCase s(b)>s(g), L2>M+1\nT_{M+1}<g_{M+1}\nT_{L2}<g_{L2}\n";
            std::cerr<<"tail test failed!\nCase s(b)>s(g), L2>M+1\nT_{M+1}<g_{M+1}\nT_{L2}<g_{L2}\n";
            throw std::runtime_error("tail test failed! T_{L2}<g_{L2}\n");
          }
        } else { //L2<=M+1
          if(s(m_W_2) > s(gt)) {
            tailDebug<<"Tail test failed!\nCase s(b)<s(x_0)\nUnexpected values of L2=infty and s(g)<s(W_2).\n";
            std::cerr<<"Tail test failed!\nCase s(b)<s(x_0)\nUnexpected values of L2=infty and s(g)<s(W_2).\n";
            throw std::runtime_error("Tail test failed ! W_2[0,h] subset W_2 \n.");
          }
          if(!(C(m_W_2) >= C(gt) * power(m_diffIncl.M, s(m_W_2)-s(gt)))) {
            tailDebug<<"tail test failed!\nCase s(b)<s(x_0), L2<=M+1\nT_{M+1}<g_{M+1}\n";
            std::cerr<<"tail test failed!\nCase s(b)<s(x_0), L2<=M+1\nT_{M+1}<g_{M+1}\n";
            throw std::runtime_error("tail test failed!\nCase s(b)<s(x_0), L2<=M+1\nT_{M+1}<g_{M+1}\n");
          }
        }
      } else { //L2==\infty
        if(s(x_0) >= s(m_W_2)) {
          if(!(C(m_W_2) >= C(x_0) * power(m_diffIncl.M, s(m_W_2)-s(x_0)))) {
            tailDebug<<"tail test failed !\nCase s(b)<s(x_0)\nT_{M+1}<T0_{M+1}\n";
            std::cerr<<"tail test failed !\nCase s(b)<s(x_0)\nT_{M+1}<T0_{M+1}\n";
            throw std::runtime_error("Tail test failed ! T_{M+1}<T0_{M+1}\n");
          }
        } else {
          tailDebug<<"Tail test failed!\nCase s(b)<s(x_0)\nUnexpected values of L2=infty and s(x_0)<s(W_2).\n";
          std::cerr<<"Tail test failed!\nCase s(b)<s(x_0)\nUnexpected values of L2=infty and s(x_0)<s(W_2).\n";
          throw std::runtime_error("Tail test failed ! W_2[0,h] subset W_2 \n.");
        }
      }
    }
  #endif

  }


  ///07.08.2013 temporary printing C++ code that creates this matrix for the purpose of test program
  void printMatrix(const MatrixType& J, capd::auxil::OutputStream& out) const{
    int rows = J.numberOfRows(),
        cols = J.numberOfColumns();

    for(int i = 0; i < rows; i++){
      for(int j = 0; j < cols; j++){
        bool dot = (J[i][j] == 0. ? true : false);
        out << "J[" << i << "][" << j << "]=" << "Scalar(" << J[i][j].leftBound() << (dot ? "." : "") << ", " << J[i][j].rightBound() << (dot ? "." : "") << ");";
      }
      out << "\n";
    }

  }

  ///07.08.2013 temporary printing C++ code that creates this vector for the purpose of test program
  void printVector(const VectorType& C, capd::auxil::OutputStream& out) const{
    int size = C.size();
    for(int i = 0; i < size; i++){
      bool dot = (C[i] == 0. ? true : false);
      out << "C[" << i << "]=" << "Scalar(" << C[i].leftBound() << (dot ? "." : "") << ", " << C[i].rightBound() << (dot ? "." : "") << ");\n";
    }
  }

  ///calculates the [\Delta] in the step 8 of Algorithm 1.
  ///Input
  ///-x, as in the description of Algorithm 1,
  ///-T, the enclosure (candidate) for the T([0,h]),
  ///Output
  ///-Nt, is N_k(x+T),
  ///-W_2, the enclosure [W_2] in the step 1 of Algorithm 1,
  ///-y_c, y_c in the step 3 of Algorithm 1.
  VectorType perturbations(const VectorType & x, ParamType & x_0){
//    COUNT_OPERATIONS = true;
    copyFinitePart(x, x_0);
    generalDebug << "current x_0:\n" << x_0;
//    clock_t start, end;
    VectorType y_c;

//    start = clock();
    ///in blocks {} we collect all the steps of Algorithm 1, with reference to corresponding step number
    getDynamicalSystem().eraseYc();


    {///step 1
      enclosureWithTail(x_0); ///returns enclosures [W_2], T: T([0,h])\subset T, N: N_k(x+T)
    }///end of step 1

    generalDebug<<"current T:\n"; generalDebug << m_W_2;
//    end = clock();
//    std::cout << "enclosureWithTail time: " << end-start << "\n";

    VectorType galerkinProjectionError, delta, C, D, result;
    {///step 2
      m_diffIncl.perturbations(m_W_2, galerkinProjectionError);
    }///end of step 2

    generalDebug << "found perturbations:\n"; m_diffIncl.printModes(galerkinProjectionError, generalDebug);

    {///step 3
    y_c=midVector(galerkinProjectionError);
    }///end of step 3

    {///step 5
    delta=y_c-galerkinProjectionError;
    }///end of step 5
    generalDebug<<"found mid(delta):\n"; m_diffIncl.printModes(y_c, generalDebug);
    generalDebug<<"found delta:\n"; m_diffIncl.printModes(delta, generalDebug);

    ///Component-wise estimates step 6-8
    {///step 6
    C = rightVector(abs(delta));
    }///end of step 6    

    inclDebug << "found C:\n" << C << "\n";
    printVector(C, inclDebug);

 //   COUNT_OPERATIONS = true;

    int i, j;
    MatrixType J;
    {///step 7
    J = (getDynamicalSystem().getJetDynSys()).jacobian(m_W_2); /// calculating Jacobian on W_2\times y_c
    for( i=0; i < J.numberOfRows(); ++i){
      for( j = 0; j < J.numberOfColumns(); ++j){
        if(i != j)
          J[i][j] = (abs(J[i][j])).right();
        else
          J[i][j] = (J[i][j]).right();
      }
    }
    }///end of step 7
    inclDebug << "found J:\n" << J << "\n";
    printMatrix(J, inclDebug);
    getDynamicalSystem().setYc(y_c); ///setting y_c for the next step, when the solution of the Cauchy problem with y_c
                                           ///is calculated.    
    {///step 8
    MatrixType At = getStep() * J;
    MatrixType A = MatrixType::Identity(J.numberOfRows());    
    MatrixType Sum = A;
    generalDebug << "At:"<<At<<"\nA:"<<A<<"\n";
    generalDebug << "norm:"<<right((*m_norm)(y_c))<<"\n";
    ScalarType AtNorm = right((*m_norm)(At)),
        AnNorm = right((*m_norm)(A));
    ScalarType n = 2.0;  // n = i + 2
    // remainder = |A| * |At/(N+2)| / (1 - |At/(N+2)|)   (the sum of geometric series from N to infinity)
    ScalarType q = AtNorm/n;
    int Aiterations = 0;
    while(true){
      Aiterations ++;
      A = A * At / n;
      Sum += A;
      AnNorm *= q;
      n += ScalarType(1.0);
      q = AtNorm / n;
      if(q < 1){
        // remainder = |A| * |At/(N+2)| / (1 - |At/(N+2)|)   (the sum of geometric series from N to infinity)
        ScalarType remainder = right(AnNorm * AtNorm / (n - AtNorm ));
        if(remainder < this->m_errorTolerance)
          break;
      }
    }
    std::cout << "Aiterations=" << Aiterations << "\n";
    // we recompute remainder because norm of A can be smaller than approximation : AnNorm
    ScalarType remainder = right((*m_norm)(A) * AtNorm / (n  - AtNorm ));    
    for(i=0; i < J.numberOfRows(); ++i)
      for(int j = 0; j < J.numberOfColumns(); ++j)
        Sum[i][j] += remainder * ScalarType(-1.0, 1.0);
    VectorType D = getStep() * (Sum * C);
    result=VectorType(D.dimension());    
    for(int i=0; i< D.dimension(); ++i)
      result[i] = ScalarType(-D[i].rightBound(), D[i].rightBound());

    }///end of step 8
    //COUNT_OPERATIONS = false;

    generalDebug<<"found Deltha:\n";
    inclDebug << "found Deltha:\n" << result << "\n";
    m_diffIncl.printModes(result, generalDebug);
//    COUNT_OPERATIONS = false;
    return result;
  }

  ///calculates the enclosure [W_1], \varphi([0,h], [x_0], y_c)\subset [W_1].
  VectorType dynamicalSystemEnclosure(const VectorType& x);

  ///calculates enclosure, [W_2] in the description of Algorithm 1, using provided enclosure (candidate)
  ///for the tail T.
  ParamType diffInclusionEnclosure(const ParamType& x) const{
    return capd::jaco::enclosure(m_diffIncl, x, this->getStep());
  }

  ///calculates T(h), last step in Algorithm 1.
  void moveParams(const VectorType& x, ParamType & x_0){
    ScalarType step = getStep();    
    m_diffIncl.gt(step, x_0, m_N, m_g);

    copyFinitePart(x, m_g);//here we copy the projection part to m_g object, which stores the current set
    x_0 = m_g;    
  }

  using BaseClass::enclosure;
  using BaseClass::getStep;
  using BaseClass::getDynamicalSystem;

protected:
  using BaseClass::m_norm;
  using BaseClass::m_diffIncl;

};

template <typename MapT, typename DynSysT>
DPDEInclusionCW<MapT, DynSysT>::DPDEInclusionCW(MultiMapType &     diffInclusion, int order, ScalarType const & step,            ///< time step for integration
    NormType const &   norm, ScalarType const & expErrorTolerance) :
    BaseClass(diffInclusion, order, step, norm){
  m_errorTolerance = expErrorTolerance;
}

template <typename MapT, typename DynSysT>
DPDEInclusionCW<MapT, DynSysT>::DPDEInclusionCW(int m, int dftPts, int M, int dftPts2, const ScalarType & pi, const ScalarType & nu, 
    int order, ScalarType const & step, NormType const &   norm, ScalarType const & expErrorTolerance) :
    BaseClass(m, dftPts, M, dftPts2, pi, nu, order, step, norm){
  m_errorTolerance = expErrorTolerance;
}

template <typename MapT, typename DynSysT>
inline typename DPDEInclusionCW<MapT, DynSysT>::VectorType DPDEInclusionCW<MapT, DynSysT>::dynamicalSystemEnclosure(const VectorType& x){
  return getDynamicalSystem().enclosure(x);
}



}}

#endif /* DPDEInclusionCW3_H_ */
