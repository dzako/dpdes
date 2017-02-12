/*
 * DPDE2DDynSys.h
 *
 *  Created on: Nov 3, 2011
 *      Author: cyranka
 */

#ifndef FFTDYNSYS_H_
#define FFTDYNSYS_H_

#include "capd/vectalg/Vector.h"
#include "capd/vectalg/Matrix.h"
#include "capd/dynsys/DynSys.hpp"

#include "dissipative_enclosure.h"
#include "time.h"

namespace capd{
namespace jaco{

class Dummy{
public:
  Dummy(){
    //std::cout << "HRERERERE\n";
  }
};

/**This enum is for selecting the algorithm type that is used. The following three types are avaialble:
 * -Direct evaluation of convolutions,
 * -Evaluation of convolutions using the FFT,
 * -Evaluation of convolutions using the FFT, but the first order normalized derivative is evaluated directly in order to avoid 
 *  the blow-ups (see the paper for description).
 */
enum AlgorithmType{ direct, FFT, FFTButFirstOrderDirect};

/**Calculates only solutions value (or bounds, when intervals are used), without the partial derivative with respect to
 * initial condition.
 */
template <class EquationTaylorT, int ORD>
class FFTBasicDynSys : public EquationTaylorT{
public:
  typedef EquationTaylorT EquationTaylorType;
  typedef typename EquationTaylorType::GridsContainerType GridsContainerType;
  typedef typename EquationTaylorType::ModesContainerContainerType ModesContainerContainerType;
  typedef typename EquationTaylorType::FFTType FFTType;
  typedef typename EquationTaylorType::ScalarType ScalarType; ///<this is interval, this typename is used for enclosure function purpose
  typedef typename EquationTaylorType::ScalarT ScalarT; ///<this is the actual scalartype (jet or complex number)
  typedef typename EquationTaylorType::ComplexScalarType ComplexScalarType;
  typedef typename ComplexScalarType::ScalarType RealType;
  typedef typename FFTType::IndexRangeType IndexRangeType;
  typedef typename FFTType::ModesContainerType ModesContainerType;
  typedef typename ModesContainerType::DPDEContainerType DPDEContainerType;
  typedef typename FFTType::DFTGridType DFTGridType;
  typedef typename FFTType::IndexType IndexType;

  int m;
  int dftPts;
  RealType step;
  int order;
  int algorithmType;
  bool changeStepAdaptively;

  ModesContainerType rhsSeries;///< temporary variable
  DFTGridType rhsFunctionSpace;
  ModesContainerType td;
  ModesContainerType tl;
  ModesContainerContainerType modes;
  GridsContainerType grids; ///<DFT grids calculated for each order
  Dummy dummy1, dummy2;

  /**IMPORTANT: This constructor is used to create a FINITE dimensional integrator (only projection is taken into account)   
   */
  FFTBasicDynSys(int m_, int dftPts_, RealType step_, int order_, RealType pi_, RealType nu_) :
    EquationTaylorType(m_, dftPts_, nu_, pi_, order_), m(m_), dftPts(dftPts_), step(step_), order(order_), changeStepAdaptively(false),
    rhsSeries(m), rhsFunctionSpace(dftPts_), td(m), tl(m), modes(order_ + 1), grids(order_ + 1), dummy1(), dummy2(){

    int i;
    for(i=0; i <= order; ++i){
      modes[i] = ModesContainerType(m);
      grids[i] = DFTGridType(dftPts_);
    }    
  }
  
  /**IMPORTANT: This constructor is used to create a INFINITE dimensional integrator (only projection is taken into account)   
   */
  FFTBasicDynSys(int m_, int dftPts_, int M_, int dftPts2_, RealType step_, int order_, RealType pi_, RealType nu_, bool initializeHigherDFT = true) :
  EquationTaylorType(m_, M_, dftPts_, dftPts2_, nu_, pi_, order_, initializeHigherDFT), m(m_), dftPts(dftPts_), step(step_), order(order_), changeStepAdaptively(false),
  rhsSeries(m), rhsFunctionSpace(dftPts_), td(m), tl(m), modes(order_+1, false), grids(order_+1, false){
    int i;
    for(i=0; i <= order; ++i){
      modes[i] = ModesContainerType(m);
      grids[i] = DFTGridType(dftPts_);
    }
  }

  void setStep(RealType newStep){
    step = newStep;
    changeStepAdaptively = false;
  }

  inline void setInitialCondition(const ModesContainerType& ic) {
    modes[0] = ic;
  }

  template<class VectorType>
  inline void setInitialCondition(VectorType& v){
    modes[0] = v;    
    ScalarT::setContainer(modes[0]);    
  }

  inline void calculateGrid2(int i){
    fft1.fastTransform(modes[i], grids[i]);
  }

  template<class VectorType>
  inline void current(VectorType& r){
    modes[0].projectOntoSubspace();
    r = modes[0];
  }

  template<class VectorType>
  inline void remainder(VectorType& r){
    modes[order].projectOntoSubspace();
    r = modes[order];
    RealType fac = power(step, order);
    r *= fac;
  }

  /** without wrapping effect control*/
  template<class VectorType>
  inline void move( VectorType& in, VectorType& out){
    setInitialCondition(in);    
    doStep();
    current(out);
  }

  RealType computeNextTimeStep(){
    RealType epsilon = power(10, -TypeTraits<RealType>::numberOfDigits()-2),
             minStep = 1e-15,
             maxStep = 0.1,
             optStep = maxStep;

    double coeffNorm = toDouble( rightBound(modes[order].maxOfNorms()) );
//  capd::rounding::DoubleRounding::roundNearest();
    RealType Cstep = exp(log(epsilon / coeffNorm) / (order));
    Cstep = capd::max(Cstep, minStep) ;

    return capd::min(maxStep, Cstep);
  }

  /**Direct FFT method.
   */
  inline void doStepFFT(){
    int i;
    ScalarT::switchToComplexValued();
    ScalarT::switchToGlobalOptimization();              
    calculateGrid2(0);       
    generalDebug2 << "grid["<<0<<"] ("<<", "<<")="<<grids[0]<<"\n";
    generalDebug2 << "modes["<<0<<"] ("<<modes[0].subspaceType<<", "<<modes[0].solutionType<<")="<<modes[0]<<"\n";
    for(i=0; i < order; ++i){       
      rightHandSide( i, grids, modes, rhsSeries, rhsFunctionSpace, (i < order - 1 ? true : false));      
      //this is not needed, because we divide rhsSeries inside rhs function, and then transform to rhsFunctionSpace
      //rhsSeries *= RealType(1) / RealType(i+1);
      //rhsFunctionSpace *= RealType(1) / RealType(i+1);
      modes[i+1] = rhsSeries;
      generalDebug2 << "modes["<<i+1<<"] ("<<modes[i+1].subspaceType<<", "<<modes[i+1].solutionType<<")="<<modes[i+1]<<"\n";
      grids[i+1] = rhsFunctionSpace;
      generalDebug2 << "grid["<<i+1<<"] ("<<")="<<grids[i+1]<<"\n";
    }    
    if( changeStepAdaptively ){
      step = computeNextTimeStep();
      //std::cout << "current step=" << step << "\n";
    }
    evaluatePolynomials(); ///<evaluates polynomials and calculates new modes[0]    
  }

  /**Direct FFT method. First order normalized derivative is calculated directly.
   */
  inline void doStepFFTButFirstOrderDirect(){
    int i;
    ScalarT::switchToComplexValued();    
    ScalarT::switchToLocalOptimization();
    rightHandSide(0, modes, rhsSeries);
    generalDebug2 << "modes["<<0<<"] ("<<modes[0].subspaceType<<", "<<modes[0].solutionType<<")="<<modes[0]<<"\n";
    ScalarT::switchToGlobalOptimization();
    modes[1] = rhsSeries;
    calculateGrids(1, modes, grids);
    generalDebug2 << "modes["<<1<<"] ("<<modes[1].subspaceType<<", "<<modes[1].solutionType<<")="<<modes[1]<<"\n";
    generalDebug2 << "grid["<<0<<"] ("<<grids[0].subspaceType<<", "<<grids[0].solutionType<<")="<<grids[0]<<"\n";
    generalDebug2 << "grid["<<1<<"] ("<<grids[1].subspaceType<<", "<<grids[1].solutionType<<")="<<grids[1]<<"\n";
    for(i=1; i < order; ++i){
      rightHandSide( i, grids, modes, rhsSeries, rhsFunctionSpace, (i < order - 1 ? true : false));
      modes[i+1] = rhsSeries;
      generalDebug2 << "modes["<<i+1<<"] ("<<modes[i+1].subspaceType<<", "<<modes[i+1].solutionType<<")="<<modes[i+1]<<"\n";
      grids[i+1] = rhsFunctionSpace;
      generalDebug2 << "grid["<<i+1<<"] ("<<grids[i+1].subspaceType<<", "<<grids[i+1].solutionType<<")="<<grids[i+1]<<"\n";
    }
    if( changeStepAdaptively ){
      step = computeNextTimeStep();
      //std::cout << "current step=" << step << "\n";
    }
    evaluatePolynomials(); ///<evaluates polynomials and calculates new modes[0]
  }

  /**Calculates convolutions directly, without FFT
   */
  inline void doStepDirect(){
    int i;
    generalDebug2 << "modes["<<0<<"] ("<<modes[0].subspaceType<<", "<<modes[0].solutionType<<")="<<modes[0]<<"\n";
    for(i=0; i < order; ++i){
      ScalarT::switchToLocalOptimization();
      rightHandSide(i, modes, rhsSeries);
      ScalarT::switchToGlobalOptimization();
      //this is not needed, because we divide rhsSeries inside rhs function
      //rhsSeries *= RealType(1) / RealType(i+1);
      modes[i+1] = rhsSeries;
      generalDebug2 << "modes["<<i+1<<"] ("<<modes[i+1].subspaceType<<", "<<modes[i+1].solutionType<<")="<<modes[i+1]<<"\n";
    }
    if( changeStepAdaptively ){
      step = computeNextTimeStep();
      //std::cout << "current step=" << step << "\n";
    }
    evaluatePolynomials();
  }

  inline void doStep(){    
    switch(algorithmType){
      case direct : doStepDirect(); break;
      case FFT : doStepFFT(); break;
      case FFTButFirstOrderDirect : doStepFFTButFirstOrderDirect(); break;
      default : doStepDirect();
    }
  }

  inline void evaluatePolynomials(){
    int i;
    if(modes[0].isRealValued()) ///optimizing thing
      ScalarT::switchToRealValued();
    for(i = order; i > 0; --i){
      if(i == order)
        rhsSeries = modes[i];
      else
        rhsSeries += modes[i];
      rhsSeries *= step;
    }
    rhsSeries += modes[0];
    modes[0] = rhsSeries;
    generalDebug2 << "rhsSeries:\n" << rhsSeries << "\n";

  }

  using EquationTaylorType::calculateGrids;
  using EquationTaylorType::useFFT;
  using EquationTaylorType::L;
  using EquationTaylorType::fft1;
  using EquationTaylorType::mode;
  using EquationTaylorType::modes2arraySize;
  using EquationTaylorType::mode2array;
  using EquationTaylorType::rightHandSide;

};

/**Used for rigorous calculations, calculates solutions value and the partial derivative with respect to initial condition.
 */
template <class TaylorT, class JetTaylorT, int ORD>
class FFTTaylorDynSys : public capd::dynsys::DynSys<typename TaylorT::MatrixType>{
public:
  typedef TaylorT TaylorType;
  typedef typename TaylorType::VectorType VectorType;
  typedef typename TaylorType::PolyBdType PolyBdType;

  typedef capd::dynsys::DynSys<typename TaylorT::MatrixType > DynSysType;
  typedef typename DynSysType::VectorType DynSysVectorType;
  typedef typename TaylorType::FFTType FFTType;
  typedef typename TaylorType::RealType RealType;
  typedef typename FFTType::IndexRangeType IndexRangeType;
  typedef typename FFTType::IndexType IndexType;

  typedef typename DynSysType::MatrixType MatrixType;

  typedef capd::jaco::FFTBasicDynSys<TaylorT, ORD> BasicDynSys;
  typedef capd::jaco::FFTBasicDynSys<JetTaylorT, ORD> JetDynSys;
  typedef typename JetDynSys::ModesContainerType ModesContainerType;
  typedef typename ModesContainerType::DPDEContainerType DPDEContainerType;


  BasicDynSys basicDynSys;
  JetDynSys jetDynSys;

  PolyBdType x;
  VectorType xx;
  clock_t start, end;
  DynSysVectorType W;
  DynSysVectorType forcing;
  DynSysVectorType yc;

  
  /**IMPORTANT: This constructor is used to create a finite dimensional integrator (only projection is taken into account)   
   */
  FFTTaylorDynSys(int m_, int dftPts_, RealType step_, int order_, RealType pi_, RealType nu_) :
    basicDynSys(m_, dftPts_, step_, order_ + 1, pi_, nu_), jetDynSys(m_, dftPts_, step_, order_, pi_, nu_), x(m_), xx(m_),
    forcing(PolyBdType::modes2realArraySizeStatic(m_)), yc(PolyBdType::modes2realArraySizeStatic(m_)){
    jetDynSys.algorithmType = direct;
    basicDynSys.algorithmType = direct;
  }
  
  /**IMPORTANT: This constructor is for the differential inclusion. The infinite part is taken into account (that's why M and dftPts2 are needed)
   */
  FFTTaylorDynSys(int m_, int dftPts_, int M_, int dftPts2_, RealType step_, int order_, RealType pi_, RealType nu_) :
    basicDynSys(m_, dftPts_, M_, dftPts2_, step_, order_ + 1, pi_, nu_, true), jetDynSys(m_, dftPts_, step_, order_, pi_, nu_), x(m_), 
    xx(m_), forcing(PolyBdType::modes2realArraySizeStatic(m_)), yc(PolyBdType::modes2realArraySizeStatic(m_)){
    jetDynSys.algorithmType = direct;
    basicDynSys.algorithmType = direct;
  }
  
  void setJetDynSysAlgorithmType(int algorithmType){
    jetDynSys.algorithmType = algorithmType;
  }
  
  void setBasicDynSysAlgorithmType(int algorithmType){    
    basicDynSys.algorithmType = algorithmType;
    if(algorithmType == FFT)
      basicDynSys.useFFT = true;
    else
      basicDynSys.useFFT = false;
  }
  
  inline void setSecondInitialCondition(const DynSysVectorType& v){
    jetDynSys.modes[0].setSecondFreeCoeffs(v);
  }

  void setInitialConditionSubspace(DPDEContainerType& c){
    basicDynSys.modes[0].setSubspaceType(c);
    jetDynSys.modes[0].setSubspaceType(c);
  }

  DynSysVectorType Phi(const DynSysVectorType &iv){
    DynSysVectorType r(iv.size());
    jetDynSys.current(r);
    return r;
  }
  MatrixType JacPhi(const DynSysVectorType &iv){
    MatrixType r(iv.size(), iv.size());
    currentD(r);
    return r;
  }
  DynSysVectorType Remainder(const DynSysVectorType &iv){
    DynSysVectorType r(iv.size());
    basicDynSys.remainder(r);
    return r;
  }

  
  inline void encloseMap(
      const DynSysVectorType& x,
      const DynSysVectorType& xx,
      DynSysVectorType& o_phi,
      MatrixType& o_jacPhi,
      DynSysVectorType& o_rem
  ){
    //IMPORTANT: CAPD compatible part START
    
    //attention: stupid notation
    //vector xx is in fact the current [x] (values of the solutions of the finite projection) - interval set
    jetDynSys.setInitialCondition(xx);        
    //vector x is in fact the middle point of the current [x]
    //this middle point is assigned by using the "setSecondInitialCondition()" because in the jets 
    //setInitialCondition() sets the value of current interval set [x] - this value 
    //is multiplied by the partial derivatives, when jets are multiplied during calculation of \frac{\partial \Phi}{\partial x}. 
    
    //setSecondInitialCondition() sets the value of the middle point m([x]) of the current [x], values of the middle    
    //point m([x]) are used to calculate \Phi([x]) and are not multiplied by the partial derivatives, 
    //when jets are multiplied during calculation of \frac{\partial \Phi}{\partial x}. 
    setSecondInitialCondition(x);     
    
//    clock_t start = clock();    
    W = enclosure(xx);
//    clock_t end = clock();
//    generalDebug << "enc time: "<<(end-start)<<"\n";    
        
    basicDynSys.setInitialCondition(W);    
//    generalDebug << "enc: " << W << "\n";        
    
    jetDynSys.doStep(); ///calculates Phi and JacPhi (derivative with respect to initial conditions)    
    //below the Taylor remainder is calculated by using basicDynSys (W the enclosure is the initial condition)
    basicDynSys.doStep();
    //below set the internal CAPD representation of the remainder
    basicDynSys.remainder(o_rem);
//    generalDebug << "rem: " << o_rem << "\n";
    //below set the internal CAPD representation of \Phi( m([x]) )
    jetDynSys.current(o_phi);
    //below set the internal CAPD representation of the partial derivative matrix \frac{\partial \Phi([x])}{\partial x} 
    currentD(o_jacPhi);
    
    // CAPD compatible part END
  }

  
  const DynSysVectorType& enclosure(const DynSysVectorType& in){
    x = in;
    xx = capd::jaco::enclosure(basicDynSys, x, jetDynSys.step);
    return xx;
  }

  template<class MapT>
  typename MapT::VectorType enclosure(MapT& map, typename MapT::VectorType& x){
    return capd::jaco::enclosure(map, x, jetDynSys.step);
  }

  inline void currentD(MatrixType& m){
    jetDynSys.modes[0].monodromyMatrix(m);
  }

  inline BasicDynSys& getVectorField(){
    return basicDynSys;
  }
  
  inline JetDynSys& getJetDynSys(){
    return jetDynSys;
  }

  inline const RealType& getStep() const{
    return jetDynSys.step;
  }

  inline void eraseYc(){
    basicDynSys.yc = 0;
    jetDynSys.yc = 0;
  }

  inline void setYc(const VectorType& yc){
    basicDynSys.yc = yc;
    jetDynSys.yc = yc;
  }

  inline void setForcing(const VectorType& f){
    basicDynSys.forcing = f;
    jetDynSys.forcing = f;
  }
  
  void printGridGnuplotFormat(capd::auxil::OutputStream& out){
    
    int i, j;
    for(i=0; i < jetDynSys.grids[0].m; ++i){
      for(j=0; j < jetDynSys.grids[0].m; ++j)
        out << i*2*3.14159265358979323846264338 / double(jetDynSys.grids[0].m) << " " << j*2*3.14159265358979323846264338 / double(jetDynSys.grids[0].m)
        << " " << leftBound((jetDynSys.grids[0])[i][j].value().re) << " " << 
        rightBound((jetDynSys.grids[0])[i][j].value().re) << "\n";/* leftBound(c[i][j].re) << " " << rightBound(c[i][j].re) << "\n";*/
    }
    out << "\n";
  }
  
};

}
}

#endif /* FFTDYNSYS_H_ */
