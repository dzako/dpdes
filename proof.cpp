///this file contains the main program

#include "capd/intervals/DoubleInterval.h"
//#include "capd/filib/Interval.h"
#include "capd/dynsys/BasicFadTaylor.h"
#include "capd/dynsys/FadTaylor.h"
#include "DPDETaylor.h"
#include "DPDEPerturbation.h"
#include "DPDEInclusionCW.h"
#include "DPDEMultiMap.h"
#include "Algorithms2.h"
#include "Algorithms2TailManager.h"
#include "Burgers.h"
#include "BasinOfAttraction.h"

#include "Integrator.h"//temp
//new versions
#include "DPDEPerturbation2.h"
#include "DPDEInclusionCW2.h"
#include "DPDEMultiMap2.h"
#include "Burgers2.h"
#include "Real.h"
#include "Index.h"
#include "Algorithms2_2.h"
#include "Algorithms2TailManager2.h"

#include "IndicesTester.h"

#define DOUBLE double
#if __FILIB__
  typedef capd::filib::Interval<DOUBLE> Interval;
#else
  #define Interval interval
#endif

//===================dimension and conditions dependent classes and tails===========================
typedef capd::vectalg::Matrix<Interval,_D,_D> IntervalMatrix;
typedef capd::vectalg::Vector<Interval, _D> IntervalVector;
typedef capd::vectalg::EuclNorm<IntervalVector, IntervalMatrix> EuclNorm;

typedef capd::jaco::D1Real<Interval> Jaco1Real;
typedef capd::jaco::Index1D Index1D;
typedef capd::jaco::Real<Interval, Index1D, EuclNorm> Real;

typedef capd::jaco::PolyBd<Interval, Jaco1Real, _D> PolyBd1Real;
typedef capd::jaco::PolyBd2<Interval, Real, _D> PolyBd1Real2;

typedef capd::jaco::DPDEInclRect2Set<IntervalMatrix, PolyBd1Real> DPDEInclRect2Set;
//define specific dPDE
typedef capd::jaco::Burgers<Jaco1Real, _D> Burgers1PDE;

typedef capd::jaco::Burgers2<Real, _D> Burgers1PDE2;

//chose one dPDE and provide it as a template parameter to an algorithm as a template parameter
typedef capd::jaco::Algorithms2TailManager<PolyBd1Real, Burgers1PDE> Algorithms2TailManager;

typedef capd::jaco::Algorithms2TailManager2<PolyBd1Real2, Burgers1PDE2> Algorithms2TailManager2;

typedef capd::jaco::DPDE<Interval, Algorithms2TailManager, _D> Burgers1Real;

typedef capd::jaco::DPDE<Interval, Algorithms2TailManager2, _D> Burgers1Real2;

typedef capd::jaco::DPDEPerturbation<Interval, Algorithms2TailManager, _D> DPDEPerturbation;

typedef capd::jaco::DPDEPerturbation<Interval, Algorithms2TailManager2, _D> DPDEPerturbation2;

typedef capd::jaco::DPDETaylor<Burgers1Real> DPDETaylor;

typedef capd::jaco::DPDETaylor<Burgers1Real2> DPDETaylor2;

typedef capd::jaco::DPDEMultiMap<Burgers1Real, DPDEPerturbation > DPDEMultiMap;
typedef capd::jaco::DPDEInclusionCW<DPDEMultiMap, PolyBd1Real, DPDETaylor> DPDEInclusionCW;
typedef capd::jaco::BasinOfAttraction<DPDEMultiMap, DPDEInclusionCW, DPDEInclRect2Set, EuclNorm, _D> BasinOfAttraction;

typedef capd::jaco::DPDEMultiMap2<Burgers1Real2, DPDEPerturbation2 > DPDEMultiMap2;
typedef capd::jaco::DPDEInclusionCW2<DPDEMultiMap2, DPDETaylor2> DPDEInclusionCW2;
typedef capd::diffIncl2::InclRect2Set<IntervalMatrix, PolyBd1Real2> InclRect2Set;

//typedef capd::jaco::Sine<DPDEMultiMap2, DPDEInclusionCW2, InclRect2Set, EuclNorm> Sine;
//typedef capd::jaco::TestCase<DPDEMultiMap2, DPDEInclusionCW2, InclRect2Set, EuclNorm> TestCase;

typedef capd::jaco::Sine<DPDEMultiMap, DPDEInclusionCW, DPDEInclRect2Set, EuclNorm> Sine;
typedef capd::jaco::TestCase<DPDEMultiMap, DPDEInclusionCW, DPDEInclRect2Set, EuclNorm> TestCase;

void printInfo(){
  std::cout<<"Wrong number of parameters.\nUsage:\n capd file_name, file_name indicates name of the file with input data. "<<
      "Refer the documentation for details.\n";
}

void burgers1Experiment6(){
  Interval i(-2, 1);
  std::cout<<capd::abs(i)<<"\n";

  //TestCase testCase(25, 102, 0.2, 0, 1);
  Sine testCase(15, 41, 1, 0, 0.05);
  //capd::jaco::SineCh testCase(_m, _M, _ni, 0, __END_TIME_INTERVAL__);
  //capd::jaco::POForKS testCase(14, 42, 0.1212, 0, __END_TIME_INTERVAL__, 6., 0.0001);
  //capd::jaco::SetForBurgers testCase(30, 61, 0.01, 0, __END_TIME_INTERVAL__, 6., 0.001);
  //capd::jaco::FPForBurgers2 testCase(30, 61, 0.1, 0, __END_TIME_INTERVAL__, 6., 0.001);

  const int steps=testCase.getSteps();
  std::ofstream f;
  f.open("test.txt");
  f<<"steps="<<steps<<"\n";
  clock_t start = clock();
  int j=0;
  while(j<steps)
  {
    testCase.doStep(j, f);
    ++j;
  }
  clock_t end=clock();
  testCase.finish((double)(end-start)/CLOCKS_PER_SEC);
}


int main(int argc, char * argv[])
{
  setLoggers();
  char *fileName, *basinInFileName;
   try
   {
       if(argc!=2 && argc!=1){
         printInfo();
       }else{
         if(argc==1){
           burgers1Experiment6();
         }
         if(argc==2){
           fileName=argv[1];
           BasinOfAttraction b(fileName);
           b.mainAlgorithm();
           std::cout<<"Program completed successfully, see output files.\n";
         }
         if(argc==3){
           fileName=argv[1];
           basinInFileName=argv[2];
           BasinOfAttraction b(fileName, basinInFileName);
           b.mainAlgorithm();
         }
     }

   }catch(std::exception& e)
   {
      std::ofstream plik;
      plik.open("report");
      plik << e.what();
      plik.close();
      std::cout << "\n\nException caught! See 'report' file for details.";
   }
   return 0;
}
