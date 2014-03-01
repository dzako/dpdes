/*
 * tests.h
 *
 *  Created on: Aug 18, 2011
 *      Author: cyranka
 */

#ifndef TESTS_H_
#define TESTS_H_

#include "config.h"
#include "capd/vectalg/Vector.h"
//#include "capd/filib/Interval.h"
#include "capd/intervals/Interval.hpp"
#include "Index.h"
#include "Real.h"
#include "norms.h"

class IndicesTester{
public:
  typedef capd::intervals::Interval<double> Interval;
  typedef capd::jaco::Index2DTwoComponents Index2DTwoComponents;
  typedef capd::jaco::Index2D Index2D;
  typedef capd::jaco::Index1D Index1D;
  typedef capd::jaco::MaximumNorm<Index1D> MaximumNorm1D;
  typedef capd::jaco::MaximumNorm<Index2D> MaximumNorm2D;
  typedef capd::jaco::MaximumNorm<Index2DTwoComponents> MaximumNorm2DTwoComponents;
  typedef capd::jaco::FourierConvolutionIndexRange<Index1D, typename capd::jaco::MaximumNorm<Index2D> > IndexRange1D;
  typedef capd::jaco::FourierConvolutionIndexRange<Index2D, typename capd::jaco::MaximumNorm<Index2D> > IndexRange2D;
  typedef capd::jaco::FourierConvolutionIndexRange<Index2DTwoComponents, typename capd::jaco::MaximumNorm<Index2DTwoComponents> > IndexRange2DTwoComponents;
  typedef capd::jaco::Real<Interval, Index2DTwoComponents, MaximumNorm2DTwoComponents> Real2DTwoComponents;
  typedef capd::jaco::Real<Interval, Index2D, MaximumNorm2D> Real2D;
  typedef capd::jaco::Real<Interval, Index1D, MaximumNorm1D> Real1D;
  std::ofstream out;
  std::ofstream out2;
  std::ofstream out_2;
  std::ofstream out2_2;

  ///indices in a range
  void test1(int m){
    int i;
    out.open("test1.txt");
    Index2D index;
    IndexRange2D ir(Index2D::zero());
    ir.setRange(0, capd::jaco::weak, m, capd::jaco::weak);
    for(index[0]=-m, index[1]=0; !index.limitReached(ir); index.inc(ir)){
        std::cout<<" index="<<index<<", norm="<<sqrt(index.squareEuclNorm())<<"\n";
        out<<index.k[0]<<" "<<index.k[1]<<"\n";
    }
    out.close();
  }

  ///two ranges of indices
  void test2(int m){
    int i;
    out.open("test1.txt");
    out2.open("test2.txt");
    Index2D index;
    IndexRange2D ir(Index2D::zero()),
                 ir2(Index2D::zero());
    Real2D real(m, 2*m);
    ir.setRange(2, capd::jaco::strong, m, capd::jaco::weak);
    ir2.setRange(m, capd::jaco::strong, 2*m, capd::jaco::weak);
    for(index=real.firstModeIndex(ir); !index.limitReached(ir); index.inc(ir)){
      std::cout<<" index="<<index<<", norm="<<sqrt(index.squareEuclNorm())<<"\n";
      out<<index.k[0]<<" "<<index.k[1]<<"\n";
    }
    for(index=real.firstModeIndex(ir2); !index.limitReached(ir2); index.inc(ir2)){
      std::cout<<"2 index="<<index<<", norm="<<sqrt(index.squareEuclNorm())<<"\n";
      out2<<index.k[0]<<" "<<index.k[1]<<"\n";
    }
    ir.setRange(m, capd::jaco::strong, 2*m, capd::jaco::weak);

    out.close();
    out2.close();
  }

  ///symmetry check
  void test3(int m, int k1, int k2){
    out_2.open("test1_2.txt");
    out2_2.open("test2_2.txt");
    out.open("test1.txt");
    out2.open("test2.txt");
    Index2D k, index;
    k[0]=k1; k[1]=k2;

    IndexRange2D ir(k);

    Real2D real(m, 2*m);
    ir.setRange(0, capd::jaco::weak, m, capd::jaco::weak);

    for(index=real.firstWithinRange(ir); index<(k/2); index.inc(ir)){
      std::cout<<" index="<<index<<", norm="<<sqrt(index.squareEuclNorm())<<"\n";
      out<<index.k[0]<<" "<<index.k[1]<<"\n";
      out2<<ir.k[0]-index.k[0]<<" "<<ir.k[1]-index.k[1]<<"\n";
    }
    for(index=k/2, index.inc(ir); !index.limitReached(ir); index.inc(ir)){
      out_2<<index.k[0]<<" "<<index.k[1]<<"\n";
      out2_2<<ir.k[0]-index.k[0]<<" "<<ir.k[1]-index.k[1]<<"\n";
    }
    index=k/2; index.inc(ir);

    out_2.close();
    out2_2.close();
    out.close();
    out2.close();
  }
  ///symmetry and two ranges
  void test4(int m, int k1, int k2){
    out_2.open("test1_2.txt");
    out2_2.open("test2_2.txt");
    out.open("test1.txt");
    out2.open("test2.txt");
    Index2D k, t, index;
    k[0]=k1; k[1]=k2;
    t = k;

    t = t /2 ;
    t.cast();
    std::cout << "k/2.cast()=" << t << "\n";
    IndexRange2D ir(k), ir2(k);

    Real2D real(m, 2*m);
    ir.setRange(0, capd::jaco::weak, m, capd::jaco::strong);
//   ir2.setKmk_1Range(0, capd::jaco::weak, m, capd::jaco::strong);
    ir2.setRange(m, capd::jaco::weak, 2*m, capd::jaco::weak);
//    for(index=real.firstWithinRange(ir); index<(k/2); index.inc(ir)){
//      std::cout<<" index="<<index<<", norm="<<sqrt(index.squareNorm())<<"\n";
//      out<<index.k[0]<<" "<<index.k[1]<<"\n";
//      out2<<ir.k[0]-index.k[0]<<" "<<ir.k[1]-index.k[1]<<"\n";
//    }
    for(index=k/2, index.inc(ir); !index.limitReached(ir); index.inc(ir)){
      out<<index.k[0]<<" "<<index.k[1]<<"\n";
      out2<<ir.k[0]-index.k[0]<<" "<<ir.k[1]-index.k[1]<<"\n";
    }
    for(index=real.firstWithinRange(ir2); index<k/2; index.inc(ir2)){
      std::cout<<" index="<<index<<", norm="<<sqrt(index.squareEuclNorm())<<"\n";
      out_2<<index.k[0]<<" "<<index.k[1]<<"\n";
      out2_2<<ir2.k[0]-index.k[0]<<" "<<ir2.k[1]-index.k[1]<<"\n";
    }

    out_2.close();
    out2_2.close();
    out.close();
    out2.close();
  }

  /**IMPORTANT: testing mode2array functions - returning index of the given mode in an array
   */
  void test5(int m){
    out.open("test.txt");
    Index2D index;
    IndexRange2D ir(index);
    Real2D real(m, 2*m);
    ir.setRange(0, capd::jaco::strong, m, capd::jaco::weak);
    for(index = real.firstModeIndex(ir); !index.limitReached(ir); index.inc(ir)){
      out<<index<<" "<<index.mode2array(m, 1)<<" "<<index.mode2array(m, 0);
      out<<" "<<Index2D::array2modeIndex(m, index.mode2array(m, 1))<< 
      " "<<Index2D::array2modeIndex(m, index.mode2array(m, 0)) << "\n";
    }

    out.close();
  }

  ///Index2DTwoComponents test
  void test6(int m){
    out.open("test.txt");

    Index2DTwoComponents index;
    IndexRange2DTwoComponents ir(index);
    Real2DTwoComponents real(m, 2*m);
    ir.setRange(0, capd::jaco::strong, m, capd::jaco::weak);
    for(index = real.firstModeIndex(ir); !index.limitReached(ir); index.inc(ir)){
      out<<index<<" "<<index.mode2array(m, 1)<<" "<<index.mode2array(m, 0);
      out<<" "<<Index2DTwoComponents::array2modeIndex(m, index.mode2array(m, 1))<<
      " "<<Index2DTwoComponents::array2modeIndex(m, index.mode2array(m, 0)) << "\n";
    }

    out.close();
  }



};

int main(int argc, char * argv[])
{
  IndicesTester tester;
//  if(argc!=4)
//    std::cerr<<"Wrong nr of params!\n";
//  else
//    tester.test4(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
  
  if(argc!=2)
    std::cerr<<"Wrong nr of params!\n";
  else
    tester.test6(atoi(argv[1]));
  
//  if(argc != 2)
//    std::cerr << "Wrong nr of params!\n";
//  else
//    tester.test5(atoi(argv[1]));
}


#endif /* TESTS_H_ */
