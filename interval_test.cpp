#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <cstring>

//#include "capd/filib/Interval.h"
#include "capd/intervals/Interval.hpp"

#define DOUBLE double
  
//typedef capd::filib::Interval<DOUBLE> Interval;
typedef capd::intervals::Interval<DOUBLE> Interval;
#define PI Interval::pi()

  

int main(int argc, char * argv[]){
  Interval int1 = Interval(1);
  Interval int3 = Interval(3);
  
  Interval r = int1 / int3;
  
  std::cout << "r: " << r << ", inf==sup: " << (leftBound(r) == rightBound(r)) << "\n"; 
  
}