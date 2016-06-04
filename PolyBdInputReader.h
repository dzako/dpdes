#ifndef _CAPD_JACO_POLYBDINPUTREADER_H_
#define _CAPD_JACO_POLYBDINPUTREADER_H_

#include <stdio.h>
#include "DPDEContainer.h"

/* Reads a polynomial bound from input file.
 * Works for rigorous computations only (input modes are intervals)
 */
template< class PolyBdT >
class PolyBdInputReader {
public:
  typedef PolyBdT PolynomialBoundType;
  typedef typename PolynomialBoundType::IndexType IndexType;
  typedef typename PolynomialBoundType::ScalarType ScalarType;
  typedef double BoundType;
  typedef typename PolynomialBoundType::ComplexScalarType ComplexScalarType;

  int m;
  int M;
  capd::jaco::DPDEContainer container;

  PolynomialBoundType polyBd;
  char str[50];

  PolyBdInputReader(const char* fileName){
    FILE* file;
    char str[50];
    if(!(file = fopen(fileName, "r"))) {
      std::cerr << "The specified file does not exists.\n";
      throw std::runtime_error("The specified file does not exists.\n");
    }
    std::cout << "Using the following initial data (filename " << fileName << "): \n";
    if( fscanf(file, "%s\n", str) <= 0 ) {
      std::cerr << "error in input file string\n";
    }
    if( fscanf(file, "m=%d\n", &m) <= 0 ) {
      std::cerr << "error in input file m\n";
    }
    if( fscanf(file, "M=%d\n", &M) <= 0 ) {
      std::cerr << "error in input file M\n";
    }
    fscanf(file, "%s\n", &str);
    if( strstr(str, "DPDEContainer=\n") != 0 ){
      std::cerr << "error in input file\n";
    }
    int solutionType, solutionType2, subspaceType, baseImZero, baseReZero, partialReReZero, partialReImZero, partialImReZero, partialImImZero;

    if( fscanf(file, "[%d, %d, %d]\n", &subspaceType, &solutionType, &solutionType2) <= 0 ){
      std::cerr << "error in input file[%d, %d, %d]\n";
    }
    if( fscanf(file, "(%d, %d)\n", &baseReZero, &baseImZero) <= 0 ){
      std::cerr << "error in input file(%d, %d)\n";
    }
    if( fscanf(file, "%d %d\n", &partialReReZero, &partialReImZero) <= 0 ){
      std::cerr << "error in input file%d %d\n";
    }
    if( fscanf(file, "%d %d\n", &partialImReZero, &partialImImZero) <= 0 ){
      std::cerr << "error in input file%d %d\n";
    }
    container = capd::jaco::DPDEContainer();
    container.subspaceType = subspaceType;
    container.solutionType = solutionType;
    container.solutionType2 = solutionType2;
    container.baseReZero = baseReZero;
    container.baseImZero = baseImZero;
    container.partialReReZero = partialReReZero;
    container.partialReImZero = partialReImZero;
    container.partialImReZero = partialImReZero;
    container.partialImImZero = partialImImZero;

    polyBd = PolynomialBoundType(m, M, container);
    BoundType leftRe, rightRe, leftIm, rightIm;
    for(int i=0; i <= M; i++){
      if( fscanf(file, "([%lg,%lg],[%lg,%lg])\n", &leftRe, &rightRe, &leftIm, &rightIm) <= 0 ){
        std::cerr << i << "error in input file ([%le,%le],[%le,%le])\n";
      }
      polyBd[IndexType(i)] = ComplexScalarType( ScalarType(leftRe, rightRe), ScalarType(leftIm, rightIm) );
    }
    //fscanf(file, "%s\n", &str);
    //if( strstr(str, "far_tail=\n") != 0 ){
    //  std::cerr << "error in input file far_tail=\n";
    //}
    BoundType Cleft, Cright;
    int s;
    if( fscanf(file, "[%le,%le]\n", &Cleft, &Cright) <= 0 ){
      std::cerr << "error in input file\n";
    }
    if( fscanf(file, "%d\n", &s) <= 0 ){
      std::cerr << "error in input file\n";
    }
    setC( polyBd, ScalarType(Cleft, Cright) );
    setS( polyBd, s );
  }


};


#endif
