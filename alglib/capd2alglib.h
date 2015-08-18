#ifndef _CAPD_ALGLIB_CAPD2ALGLIB_H_
#define _CAPD_ALGLIB_CAPD2ALGLIB_H_

#include <stdexcept>
#include <sstream>
#include "ap.h"
#include "ap.cpp"
#include "apvt.h"
#include "blas.h"
#include "blas.cpp"
#include "nsevd.h"
#include "nsevd.cpp"
#include "hessenberg.h"
#include "hessenberg.cpp"
#include "hsschur.h"
#include "hsschur.cpp"
#include "reflections.h"
#include "reflections.cpp"
#include "rotations.h"
#include "rotations.cpp"
  

/**
 * Function computes Eigenvalues of a general matrix 
 * 
 * The algorithm finds eigenvalues of a general matrix 
 * by using the QR algorithm with multiple shifts.
 * @remark Nonrigorous! Works only for double matrices!
 * 
 * Input:
 * @param[in] A  square matrix
 * 
 * Output:
 * @param[out] eigenRealPart  vector of real parts of eigenvalues.
 * @param[out] eigenImPart    vector of imaginary parts of eigenvalues.
 */ 
template <typename MatrixT, typename VectorT>
void computeEigenvalues(const MatrixT & A,
                          VectorT &eigenRealPart, VectorT & eigenImPart){

  typedef int size_type;

  if(A.numberOfRows() != A.numberOfColumns())
     throw std::invalid_argument("calculateEigenvalues works only for square matrices"); 
   ap::real_2d_array a;
   size_type n = A.numberOfRows();
   a.setbounds(0,n-1, 0, n-1) ;
   for(size_type i =0; i < n; i++)
     for(size_type j=0; j< n; j++){
   	   a(i,j) = A[i][j];                   	
     }
     
     int vneeded = 0;
     ap::real_1d_array wr;
     ap::real_1d_array wi;
     ap::real_2d_array vl;
     ap::real_2d_array vr;
   if(! rmatrixevd(a, n, vneeded, wr, wi, vl, vr))
     throw std::runtime_error("algorithm for eigenvalues computation did not converge!");
   if((eigenRealPart.dimension() < n) or (eigenImPart.dimension() < n))
     throw std::invalid_argument("computeEigenvalues: cannot store eigenvalues (vectors are too small)");
   for(size_type i =0; i < n; i++){
     eigenRealPart[i] = wr(i);
     eigenImPart[i] = wi(i);
   } 
}

/**
 * Computes Eigenvalues and corresponding right Eigenvectors
 * 
 * The algorithm finds eigenvalues and eigenvectors of 
 * a general matrix by using the QR algorithm with multiple shifts.
 * @remark Works only for double matrices nad vectors!
 * 
 * Input:
 * @param[in] A  square matrix
 *  
 * Output:
 * @param[out] eigenRealPart  vector of real parts of eigenvalues.
 * @param[out] eigenImPart    vector of imaginary parts of eigenvalues.
 * @param[out] vectorRealPart matrix whose i-th column is real part of 
 *                            eigenvector corresponding to i-th eigenvalue
 * @param[out] vectorImPart   matrix whose i-th column is imaginary part of 
 *                            eigenvector corresponding to i-th eigenvalue
 */  
template <typename MatrixT, typename VectorT>
void computeEigenvaluesAndEigenvectors(
  const MatrixT & A, 
  VectorT &eigenRealPart,  VectorT &eigenImPart,
  MatrixT &vectorRealPart, MatrixT &vectorImPart){
 
  typedef typename VectorT::ScalarType VectorScalar;
  typedef typename MatrixT::ScalarType MatrixScalar;
  typedef int size_type;


  if(A.numberOfRows() != A.numberOfColumns())
    throw std::invalid_argument("calculateEigenvalues works only for square matrices"); 
  ap::real_2d_array a;
  size_type n = A.numberOfRows();
  a.setbounds(0,n-1, 0, n-1) ;
  for(size_type i =0; i < n; i++)
    for(size_type j=0; j< n; j++){
      a(i,j) = toDouble(A[i][j]);
    }
  int vneeded = 3;
  ap::real_1d_array wr;
  ap::real_1d_array wi;
  ap::real_2d_array vl;
  ap::real_2d_array vr;
  if(! rmatrixevd(a, n, vneeded, wr, wi, vl, vr))
    throw std::runtime_error("algorithm for eigenvalues computation did not converge!");
  if((eigenRealPart.dimension() < n) or (eigenImPart.dimension() < n))
     throw std::invalid_argument("computeEigenvalues: cannot store eigenvalues (vectors are too small)");
  for(size_type i =0; i < n; i++){
    eigenRealPart[i] = VectorScalar(wr(i));
    eigenImPart[i] = VectorScalar(wi(i));
    
    if(wi(i) == 0){   // real eigenvalue
      for(size_type j=0; j< n; j++){
        vectorRealPart[j][i] = MatrixScalar(vr(j,i));
        vectorImPart[j][i] = MatrixScalar(0.0);
      }
    }else if(wi(i)>0){ // complex eigenvalue 
      for(size_type j=0; j< n; j++){
        vectorRealPart[j][i] = MatrixScalar(vr(j,i));
        vectorImPart[j][i] = MatrixScalar(vr(j,i+1));
      }
    }else{             // complex conjugate eigenvalue
      for(size_type j=0; j< n; j++){
        vectorRealPart[j][i] = MatrixScalar(vr(j,i-1));
        vectorImPart[j][i] = MatrixScalar(-vr(j,i));
      }
    }
  } 
}

/** 
 *  Writes complex number to string in human readable form :)
 */ 
template <typename T>
std::string complexToString(const T & re, const T &im){
  std::ostringstream str;
  if(im !=0){
    if(re != 0)
      str << re;
    if(im > 0)
      str <<"+";
    str << im << "i";
  }else{
    str << re;
  } 
  return str.str();
}

/**
 * Converts vector of eigenvalues to text (string)
 * 
 * @param[in] re        vector of real parts of eigenvalues
 * @param[in] im        vector of imaginary parts of eigenvalues
 * @param[in] separator string that separates eigenvalues
 * 
 * @return    eigenvalues in text format
 */ 
template <typename VectorType>
std::string eigenvaluesToString(const VectorType & re, const VectorType &im,
                                std::string separator = ", "){
  //typedef typename VectorType::size_type size_type;
  const int dim = re.dimension();
  std::ostringstream str;
  str << "(";
  for(int i=0; i< dim; i++){
    str << complexToString(re[i], im[i]);
    if(i != dim-1)
      str << separator;
  } 
  str << ")";
  return str.str();
}
/**
 * Converts eigenvectors stored as columns of matrices to text (string)
 *
 * @param[in] re        matrix of real parts of eigenvectors (as its columns)
 * @param[in] im        matrix of imaginary parts of eigenvectors
 * @param[in] separator string that separates vectors coeficients
 * @param[in] vectorSeparator  string that seperates vectors (default is ,\n)
 *
 * @return    eigenvectors in text format
 */
template <typename MatrixType>
std::string eigenvectorsToString(const MatrixType &re, const MatrixType &im,
                                 std::string separator =", ", std::string vectorSeparator=",\n"){
  typedef typename MatrixType::size_type size_type;
  const size_type dim = re.numberOfColumns();
  std::ostringstream str;
  str << "{";

  for(size_type i=0; i< dim; i++){
    str << eigenvaluesToString(re.column(i), im.column(i), separator);
    if(i != dim-1)
      str << vectorSeparator;
  }
  str << "}";
  return str.str();
}


#endif /*CAPD2ALGLIB_H_*/
