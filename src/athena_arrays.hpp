#ifndef ATHENA_ARRAYS_HPP
#define ATHENA_ARRAYS_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file athena_arrays.hpp
//  \brief provides array classes valid in 1D to 5D.
//
//  The operator() is overloaded, e.g. elements of a 4D array of size [N4xN3xN2xN1]
//  are accessed as:  A(n,k,j,i) = A[i + N1*(j + N2*(k + N3*n))]
//  NOTE THE TRAILING INDEX INSIDE THE PARENTHESES IS INDEXED FASTEST

// C++ headers
#include <cstddef>  // size_t
#include <iostream>
#include <fstream>

template<typename T>
class AthenaArray {
public:
  AthenaArray();
  AthenaArray(std::string fname, char c = ' ');
  ~AthenaArray();
  // define copy constructor and overload assignment operator so both do deep copies.
  AthenaArray(const AthenaArray<T>& t);
  AthenaArray<T> &operator= (const AthenaArray<T> &t);

  // public functions to allocate/deallocate memory for 1D-5D data
  void NewAthenaArray(int nx1);
  void NewAthenaArray(int nx2, int nx1);
  void NewAthenaArray(int nx3, int nx2, int nx1);
  void NewAthenaArray(int nx4, int nx3, int nx2, int nx1);
  void NewAthenaArray(int nx5, int nx4, int nx3, int nx2, int nx1);
  void DeleteAthenaArray();

  // functions to get array dimensions 
  int GetDim1() const { return nx1_; }
  int GetDim2() const { return nx2_; }
  int GetDim3() const { return nx3_; }
  int GetDim4() const { return nx4_; }
  int GetDim5() const { return nx5_; }

  // a function to get the total size of the array
  int GetSize() const { return nx1_*nx2_*nx3_*nx4_*nx5_; }
  size_t GetSizeInBytes() const {return nx1_*nx2_*nx3_*nx4_*nx5_*sizeof(T); }

  bool IsShallowCopy() { return (scopy_ == true); }
  T *data() { return pdata_; }
  const T *data() const	{ return pdata_; }

  // overload operator() to access 1d-5d data
  T &operator() (const int n) { 
    return pdata_[n]; }
  T operator() (const int n) const { 
    return pdata_[n]; }

  T &operator() (const int n, const int i) { 
    return pdata_[i + nx1_*n]; }
  T operator() (const int n, const int i) const { 
    return pdata_[i + nx1_*n]; }

  T &operator() (const int n, const int j, const int i) { 
    return pdata_[i + nx1_*(j + nx2_*n)]; }
  T operator() (const int n, const int j, const int i) const { 
    return pdata_[i + nx1_*(j + nx2_*n)]; }

  T &operator() (const int n, const int k, const int j, const int i) { 
    return pdata_[i + nx1_*(j + nx2_*(k + nx3_*n))]; }
  T operator() (const int n, const int k, const int j, const int i) const { 
    return pdata_[i + nx1_*(j + nx2_*(k + nx3_*n))]; }

  T &operator() (const int m, const int n, const int k, const int j, const int i) { 
    return pdata_[i + nx1_*(j + nx2_*(k + nx3_*(n + nx4_*m)))]; }
  T operator() (const int m, const int n, const int k, const int j, const int i) const {
    return pdata_[i + nx1_*(j + nx2_*(k + nx3_*(n + nx4_*m)))]; }

  // functions that initialize an array with shallow copy or slice from another array
  void InitWithShallowCopy(AthenaArray<T> const& src);
  void InitWithShallowSlice(AthenaArray<T> const& src, const int dim, const int indx,
    const int nvar);

private:
  T *pdata_;
  int nx1_, nx2_, nx3_, nx4_, nx5_;
  bool scopy_;  // true if shallow copy (prevents source from being deleted)
};

//constructor

template<typename T>
AthenaArray<T>::AthenaArray()
  : pdata_(0), nx1_(0), nx2_(0), nx3_(0), nx4_(0), nx5_(0), scopy_(true)
{
}

// destructor

template<typename T>
AthenaArray<T>::~AthenaArray()
{
  DeleteAthenaArray();
}

// copy constructor (does a deep copy)

template<typename T>
AthenaArray<T>::AthenaArray(const AthenaArray<T>& src) {
  nx1_ = src.nx1_;
  nx2_ = src.nx2_;
  nx3_ = src.nx3_;
  nx4_ = src.nx4_;
  nx5_ = src.nx5_;
  if (src.pdata_) {
    std::size_t size = (src.nx1_)*(src.nx2_)*(src.nx3_)*(src.nx4_)*(src.nx5_);
    pdata_ = new T[size]; // allocate memory for array data
    for (std::size_t i=0; i<size; ++i) {
      pdata_[i] = src.pdata_[i]; // copy data (not just addresses!) into new memory
    } 
  }
}

// assignment operator (does a deep copy).  Does not allocate memory for destination.
// THIS REQUIRES DESTINATION ARRAY BE ALREADY ALLOCATED AND SAME SIZE AS SOURCE

template<typename T>
AthenaArray<T> &AthenaArray<T>::operator= (const AthenaArray<T> &src) {
  if (this != &src){
    std::size_t size = (src.nx1_)*(src.nx2_)*(src.nx3_)*(src.nx4_)*(src.nx5_);
    for (std::size_t i=0; i<size; ++i) {
      this->pdata_[i] = src.pdata_[i]; // copy data (not just addresses!)
    } 
  }
  return *this;
}

//----------------------------------------------------------------------------------------
//! \fn AthenaArray::InitWithShallowCopy()
//  \brief shallow copy of array (copies ptrs, but not data)

template<typename T>
void AthenaArray<T>::InitWithShallowCopy(AthenaArray<T> const& src) {
  nx1_=src.nx1_;
  nx2_=src.nx2_;
  nx3_=src.nx3_;
  nx4_=src.nx4_;
  nx5_=src.nx5_;
  pdata_ = src.pdata_;
  scopy_ = true;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn AthenaArray::InitWithShallowSlice()
//  \brief shallow copy of nvar elements in dimension dim of an array, starting at
//  index=indx.  Copies pointers to data, but not data itself.

template<typename T>
void AthenaArray<T>::InitWithShallowSlice(AthenaArray<T> const& src, const int dim,
  const int indx, const int nvar)
{
  pdata_ = src.pdata_;

  if (dim == 5) {
    nx5_=nvar;
    nx4_=src.nx4_;
    nx3_=src.nx3_;
    nx2_=src.nx2_;
    nx1_=src.nx1_;
    pdata_ += indx*(nx1_*nx2_*nx3_*nx4_);
  } else if (dim == 4) {
    nx5_=1;
    nx4_=nvar;
    nx3_=src.nx3_;
    nx2_=src.nx2_;
    nx1_=src.nx1_;
    pdata_ += indx*(nx1_*nx2_*nx3_);
  } else if (dim == 3) {
    nx5_=1;
    nx4_=1;
    nx3_=nvar;
    nx2_=src.nx2_;
    nx1_=src.nx1_;
    pdata_ += indx*(nx1_*nx2_);
  } else if (dim == 2) {
    nx5_=1;
    nx4_=1;
    nx3_=1;
    nx2_=nvar;
    nx1_=src.nx1_;
    pdata_ += indx*(nx1_);
  } else if (dim == 1) {
    nx5_=1;
    nx4_=1;
    nx3_=1;
    nx2_=1;
    nx1_=nvar;
    pdata_ += indx;
  }
  scopy_ = true;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn
//  \brief allocate new 1D array with elements initialized to zero.

template<typename T>
void AthenaArray<T>::NewAthenaArray(int nx1)
{
  scopy_ = false;
  nx1_ = nx1;
  nx2_ = 1;
  nx3_ = 1;
  nx4_ = 1;
  nx5_ = 1;
  pdata_ = new T[nx1](); // allocate memory and initialize to zero
}
 
//----------------------------------------------------------------------------------------
//! \fn
//  \brief 2d data allocation

template<typename T>
void AthenaArray<T>::NewAthenaArray(int nx2, int nx1)
{
  scopy_ = false;
  nx1_ = nx1;
  nx2_ = nx2;
  nx3_ = 1;
  nx4_ = 1;
  nx5_ = 1;
  pdata_ = new T[nx1*nx2](); // allocate memory and initialize to zero
}
 
//----------------------------------------------------------------------------------------
//! \fn
//  \brief 3d data allocation

template<typename T>
void AthenaArray<T>::NewAthenaArray(int nx3, int nx2, int nx1)
{
  scopy_ = false;
  nx1_ = nx1;
  nx2_ = nx2;
  nx3_ = nx3;
  nx4_ = 1;
  nx5_ = 1;
  pdata_ = new T[nx1*nx2*nx3](); // allocate memory and initialize to zero
}
 
//----------------------------------------------------------------------------------------
//! \fn
//  \brief 4d data allocation

template<typename T>
void AthenaArray<T>::NewAthenaArray(int nx4, int nx3, int nx2, int nx1)
{
  scopy_ = false;
  nx1_ = nx1;
  nx2_ = nx2;
  nx3_ = nx3;
  nx4_ = nx4;
  nx5_ = 1;
  pdata_ = new T[nx1*nx2*nx3*nx4](); // allocate memory and initialize to zero
}

//----------------------------------------------------------------------------------------
//! \fn
//  \brief 5d data allocation

template<typename T>
void AthenaArray<T>::NewAthenaArray(int nx5, int nx4, int nx3, int nx2, int nx1)
{
  scopy_ = false;
  nx1_ = nx1;
  nx2_ = nx2;
  nx3_ = nx3;
  nx4_ = nx4;
  nx5_ = nx5;
  pdata_ = new T[nx1*nx2*nx3*nx4*nx5](); // allocate memory and initialize to zero
}

//----------------------------------------------------------------------------------------
//! \fn
//  \brief  free memory allocated for data array

template<typename T>
void AthenaArray<T>::DeleteAthenaArray()
{
  if (scopy_) {
    pdata_ = NULL;
  } else {
    delete[] pdata_;
    scopy_ = true;
  }
} 

template<typename T>
AthenaArray<T>::AthenaArray(std::string fname, char c)
{
  // read first time to determine dimension
  std::ifstream inp(fname.c_str(), std::ios::in);
  std::string line;
  std::getline(inp, line);
  int rows = 0, cols = 0;
  if (!line.empty()) {
    rows = 1;
    cols = line[0] == c ? 0 : 1;
    for (int i = 1; i < line.length(); ++i)
      if (line[i-1] == c && line[i] != c) cols++;
  }
  while (std::getline(inp, line)) ++rows;
  inp.close();

  // allocate space
  scopy_ = false;
  nx1_ = cols;
  nx2_ = rows;
  nx3_ = 1;
  nx4_ = 1;
  nx5_ = 1;
  pdata_ = new T[nx1_*nx2_](); // allocate memory and initialize to zero

  // read second time
  inp.open(fname.c_str(), std::ios::in);
  for (int j = 0; j < rows; ++j)
    for (int i = 0; i < cols; ++i)
      inp >> pdata_[j*cols + i];
}

#endif // ATHENA_ARRAYS_HPP
