#ifndef LITTLEMATH_H
#define LITTLEMATH_H

#include <iostream>
#include <array>
#include <cassert>
#include <cmath>

template <typename T, size_t N>
struct VectorT {
  std::array<T,N> _data;
  VectorT() {
    for (int i=0; i<N; i++)
      (*this)(i)=T(0);
  }
  VectorT(const VectorT& v) {
    for (int i=0; i<N; i++)
      (*this)(i)=v(i);
  }
  const T& operator()(int i) const {assert(i<N); return _data[i];}
  T& operator()(int i) {assert(i<N); return _data[i];}
  VectorT operator+(const VectorT& v) const {
    VectorT vo;
    for (int i=0; i<N; i++)
      vo(i) = (*this)(i)+v(i);
    return vo;
  }
  VectorT& operator+=(const VectorT& v) {
    for (int i=0; i<N; i++)
      (*this)(i)+=v(i);
    return (*this);
  }
  VectorT operator*(const double a) const {
    VectorT vo;
    for (int i=0; i<N; i++)
      vo(i) = a*(*this)(i);
    return vo;
  }
  VectorT operator/(const double a) const {
    VectorT vo;
    for (int i=0; i<N; i++)
      vo(i) = (*this)(i)/a;
    return vo;
  }
  T operator*(const VectorT& v) const {
    T p = 0;
    for (int i=0; i<N; i++)
      p += (*this)(i)*v(i);
    return p;
  }
  T dot(const VectorT& v) const { return (*this)*v; }
  T norm() const { return std::sqrt((*this)*(*this)); }
};

template <typename T,size_t N>
std::ostream& operator<<(std::ostream& os, const VectorT<T,N>& v) {
  for (int i=0; i<N; i++)
    os << v(i) << " ";
  return os;
}

template <typename T, size_t N, size_t M>
struct MatrixT {
  std::array<std::array<T,M>,N> _data;
  MatrixT() {
    for (int i=0; i<N; i++)
      for (int j=0; j<M; j++)
        (*this)(i,j)=T(0);
  }
  MatrixT(const MatrixT& m) {
    for (int i=0; i<N; i++)
      for (int j=0; j<M; j++)
        (*this)(i,j)=m(i,j);
  }
  MatrixT(const VectorT<T,N>& v) {
    for (int i=0; i<N; i++)
      (*this)(i,0)=v(i);
  }
  static MatrixT identityMatrix() {
    assert(N==M); MatrixT m;
    for (int i=0;i<N;i++)
      m(i,i)=1.0;
    return m;
  }
  const T& operator()(int i, int j) const { assert(i<N && j<M);
    return _data[i][j];}
  T& operator()(int i, int j) { assert(i<N && j<M);
    return _data[i][j];}
  MatrixT operator+(const MatrixT& m) const {
    MatrixT mo;
    for (int i=0; i<N; i++)
      for (int j=0; j<M; j++)
	mo(i,j) = (*this)(i,j)+m(i,j);
    return mo;
  }
  MatrixT operator*(const double a) const {
    MatrixT mo;
    for (int i=0; i<N; i++)
      for (int j=0; j<M; j++)
	mo(i,j) = (*this)(i,j)*a;
    return mo;
  }
  MatrixT operator/(const double a) const {
    MatrixT mo;
    for (int i=0; i<N; i++)
      for (int j=0; j<M; j++)
	mo(i,j) = (*this)(i,j)/a;
    return mo;
  }
  template<template<typename,size_t,size_t> class MatT, size_t P>
  MatT<T,N,P> operator*(const MatT<T,M,P>& m) const {
    MatT<T,N,P> mo;
    for (int i=0; i<N; i++)
      for (int j=0; j<P; j++)
	for (int k=0; k<M; k++)
	  mo(i,j) += (*this)(i,k)*m(k,j);
    return mo;
  }
  MatrixT<T,M,N> transpose() const {
    MatrixT<T,M,N> mo;
    for (int i=0; i<N; i++)
      for (int j=0; j<M; j++)
	  mo(j,i) = (*this)(i,j);
    return mo;
  }
  template<template<typename,size_t,size_t> class MatT, size_t P>
    MatT<T,N,M+P> cath(const MatT<T,N,P>& m) const {
    MatT<T,N,M+P> mo;
    for (int i=0; i<N; i++) {
      for (int j=0; j<M; j++)
	mo(i,j) += (*this)(i,j);
      for (int j=0; j<P; j++)
	mo(i,j+M) += m(i,j);
    }
    return mo;
  }
  template<template<typename,size_t,size_t> class MatT, size_t P>
    MatT<T,N+P,M> catv(const MatT<T,P,M>& m) const {
    MatT<T,N+P,M> mo;
    for (int j=0; j<M; j++) {
      for (int i=0; i<N; i++)
	mo(i,j) += (*this)(i,j);
      for (int i=0; i<P; i++)
	mo(i+N,j) += m(i,j);
    }
    return mo;
  }
};

template <typename T,size_t N,size_t M>
std::ostream& operator<<(std::ostream& os,
			 const MatrixT<T,N,M>& m) {
  for (int i=0; i<N; i++, os << std::endl)
    for (int j=0; j<M; j++)
      os << m(i,j) << " ";
  return os;
}

typedef VectorT<double, 2> Vector2d; 
typedef MatrixT<double, 2, 2> Matrix2d; 

struct SpatialVector2D: VectorT<double,3> {
  SpatialVector2D() {
    (*this)(0)=(*this)(1)=(*this)(2)=0.0;
  }
  SpatialVector2D( double th, Vector2d v ) {
    (*this)(0)=th;(*this)(1)=v(1);(*this)(2)=v(2);
  }
  SpatialVector2D( double v0, double v1, double v2 ) {
    (*this)(0)=v0;(*this)(1)=v1;(*this)(2)=v2;
  }
  SpatialVector2D( const SpatialVector2D& V ) {
    for (int i=0; i<3; i++, _data[i]=V(i));
  }
};

struct SpatialTransform2D {
  Matrix2d E;
  Vector2d r;
  SpatialTransform2D()
  : E(Matrix2d::identityMatrix()), r() {}
  SpatialTransform2D(double th, Vector2d& v)
  : E(Matrix2d::identityMatrix()), r(v) {
    setTheta(th);
  }
  void setTheta(double th) {
    E(0,0)=E(1,1)=std::cos(th);
    E(1,0)=-std::sin(th); E(0,1)=std::sin(th);
  }
  double getTheta() const {
    return std::atan2(E(0,1),E(0,0));
  }
};

std::ostream& operator<<(std::ostream& os,
			 const SpatialTransform2D& V) {
  os << "r:\n" << V.r
     << "\nE:\n" << V.E;
  return os;
}

#endif // #ifndef LITTLEMATH_H
