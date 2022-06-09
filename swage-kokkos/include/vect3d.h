#pragma once

template<typename T>
struct Vector3D {
  T x, y, z;

public:
  Vector3D(T _x=0, T _y=0, T _z=0)
    : x(_x), y(_y), z(_z) { }
    
public:
  T & operator[](const int idx) {
    assert(idx<3);
    return (&x)[idx];
  }
  const T & operator[](const int idx) const {
    assert(idx<3);
    return (&x)[idx];
  }
  T & operator()(const int idx) {
    assert(idx<3);
    return (&x)[idx];
  }
  const T & operator()(const int idx) const {
    assert(idx<3);
    return (&x)[idx];
  }

public:
  template<typename S>
  Vector3D<T> & operator=(const Vector3D<S> & rhs) {
    x = rhs.x; y = rhs.y; z = rhs.z;
    return *this;
  }
  template<typename S>
  Vector3D<T> & operator+=(const Vector3D<S> & rhs) {
    x += rhs.x; y += rhs.y; z += rhs.z;
    return *this;
  }
  template<typename S>
  Vector3D<T> & operator-=(const Vector3D<S> & rhs) {
    x -= rhs.x; y -= rhs.y; z -= rhs.z;
    return *this;
  }
  template<typename S>
  Vector3D<T> & operator*=(const S & rhs) {
    x *= rhs; y *= rhs; z *= rhs;
    return *this;
  }
  template<typename S>
  Vector3D<T> & operator/=(const S & rhs) {
    const T den(1.0/rhs);
    x *= den; y *= den; z *= den;
    return *this;
  }
};

typedef Vector3D<double> vect3d;
typedef Vector3D<float> vect3f;
