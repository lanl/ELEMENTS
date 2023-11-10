#pragma once

template<typename T>
struct Vector4D {
  T x, y, z, w;

public:
  Vector4D(T _x=0, T _y=0, T _z=0, T _w=0)
    : x(_x), y(_y), z(_z), w(_w) { }
    
public:
  T & operator[](const int idx) {
    assert(idx<4);
    return (&x)[idx];
  }
  const T & operator[](const int idx) const {
    assert(idx<4);
    return (&x)[idx];
  }
  T & operator()(const int idx) {
    assert(idx<4);
    return (&x)[idx];
  }
  const T & operator()(const int idx) const {
    assert(idx<4);
    return (&x)[idx];
  }
};


typedef Vector4D<double> vect4d;
typedef Vector4D<float> vect4f;
