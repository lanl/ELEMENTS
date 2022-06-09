#pragma once

template<typename T>
struct Vector2D {
  T x, y;

public:
  Vector2D(T _x=0, T _y=0)
    : x(_x), y(_y) { }
    
public:
  T & operator[](const int idx) {
    assert(idx<2);
    return (&x)[idx];
  }
  const T & operator[](const int idx) const {
    assert(idx<2);
    return (&x)[idx];
  }
  T & operator()(const int idx) {
    assert(idx<2);
    return (&x)[idx];
  }
  const T & operator()(const int idx) const {
    assert(idx<2);
    return (&x)[idx];
  }
};

typedef Vector2D<double> vect2d;
typedef Vector2D<float> vect2f;
