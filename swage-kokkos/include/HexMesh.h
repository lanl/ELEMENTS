#pragma once

#include <vect3d.h>

struct HexMesh {
  enum { ND = 3 };
  
  CArrayKokkos<vect3d> cpts;
  
  HexMesh() {
  }
  
};
