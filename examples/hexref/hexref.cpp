#include "HexRef.h"

#include "Kokkos_Core.hpp"

#include <iostream>


int main(int argc, char *argv[]) {
  Kokkos::initialize(argc, argv);
  int order_basis = 2;
  BasisType type_basis = BasisType::LAGRANGE;
  HexRef elem(order_basis, type_basis);

  // Print out the quadrature points pertaining the the cells (subelements)
  for (int i = 0; i < elem.num_ref_cells_in_elem(); ++i) {
    std::cout << elem.ref_cell_positions(i, 0) << " "
              << elem.ref_cell_positions(i, 1) << " "
              << elem.ref_cell_positions(i, 2) << " "
              << std::endl;
  }

  // Print out the basis function evaluations on the cell (subelements)
  for (int i = 0; i < elem.num_ref_cells_in_elem(); ++i) {
    for (int j = 0; j < elem.num_basis(); ++j)
      std::cout << elem.ref_cell_basis(i, j) << " ";
    std::cout << std::endl;
  }

  return 0;
}
