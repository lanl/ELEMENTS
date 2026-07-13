#ifndef GEOMETRY_H
#define GEOMETRY_H

// This pulls in kokkos, matar, mesh, ref_elem stuff, and PT-Scotch
#include "ELEMENTS.h"

#include "matar.h"
#include "shapes.h"

using namespace mtr;
using namespace swage; // unstructured mesh and hash
using namespace elements;

/////////////////////////////////////////////////////////////////////////////
///
/// \fn jacobian
///
/// \brief Calculates the jacobian matrix in 2D, 3D, 4D, ...
///
/// \param jacobian     The jacobian matrix calculated in the routine, J[dims,dims]
/// \param node_coords  An array containing all node coords on the mesh
/// \param nodes_in_an_elem  The node indices in a single elem, its a 1D array
/// \param a_grad_basis The gradient of the basis at a single point, Grad[DOFs,dims]
///
/////////////////////////////////////////////////////////////////////////////
template <typename T1, typename T2, typename T3, typename T4>
KOKKOS_INLINE_FUNCTION
void jacobian(
    const T1 &jacobian,         // e.g., ViewCArrayKokkos <double>
    const T2 &node_coords,      // e.g., DCArrayKokkos    <double>
    const T3 &nodes_in_an_elem, // e.g., ViewCArrayKokkos <size_t>
    const T4 &a_grad_basis){    // e.g., ViewCArrayKokkos <double>

    const size_t dims = a_grad_basis.dims(1);
    const size_t num_dofs_in_elem = nodes_in_an_elem.size();

    // setting jacobian matrix to all zeros
    for(size_t i = 0; i < dims; i++)  // looping over dimension
    for(size_t j = 0; j < dims; j++){ // looping over dimension
        jacobian(i, j) = 0.0;
    } // end for 

    // Calculate Jacobian: J[i,j] = partial x_i/partial \xi_j
    for(size_t i = 0; i < dims; i++) // looping over dimension (partial)
    for(size_t j = 0; j < dims; j++) // looping over dimension (node position)
    for(size_t node_lid = 0; node_lid < num_dofs_in_elem; node_lid++){ 
        const size_t node_gid = nodes_in_an_elem(node_lid);  
        jacobian(i, j) += node_coords(node_gid, i)*a_grad_basis(node_lid, j);           
    } // end for  


} // end of jacobian function



#endif // GEOMETRY_H