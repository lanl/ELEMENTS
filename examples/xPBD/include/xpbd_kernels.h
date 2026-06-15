#ifndef XPBD_KERNELS_H
#define XPBD_KERNELS_H

#include "matar.h"
#include "ELEMENTS.h"

using namespace mtr;

namespace xpbd {


KOKKOS_INLINE_FUNCTION
void compute_jacobian(
    const swage::Mesh& mesh,
    const DCArrayKokkos<double>& jacobian_matrix,
    const MPICArrayKokkos<double>& node_coords,
    const elements::fe_ref_elem_t& ref_elem,
    const int num_dims,
    const int num_nodes_in_element,
    const int elem_gid,
    const int gauss_pt_gid,
    const int gauss_pt_lid)
{
    
    for(int dim_i = 0; dim_i < num_dims; dim_i++){
        for(int dim_j = 0; dim_j < num_dims; dim_j++){
            jacobian_matrix(gauss_pt_gid, dim_i, dim_j) = 0.0;
            
            for(int node_lid = 0; node_lid < num_nodes_in_element; node_lid++){
                int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
                jacobian_matrix( gauss_pt_gid, dim_i, dim_j) +=
                    node_coords(node_gid, dim_i) 
                    *ref_elem.gauss_point_grad_basis(gauss_pt_lid, node_lid, dim_j);
                    
            } // end loop over nodes
        } // end loop over dimensions
    } // end loop over dimensions
}

KOKKOS_INLINE_FUNCTION
void compute_determinant_and_inverse(
    const DCArrayKokkos<double>& jacobian_matrix,
    ViewCArrayDevice<double>& inverse_matrix,
    double& determinant,
    const int num_dims,
    const int gauss_gid)
{


    // Copy the jacobian into a fixed-size array for inversion and determinant
    double jac[3][3] = {0.0};
    for(int i = 0; i < num_dims; i++) {
        for(int j = 0; j < num_dims; j++) {
            jac[i][j] = jacobian_matrix(gauss_gid, i, j);
        }
    }

    // Compute determinant and inverse for 2D or 3D only
    if (num_dims == 2) {
        determinant = jac[0][0]*jac[1][1] - jac[0][1]*jac[1][0];
        double inv_det = 1.0 / determinant;
        inverse_matrix(0,0) =  jac[1][1] * inv_det;
        inverse_matrix(0,1) = -jac[0][1] * inv_det;
        inverse_matrix(1,0) = -jac[1][0] * inv_det;
        inverse_matrix(1,1) =  jac[0][0] * inv_det;
    } else if (num_dims == 3) {
        determinant =
            jac[0][0]*(jac[1][1]*jac[2][2] - jac[1][2]*jac[2][1]) -
            jac[0][1]*(jac[1][0]*jac[2][2] - jac[1][2]*jac[2][0]) +
            jac[0][2]*(jac[1][0]*jac[2][1] - jac[1][1]*jac[2][0]);
        double inv_det = 1.0 / determinant;

        inverse_matrix(0,0) =  (jac[1][1]*jac[2][2] - jac[1][2]*jac[2][1]) * inv_det;
        inverse_matrix(0,1) = -(jac[0][1]*jac[2][2] - jac[0][2]*jac[2][1]) * inv_det;
        inverse_matrix(0,2) =  (jac[0][1]*jac[1][2] - jac[0][2]*jac[1][1]) * inv_det;
        inverse_matrix(1,0) = -(jac[1][0]*jac[2][2] - jac[1][2]*jac[2][0]) * inv_det;
        inverse_matrix(1,1) =  (jac[0][0]*jac[2][2] - jac[0][2]*jac[2][0]) * inv_det;
        inverse_matrix(1,2) = -(jac[0][0]*jac[1][2] - jac[0][2]*jac[1][0]) * inv_det;
        inverse_matrix(2,0) =  (jac[1][0]*jac[2][1] - jac[1][1]*jac[2][0]) * inv_det;
        inverse_matrix(2,1) = -(jac[0][0]*jac[2][1] - jac[0][1]*jac[2][0]) * inv_det;
        inverse_matrix(2,2) =  (jac[0][0]*jac[1][1] - jac[0][1]*jac[1][0]) * inv_det;
    } else if (num_dims == 1) {
        determinant = jac[0][0];
        inverse_matrix(0,0) = 1.0 / jac[0][0];
    } else {
        // Dimensions not supported
        determinant = 0.0;
    }

}

KOKKOS_INLINE_FUNCTION
void compute_F(
    const swage::Mesh& mesh,
    const MPICArrayKokkos<double>& coords_predicted,
    const MPICArrayKokkos<double>& J_inv_initial,
    const elements::fe_ref_elem_t& ref_elem,
    const DCArrayKokkos<double>& jacobian_matrix,
    const ViewCArrayDevice<double>& F,
    const int num_dims,
    const int elem_gid,
    const int gauss_gid,
    const int gauss_lid )
{

    for(int i = 0; i < num_dims; i++){
        for(int j = 0; j < num_dims; j++){
            F(i, j) = 0.0;
            for(int k = 0; k < num_dims; k++){
                F(i, j) += jacobian_matrix(gauss_gid, i, k) * J_inv_initial(gauss_gid, k, j);
            }
        }
    }
}

KOKKOS_INLINE_FUNCTION
void compute_cofactor_and_invariants(
    const ViewCArrayDevice<double>& F,
    const ViewCArrayDevice<double>& cofactor_F,
    const int num_dims,
    double& I2,
    double& I3) 
{
    for(int i=0; i<3; i++) for(int j=0; j<3; j++) cofactor_F(i, j) = 0.0;
    
    if (num_dims == 3) {
        cofactor_F(0, 0) = F(1, 1)*F(2, 2) - F(1, 2)*F(2, 1);
        cofactor_F(0, 1) = F(1, 2)*F(2, 0) - F(1, 0)*F(2, 2);
        cofactor_F(0, 2) = F(1, 0)*F(2, 1) - F(1, 1)*F(2, 0);
        cofactor_F(1, 0) = F(0, 2)*F(2, 1) - F(0, 1)*F(2, 2);
        cofactor_F(1, 1) = F(0, 0)*F(2, 2) - F(0, 2)*F(2, 0);
        cofactor_F(1, 2) = F(0, 1)*F(2, 0) - F(0, 0)*F(2, 1);
        cofactor_F(2, 0) = F(0, 1)*F(1, 2) - F(0, 2)*F(1, 1);
        cofactor_F(2, 1) = F(0, 2)*F(1, 0) - F(0, 0)*F(1, 2);
        cofactor_F(2, 2) = F(0, 0)*F(1, 1) - F(0, 1)*F(1, 0);
        
        I3 = F(0, 0)*cofactor_F(0, 0) + F(0, 1)*cofactor_F(1, 0) + F(0, 2)*cofactor_F(2, 0);
    } else if (num_dims == 2) {
        cofactor_F(0, 0) =  F(1, 1);
        cofactor_F(0, 1) = -F(1, 0);
        cofactor_F(1, 0) = -F(0, 1);
        cofactor_F(1, 1) =  F(0, 0);
        
        I3 = F(0, 0)*F(1, 1) - F(0, 1)*F(1, 0);
    } else {
        I3 = F(0, 0);
        cofactor_F(0, 0) = 1.0;
    }

    I2 = 0.0;
    for(int i=0; i<num_dims; i++) {
        for(int j=0; j<num_dims; j++) {
            I2 += F(i, j)*F(i, j);
        }
    }
}

KOKKOS_INLINE_FUNCTION
void compute_dN_dX(
    const MPICArrayKokkos<double>& J_inv_initial,
    const elements::fe_ref_elem_t& ref_elem,
    const int num_dims,
    const int gauss_gid,
    const int gauss_lid,
    const int node_lid,
    double dN_dX[3]) 
{
    for(int i = 0; i < num_dims; i++){
        dN_dX[i] = 0.0;
        for(int j = 0; j < num_dims; j++){
            dN_dX[i] += ref_elem.gauss_point_grad_basis(gauss_lid, node_lid, j) * J_inv_initial(gauss_gid, j, i);
        }
    }
}

} // namespace xpbd

#endif
