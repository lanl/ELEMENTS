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

// KOKKOS_INLINE_FUNCTION
// void invert_3x3_matrix(
//     const double* input_matrix,
//     double* inverse_matrix,
//     double& matrix_determinant)
// {
//     const double a00 = input_matrix[0];
//     const double a01 = input_matrix[1];
//     const double a02 = input_matrix[2];
//     const double a10 = input_matrix[3];
//     const double a11 = input_matrix[4];
//     const double a12 = input_matrix[5];
//     const double a20 = input_matrix[6];
//     const double a21 = input_matrix[7];
//     const double a22 = input_matrix[8];

//     const double cofactor_00 = a11 * a22 - a12 * a21;
//     const double cofactor_01 = a12 * a20 - a10 * a22;
//     const double cofactor_02 = a10 * a21 - a11 * a20;
//     const double cofactor_10 = a02 * a21 - a01 * a22;
//     const double cofactor_11 = a00 * a22 - a02 * a20;
//     const double cofactor_12 = a01 * a20 - a00 * a21;
//     const double cofactor_20 = a01 * a12 - a02 * a11;
//     const double cofactor_21 = a02 * a10 - a00 * a12;
//     const double cofactor_22 = a00 * a11 - a01 * a10;

//     matrix_determinant = a00 * cofactor_00 + a01 * cofactor_10 + a02 * cofactor_20;

//     if (std::fabs(matrix_determinant) < 1.0e-14) {
//         matrix_determinant = 0.0;
//         for (int index = 0; index < 9; ++index) {
//             inverse_matrix[index] = 0.0;
//         }
//         return;
//     }

//     const double inverse_det = 1.0 / matrix_determinant;
//     inverse_matrix[0] = cofactor_00 * inverse_det;
//     inverse_matrix[1] = cofactor_10 * inverse_det;
//     inverse_matrix[2] = cofactor_20 * inverse_det;
//     inverse_matrix[3] = cofactor_01 * inverse_det;
//     inverse_matrix[4] = cofactor_11 * inverse_det;
//     inverse_matrix[5] = cofactor_21 * inverse_det;
//     inverse_matrix[6] = cofactor_02 * inverse_det;
//     inverse_matrix[7] = cofactor_12 * inverse_det;
//     inverse_matrix[8] = cofactor_22 * inverse_det;
// }






// KOKKOS_INLINE_FUNCTION
// void compute_deformation_gradient(
//     const double* element_node_positions,
//     const double* element_reference_positions,
//     const size_t num_nodes_in_element,
//     const double* physical_basis_gradient,
//     double* deformation_gradient)
// {
//     for (int row = 0; row < 9; ++row) {
//         deformation_gradient[row] = 0.0;
//     }

//     deformation_gradient[0] = 1.0;
//     deformation_gradient[4] = 1.0;
//     deformation_gradient[8] = 1.0;

//     for (size_t node_lid = 0; node_lid < num_nodes_in_element; ++node_lid) {
//         for (int spatial_dim = 0; spatial_dim < 3; ++spatial_dim) {
//             const double displacement = element_node_positions[3 * node_lid + spatial_dim]
//                 - element_reference_positions[3 * node_lid + spatial_dim];
//             for (int reference_dim = 0; reference_dim < 3; ++reference_dim) {
//                 deformation_gradient[3 * spatial_dim + reference_dim] +=
//                     displacement * physical_basis_gradient[3 * node_lid + reference_dim];
//             }
//         }
//     }
// }

// KOKKOS_INLINE_FUNCTION
// void compute_stvk_strain_energy_density(
//     const double* deformation_gradient,
//     const double lame_lambda,
//     const double lame_mu,
//     double& strain_energy_density)
// {
//     const double f00 = deformation_gradient[0];
//     const double f01 = deformation_gradient[1];
//     const double f02 = deformation_gradient[2];
//     const double f10 = deformation_gradient[3];
//     const double f11 = deformation_gradient[4];
//     const double f12 = deformation_gradient[5];
//     const double f20 = deformation_gradient[6];
//     const double f21 = deformation_gradient[7];
//     const double f22 = deformation_gradient[8];

//     const double c00 = f00 * f00 + f10 * f10 + f20 * f20;
//     const double c11 = f01 * f01 + f11 * f11 + f21 * f21;
//     const double c22 = f02 * f02 + f12 * f12 + f22 * f22;
//     const double c01 = f00 * f01 + f10 * f11 + f20 * f21;
//     const double c02 = f00 * f02 + f10 * f12 + f20 * f22;
//     const double c12 = f01 * f02 + f11 * f12 + f21 * f22;

//     const double e00 = 0.5 * (c00 - 1.0);
//     const double e11 = 0.5 * (c11 - 1.0);
//     const double e22 = 0.5 * (c22 - 1.0);
//     const double e01 = 0.5 * c01;
//     const double e02 = 0.5 * c02;
//     const double e12 = 0.5 * c12;

//     const double trace_strain = e00 + e11 + e22;
//     const double trace_strain_squared = trace_strain * trace_strain;
//     const double strain_double_contraction =
//         e00 * e00 + e11 * e11 + e22 * e22 + 2.0 * (e01 * e01 + e02 * e02 + e12 * e12);

//     strain_energy_density =
//         lame_mu * trace_strain_squared
//         + lame_lambda * trace_strain_squared
//         + lame_mu * strain_double_contraction;
// }

// KOKKOS_INLINE_FUNCTION
// void compute_stvk_first_piola_stress(
//     const double* deformation_gradient,
//     const double lame_lambda,
//     const double lame_mu,
//     double* first_piola_kirchhoff_stress)
// {
//     const double f00 = deformation_gradient[0];
//     const double f01 = deformation_gradient[1];
//     const double f02 = deformation_gradient[2];
//     const double f10 = deformation_gradient[3];
//     const double f11 = deformation_gradient[4];
//     const double f12 = deformation_gradient[5];
//     const double f20 = deformation_gradient[6];
//     const double f21 = deformation_gradient[7];
//     const double f22 = deformation_gradient[8];

//     const double c00 = f00 * f00 + f10 * f10 + f20 * f20;
//     const double c11 = f01 * f01 + f11 * f11 + f21 * f21;
//     const double c22 = f02 * f02 + f12 * f12 + f22 * f22;
//     const double c01 = f00 * f01 + f10 * f11 + f20 * f21;
//     const double c02 = f00 * f02 + f10 * f12 + f20 * f22;
//     const double c12 = f01 * f02 + f11 * f12 + f21 * f22;

//     const double e00 = 0.5 * (c00 - 1.0);
//     const double e11 = 0.5 * (c11 - 1.0);
//     const double e22 = 0.5 * (c22 - 1.0);
//     const double e01 = 0.5 * c01;
//     const double e02 = 0.5 * c02;
//     const double e12 = 0.5 * c12;
//     const double trace_strain = e00 + e11 + e22;

//     const double s00 = lame_lambda * trace_strain + 2.0 * lame_mu * e00;
//     const double s11 = lame_lambda * trace_strain + 2.0 * lame_mu * e11;
//     const double s22 = lame_lambda * trace_strain + 2.0 * lame_mu * e22;
//     const double s01 = 2.0 * lame_mu * e01;
//     const double s02 = 2.0 * lame_mu * e02;
//     const double s12 = 2.0 * lame_mu * e12;

//     first_piola_kirchhoff_stress[0] = f00 * s00 + f01 * s01 + f02 * s02;
//     first_piola_kirchhoff_stress[1] = f00 * s01 + f01 * s11 + f02 * s12;
//     first_piola_kirchhoff_stress[2] = f00 * s02 + f01 * s12 + f02 * s22;
//     first_piola_kirchhoff_stress[3] = f10 * s00 + f11 * s01 + f12 * s02;
//     first_piola_kirchhoff_stress[4] = f10 * s01 + f11 * s11 + f12 * s12;
//     first_piola_kirchhoff_stress[5] = f10 * s02 + f11 * s12 + f12 * s22;
//     first_piola_kirchhoff_stress[6] = f20 * s00 + f21 * s01 + f22 * s02;
//     first_piola_kirchhoff_stress[7] = f20 * s01 + f21 * s11 + f22 * s12;
//     first_piola_kirchhoff_stress[8] = f20 * s02 + f21 * s12 + f22 * s22;
// }




} // namespace xpbd

#endif
