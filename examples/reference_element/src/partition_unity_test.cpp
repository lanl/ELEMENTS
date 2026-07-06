/**********************************************************************************************
© 2020. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and
to permit others to do so.
This program is open source under the BSD-3 License.
Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:
1.  Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer.
2.  Redistributions in binary form must reproduce the above copyright notice, this list of
conditions and the following disclaimer in the documentation and/or other materials
provided with the distribution.
3.  Neither the name of the copyright holder nor the names of its contributors may be used
to endorse or promote products derived from this software without specific prior
written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**********************************************************************************************/
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

// This pulls in kokkos, matar, mesh, point cloud, ref_elem stuff, and PT-Scotch
#include "ELEMENTS.h"

//#undef NDEBUG     // Ensures NDEBUG is turned off

using namespace mtr;
using namespace swage; // unstructured mesh and point cloud
using namespace elements;

bool Verbose = false;
size_t max_num   = 19; // max number of quadrature points to test up to, the limit is 19
size_t max_order = 19; // max polynomial order to test, limit is 19th-order with Legendre


void verify_partition_of_unity(const Quadrature_t& Quad,
                               const ReferenceElement_t& RefElem) {

    // device parallel loop
    FOR_ALL(qpt, 0, Quad.num_qpts_in_elem, {
        double sum = 0.0;
        
        for (size_t basis = 0; basis < RefElem.num_dofs_in_elem; basis++) {
            sum += RefElem.qpt_basis(qpt, basis);
        }
        if (fabs(sum - 1.0) > 1.e-12) {
            printf("Error: partion of unity failed, sum of basis = %f at qpt id = %d with order = %zu \n", sum, qpt, RefElem.num_dofs_1d-1);
            Kokkos::abort("Partition of unity failed at quadrature point ");
        }
        if(Verbose)printf("sum = %f \n", sum);
    }); // end loop over qpts

} // end function

void verify_gradient(const Quadrature_t& Quad,
                     const ReferenceElement_t& RefElem) {
    
    // device parallel loop
    FOR_ALL(qpt, 0, Quad.num_qpts_in_elem, {
        double sum[3];
        sum[0] = 0.0;
        sum[1] = 0.0;
        sum[2] = 0.0;
        
        for(size_t dim=0; dim<RefElem.elem_dims; dim++){
            for (size_t basis = 0; basis < RefElem.num_dofs_in_elem; basis++) {
                sum[dim] += RefElem.qpt_grad_basis(qpt, basis, dim);
            }
            if (fabs(sum[dim]) > 1.e-12) {
                printf("Error: gradient failed, sum of gradient basis = %f at rid = %d, for dim = %zu, with order = %zu \n", sum[dim], qpt, dim, RefElem.num_dofs_1d);
                Kokkos::abort("Gradient of basis failed at quadrature point ");
            }
            if(Verbose)printf("dim = %zu, sum = %f \n", dim, sum[dim]);
        } // for dim

    }); // end loop over qpts

} // end function

int main(int argc, char** argv) {

MATAR_INITIALIZE(argc, argv);
{ // MATAR scope

    printf("\n--- partition unity and gradient test ---\n");

    printf("\n--- DG element with Legendre Quadrature & Legendre DOFs ---\n");
    for(size_t num_qpts_1D = 1; num_qpts_1D<=max_num; num_qpts_1D++){

        if(Verbose)printf("num quadrature points in 1D = %zu \n", num_qpts_1D);

        Quadrature_t Quad;
        
        // elem_dims=1,2,3
        for(size_t elem_dims_test = 1; elem_dims_test<=3; elem_dims_test++){   
            Quad.initialize_quadrature(reference_space::GaussLegendre,
                                       num_qpts_1D,
                                       elem_dims_test);

            // build reference elements of varing orders
            for (size_t p_order = 0; p_order<max_order; p_order++){
                if(Verbose)printf("p_order = %zu: \n", p_order);
                ReferenceElement_t FERefElem;

                // p_order is the basis order for Lagrange polynomial
                FERefElem.initialize_ref_elem(reference_space::arbitraryOrderElement,
                                              reference_space::LagrangeLegendre,
                                              Quad,
                                              p_order);

                if(Verbose)printf("basis check: \n");
                verify_partition_of_unity(Quad, FERefElem);
                if(Verbose)printf("\n");

                if(Verbose)printf("gradient basis check: \n");
                verify_gradient(Quad, FERefElem);
                if(Verbose)printf("\n");
            } // end p_order loop
        } // elem

        if(Verbose)printf("\n");
    } // end loop of num qpts 


    printf("\n--- DG element with Lobatto Quadrature & Legendre DOFs ---\n");
    for(size_t num_qpts_1D = 1; num_qpts_1D<=max_num; num_qpts_1D++){

        if(Verbose)printf("num quadrature points in 1D = %zu \n", num_qpts_1D);

        Quadrature_t Quad;
        
        // elem_dims=1,2,3
        for(size_t elem_dims_test = 1; elem_dims_test<=3; elem_dims_test++){   
            Quad.initialize_quadrature(reference_space::GaussLobatto,
                                       num_qpts_1D,
                                       elem_dims_test);

            // build reference elements of varing orders
            for (size_t p_order = 0; p_order<max_order; p_order++){
                if(Verbose)printf("p_order = %zu: \n", p_order);
                ReferenceElement_t FERefElem;

                // p_order is the basis order for Lagrange polynomial
                FERefElem.initialize_ref_elem(reference_space::arbitraryOrderElement,
                                              reference_space::LagrangeLegendre,
                                              Quad,
                                              p_order);

                if(Verbose)printf("basis check: \n");
                verify_partition_of_unity(Quad, FERefElem);
                if(Verbose)printf("\n");

                if(Verbose)printf("gradient basis check: \n");
                verify_gradient(Quad, FERefElem);
                if(Verbose)printf("\n");
            } // end p_order loop
        } // elem

        if(Verbose)printf("\n");
    } // end loop of num qpts 


    printf("\n--- FE element with Lobatto Quadrature & Lobatto DOFs ---\n");
    for(size_t num_qpts_1D = 1; num_qpts_1D<=max_num; num_qpts_1D++){

        if(Verbose)printf("num quadrature points in 1D = %zu \n", num_qpts_1D);

        Quadrature_t Quad;
        
        // elem_dims=1,2,3
        for(size_t elem_dims_test = 1; elem_dims_test<=3; elem_dims_test++){   
            Quad.initialize_quadrature(reference_space::GaussLobatto,
                                       num_qpts_1D,
                                       elem_dims_test);

            // build reference elements of varing orders
            for (size_t p_order = 1; p_order<max_order; p_order++){
                if(Verbose)printf("p_order = %zu: \n", p_order);
                ReferenceElement_t FERefElem;

                // p_order is the basis order for Lagrange polynomial
                FERefElem.initialize_ref_elem(reference_space::arbitraryOrderElement,
                                              reference_space::LagrangeLobatto,
                                              Quad,
                                              p_order);

                if(Verbose)printf("basis check: \n");
                verify_partition_of_unity(Quad, FERefElem);
                if(Verbose)printf("\n");

                if(Verbose)printf("gradient basis check: \n");
                verify_gradient(Quad, FERefElem);
                if(Verbose)printf("\n");

            } // end p_order loop
        } // elem

        if(Verbose)printf("\n");
    } // end loop of num qpts 

    printf("\n--- FE element with Legendre Quadrature & Lobatto DOFs ---\n");
    for(size_t num_qpts_1D = 1; num_qpts_1D<=max_num; num_qpts_1D++){

        if(Verbose)printf("num quadrature points in 1D = %zu \n", num_qpts_1D);

        Quadrature_t Quad;
        
        // elem_dims=1,2,3
        for(size_t elem_dims_test = 1; elem_dims_test<=3; elem_dims_test++){   
            Quad.initialize_quadrature(reference_space::GaussLegendre,
                                       num_qpts_1D,
                                       elem_dims_test);

            // build reference elements of varing orders
            for (size_t p_order = 1; p_order<max_order; p_order++){
                if(Verbose)printf("p_order = %zu: \n", p_order);
                ReferenceElement_t FERefElem;

                // p_order is the basis order for Lagrange polynomial
                FERefElem.initialize_ref_elem(reference_space::arbitraryOrderElement,
                                              reference_space::LagrangeLobatto,
                                              Quad,
                                              p_order);

                if(Verbose)printf("basis check: \n");
                verify_partition_of_unity(Quad, FERefElem);
                if(Verbose)printf("\n");

                if(Verbose)printf("gradient basis check: \n");
                verify_gradient(Quad, FERefElem);
                if(Verbose)printf("\n");
            } // end p_order loop
        } // elem

        if(Verbose)printf("\n");
    } // end loop of num qpts 

    printf("\nAll partition unity and gradient checks passed.\n");

}
MATAR_FINALIZE();

} // end main

