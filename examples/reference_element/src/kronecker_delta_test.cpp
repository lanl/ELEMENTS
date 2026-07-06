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

// This pulls in kokkos, matar, mesh, hash, ref_elem stuff, and PT-Scotch
#include "ELEMENTS.h"

//#undef NDEBUG     // Ensures NDEBUG is turned off

using namespace mtr;
using namespace swage; // unstructured mesh and hash
using namespace elements;

bool Verbose = false;
size_t max_num   = 19; // max number of quadrature points to test up to, the limit is 19


int main(int argc, char** argv) {

MATAR_INITIALIZE(argc, argv);
{ // MATAR scope

    printf("\n--- kronecker delta test ---\n");

    // Note: quadrature points will be collocated with kinematic DOFs
    for(size_t num_qpts_1D = 2; num_qpts_1D<=max_num; num_qpts_1D++){

        Quadrature_t Quad;
        
        // elem_dims=1,2,3
        for(size_t elem_dims_test = 1; elem_dims_test<=3; elem_dims_test++){   
            Quad.initialize_quadrature(reference_space::GaussLobatto,
                                       num_qpts_1D,
                                       elem_dims_test);

            // build reference element with collocated DOFs
            const size_t p_order = num_qpts_1D - 1;

            if(Verbose)printf("p_order = %zu: \n", p_order);
            ReferenceElement_t FERefElem;

            // p_order is the basis order for Lagrange polynomial
            FERefElem.initialize_ref_elem(reference_space::arbitraryOrderElement,
                                          reference_space::LagrangeLobatto,
                                          Quad,
                                          p_order);


            // Check: should be 1 at i, 0 elsewhere
            FOR_ALL(basis, 0, FERefElem.num_dofs_in_elem, {
                
                if(Verbose) printf("basis DOF id = %d \n", basis);
                for(int dof_pt=0; dof_pt<FERefElem.num_dofs_in_elem; dof_pt++){

                    // remember: dof_pt = qpt, it is a collocated grid

                    double expected = (basis == dof_pt) ? 1.0 : 0.0;

                    if(Verbose){
                        printf("   checking at DOF basis = %d \n", basis);
                        printf("   expected val = %f, basis val = %f \n", expected, FERefElem.qpt_basis(dof_pt, basis));
                        for(size_t dim=0; dim<elem_dims_test; dim++){
                            printf("   quad pos     = %f, basis pos = %f \n", Quad.qpt_positions(dof_pt, dim), FERefElem.dof_positions(dof_pt, dim));
                        }
                    }
                    if (fabs(FERefElem.qpt_basis(dof_pt, basis) - expected) > 1.0e-12) {
                        Kokkos::abort("ERROR: Kronecker delta property violated");
                    }
                } // end loop over other nodes

            });

        } // end elem_dims
    } // end number of quadrature points

    printf("\nAll kronecker delta checks passed.\n");

}
MATAR_FINALIZE();

} // end main

