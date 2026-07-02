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

void verify_partition_of_unity(const Quadrature_t& Quad,
                               const ReferenceElement_t& RefElem) {
    for (size_t qpt = 0; qpt < Quad.num_qpts_in_elem; qpt++) {
        double sum = 0.0;
        
        for (size_t basis = 0; basis < RefElem.num_dofs_in_elem; basis++) {
            sum += RefElem.qpt_basis(qpt, basis);
        }
        if (fabs(sum - 1.0) > 1e-13) {
            printf("Error: partion of unity failed, sum of basis = %f at rid = %zu \n", sum, qpt);
            Kokkos::abort("Partition of unity failed at quadrature point ");
        }
        printf("sum = %f \n", sum);
    } // end loop over qpts
} // end function

int main(int argc, char** argv) {

MATAR_INITIALIZE(argc, argv);
{ // MATAR scope

    printf("\n--- partition unity test ---\n");

    printf("\n--- 1D element ---\n");
    for(size_t num_qpts_1D = 1; num_qpts_1D<21; num_qpts_1D++){

        printf("num quadrature points in 1D = %zu \n", num_qpts_1D);

        Quadrature_t Quad1D;
        ReferenceElement_t FERefElem1D;

        const size_t elem_dims_test = 1;  
        Quad1D.initialize_quadrature(reference_space::GaussLegendre,
                                    num_qpts_1D,
                                    elem_dims_test);

        const size_t p_order = 1; // basis order for Lagrange polynomial
        FERefElem1D.initialize_ref_elem(reference_space::arbitraryOrderElement,
                                        reference_space::LagrangeLobatto,
                                        Quad1D,
                                        p_order);

        verify_partition_of_unity(Quad1D,
                                  FERefElem1D);

        printf("\n");
    } // end loop of num qpts 

    printf("\nAll partition unity checks passed.\n");

}
MATAR_FINALIZE();

} // end main

