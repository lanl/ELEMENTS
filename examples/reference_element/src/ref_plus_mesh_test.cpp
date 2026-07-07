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

// This pulls in kokkos, matar, mesh, ref_elem stuff, and PT-Scotch
#include "ELEMENTS.h"


using namespace mtr;
using namespace swage;    // unstructured mesh and point cloud
using namespace elements; // reference element space

int main(int argc, char** argv) {

MATAR_INITIALIZE(argc, argv);
{ // MATAR scope
    std::cout<<"Reference Element plus Mesh Example!"<<std::endl;
    
    Quadrature_t Quad;
    ReferenceElement_t FERefElem; // kinematic space
    ReferenceElement_t DGRefElem; // thermal space, it is discontinous

    Mesh_t Mesh; // unstructured mesh

    const size_t elem_dims = 3;
    const size_t elem_order = 3;  // cubic element
    const size_t num_elems_1D = 2;

    // the minimum quadrature for FE hydrodynamics based on elem order
    const size_t num_DOFs_1d = elem_order + 1; // 
    const size_t num_qpts_1d = 2*elem_order;   // using Legendre, but if using Lobatto, it requires 2*Order+1

    // ============================================
    // Create quadrature and reference element

    // create quadrature
    Quad.initialize_quadrature(reference_space::GaussLegendre,
                               num_qpts_1d,
                               elem_dims);

    // p_order is the basis order for the Lagrange polynomial defining the element
    FERefElem.initialize_ref_elem(reference_space::arbitraryOrderElement,
                                  reference_space::LagrangeLobatto,
                                  Quad,
                                  elem_order);    


    // for thermal DOFs, we use p_order-1 for the basis order of the Lagrange polynomial, it is DG
    DGRefElem.initialize_ref_elem(reference_space::arbitraryOrderElement,
                                  reference_space::LagrangeLegendre,
                                  Quad,
                                  elem_order-1); 


    // ==========================================
    // Build the unstructured mesh structure

    const size_t num_elems    = num_elems_1D*num_elems_1D*num_elems_1D;
    const size_t num_nodes_1D = elem_order*num_elems_1D + 1;  // number of nodes
    const size_t num_nodes    = num_nodes_1D*num_nodes_1D*num_nodes_1D;

    Mesh.initialize_dims(elem_dims);
    Mesh.initialize_elems_Pn(num_elems, elem_order, Quad.num_qpts_1d);
    Mesh.initialize_nodes(num_nodes);
    
    DCArrayKokkos <double> node_coords(Mesh.num_nodes, Mesh.num_dims);
    const double h = 1.0/((double)num_nodes_1D);

    // create indexing for a Pn order mesh
    FOR_ALL(i,0,num_elems_1D,
            j,0,num_elems_1D,
            k,0,num_elems_1D,{

        size_t elem_gid = i + (j+k*num_elems_1D)*num_elems_1D;

        size_t node_lid = 0;
        for(size_t ic=i; ic<=i+elem_order; ic++)
        for(size_t jc=j; jc<=j+elem_order; jc++)
        for(size_t kc=k; kc<=k+elem_order; kc++){
            size_t node_gid = ic + (jc+kc*num_nodes_1D)*num_nodes_1D;
            Mesh.nodes_in_elem(elem_gid,node_lid) = node_gid;
            node_lid++;

            node_coords(node_gid,0) = (double)ic*h;
            node_coords(node_gid,1) = (double)jc*h;
            node_coords(node_gid,2) = (double)kc*h;
        } // end for 

    }); // end parallel for

    Mesh.build_corner_connectivity();
    Mesh.build_elem_elem_connectivity();
    Mesh.build_patch_connectivity();

    // check mesh index sizes
    if(Mesh.num_nodes!=num_nodes){
        printf("num nodes = %zu and mesh.num_nodes = %zu", num_nodes, Mesh.num_nodes);
        Kokkos::abort("ERROR: wrong number of mesh nodes");
    }
    if(Mesh.num_gauss_in_elem!=Quad.num_qpts_in_elem){
        Kokkos::abort("ERROR: wrong number of Gauss points in elem");
    }

    // ==========================================
    // Create state on the unstructured mesh structure

    DCArrayKokkos <double> node_scalar(Mesh.num_nodes);
    DCArrayKokkos <double> gauss_scalar(Mesh.num_gauss_in_elem);


    printf("\nReference plus mesh test finished.\n");


} // end MATAR scope
MATAR_FINALIZE();

return 0;
}
