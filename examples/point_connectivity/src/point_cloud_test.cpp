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

bool VERBOSE = false; // flag for printing connectivity

int main(int argc, char** argv) {

MATAR_INITIALIZE(argc, argv);
{ // MATAR scope
    std::cout<<"--- Point cloud connectivity test ---- \n";


    // --------------------------------------------
    // Create point cloud
    // --------------------------------------------

    // number of points in each direction, creating structured point cloud
    const size_t N = 3; 
    const double h = 2.0/((double)N - 1);

    const size_t num_points = N*N*N;

    DCArrayKokkos <double> point_positions(num_points, 3, "point_positions");

    // save point locations
    size_t point_gid = 0;  
    for(size_t k=0; k<N; k++)
    for(size_t j=0; j<N; j++)
    for(size_t i=0; i<N; i++){
        point_positions.host(point_gid, 0) = ((double)i)*h;
        point_positions.host(point_gid, 1) = ((double)j)*h;
        point_positions.host(point_gid, 2) = ((double)k)*h;
        point_gid++;
    } // end for

    point_positions.update_device();


    SpatialConnectivity_t sc;

    sc.build_bin_mesh(-0.5,-0.5,-0.5, 2.5,2.5,2.5, 6,6,6);

    const size_t min_num_points = 1; 
    sc.initialize_point_cloud_vars(1.5*h, min_num_points);


    sc.build_point_cloud_connectivity(point_positions);


    printf("\n--- point-point connectivity symmetry test ---\n");

    RUN({
        for (size_t i = 0; i < num_points; i++)
        for (size_t lid = 0; lid < sc.points_num_neighbors(i); lid++) {
            size_t j = sc.points_in_point(i, lid);
            bool found = false;
            for (size_t jlid = 0; jlid < sc.points_num_neighbors(j); jlid++)
                if (sc.points_in_point(j, jlid) == i) { found = true; break; }

            if (found==false) Kokkos::abort("asymmetry bug \n");  
        } // end for
    });

    printf("\n--- point-point connectivity reverse map test ---\n");

    RUN({
        for (size_t i = 0; i < num_points; i++)
        for (size_t lid = 0; lid <sc.points_num_neighbors(i); lid++) {
            size_t j    = sc.points_in_point(i, lid);
            size_t jlid = sc.reverse_neighbor_lid(i, lid);
            if (sc.points_in_point(j, jlid) != i) Kokkos::abort("reverse map bug \n");  // reverse map points back to i
        }
    });

    printf("\n--- point-point connectivity no duplicates test ---\n");

    RUN({
        for (size_t i = 0; i < num_points; i++) {
            
            size_t num_unique_nodes = 0;
            size_t unique_nodes[N*N*N];  // holds points seen

            for (size_t lid = 0; lid < sc.points_num_neighbors(i); lid++) {

                size_t j = sc.points_in_point(i, lid);
                if(VERBOSE) printf("point %zu neighbors point %zu \n", j, i);
                
                if(j == i) Kokkos::abort("self seen in neighbor list \n");  // self should never appear

                // check for duplicates, if true then its a duplicate
                bool found = false;

                // loop over all neighboring nodes to see if the node already appeared
                for(size_t node_lid = 0; node_lid < num_unique_nodes; node_lid++){
                    if (j==unique_nodes[node_lid]) found = true;
                } // end for

                if(found==false){
                    unique_nodes[num_unique_nodes]=j;
                    num_unique_nodes++;
                } 

                if(found==true) Kokkos::abort("duplicate neighbor \n"); 
                
            } // end lid
            if(VERBOSE) printf("\n");
        } // end i
    });


    printf("\n--- point-point connectivity number of neighbors test ---\n");

    RUN({

        // center point of 3x3x3 grid is index 13
        if (sc.points_num_neighbors(13) != 26) Kokkos::abort("wrong number of neighbors at center \n");

        // corner point (index 0) should have 2x2x2-1 neighbors
        if (sc.points_num_neighbors(0) != 7) Kokkos::abort("wrong number of neighbors at corner \n");

        for (size_t i = 0; i < num_points; i++) {
            printf("number of neighbors around point %zu = %zu\n", i, sc.points_num_neighbors(i));
        } // end i
    });
    
    printf("\nAll point cloud checks passed.\n");
    
} // end MATAR scope
MATAR_FINALIZE();

return 0;
}