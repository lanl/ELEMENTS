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

// This pulls in kokkos, matar, mesh, PointCloud, ref_elem stuff, and PT-Scotch
#include "ELEMENTS.h"

//#undef NDEBUG     // Ensures NDEBUG is turned off

using namespace mtr;
using namespace swage; // unstructured mesh and PointCloud



int main(int argc, char** argv) {

MATAR_INITIALIZE(argc, argv);
{ // MATAR scope
    std::cout<<"--- Points in box test ---- \n";

    // --------------------------------------------
    // Create point cloud
    // --------------------------------------------

    // number of points in each direction, creating structured point cloud
    const size_t N = 4; 
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


    // --------------------------------------------
    // Create PointCloud class
    // --------------------------------------------
    PointCloud_t pc;

    // build the bin mesh, ensuring it covers point cloud domain
    double xmin; double ymin; double zmin;
    double xmax; double ymax; double zmax;
    pc.get_bounds_point_cloud(xmin, ymin, zmin,
                              xmax, ymax, zmax,
                              point_positions);
    pc.build_bin_mesh(xmin-0.5,ymin-0.5,zmin-0.5, xmax+0.5,ymax+0.5,zmax+0.5, 6,6,6);

    printf("\n---bounds test ---\n");

    printf("min = (%f,%f,%f) \n", xmin,ymin,zmin);
    printf("max = (%f,%f,%f) \n", xmax,ymax,zmax);
    if (fabs(xmin+0.00001) <= 1e-12) throw std::runtime_error("ERROR: xmin bound is wrong\n");
    if (fabs(ymin+0.00001) <= 1e-12) throw std::runtime_error("ERROR: ymin bound is wrong\n");
    if (fabs(ymin+0.00001) <= 1e-12) throw std::runtime_error("ERROR: ymin bound is wrong\n");

    if (fabs(xmax-2.00001) <= 1e-12) throw std::runtime_error("ERROR: xmax bound is wrong\n");
    if (fabs(ymax-2.00001) <= 1e-12) throw std::runtime_error("ERROR: ymax bound is wrong\n");
    if (fabs(ymax-2.00001) <= 1e-12) throw std::runtime_error("ERROR: ymax bound is wrong\n");
    

    // declare array holding points in box, it is allocated & populated in PointCloud_t
    RaggedRightArrayKokkos <size_t> points_in_box;    
    

    // define boxes to capture points
    const size_t num_boxes = 2;
    CArrayKokkos <double> bounding_box(num_boxes,2,3); // always these dims
    DCArrayKokkos <size_t> num_points_in_box(num_boxes); //always this dim   

    RUN({
        // first box corners = (0,0,0) & (1,1,1)
        bounding_box(0,0,0) = 0.0; // xmin
        bounding_box(0,0,1) = 0.0; // ymin
        bounding_box(0,0,2) = 0.0; // zmin

        bounding_box(0,1,0) = 1.0; // xmax
        bounding_box(0,1,1) = 1.0; // ymax
        bounding_box(0,1,2) = 1.0; // zmax

        // second box corners = (0.5,0.5,0.5) & (1.5,1.5,1.5)
        bounding_box(1,0,0) = 0.5; // xmin
        bounding_box(1,0,1) = 0.5; // ymin
        bounding_box(1,0,2) = 0.5; // zmin

        bounding_box(1,1,0) = 1.5; // xmax
        bounding_box(1,1,1) = 1.5; // ymax
        bounding_box(1,1,2) = 1.5; // zmax
    });


    pc.get_points_in_box(point_positions, 
                         points_in_box, 
                         num_points_in_box,
                         bounding_box);



    printf("\n--- points in box test ---\n");

    RUN({
        for (size_t box_gid = 0; box_gid < num_boxes; box_gid++){

            printf("box = %zu\n", box_gid);
            printf("num_points_in_box = %zu \n", num_points_in_box(box_gid));

            for (size_t lid = 0; lid < num_points_in_box(box_gid); lid++) {
                size_t j = points_in_box(box_gid, lid);
                printf("point = %zu \n", j); 
            } // end for
            printf("\n");
        }

       if (points_in_box(0, 0) != 0)  Kokkos::abort("first point in first box is wrong \n"); 
       if (points_in_box(0, 7) != 21) Kokkos::abort("last point in first box is wrong \n"); 
       if (points_in_box(1, 0) != 21) Kokkos::abort("first point in second box is wrong \n"); 
       if (points_in_box(1, 7) != 42) Kokkos::abort("last point in second box is wrong \n");
    });
    
    printf("\nAll points in box checks passed.\n");
    
} // end MATAR scope
MATAR_FINALIZE();

return 0;
}