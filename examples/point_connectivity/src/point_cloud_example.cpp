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

using namespace mtr;
using namespace swage; // unstructured mesh and hash


int main(int argc, char** argv) {

MATAR_INITIALIZE(argc, argv);
{ // MATAR scope
    std::cout<<"Hello, running point-point connectivity examples! \n";


    // --------------------------------------------
    // Create point cloud
    // --------------------------------------------

    // Domain of point cloud
    const double X0 = 0.0;   // origin
    const double Y0 = 0.0;
    const double Z0 = 0.0;

    const double LX = 2.0;   // length 
    const double LY = 2.0;
    const double LZ = 2.0;

    // number of points in each direction, creating structured point cloud
    const size_t num_1d_x = 4;
    const size_t num_1d_y = 4;
    const size_t num_1d_z = 4;

    const size_t num_points = num_1d_x*num_1d_y*num_1d_z;


    const double dx = (LX-X0)/((double)num_1d_x - 1.0);
    const double dy = (LY-Y0)/((double)num_1d_y - 1.0);
    const double dz = (LZ-Z0)/((double)num_1d_z - 1.0);


    DCArrayKokkos <double> point_positions(num_points, 3, "point_positions");

    // save point locations
    size_t point_gid = 0;  
    for(size_t k=0; k<num_1d_z; k++){
        for(size_t j=0; j<num_1d_y; j++){
            for(size_t i=0; i<num_1d_x; i++){
                // x = x0 + i*dx
                point_positions.host(point_gid, 0) = X0 + ((double)i)*dx;
                point_positions.host(point_gid, 1) = Y0 + ((double)j)*dy;
                point_positions.host(point_gid, 2) = Z0 + ((double)k)*dz;
                point_gid++;
            } // end i
        } // end j
    } // end k
    point_positions.update_device();



    // ------------------------------------------------
    // create the spatial connectivity data structure
    // ------------------------------------------------
    SpatialConnectivity_t sc;
    

    // fitting radius of the cloud
    const double h_kernel = fmax(fmax(dx,dy),dz); 
    const double cutoff_coeff = 1.5; // coeff*h_kernel = search radius
    const double min_num_points_fit = 7; 

    // set the data structure variables
    sc.initialize_point_cloud_vars(h_kernel*cutoff_coeff,min_num_points_fit);

    
    // --------------------------------------------
    // create bin mesh for building connectivity
    // --------------------------------------------
    const size_t num_bins_x_in = 6;
    const size_t num_bins_y_in = 6;
    const size_t num_bins_z_in = 6;
    double dx_cloud = LX/((double)num_bins_x_in); // or use search radius
    sc.build_bin_mesh(X0-dx_cloud, 
                      Y0-dx_cloud, 
                      Z0-dx_cloud,
                      LX+X0+dx_cloud, 
                      LY+Y0+dx_cloud, 
                      LZ+Z0+dx_cloud,
                      num_bins_x_in,
                      num_bins_y_in,
                      num_bins_z_in);



    std::cout<<"Building connectivity in a point-cloud \n";




    // build the point cloud connectivity
    sc.build_point_cloud_connectivity(point_positions);

    // verify correctness of the connectivity data structure looping interior points
    FOR_ALL(i, 1, num_1d_x-1,
            j, 1, num_1d_y-1,
            k, 1, num_1d_z-1, {

        size_t point_gid = get_id_of_ijk(i, j, k, num_1d_x, num_1d_y);
        
        printf("Number of neighbors = %zu, correct answer is 26 \n", sc.points_num_neighbors(point_gid));

    }); //end for all

    printf("\n");

    // verify correctness of the connectivity data structure showing points
    RUN({
        for (size_t i = 0; i < num_points; i++) {
            for (size_t lid = 0; lid < sc.points_num_neighbors(i); lid++) {

                size_t j = sc.points_in_point(i, lid);
                printf("point %zu neighbors point %zu \n", j, i);
                
            } // end lid
            printf("\n");
        } // end i
    });


    std::cout<<"Finished"<<std::endl;
    
} // end MATAR scope
MATAR_FINALIZE();

return 0;
}