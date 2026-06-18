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
#ifndef SEARCH_UTILS_H
#define SEARCH_UTILS_H

#include <cmath>
#include "matar.h"

using namespace mtr;


#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))


struct bin_keys_t{
    size_t i,j,k;
};


// This function is for seriallizing the integation locations in the reference element
KOKKOS_INLINE_FUNCTION
size_t get_id_of_ijk(size_t i, size_t j, size_t k, size_t num_x, size_t num_y){
    return i + (j + k*num_y)*num_x;
}


KOKKOS_INLINE_FUNCTION
bin_keys_t get_bin_keys(const double x_pt, 
                        const double y_pt, 
                        const double z_pt){
            

    double i_dbl = fmax(0, round((x_pt - X0 - bin_dx*0.5)/bin_dx - 1.0e-10)); // x = ih + X0 + dx_bin*0.5
    double j_dbl = fmax(0, round((y_pt - Y0 - bin_dy*0.5)/bin_dy - 1.0e-10));
    double k_dbl = fmax(0, round((z_pt - Z0 - bin_dz*0.5)/bin_dz - 1.0e-10));

    bin_keys_t bin_keys; // save i,j,k to the bin keys

    // get the integer for the bins
    bin_keys.i = (size_t)i_dbl;
    bin_keys.j = (size_t)j_dbl;
    bin_keys.k = (size_t)k_dbl;

    return bin_keys;

} // end function


KOKKOS_INLINE_FUNCTION
size_t get_bin_gid(const double x_pt, 
                   const double y_pt, 
                   const double z_pt, 
                   const double xmin,
                   const double ymin,
                   const double zmin,
                   const double bin_dx,
                   const double bin_dy,
                   const double bin_dz,
                   const size_t num_bins_x,
                   const size_t num_bins_y,
                   const size_t num_bins_z){
            

    double i_dbl = fmin(num_bins_x-1, fmax(0.0, round((x_pt - xmin)/bin_dx - 1.0e-8))); // x = ih + Xmin
    double j_dbl = fmin(num_bins_y-1, fmax(0.0, round((y_pt - ymin)/bin_dy - 1.0e-8)));
    double k_dbl = fmin(num_bins_z-1, fmax(0.0, round((z_pt - zmin)/bin_dz - 1.0e-8)));

    // get the integers for the bins
    size_t i = (size_t)i_dbl;
    size_t j = (size_t)j_dbl;
    size_t k = (size_t)k_dbl;
    
    // get the 1D index for this bin                               
    return get_id_of_ijk(i, j, k, num_bins_x, num_bins_y);

} // end function




// 
//
//
void get_bounds_STL(double &xmin, double &ymin, double &zmin,
                    double &xmax, double &ymax, double &zmax,
                    const DViewCArrayKokkos <double> &v0X,
                    const DViewCArrayKokkos <double> &v0Y,
                    const DViewCArrayKokkos <double> &v0Z,
                    const DViewCArrayKokkos <double> &v1X,
                    const DViewCArrayKokkos <double> &v1Y,
                    const DViewCArrayKokkos <double> &v1Z,
                    const DViewCArrayKokkos <double> &v2X,
                    const DViewCArrayKokkos <double> &v2Y,
                    const DViewCArrayKokkos <double> &v2Z){

    // ----
    // find (xmin, ymin, zmin) and (xmax, ymax, zmax) for building bin mesh
    double xmin_lcl, ymin_lcl, zmin_lcl;
    double xmin, ymin, zmin;

    // no MATAR reductions for multiple values yet, thus hard coding the kokkos loop
    Kokkos::parallel_reduce(
        "stl_min_domain_extents",
        num_inp_triangles,
    
        // this is the for loop coding
        KOKKOS_LAMBDA(const int tri,         
                      double& xmin_lcl,
                      double& ymin_lcl,
                      double& zmin_lcl) 
        {
            // x
            xmin_lcl = fmin(v2X(tri),
                       fmin(v1X(tri),
                       fmin(v0X(tri), xmin_lcl)));

            // y
            ymin_lcl = fmin(v2Y(tri),
                       fmin(v1Y(tri),
                       fmin(v0Y(tri), ymin_lcl)));

            // z
            zmin_lcl = fmin(v2Z(tri),
                       fmin(v1Z(tri),
                       fmin(v0Z(tri), zmin_lcl)));
        },
    Kokkos::Min<double>(xmin), Kokkos::Min<double>(ymin), Kokkos::Min<double>(zmin)); 
    // end parallel reduction over all STL triangles

    xmin -= 1.e-6; // decrease by a small fraction
    ymin -= 1.e-6; // decrease by a small fraction
    zmin -= 1.e-6; // decrease by a small fraction


    double xmax_lcl, ymax_lcl, zmax_lcl;
    double xmax, ymax, zmax;

     // no MATAR reductions for multiple values yet, thus hard coding the kokkos loop
    Kokkos::parallel_reduce(
        "stl_max_domain_extents",
        num_inp_triangles,
    
        // this is the for loop coding
        KOKKOS_LAMBDA(const int tri,         
                      double& xmax_lcl,
                      double& ymax_lcl,
                      double& zmax_lcl) 
        {
            // x                
            xmax_lcl = fmax(v2X(tri),
                       fmax(v1X(tri),
                       fmax(v0X(tri), xmax_lcl)));

            // y               
            ymax_lcl = fmax(v2Y(tri),
                       fmax(v1Y(tri),
                       fmax(v0Y(tri), ymax_lcl)));
                    
            // z                
            zmax_lcl = fmax(v2Z(tri),
                       fmax(v1Z(tri),
                       fmax(v0Z(tri), zmax_lcl)));
  
        },
    Kokkos::Max<double>(xmax), Kokkos::Max<double>(ymax), Kokkos::Max<double>(zmax)); 
    // end parallel reduction over all STL triangles

    xmax += 1.e-6; // increase by a small fraction
    ymax += 1.e-6; // increase by a small fraction
    zmax += 1.e-6; // increase by a small fraction

    return;

} // end function




struct Hash_t{

    // bin mesh memory
    size_t num_bins_x;
    size_t num_bins_y;
    size_t num_bins_z;
    size_t num_bins;

    double xmin;
    double ymin;
    double zmin;

    double xmax;
    double ymax;
    double zmax;

    // 
    double search_radius;

    // bins and their connectivity to each other and points
    DCArrayKokkos <bin_keys_t> keys_in_bin;
    DCArrayKokkos <size_t> num_points_in_bin;
    DRaggedRightArrayKokkos <size_t> points_in_bin;

    // connectivity from points to bins
    DCArrayKokkos <size_t> points_bin_gid;
    CArrayKokkos <size_t>  points_bin_lid_storage;  // only used to create storage
    DCArrayKokkos <int> points_bin_stencil;   // how imin,imax,jmin,jmax,kmin,kmax range for bins in stencil
    DCArrayKokkos <size_t> points_num_neighbors;


    // connectivity between tri corners to continious nodal mesh


    // make a bin mesh to find pairs or build connectivity
    void build_bin_mesh(const double xmin, const double ymin, const double zmin,
                        const double xmax, const double ymax, const double zmax,
                        const size_t num_bins_x_in,
                        const size_t num_bins_y_in,
                        const size_t num_bins_z_in){


        this->num_bins_x = num_bins_x_in;   
        this->num_bins_y = num_bins_y_in;   
        this->num_bins_z = num_bins_z_in;   
        this->num_bins = num_bins_x_in*num_bins_y_in*num_bins_z_in;

        this->bin_dx = (xmax - xmin)/(num_bins_x_in - 1.0);
        this->bin_dy = (ymax - ymin)/(num_bins_y_in - 1.0);
        this->bin_dz = (zmax - zmin)/(num_bins_z_in - 1.0);

        // allocate bin key memory
        keys_in_bin = DCArrayKokkos <bin_keys_t> (num_bins, "keys_in_bin"); // mapping from gid to (i,j,k)
        num_points_in_bin = DCArrayKokkos <size_t> (num_bins, "num_bins");
        num_points_in_bin.set_values(0);

        // build reverse mapping between gid and i,j,k
        FOR_ALL(i, 0, num_bins_x,
                j, 0, num_bins_y,
                k, 0, num_bins_z, {
            

            // get bin gid for this i,j,k
            size_t bin_gid = get_gid(i, j, k, num_bins_x, num_bins_y);

            // the i,j,k for this bin
            bin_keys_t bin_keys;
            bin_keys.i = i;
            bin_keys.j = j;
            bin_keys.k = k;

            // save mapping from bin_gid to bin_keys i,j,k
            keys_in_bin(bin_gid) = bin_keys;

        });
        Kokkos::fence();
        keys_in_bin.update_host();

    } // end function build bin mesh




    void build_multi_node_connectivity(const CArrayKokkos <double>  &corner_point_positions,
                                       const DCArrayKokkos <size_t> node_in_corner_point,
                                       size_t &num_nodes){

        const size_t num_points = corner_point_positions.dims[0];
        
        // allocate arrays for connectivity from points to bins
        points_bin_gid         = DCArrayKokkos <size_t> (num_points, "points_in_gid");
        points_bin_lid_storage = CArrayKokkos  <size_t> (num_points, "bin_lid_storage");  // only used to create storage

        // bins and their connectivity to each other and points
        // Member variable: keys_in_bin(num_bins, "keys_in_bin") is for mapping from gid to (i,j,k) --> now a member var
        // Member variable: num_points_in_bin(num_bins, "num_bins")
        // Member variable: points_in_bin;
        num_points_in_bin.set_values(0); // initializing the member variable to 0

        // save bin id to points
        FOR_ALL(point_gid, 0, num_points, {

            // get the 1D index for this bin
            size_t bin_gid = get_bin_gid(corner_point_positions(point_gid,0), 
                                         corner_point_positions(point_gid,1), 
                                         corner_point_positions(point_gid,2),
                                         num_bins_x, 
                                         num_bins_y,
                                         num_bins_z);

            size_t storage_lid = Kokkos::atomic_fetch_add(&num_points_in_bin(bin_gid), 1);
            points_bin_gid(point_gid) = bin_gid; // the id of the bin
            points_bin_lid_storage(point_gid) = storage_lid; // the storage place in the bin

        }); // end for all
        Kokkos::fence();
        points_bin_gid.update_host();
        num_points_in_bin.update_host();


        // allocate points in bin connectivity
        points_in_bin = DRaggedRightArrayKokkos <size_t> (num_points_in_bin, "num_points_in_bin");

        // save points in bin
        FOR_ALL(point_gid, 0, num_points, {

            // get bin gid
            size_t bin_gid = points_bin_gid(point_gid);

            // get it's storage location in the ragged right compressed storage
            size_t storage_lid = points_bin_lid_storage(point_gid);

            // save the point to this bin
            points_in_bin(bin_gid, storage_lid) = point_gid;

        }); // end for all



        // build the map from corners to node
        node_in_corner_point = CArrayKokkos <size_t> (num_points);     // saves the node_id to point
        smallest_point_id_in_set = CArrayKokkos <size_t> (num_points); // saves the smallest point id in the point-point set

        // ---------------------------------------------------------
        // step 1: get smallest point index in the point-point set
        FOR_ALL(point_gid, 0, num_points, {
           
            size_t bin_gid = points_bin_gid(point_gid);
            bin_keys_t keys = keys_in_bin(bin_gid); // get i,j,k for this bin

            // loop over (i-1:i+1), (j-1:j+1), (k-1:k+1) if node point falls on a bin mesh edge
            const int imin = MAX(0, keys.i-1); const int imax = MIN(num_bins_x-1, keys.i+1);
            const int jmin = MAX(0, keys.j-1); const int jmax = MIN(num_bins_y-1, keys.j+1);
            const int kmin = MAX(0, keys.k-1); const int kmax = MIN(num_bins_z-1, keys.k+1);

            size_t smallest_point_id = point_gid;

            // serial loop over all bins around me, tagging all overlapping points
            for(int i=imin; i<=imax; i++){
                for(int j=jmin; j<=jmax; j++){
                    for(int k=kmin; k<=kmax; k++){
                        
                        // get this bin gid for (i,j,k)
                        size_t a_bin_gid = get_gid(i, j, k, num_bins_x, num_bins_y);

                        // loop over the corner points in this bin
                        for(size_t storage_lid=0; storage_lid<num_points_in_bin(a_bin_gid); storage_lid++){
                            
                            // get the nieghboring point gid 
                            size_t neighbor_point_gid = points_in_bin(a_bin_gid, storage_lid);

                            if(neighbor_point_gid >= smallest_point_id) continue; // only look at lower id's

                            const double pnt_dx = (corner_point_positions(point_gid,0)-corner_point_positions(neighbor_point_gid,0));
                            const double pnt_dy = (corner_point_positions(point_gid,1)-corner_point_positions(neighbor_point_gid,1));
                            const double pnt_dz = (corner_point_positions(point_gid,2)-corner_point_positions(neighbor_point_gid,2));

                            const double distance_sqrd = pnt_dx*pnt_dx + pnt_dy*pnt_dy + pnt_dz*pnt_dz;

                            if(distance_sqrd < 1.e-25)
                            {
                                // this is a matching node
                                if(neighbor_point_gid < smallest_point_id){
                                    smallest_point_id = neighbor_point_gid;
                                } // end if smallest point_id
                            } // end if

                        } // end for over stored points in bin

                    } // end for k
                } // end for j
            } // end for i

            smallest_point_id_in_set(point_gid) = smallest_point_id; // write once, no atomics needed

        }); // end parallel loop over points
        Kokkos::fence();

        // ----------------------------------------------------------------------
        // step 2: count the number of point-point set using smallest point save

        // array holding node_gid on the device and host side
        DCArrayKokkos <size_t> node_gid_counter(1);
        RUN({
            node_gid_counter(0) = 0; // initializing on device to be = 0
        });

        // tally the unique point-point sets using the smallest_point_id as a flag
        FOR_ALL(point_gid, 0, num_points, {

            if (smallest_point_id_in_set(point_gid) == point_gid){
                // this point_id the smallest in the point-point set, so it is a unique node
                const size_t node_gid = Kokkos::atomic_fetch_add(&node_gid_counter(0), 1); 

                node_in_corner_point(point_gid) = node_gid;
            } // end if

        }); // end for all
        Kokkos::fence();
        node_gid_counter.update_host();

        // save the number of nodes, these overlap the unique point-point sets
        num_nodes = node_gid_counter.host(0);
        

        // ---------------------------------------------------------------------------
        // step 3: save the node_gid to the rest of the points in the point-point set

        // save the node_gid to all points in the point-point set
        FOR_ALL(point_gid, 0, num_points, {

            // if not the smallest, then save the node_gid saved to the smallest point
            if (smallest_point_id_in_set(point_gid) != point_gid){

                const size_t smallest_point_gid = smallest_point_id_in_set(point_gid);
                node_in_corner_point(point_gid) = node_in_corner_point( smallest_point_gid);
            } // end if

        }); // end for all


    }// end function




    void build_point_cloud_connectivity(const DCArrayKokkos <double> &point_positions,
                                        DRaggedRightArrayKokkos <size_t> &points_in_point,
                                        DCArrayKokkos <size_t> &points_num_neighbors,
                                        DRaggedRightArrayKokkos <size_t> &reverse_neighbor_lid
                                    ){

        const size_t num_points = corner_point_positions.dims[0];

        // allocate arrays for connectivity from points to bins
        points_bin_gid         = DCArrayKokkos <size_t> (num_points, "points_in_gid");
        points_bin_lid_storage = CArrayKokkos  <size_t> (num_points, "bin_lid_storage");  // only used to create storage
        points_bin_stencil     = DCArrayKokkos <int>    (num_points, 6, "bin_stencil");   // how imin,imax,jmin,jmax,kmin,kmax range for bins in stencil
        points_num_neighbors   = DCArrayKokkos <size_t> (num_points, "num_neighbors");

        // bins and their connectivity to each other and points
        // Member variable: keys_in_bin(num_bins, "keys_in_bin") is for mapping from gid to (i,j,k)
        // Member variable: num_points_in_bin(num_bins, "num_bins")
        // Member variable: points_in_bin;
        num_points_in_bin.set_values(0); // initializing the member variable to 0
        

        // save bin id to points
        FOR_ALL(point_gid, 0, num_points, {

            // get the 1D index for this bin
            size_t bin_gid = get_bin_gid(point_positions(point_gid,0), 
                                         point_positions(point_gid,1), 
                                         point_positions(point_gid,2),
                                         num_bins_x, 
                                         num_bins_y,
                                         num_bins_z);

            size_t storage_lid = Kokkos::atomic_fetch_add(&num_points_in_bin(bin_gid), 1);
            points_bin_gid(point_gid) = bin_gid; // the id of the bin
            points_bin_lid_storage(point_gid) = storage_lid; // the storage place in the bin

        }); // end for all
        Kokkos::fence();
        points_bin_gid.update_host();
        num_points_in_bin.update_host();


        // allocate points in bin connectivity
        points_in_bin = DRaggedRightArrayKokkos <size_t> (num_points_in_bin, "num_points_in_bin");

        // save points in bin
        FOR_ALL(point_gid, 0, num_points, {

            // get bin gid
            size_t bin_gid = points_bin_gid(point_gid);

            // get it's storage location in the ragged right compressed storage
            size_t storage_lid = points_bin_lid_storage(point_gid);

            // save the point to this bin
            points_in_bin(bin_gid, storage_lid) = point_gid;

        }); // end for all


        // ------------------------------------------------
        // Find the neighbors around each point using bins
        // ------------------------------------------------
        
        FOR_ALL(point_gid, 0, num_points, {

            // get bin gid
            size_t bin_gid = points_bin_gid(point_gid);
            
            // get i,j,k for this bin
            bin_keys_t bin_keys = keys_in_bin(bin_gid);

            // loop over neighboring bins
            size_t num_points_found;

            // establish the stencil size to get enough particles
            for(int stencil=1; stencil<100000; stencil++){

                num_points_found = 0;

                const int i = bin_keys.i;
                const int j = bin_keys.j;
                const int k = bin_keys.k;

                const int imin = MAX(0, i-stencil);
                const int imax = MIN(num_bins_x-1, i+stencil);

                const int jmin = MAX(0, j-stencil);
                const int jmax = MIN(num_bins_y-1, j+stencil);

                const int kmin = MAX(0, k-stencil);
                const int kmax = MIN(num_bins_z-1, k+stencil);
                    
                for (int kcount=kmin; kcount<=kmax; kcount++){
                    for (int jcount=jmin; jcount<=jmax; jcount++) {
                        for (int icount=imin; icount<=imax; icount++){

                            // get bin neighbor gid 
                            size_t neighbor_bin_gid = get_gid(icount, jcount, kcount, num_bins_x, num_bins_y);
                            num_points_found += num_points_in_bin(neighbor_bin_gid);

                        } // end for kcount
                    } // end for jcount
                } // end for icount

                // the min number of points required to solve the system is num_poly_basis+1, was 2*num_poly_basis
                if (num_points_found >= num_points_fit  || num_points_found==num_points){

                    const double x_pt_middle = bin_dx*((double)i) + X0; 
                    const double y_pt_middle = bin_dy*((double)j) + Y0; 
                    const double z_pt_middle = bin_dz*((double)k) + Z0; 

                    const double x_pt_minus = bin_dx*((double)imin) + X0; 
                    const double y_pt_minus = bin_dy*((double)jmin) + Y0; 
                    const double z_pt_minus = bin_dz*((double)kmin) + Z0; 
                    
                    const double x_pt_plus = bin_dx*((double)imax) + X0; 
                    const double y_pt_plus = bin_dy*((double)jmax) + Y0; 
                    const double z_pt_plus = bin_dz*((double)kmax) + Z0; 

                    const double dist_minus = sqrt( (x_pt_minus - x_pt_middle)*(x_pt_minus - x_pt_middle) +
                                                    (y_pt_minus - y_pt_middle)*(y_pt_minus - y_pt_middle) +
                                                    (z_pt_minus - z_pt_middle)*(z_pt_minus - z_pt_middle) );

                    const double dist_plus = sqrt( (x_pt_plus - x_pt_middle)*(x_pt_plus - x_pt_middle) +
                                                   (y_pt_plus - y_pt_middle)*(y_pt_plus - y_pt_middle) +
                                                   (z_pt_plus - z_pt_middle)*(z_pt_plus - z_pt_middle) );

                    //printf("h = %f, dist_m = %f, dist_p = %f, num_points=%zu, imin = %d, imax = %d, jmin = %d,  jmax = %d, kmin = %d,  kmax = %d, \n", 
                    //         h_kernel, dist_minus, dist_plus, num_points_found, imin, imax, jmin, jmax, kmin, kmax);

                    // only exit when we exceed kernel distance
                    if (dist_minus >= h_kernel || dist_plus >= h_kernel || num_points_found==num_points){

                        //printf("exiting \n\n");

                        points_bin_stencil(point_gid,0) = imin;
                        points_bin_stencil(point_gid,1) = imax;
                        points_bin_stencil(point_gid,2) = jmin;
                        points_bin_stencil(point_gid,3) = jmax;
                        points_bin_stencil(point_gid,4) = kmin;
                        points_bin_stencil(point_gid,5) = kmax;

                        points_num_neighbors(point_gid) = num_points_found; // including node_i in the list of neighbors
                        points_num_neighbors(point_gid) = num_points_found - 1; // the -1 is because we counted point i as a neighbor

                        break;
                    }
                    // else increase stencil size


                } // end of check
                
            } // end for stencil


        }); // end for all
        Kokkos::fence();
        points_bin_stencil.update_host();



        // account for stencels not overlapping, fixing assymetry in points connectivity
        FOR_ALL(point_gid, 0, num_points, {

            // get bin gid for this point
            size_t bin_gid = points_bin_gid(point_gid);
                    
            // get i,j,k for this bin
            bin_keys_t bin_keys = keys_in_bin(bin_gid);

            const int i = bin_keys.i;
            const int j = bin_keys.j;
            const int k = bin_keys.k;

            // walk over the stencil to get neighbors of this bin
            const int imin = points_bin_stencil(point_gid,0);
            const int imax = points_bin_stencil(point_gid,1);
            const int jmin = points_bin_stencil(point_gid,2);
            const int jmax = points_bin_stencil(point_gid,3);
            const int kmin = points_bin_stencil(point_gid,4);
            const int kmax = points_bin_stencil(point_gid,5);

            // loop over my bin stencil
            for (int kcount=kmin; kcount<=kmax; kcount++){
                for (int jcount=jmin; jcount<=jmax; jcount++) {
                    for (int icount=imin; icount<=imax; icount++){

                        // get bin neighbor gid 
                        size_t neighbor_bin_gid = get_gid(icount, jcount, kcount, num_bins_x, num_bins_y);

                        // save the points in this bin
                        for(size_t neighbor_pt_lid=0; neighbor_pt_lid<num_points_in_bin(neighbor_bin_gid); neighbor_pt_lid++){

                            size_t neighbor_point_gid = points_in_bin(neighbor_bin_gid, neighbor_pt_lid);

                            // check if the point-point pairs have identical, overlapping stencils, if not, increment the number of neighbors
                            const int neighbor_imin = points_bin_stencil(neighbor_point_gid,0);
                            const int neighbor_imax = points_bin_stencil(neighbor_point_gid,1);
                            const int neighbor_jmin = points_bin_stencil(neighbor_point_gid,2);
                            const int neighbor_jmax = points_bin_stencil(neighbor_point_gid,3);
                            const int neighbor_kmin = points_bin_stencil(neighbor_point_gid,4);
                            const int neighbor_kmax = points_bin_stencil(neighbor_point_gid,5);
                            
                            // i,j,k is the bin where point_gid lives
                            bool inside =
                                (i >= neighbor_imin && i <= neighbor_imax) &&
                                (j >= neighbor_jmin && j <= neighbor_jmax) &&
                                (k >= neighbor_kmin && k <= neighbor_kmax);

                            if(!inside){
                                Kokkos::atomic_increment(&points_num_neighbors(neighbor_point_gid)); 
                                // the other stencil didn't see my point because it was smaller, now it does see it
                            }

                        } // neighbor_point_lid

                    } // end for kcount
                } // end for jcount
            } // end for icount        

        }); // end for all
        Kokkos::fence();
        points_num_neighbors.update_host();


        for(size_t point_gid=0; point_gid<num_points; point_gid++){
            printf("point_gid = %zu, num_neighbors = %zu, num_points = %zu \n", point_gid, points_num_neighbors.host(point_gid), num_points);
        }
        
        // allocate memory for points in point
        points_in_point = DRaggedRightArrayKokkos <size_t> (points_num_neighbors, "points_in_point");
        points_num_neighbors.set_values(0);  // this is a num saved counter now



        // ---------------------
        // Save the neighbors
        // ---------------------

        // find neighbors using bins
        FOR_ALL(point_gid, 0, num_points, {

            // get bin gid for this point
            size_t bin_gid = points_bin_gid(point_gid);
                    
            // get i,j,k for this bin
            bin_keys_t bin_keys = keys_in_bin(bin_gid);

            const int i = bin_keys.i;
            const int j = bin_keys.j;
            const int k = bin_keys.k;

            // walk over the stencil to get neighbors
            int imin = points_bin_stencil(point_gid,0);
            int imax = points_bin_stencil(point_gid,1);
            int jmin = points_bin_stencil(point_gid,2);
            int jmax = points_bin_stencil(point_gid,3);
            int kmin = points_bin_stencil(point_gid,4);
            int kmax = points_bin_stencil(point_gid,5);


            for (int kcount=kmin; kcount<=kmax; kcount++){
                for (int jcount=jmin; jcount<=jmax; jcount++) {
                    for (int icount=imin; icount<=imax; icount++){

                        // get bin neighbor gid 
                        size_t neighbor_bin_gid = get_gid(icount, jcount, kcount, num_bins_x, num_bins_y);

                        // save the points in this bin
                        for(size_t neighbor_pt_lid=0; neighbor_pt_lid<num_points_in_bin(neighbor_bin_gid); neighbor_pt_lid++){

                            size_t neighbor_point_gid = points_in_bin(neighbor_bin_gid, neighbor_pt_lid);
                            
                            // make sure its a neighbor
                            if(neighbor_point_gid != point_gid){

                                // save the neighbor
                                size_t num_saved = Kokkos::atomic_fetch_add(&points_num_neighbors(point_gid), 1);
                                points_in_point(point_gid, num_saved) = neighbor_point_gid;
                                
                                
                                // if point j's stencil did not see point i, then save i to j's list
                                const int neighbor_imin = points_bin_stencil(neighbor_point_gid,0);
                                const int neighbor_imax = points_bin_stencil(neighbor_point_gid,1);
                                const int neighbor_jmin = points_bin_stencil(neighbor_point_gid,2);
                                const int neighbor_jmax = points_bin_stencil(neighbor_point_gid,3);
                                const int neighbor_kmin = points_bin_stencil(neighbor_point_gid,4);
                                const int neighbor_kmax = points_bin_stencil(neighbor_point_gid,5);

                                // i,j,k is the bin where point_gid lives
                                bool inside =
                                    (i >= neighbor_imin && i <= neighbor_imax) &&
                                    (j >= neighbor_jmin && j <= neighbor_jmax) &&
                                    (k >= neighbor_kmin && k <= neighbor_kmax);

                                if(!inside){

                                    size_t num_saved_neighbor = Kokkos::atomic_fetch_add(&points_num_neighbors(neighbor_point_gid), 1);
                                    points_in_point(neighbor_point_gid, num_saved_neighbor) = point_gid;
                                    // the other stencil didn't see my point because it was smaller, now it does see it

                                } // end if

                            } // end if neighbor != point_gid

                        } // neighbor_point_lid

                    } // end for kcount
                } // end for jcount
            } // end for icount        

        }); // end for all
        Kokkos::fence();
        points_in_point.update_host();


        // build the reverse map
        reverse_neighbor_lid = DRaggedRightArrayKokkos <size_t> (points_num_neighbors); 

        FOR_ALL(point_gid, 0, num_points, {
                
            for(int neighbor_point_lid = 0; neighbor_point_lid<points_num_neighbors(point_gid); neighbor_point_lid++){
                
                // get the point gid for this neighbor
                int neighbor_point_gid = points_in_point(point_gid, neighbor_point_lid);
                
                // loop over the neighbors of my neighbor
                size_t found = 0;
                for(int j_lid = 0; j_lid<points_num_neighbors(neighbor_point_gid); j_lid++){

                    // get the neighboring point gid of my neighbor
                    int j_point_gid = points_in_point(neighbor_point_gid, j_lid);
                    if (point_gid == j_point_gid){
                        reverse_neighbor_lid(point_gid, neighbor_point_lid) = j_lid;
                        found = 1;
                        //printf("found \n");
                        break;
                    }
                } // end loop over j's neighboring points
                if(found==0)printf("reverse map for i=%d and j=%d pair not found \n", point_gid, neighbor_point_gid);
            } // end loop over i's neighboring points
                
        });


    } // end function build point connectivity





} // Hash_t


void build_bin_mesh(){

    // the number of nodes in the bin mesh
    const double num_avg_inp_triangles_1D = pow((double)num_inp_triangles, 0.33333333);

    const size_t num_bins_x = (size_t)num_avg_inp_triangles_1D;   
    const size_t num_bins_y = (size_t)num_avg_inp_triangles_1D;   
    const size_t num_bins_z = (size_t)num_avg_inp_triangles_1D;   

    const double bin_dx = (xmax - xmin)/(num_avg_inp_triangles_1D - 1.0);
    const double bin_dy = (ymax - ymin)/(num_avg_inp_triangles_1D - 1.0);
    const double bin_dz = (zmax - zmin)/(num_avg_inp_triangles_1D - 1.0);

    const size_t num_bins = num_bins_x*num_bins_y*num_bins_z; // for ease of use
    DCArrayKokkos <size_t> num_tris_in_bin(num_bins);
    num_tris_in_bin.set_values(0.0);

    
    FOR_ALL(tri_gid, 0, num_inp_triangles, {

        // get the 1D index for this bin
        size_t bin_gid = get_bin_gid(tri_coords(tri_gid,0), 
                                     tri_coords(tri_gid,1), 
                                     tri_coords(tri_gid,2),
                                     xmin,
                                     ymin,
                                     zmin,
                                     bin_dx,
                                     bin_dy,
                                     bin_dz,
                                     num_bins_x, 
                                     num_bins_y,
                                     num_bins_z);

    });


}



#endif