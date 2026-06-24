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
#ifndef POINT_CLOUD_H
#define POINT_CLOUD_H

#include <cmath>
#include "matar.h"

using namespace mtr;


#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

namespace swage
{

/////////////////////////////////////////////////////////////////////////////
///
/// \struct bin_keys_t
///
/// \brief Stores hash keys i,j,k
///
/////////////////////////////////////////////////////////////////////////////
struct bin_keys_t{
    int i,j,k;
};


/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_id_of_ijk
///
/// \brief This function is for seriallizing locations in a reference 
///        element or a structured mesh
///
/// \param i is the x direction index
/// \param j is the y direction index
/// \param k is the z direction index
/// \param num_x is the number of locations in the x direction
/// \param num_y is the number of locations in the y direction
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_INLINE_FUNCTION
size_t get_id_of_ijk(size_t i, size_t j, size_t k, size_t num_x, size_t num_y) {
    return i + (j + k*num_y)*num_x;
};


/////////////////////////////////////////////////////////////////////////////
///
/// \struct PointCloud_t
///
/// \brief Builds and stores spatial connectivity for arbitrary point clouds.
///        Provides overlapping node detection and point-to-point neighbor
///        lists for meshfree methods such as SPH and RKPM.  A structured
///        bin mesh is used internally to accelerate all proximity queries.
///
/////////////////////////////////////////////////////////////////////////////
struct PointCloud_t{

    // --- bin mesh ---
    size_t num_bins_x = 0;
    size_t num_bins_y = 0;
    size_t num_bins_z = 0;
    size_t num_bins = 0;

    double bin_dx = 0.0;
    double bin_dy = 0.0;
    double bin_dz = 0.0;

    double xmin = 0.0;
    double ymin = 0.0;
    double zmin = 0.0;

    double xmax = 0.0;
    double ymax = 0.0;
    double zmax = 0.0;

    // --- point cloud parameters ---
    double search_radius = 0.0;     // kernel support radius
    size_t min_num_points_fit = 0;  // minimum neighbors required for basis fit

    // --- bin mesh connectivity ---
    DCArrayKokkos <bin_keys_t> keys_in_bin;        // bin gid -> (i,j,k)
    DCArrayKokkos <size_t>     num_points_in_bin;  // number of points per bin
    DRaggedRightArrayKokkos <size_t> points_in_bin; // points stored in each bin

    // --- point-to-bin connectivity (arrays) ---
    DCArrayKokkos <size_t> points_bin_gid;          // bin gid for each point
    CArrayKokkos  <size_t> points_bin_lid_storage;  // storage slot within bin
    DCArrayKokkos <int>    points_bin_stencil;      // stencil bounds: imin,imax,jmin,jmax,kmin,kmax

    // --- point cloud connectivity (arrays) ---
    DRaggedRightArrayKokkos <size_t> points_in_point;      // point connectivity array
    DCArrayKokkos <size_t> points_num_neighbors;           // neighbor count per point
    DRaggedRightArrayKokkos <size_t> reverse_neighbor_lid; // accessing neighbors connectivity array


    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn initialize_point_cloud_vars
    ///
    /// \brief Sets parameters for point cloud neighbor search
    ///
    /// \param search_radius_in      kernel support radius
    /// \param min_num_points_fit_in minimum neighbors required for basis fit
    ///
    /////////////////////////////////////////////////////////////////////////////
    void initialize_point_cloud_vars(const double search_radius_in,
                                     const size_t min_num_points_fit_in){

        search_radius = search_radius_in;
        min_num_points_fit = min_num_points_fit_in;

        return;

    } // end initialize_point_cloud_vars


    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn build_bin_mesh
    ///
    /// \brief Builds the structured bin mesh used for all spatial queries.
    ///        Must be called before any connectivity function.
    ///
    /// \param xmin_in       minimum x coordinate of the bin mesh domain
    /// \param ymin_in       minimum y coordinate of the bin mesh domain
    /// \param zmin_in       minimum z coordinate of the bin mesh domain
    /// \param xmax_in       maximum x coordinate of the bin mesh domain
    /// \param ymax_in       maximum y coordinate of the bin mesh domain
    /// \param zmax_in       maximum z coordinate of the bin mesh domain
    /// \param num_bins_x_in number of bins in the x direction
    /// \param num_bins_y_in number of bins in the y direction
    /// \param num_bins_z_in number of bins in the z direction
    ///
    /////////////////////////////////////////////////////////////////////////////
    void build_bin_mesh(const double xmin_in, const double ymin_in, const double zmin_in,
                        const double xmax_in, const double ymax_in, const double zmax_in,
                        const size_t num_bins_x_in,
                        const size_t num_bins_y_in,
                        const size_t num_bins_z_in) {

        xmin = xmin_in;  ymin = ymin_in;  zmin = zmin_in;
        xmax = xmax_in;  ymax = ymax_in;  zmax = zmax_in;

        num_bins_x = num_bins_x_in;
        num_bins_y = num_bins_y_in;
        num_bins_z = num_bins_z_in;
        num_bins   = num_bins_x * num_bins_y * num_bins_z;

        bin_dx = fmax(1.e-13, (xmax - xmin) / (double)num_bins_x);
        bin_dy = fmax(1.e-13, (ymax - ymin) / (double)num_bins_y);
        bin_dz = fmax(1.e-13, (zmax - zmin) / (double)num_bins_z);

        keys_in_bin       = DCArrayKokkos <bin_keys_t> (num_bins, "keys_in_bin");
        num_points_in_bin = DCArrayKokkos <size_t>     (num_bins, "num_points_in_bin");
        num_points_in_bin.set_values(0);

        // build reverse mapping between gid and i,j,k
        FOR_ALL_CLASS(i, 0, num_bins_x,
                      j, 0, num_bins_y,
                      k, 0, num_bins_z, {
            
            size_t bin_gid = get_id_of_ijk(i, j, k, num_bins_x, num_bins_y);

            bin_keys_t bk;
            bk.i = i;
            bk.j = j;
            bk.k = k;

            // save mapping from bin_gid to ben keys i,j,k
            keys_in_bin(bin_gid) = bk;

        });
        Kokkos::fence();
        keys_in_bin.update_host();

    }  // end function build bin mesh


    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \struct get_bin_keys
    ///
    /// \brief returns the i,j,k indices for the bin cell center given the (x,y,z)
    ///        coordinates of a point
    ///
    /// \param x_pt is the x coordinate
    /// \param y_pt is the y coordinate
    /// \param z_pt is the z coordinate
    /////////////////////////////////////////////////////////////////////////////
    KOKKOS_INLINE_FUNCTION
    bin_keys_t get_bin_keys(const double x_pt, 
                            const double y_pt, 
                            const double z_pt) {

        bin_keys_t bk;
        bk.i = (int)fmin((double)num_bins_x-1., fmax(0., round((x_pt - xmin)/bin_dx - 1.0e-10)));
        bk.j = (int)fmin((double)num_bins_y-1., fmax(0., round((y_pt - ymin)/bin_dy - 1.0e-10)));
        bk.k = (int)fmin((double)num_bins_z-1., fmax(0., round((z_pt - zmin)/bin_dz - 1.0e-10)));

        return bk;

    } // end get_bin_keys


    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \struct get_bin_gid
    ///
    /// \brief returns the 1D global index for the bin mesh given the (x,y,z)
    ///        coordinates of a point
    ///
    /// \param x_pt is the x coordinate
    /// \param y_pt is the y coordinate
    /// \param z_pt is the z coordinate
    /////////////////////////////////////////////////////////////////////////////
    KOKKOS_INLINE_FUNCTION
    size_t get_bin_gid(const double x_pt, 
                       const double y_pt, 
                       const double z_pt) const {

        size_t i = (size_t)fmin((double)num_bins_x-1., fmax(0., round((x_pt - xmin)/bin_dx - 1.0e-10)));
        size_t j = (size_t)fmin((double)num_bins_y-1., fmax(0., round((y_pt - ymin)/bin_dy - 1.0e-10)));
        size_t k = (size_t)fmin((double)num_bins_z-1., fmax(0., round((z_pt - zmin)/bin_dz - 1.0e-10)));

        return get_id_of_ijk(i, j, k, num_bins_x, num_bins_y);

    } // end get_bin_gid


    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \struct get_bounds_STL
    ///
    /// \brief Computes the axis-aligned bounding box of an STL triangle mesh
    ///        using parallel reductions.  Results are padded by 1e-6 on each
    ///        side so the bin mesh strictly contains all geometry.
    ///
    /// \param x_min is the minimum x coordinate, the start location of bin mesh
    /// \param y_min is the minimum y coordinate, the start location of bin mesh
    /// \param z_min is the minimum z coordinate, the start location of bin mesh
    /// \param x_max is the maximum x coordinate, the end of the bin mesh
    /// \param y_max is the maximum y coordinate, the end of the bin mesh
    /// \param z_max is the maximum z coordinate, the end of the bin mesh
    /// \param v0X is the x coordinate of the first vertex on the triangle 
    /// \param v0Y is the y coordinate of the first vertex on the triangle 
    /// \param v0Z is the z coordinate of the first vertex on the triangle 
    /// \param v1X is the x coordinate of the second vertex on the triangle 
    /// \param v1Y is the y coordinate of the second vertex on the triangle 
    /// \param v1Z is the z coordinate of the second vertex on the triangle 
    /// \param v2X is the x coordinate of the third vertex on the triangle 
    /// \param v2Y is the y coordinate of the third vertex on the triangle 
    /// \param v2Z is the z coordinate of the third vertex on the triangle 
    /////////////////////////////////////////////////////////////////////////////
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
                        const DViewCArrayKokkos <double> &v2Z) {

        
        const size_t num_inp_triangles = v0X.dims(0);

        // find (xmin, ymin, zmin) for building bin mesh
        Kokkos::parallel_reduce(
            "stl_min_domain_extents",
            num_inp_triangles,
            KOKKOS_LAMBDA(const int tri,         
                        double& xmin_lcl,
                        double& ymin_lcl,
                        double& zmin_lcl) 
            {
                xmin_lcl = fmin(v2X(tri),fmin(v1X(tri),fmin(v0X(tri), xmin_lcl)));
                ymin_lcl = fmin(v2Y(tri),fmin(v1Y(tri),fmin(v0Y(tri), ymin_lcl)));
                zmin_lcl = fmin(v2Z(tri),fmin(v1Z(tri),fmin(v0Z(tri), zmin_lcl)));
            },
            Kokkos::Min<double>(xmin), 
            Kokkos::Min<double>(ymin), 
            Kokkos::Min<double>(zmin)); 
            // end parallel reduction over all STL triangles

        xmin -= 1.e-6; // decrease by a small fraction
        ymin -= 1.e-6; // decrease by a small fraction
        zmin -= 1.e-6; // decrease by a small fraction

        // find (xmax, ymax, zmax) for building bin mesh
        Kokkos::parallel_reduce(
            "stl_max_domain_extents",
            num_inp_triangles,
            // this is the for loop coding
            KOKKOS_LAMBDA(const int tri,         
                        double& xmax_lcl,
                        double& ymax_lcl,
                        double& zmax_lcl) 
            {
                xmax_lcl = fmax(v2X(tri), fmax(v1X(tri), fmax(v0X(tri), xmax_lcl)));
                ymax_lcl = fmax(v2Y(tri), fmax(v1Y(tri), fmax(v0Y(tri), ymax_lcl)));
                zmax_lcl = fmax(v2Z(tri), fmax(v1Z(tri), fmax(v0Z(tri), zmax_lcl)));
            },
            Kokkos::Max<double>(xmax), 
            Kokkos::Max<double>(ymax), 
            Kokkos::Max<double>(zmax)); 
            // end parallel reduction over all STL triangles

        xmax += 1.e-6; // increase by a small fraction
        ymax += 1.e-6; // increase by a small fraction
        zmax += 1.e-6; // increase by a small fraction

        return;

    } // end get_bounds_STL


    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \struct build_multi_node_connectivity
    ///
    /// \brief Assigns a unique node gid to each set of geometrically coincident
    ///        points (distance^2 < 1e-25).  Used to build the shared-node
    ///        connectivity of a multi-element mesh from its corner points.
    ///
    /// \param corner_point_positions (x,y,z) coordinates of every corner point 
    /// \param node_in_corner_point   output: node gid for each corner point
    /// \param num_nodes              output: total number of unique nodes
    /////////////////////////////////////////////////////////////////////////////
    void build_multi_node_connectivity(const CArrayKokkos <double>  &corner_point_positions,
                                       CArrayKokkos <size_t> &node_in_corner_point,
                                       size_t &num_nodes) {

        const size_t num_points = corner_point_positions.dims(0);
        
        // -----------------------------------------------------------------------
        // Allocate per-point and per-bin arrays
        // -----------------------------------------------------------------------
        points_bin_gid         = DCArrayKokkos <size_t> (num_points, "points_in_gid");
        points_bin_lid_storage = CArrayKokkos  <size_t> (num_points, "bin_lid_storage");  
        
        num_points_in_bin.set_values(0); 

        // -----------------------------------------------------------------------
        // Pass 1a: assign every point to a bin and record its storage slot
        // -----------------------------------------------------------------------
        FOR_ALL_CLASS(point_gid, 0, num_points, {

            size_t bin_gid = get_bin_gid(corner_point_positions(point_gid,0), 
                                         corner_point_positions(point_gid,1), 
                                         corner_point_positions(point_gid,2));

            size_t storage_lid = Kokkos::atomic_fetch_add(&num_points_in_bin(bin_gid), 1);
            points_bin_gid(point_gid)         = bin_gid; 
            points_bin_lid_storage(point_gid) = storage_lid; 

        }); // end for all
        Kokkos::fence();
        points_bin_gid.update_host();
        num_points_in_bin.update_host();


        // -----------------------------------------------------------------------
        // Pass 1b: build the ragged bin->points map
        // -----------------------------------------------------------------------
        points_in_bin = DRaggedRightArrayKokkos <size_t> (num_points_in_bin, "num_points_in_bin");

        FOR_ALL_CLASS(point_gid, 0, num_points, {

            size_t bin_gid     = points_bin_gid(point_gid);
            size_t storage_lid = points_bin_lid_storage(point_gid);
            points_in_bin(bin_gid, storage_lid) = point_gid;

        }); // end for all
        Kokkos::fence();

        // -----------------------------------------------------------------------
        // step 2: get smallest point index in the point-point set
        // -----------------------------------------------------------------------

        node_in_corner_point = CArrayKokkos <size_t> (num_points);   // saves the node_id to point
        CArrayKokkos <size_t> smallest_point_id_in_set (num_points); // saves the smallest point id in the point-point set

        FOR_ALL_CLASS(point_gid, 0, num_points, {
           
            bin_keys_t bk  = keys_in_bin(points_bin_gid(point_gid));

            // loop over (i-1:i+1), (j-1:j+1), (k-1:k+1) if node point falls on a bin mesh edge
            const int imin = MAX(0, bk.i-1); const int imax = MIN(num_bins_x-1, bk.i+1);
            const int jmin = MAX(0, bk.j-1); const int jmax = MIN(num_bins_y-1, bk.j+1);
            const int kmin = MAX(0, bk.k-1); const int kmax = MIN(num_bins_z-1, bk.k+1);

            size_t smallest_id = point_gid;

            for (int i = imin; i <= imax; i++)
            for (int j = jmin; j <= jmax; j++)
            for (int k = kmin; k <= kmax; k++) {
                        
                size_t a_bin_gid = get_id_of_ijk(i, j, k, num_bins_x, num_bins_y);

                for(size_t storage_lid=0; storage_lid<num_points_in_bin(a_bin_gid); storage_lid++){
                    
                    // get the nieghboring point gid 
                    size_t nbr_gid = points_in_bin(a_bin_gid, storage_lid);

                    if (nbr_gid >= smallest_id) continue; // only look at lower id's

                    const double pnt_dx = (corner_point_positions(point_gid,0)-corner_point_positions(nbr_gid,0));
                    const double pnt_dy = (corner_point_positions(point_gid,1)-corner_point_positions(nbr_gid,1));
                    const double pnt_dz = (corner_point_positions(point_gid,2)-corner_point_positions(nbr_gid,2));

                    if(pnt_dx*pnt_dx + pnt_dy*pnt_dy + pnt_dz*pnt_dz < 1.e-25) smallest_id = nbr_gid;

                } // end for over stored points in bin

            } // end for 

            smallest_point_id_in_set(point_gid) = smallest_id; // write once, no atomics needed

        }); // end parallel loop over points
        Kokkos::fence();

        // -----------------------------------------------------------------------
        // step 3: count the number of unique point-point sets using smallest 
        //         point flag, assign the unique node gid to each set representive
        // -----------------------------------------------------------------------

        DCArrayKokkos <size_t> node_gid_counter(1);
        RUN({
            node_gid_counter(0) = 0; // initializing on device to be = 0
        });

        // tally the unique point-point sets using the smallest_point_id as a flag
        FOR_ALL_CLASS(point_gid, 0, num_points, {

            if (smallest_point_id_in_set(point_gid) == point_gid){
                // this point_id the smallest in the point-point set, so it is a unique node
                const size_t node_gid = Kokkos::atomic_fetch_add(&node_gid_counter(0), 1); 

                node_in_corner_point(point_gid) = node_gid;
            } // end if

        }); // end for all
        Kokkos::fence();
        node_gid_counter.update_host();

        num_nodes = node_gid_counter.host(0);
        
        // -----------------------------------------------------------------------
        // step 4: propagate the node_gid to all other points in the same set.
        //         Safe because step 2 fully completed (fence above) before this 
        //         kernel reads node_in_corner_point(smallest_point_gid).
        // -----------------------------------------------------------------------
        FOR_ALL_CLASS(point_gid, 0, num_points, {

            // if not the smallest, then save the node_gid saved to the smallest point
            if (smallest_point_id_in_set(point_gid) != point_gid){

                const size_t smallest_point_gid = smallest_point_id_in_set(point_gid);
                node_in_corner_point(point_gid) = node_in_corner_point(smallest_point_gid);
            } // end if

        }); // end for all
        Kokkos::fence();

    } // end build_multi_node_connectivity



    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \struct build_point_cloud_connectivity
    ///
    /// \brief Calculates connectivity structures between points in a point cloud 
    ///
    /// \param point_positions are the (x,y,z) coordinates of every point 
    /// \param points_in_point gives the neighboring point indicies
    /// \param points_num_neighbors is the number of neighboring points
    /// \param reverse_neighbor_lid is the local index to access the point
    /////////////////////////////////////////////////////////////////////////////
    void build_point_cloud_connectivity(const DCArrayKokkos <double> &point_positions){

        const size_t num_points = point_positions.dims(0);

        // -----------------------------------------------------------------------
        // Allocate per-point and per-bin arrays
        // -----------------------------------------------------------------------
        points_bin_gid         = DCArrayKokkos <size_t> (num_points, "points_bin_gid");
        points_bin_lid_storage = CArrayKokkos  <size_t> (num_points, "bin_lid_storage");
        points_bin_stencil     = DCArrayKokkos <int>    (num_points, 6, "bin_stencil");
        points_num_neighbors   = DCArrayKokkos <size_t> (num_points, "num_neighbors");

        num_points_in_bin.set_values(0);

        // -----------------------------------------------------------------------
        // Pass 1a: assign every point to a bin and record its storage slot
        // -----------------------------------------------------------------------
        FOR_ALL_CLASS(point_gid, 0, num_points, {

            size_t bin_gid = get_bin_gid(point_positions(point_gid, 0),
                                         point_positions(point_gid, 1),
                                         point_positions(point_gid, 2));

            size_t storage_lid = Kokkos::atomic_fetch_add(&num_points_in_bin(bin_gid), 1);
            points_bin_gid(point_gid)         = bin_gid;
            points_bin_lid_storage(point_gid) = storage_lid;

        }); // end for all
        Kokkos::fence();
        points_bin_gid.update_host();
        num_points_in_bin.update_host();

        // -----------------------------------------------------------------------
        // Pass 1b: build the ragged bin->points map
        // -----------------------------------------------------------------------
        points_in_bin = DRaggedRightArrayKokkos <size_t> (num_points_in_bin, "points_in_bin");

        FOR_ALL_CLASS(point_gid, 0, num_points, {

            size_t bin_gid     = points_bin_gid(point_gid);
            size_t storage_lid = points_bin_lid_storage(point_gid);
            points_in_bin(bin_gid, storage_lid) = point_gid;

        }); // end for all
        Kokkos::fence();

        // -----------------------------------------------------------------------
        // Pass 2: determine each point's stencil and initial neighbor count.
        //
        //   Grow a cubic shell of bins outward from point_gid's bin until:
        //     (a) the stencil contains at least min_num_points_fit points, AND
        //     (b) the stencil radius meets or exceeds search_radius,
        //   OR until all points are enclosed.
        //
        //   Each thread writes only to its own row of points_bin_stencil and its
        //   own entry of points_num_neighbors, so no atomics are needed here.
        //   The initial count is (points in stencil - 1) to exclude self.
        // -----------------------------------------------------------------------
        FOR_ALL_CLASS(point_gid, 0, num_points, {

            size_t bin_gid = points_bin_gid(point_gid);
            bin_keys_t bk  = keys_in_bin(bin_gid);

            const int i = bk.i;
            const int j = bk.j;
            const int k = bk.k;

            for (int stencil = 1; stencil < 100000; stencil++){

                const int iminus = MAX(0,                  i - stencil);
                const int iplus  = MIN((int)num_bins_x-1,  i + stencil);
                const int jminus = MAX(0,                  j - stencil);
                const int jplus  = MIN((int)num_bins_y-1,  j + stencil);
                const int kminus = MAX(0,                  k - stencil);
                const int kplus  = MIN((int)num_bins_z-1,  k + stencil);

                // count all points in the candidate stencil
                size_t num_found = 0;
                for (int kc = kminus; kc <= kplus; kc++)
                for (int jc = jminus; jc <= jplus; jc++)
                for (int ic = iminus; ic <= iplus; ic++){
                    num_found += num_points_in_bin(get_id_of_ijk(ic, jc, kc, num_bins_x, num_bins_y));
                }

                // skip distance check until minimum point count is reached
                if (num_found < min_num_points_fit && num_found < num_points) continue;

                // distance from the bin center to the near and far faces of the stencil
                const double x_mid   = bin_dx*(double)i + xmin;
                const double y_mid   = bin_dy*(double)j + ymin;
                const double z_mid   = bin_dz*(double)k + zmin;

                const double dist_minus = sqrt(
                    (bin_dx*(double)iminus + xmin - x_mid)*(bin_dx*(double)iminus + xmin - x_mid) +
                    (bin_dy*(double)jminus + ymin - y_mid)*(bin_dy*(double)jminus + ymin - y_mid) +
                    (bin_dz*(double)kminus + zmin - z_mid)*(bin_dz*(double)kminus + zmin - z_mid) );

                const double dist_plus = sqrt(
                    (bin_dx*(double)iplus  + xmin - x_mid)*(bin_dx*(double)iplus  + xmin - x_mid) +
                    (bin_dy*(double)jplus  + ymin - y_mid)*(bin_dy*(double)jplus  + ymin - y_mid) +
                    (bin_dz*(double)kplus  + zmin - z_mid)*(bin_dz*(double)kplus  + zmin - z_mid) );

                // accept stencil when it spans the search radius, or all points fit
                if (dist_minus >= search_radius || dist_plus >= search_radius ||
                    num_found == num_points){

                    points_bin_stencil(point_gid, 0) = iminus;
                    points_bin_stencil(point_gid, 1) = iplus;
                    points_bin_stencil(point_gid, 2) = jminus;
                    points_bin_stencil(point_gid, 3) = jplus;
                    points_bin_stencil(point_gid, 4) = kminus;
                    points_bin_stencil(point_gid, 5) = kplus;

                    points_num_neighbors(point_gid) = num_found - 1; // exclude self
                    break;
                }

            } // end for stencil

        }); // end for all
        Kokkos::fence();
        points_bin_stencil.update_host();

        // -----------------------------------------------------------------------
        // Pass 3 (count only): fix asymmetric stencil neighbor counts.
        //
        //   When point i's stencil covers point j's bin but j's stencil does NOT
        //   cover i's bin, j would never naturally count i as a neighbor during
        //   its own stencil walk.  Thread i detects this and atomically increments
        //   j's neighbor count to reserve the extra slot.
        //
        //   Correctness of atomic_increment here:
        //     - Each thread only increments counts for neighbors that MISSED it.
        //     - Only one thread (thread i) can discover that j missed i, because
        //       the miss condition (!nbr_sees_me) is evaluated from i's perspective.
        //     - Therefore no two threads ever atomically increment the same
        //       (point_gid, nbr_gid) pair for the same reason, eliminating
        //       double-counting.
        //     - atomic_increment ensures the counter update is indivisible,
        //       preventing lost updates when multiple threads increment different
        //       missed-neighbor slots of the same point concurrently.
        //
        //   No neighbor gids are written here — counts only.
        // -----------------------------------------------------------------------
        FOR_ALL_CLASS(point_gid, 0, num_points, {

            size_t bin_gid = points_bin_gid(point_gid);
            bin_keys_t bk  = keys_in_bin(bin_gid);

            const int i = bk.i;
            const int j = bk.j;
            const int k = bk.k;

            const int iminus = points_bin_stencil(point_gid, 0);
            const int iplus  = points_bin_stencil(point_gid, 1);
            const int jminus = points_bin_stencil(point_gid, 2);
            const int jplus  = points_bin_stencil(point_gid, 3);
            const int kminus = points_bin_stencil(point_gid, 4);
            const int kplus  = points_bin_stencil(point_gid, 5);

            for (int kc = kminus; kc <= kplus; kc++)
            for (int jc = jminus; jc <= jplus; jc++)
            for (int ic = iminus; ic <= iplus; ic++){

                size_t nbr_bin = get_id_of_ijk(ic, jc, kc, num_bins_x, num_bins_y);

                for (size_t storage_lid = 0; storage_lid < num_points_in_bin(nbr_bin); storage_lid++){

                    size_t nbr_gid = points_in_bin(nbr_bin, storage_lid);
                    if (nbr_gid == point_gid) continue; // skip self

                    // does nbr's stencil cover point_gid's bin (i,j,k)?
                    const bool nbr_sees_me =
                        (i >= points_bin_stencil(nbr_gid, 0) && i <= points_bin_stencil(nbr_gid, 1)) &&
                        (j >= points_bin_stencil(nbr_gid, 2) && j <= points_bin_stencil(nbr_gid, 3)) &&
                        (k >= points_bin_stencil(nbr_gid, 4) && k <= points_bin_stencil(nbr_gid, 5));

                    if (!nbr_sees_me){
                        // nbr's own stencil walk won't find point_gid; reserve one slot now.
                        // Safe: only thread point_gid can decide nbr missed point_gid,
                        // so this increment happens exactly once per asymmetric pair.
                        Kokkos::atomic_increment(&points_num_neighbors(nbr_gid));
                    }

                } // end for storage_lid

            } // end stencil triple-loop

        }); // end for all
        Kokkos::fence();
        points_num_neighbors.update_host();

        // -----------------------------------------------------------------------
        // Allocate points_in_point using the now-exact neighbor counts, then
        // reset points_num_neighbors to zero to reuse it as a write-index counter.
        // -----------------------------------------------------------------------
        points_in_point = DRaggedRightArrayKokkos <size_t> (points_num_neighbors, "points_in_point");
        points_num_neighbors.set_values(0); // repurposed: tracks next free slot per point

        // -----------------------------------------------------------------------
        // Pass 4 (save): populate points_in_point with neighbor gids.
        //
        //   Each thread i walks its own stencil and saves every neighbor j it
        //   finds.  When j's stencil would not see i (the asymmetric case), thread
        //   i also writes itself into j's list.  This mirrors Pass 3 exactly, so
        //   the slot reserved there is always consumed here — no overflow.
        //
        //   Thread safety:
        //     - Writes to points_in_point(point_gid, ...) are serialised by the
        //       atomic fetch-add on points_num_neighbors(point_gid): each thread
        //       owns that counter exclusively via its loop variable point_gid.
        //     - Cross-writes to points_in_point(nbr_gid, ...) are serialised by
        //       the atomic fetch-add on points_num_neighbors(nbr_gid): multiple
        //       threads may cross-write different gids into the same nbr list, but
        //       each claims a unique slot atomically, so no two writes collide.
        //     - No slot is claimed twice: the cross-write fires only when
        //       !nbr_sees_me, and only thread point_gid can satisfy that condition
        //       for the (point_gid, nbr_gid) pair, so exactly one thread writes
        //       point_gid into nbr_gid's list.
        // -----------------------------------------------------------------------
        FOR_ALL_CLASS(point_gid, 0, num_points, {

            size_t bin_gid = points_bin_gid(point_gid);
            bin_keys_t bk  = keys_in_bin(bin_gid);

            const int i = bk.i;
            const int j = bk.j;
            const int k = bk.k;

            const int iminus = points_bin_stencil(point_gid, 0);
            const int iplus  = points_bin_stencil(point_gid, 1);
            const int jminus = points_bin_stencil(point_gid, 2);
            const int jplus  = points_bin_stencil(point_gid, 3);
            const int kminus = points_bin_stencil(point_gid, 4);
            const int kplus  = points_bin_stencil(point_gid, 5);

            for (int kc = kminus; kc <= kplus; kc++)
            for (int jc = jminus; jc <= jplus; jc++)
            for (int ic = iminus; ic <= iplus; ic++){

                size_t nbr_bin = get_id_of_ijk(ic, jc, kc, num_bins_x, num_bins_y);

                for (size_t storage_lid = 0; storage_lid < num_points_in_bin(nbr_bin); storage_lid++){

                    size_t nbr_gid = points_in_bin(nbr_bin, storage_lid);
                    if (nbr_gid == point_gid) continue; // skip self

                    // save nbr into point_gid's list (stencil always covers nbr's bin)
                    size_t idx = Kokkos::atomic_fetch_add(&points_num_neighbors(point_gid), 1);
                    points_in_point(point_gid, idx) = nbr_gid;

                    // if nbr's stencil does not cover point_gid's bin, also save
                    // point_gid into nbr's list to enforce symmetry
                    const bool nbr_sees_me =
                        (i >= points_bin_stencil(nbr_gid, 0) && i <= points_bin_stencil(nbr_gid, 1)) &&
                        (j >= points_bin_stencil(nbr_gid, 2) && j <= points_bin_stencil(nbr_gid, 3)) &&
                        (k >= points_bin_stencil(nbr_gid, 4) && k <= points_bin_stencil(nbr_gid, 5));

                    if (!nbr_sees_me){
                        size_t nbr_idx = Kokkos::atomic_fetch_add(&points_num_neighbors(nbr_gid), 1);
                        points_in_point(nbr_gid, nbr_idx) = point_gid;
                    }

                } // end for storage_lid

            } // end stencil triple-loop

        }); // end for all
        Kokkos::fence();
        points_in_point.update_host();
        points_num_neighbors.update_host(); // now holds the final neighbor counts

        // -----------------------------------------------------------------------
        // Pass 5: build the reverse-neighbor map.
        //
        //   reverse_neighbor_lid(i, lid) = j_lid  such that
        //   points_in_point(nbr, j_lid) == point_gid,
        //   where nbr = points_in_point(point_gid, lid).
        //
        //   Each thread only reads and writes its own row of reverse_neighbor_lid,
        //   so no atomics are needed.  The inner search is O(degree) per entry,
        //   acceptable because point degrees are small.
        // -----------------------------------------------------------------------
        reverse_neighbor_lid = DRaggedRightArrayKokkos <size_t> (points_num_neighbors, "reverse_neighbor_lid");

        FOR_ALL_CLASS(point_gid, 0, num_points, {

            for (size_t lid = 0; lid < points_num_neighbors(point_gid); lid++){

                size_t nbr_gid = points_in_point(point_gid, lid);

                for (size_t j_lid = 0; j_lid < points_num_neighbors(nbr_gid); j_lid++){

                    if (points_in_point(nbr_gid, j_lid) == point_gid){
                        reverse_neighbor_lid(point_gid, lid) = j_lid;
                        break;
                    }

                } // end for j_lid

                // uncomment to debug missing pairs:
                // bool found = false;
                // for (...) { if (...) { found = true; break; } }
                // if (!found) printf("reverse map missing: point %zu -> nbr %zu\n", point_gid, nbr_gid);

            } // end for lid

        }); // end for all
        Kokkos::fence();
        reverse_neighbor_lid.update_host();

    } // end function build point connectivity


}; // Hash_t

} // end swage namespace

#endif