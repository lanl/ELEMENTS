/**********************************************************************************************
� 2020. Triad National Security, LLC. All rights reserved.
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
#ifndef INDEXING_UTILS_H
#define INDEXING_UTILS_H

#include "matar.h"

using namespace mtr;

/*
    ==========================
    Nodal indexing convention
    ==========================

    3D:

                  K
                  ^         J
                  |        /
                  |       /
                  |      /
          6------------------7
         /|                 /|
        / |                / |
       /  |               /  |
      /   |              /   |
     /    |             /    |
    4------------------5     |
    |     |            |     | ----> I
    |     |            |     |
    |     |            |     |
    |     |            |     |
    |     2------------|-----3
    |    /             |    /
    |   /              |   /
    |  /               |  /
    | /                | /
    |/                 |/
    0------------------1

    nodes are ordered for outward normal
    patch 0: [0,4,6,2]  xi-minus dir
    patch 1: [1,3,7,5]  xi-plus  dir
    patch 2: [0,1,5,4]  eta-minus dir
    patch 3: [3,2,6,7]  eta-plus  dir
    patch 4: [0,2,3,1]  zeta-minus dir
    patch 6: [4,5,7,6]  zeta-plus  dir


    2D:

       J
       ^
       |
     3---2
     |   |  --> I
     0---1
    
    patch 0: [0, 3]  xi-minus dir  
    patch 1: [1, 2]  xi-plus  dir  
    patch 2: [0, 1]  eta-minus dir  
    patch 3: [3, 2]  eta-plus  dir 

*/

namespace swage
{

// sort in ascending order using bubble sort
KOKKOS_INLINE_FUNCTION
void bubble_sort(size_t arr[], const size_t num)
{
    for (size_t i = 0; i < (num - 1); i++) {
        for (size_t j = 0; j < (num - i - 1); j++) {
            if (arr[j] > arr[j + 1]) {
                size_t temp = arr[j];
                arr[j]     = arr[j + 1];
                arr[j + 1] = temp;
            } // end if
        } // end for j
    } // end for i
} // end function



template <typename T>
KOKKOS_INLINE_FUNCTION
void bubble_sort(T arr)
{
    const size_t num = arr.dims(0);
    for (size_t i = 0; i < (num - 1); i++) {
        for (size_t j = 0; j < (num - i - 1); j++) {
            if (arr(j) > arr(j + 1)) {
                size_t temp = arr(j);
                arr(j)     = arr(j + 1);
                arr(j + 1) = temp;
            } // end if
        } // end for j
    } // end for i
} // end function



struct zones_in_elem_t
{
    private:
        size_t num_zones_in_elem_;
    public:
        zones_in_elem_t() {
        };

        zones_in_elem_t(const size_t num_zones_in_elem_inp) {
            this->num_zones_in_elem_ = num_zones_in_elem_inp;
        };

        // return global zone index for given local zone index in an element
        size_t host(const size_t elem_gid, const size_t zone_lid) const
        {
            return elem_gid * num_zones_in_elem_ + zone_lid;
        };

        // Return the global zone ID given an element gloabl ID and a local zone ID
        KOKKOS_INLINE_FUNCTION
        size_t operator()(const size_t elem_gid, const size_t zone_lid) const
        {
            return elem_gid * num_zones_in_elem_ + zone_lid;
        };
};

// if material points are defined strictly internal to the element.
struct gauss_in_elem_t
{
    private:
        size_t num_gauss_in_elem_;
    public:
        gauss_in_elem_t() {
        };

        gauss_in_elem_t(const size_t num_gauss_in_elem_inp) {
            this->num_gauss_in_elem_ = num_gauss_in_elem_inp;
        };

        // return global gauss index for given local gauss index in an element
        size_t  host(const size_t elem_gid, const size_t leg_gauss_lid) const
        {
            return elem_gid * num_gauss_in_elem_ + leg_gauss_lid;
        };

        // Return the global gauss ID given an element gloabl ID and a local gauss ID
        KOKKOS_INLINE_FUNCTION
        size_t operator()(const size_t elem_gid, const size_t leg_gauss_lid) const
        {
            return elem_gid * num_gauss_in_elem_ + leg_gauss_lid;
        };
};


/// 
struct corners_in_elem_t
{
    private:
        size_t num_corners_in_elem_;
    public:
        corners_in_elem_t() {
        };

        corners_in_elem_t(const size_t num_corners_in_elem_inp) {
            this->num_corners_in_elem_ = num_corners_in_elem_inp;
        };

        // return global gauss index for given local gauss index in an element
        size_t host(const size_t elem_gid, const size_t corner_lid) const
        {
            return elem_gid * num_corners_in_elem_ + corner_lid;
        };

        // Return the global gauss ID given an element gloabl ID and a local gauss ID
        KOKKOS_INLINE_FUNCTION
        size_t operator()(const size_t elem_gid, const size_t corner_lid) const
        {
            return elem_gid * num_corners_in_elem_ + corner_lid;
        };
};

/// A functor to access the patches in a surface
struct patches_in_surf_t
{
    private:
        size_t num_patches_in_surf_;  
    public:
        patches_in_surf_t() {
        };

        patches_in_surf_t(const size_t num_patches_in_surf_inp) { 
            this->num_patches_in_surf_ = num_patches_in_surf_inp;
        };

        // return global patch index for given local patch index on a surface
        size_t host(const size_t surf_gid, const size_t patch_lid) const
        {
            return surf_gid * num_patches_in_surf_ + patch_lid; 
        };

        // return global patch index for given local patch index on a surface
        KOKKOS_INLINE_FUNCTION
        size_t operator()(const size_t surf_gid, const size_t patch_lid) const
        {
            return surf_gid * num_patches_in_surf_ + patch_lid;  
        };
};

////////////////////////////////////////////////////////////////////////////////////
///
/// \fn get_surf_node_lids
///
/// \brief builds the 1D local indexing to access the surface nodes in the element
///
/// The element local indexing convention follows an i,j,k access pattern. This 
/// function leverages the i,j,k convetion and saves the 1D index to access the
/// the surface nodes from the element.  The populated array is 2D, being accessed
/// using the element face local index and then a local node index. The 2D array
/// is allocated outside the function based on mesh inputs -- number of 1D nodes
/// used to build a tensor product element and the number of dimensions. 
///
/// Important: the node ordering in this 2D array is based on the i,j,k pattern
///            of the element and not the surface.
///
/// Important: the order of the faces in the element is iminus, iplus, jminus,
///            jplus, kminus, and then kplus.  This order matches the patch 
///            node lids.  The convention used to build patches must match surfaces.
///
/// The 2D array is accessed as:
///     surf_node_ordering_in_elem(face_lid, node_lid)
///
/// \param surf_node_ordering_in_elem the array populated with local node indexing
/// \param num_1D The number of nodes (DOFs) in 1D used to build the element
/// \param num_dims The dimensions
///
/// \return void
///
///////////////////////////////////////////////////////////////////////////////////
void get_surf_node_lids(DCArrayKokkos<size_t>& surf_node_ordering_in_elem, 
                        const size_t num_1D, 
                        const size_t num_dims){
    
    if (num_dims == 3) {
        // 3D arbitrary order elements
        // num_1D = Pn+1
        // Nodes indices in elem = i + j*num_1D + k*num_1D*num_1D;
        
        // iminus-dir followed by iplus-dir
        size_t face_lid = 0;
        for (size_t i = 0; i<num_1D; i+=num_1D-1){ 
            FOR_ALL(k, 0, num_1D, 
                    j, 0, num_1D, {
                size_t node_lid = j+k*num_1D;
                surf_node_ordering_in_elem(face_lid, node_lid) = i + j * num_1D + k * num_1D * num_1D;
            });
            face_lid++;
        } // end for i-dir faces

        // jminus-dir followed by jplus-dir
        for (size_t j = 0; j<num_1D; j+=num_1D-1){
            FOR_ALL(k,0,num_1D,
                    i,0,num_1D, {
                size_t node_lid = i+k*num_1D;
                surf_node_ordering_in_elem(face_lid, node_lid) = i + j * num_1D + k * num_1D * num_1D;
            });
            face_lid++;
        } // end for j-dir faces

        // kminus-dir followed by kplus-dir
        for (size_t k = 0; k < num_1D; k+=num_1D-1) {
            FOR_ALL(j,0,num_1D, 
                    i,0,num_1D, {
                size_t node_lid = i+j*num_1D;
                surf_node_ordering_in_elem(face_lid, node_lid) = i + j * num_1D + k * num_1D * num_1D;
            });
            face_lid++;
        } // end for k-dir faces

        if(face_lid!=6)Kokkos::abort("ERROR: wrong number of element faces in 3D.\n");

    } // end 3D scope
    else if (num_dims == 2){
        // 2D arbitrary order
        // num_1D = Pn+1
        // Nodes indices in elem = i + j*num_1D + k*num_1D*num_1D;
        
        // iminus-dir followed by iplus-dir
        size_t face_lid = 0;
        for (size_t i = 0; i<num_1D; i+=num_1D-1){ 
            FOR_ALL(j,0,num_1D, {
                surf_node_ordering_in_elem(face_lid, j) = i + j * num_1D;
            });
            face_lid++;
        }

        // jminus-dir followed by jplus-dir
        for (size_t j = 0; j<num_1D; j+=num_1D-1){
            FOR_ALL(i,0,num_1D,{
                surf_node_ordering_in_elem(face_lid, i) = i + j * num_1D;
            });
            face_lid++;
        }

        if(face_lid!=4)Kokkos::abort("ERROR: wrong number of element faces in 2D.\n");

    } // end 2D scope
    else
    {
        Kokkos::abort("Bad Bad Bad: Mesh class is only supported in 2D and 3D \n");
    }
    surf_node_ordering_in_elem.update_device();

    return;

} // end function


////////////////////////////////////////////////////////////////////////////////////
///
/// \fn get_patch_node_lids
///
/// \brief builds the 1D local indexing to access the patch nodes in the element
///
/// The element local indexing convention follows an i,j,k access pattern. This 
/// function leverages the i,j,k convetion and saves the 1D index to access the
/// the patch nodes from the element following an outward normal, right-hand
/// convention.  The populated array is 3D, being accessed using the element 
/// surface local index, patch local index to the surface and then a local node 
/// index. The 3D array is allocated inside the function based on supplied
/// inputs -- number of 1D nodes used to build a tensor product element and the
/// number of dimensions. 
///
/// Important: the order of the faces in the element is iminus, iplus, jminus,
///            jplus, kminus, and then kplus.  This order matches the surface 
///            node lids. The convention used to build patches must match surfaces.
///
/// The 3D array is accessed as:
///     patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid)
///
/// \param patch_node_ordering_in_elem the array populated with local node indexing
/// \param num_1D The number of nodes (DOFs) in 1D used to build the element
/// \param num_dims The dimensions
///
/// \return void
///
///////////////////////////////////////////////////////////////////////////////////
void get_patch_node_lids(DCArrayKokkos<size_t> &patch_node_ordering_in_elem,
                         const size_t num_1D, 
                         const size_t num_dims){

    size_t i_patch = 0;
    size_t j_patch = 0;
    size_t k_patch = 0;
    size_t face_lid = 0;

    if (num_dims == 3) {
        // node ordering in patches of an arbitrary-order element

        /*
    
            i,j,k layout
    
            k  j
            | /
            |/
            o-->i
       
            i=0,imax surface
      
                         o (j+1,k+1)
                       / |
              (j,k+1) o  o (j+1,k)
                      | /
                (j,k) o
       
        */

        // iminus-dir patches 
        i_patch = 0;
        FOR_ALL(k, 0, num_1D-1, 
                j, 0, num_1D-1, {

            size_t node_lid = 0;
            size_t surf_patch_lid = j+k*(num_1D-1);

            // node_lid 0 in patch
            // index = i + j*num_1D + k*num_1D*num_1D;
            patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) = 
                        i_patch + j * num_1D + k * num_1D * num_1D; // node_rid(i_patch, j, k, num_1D);
            node_lid++;

            // node_lid 1 in patch
            // index = i + j*num_1D + (k+1)*num_1D*num_1D;
            patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) = 
                        i_patch + j * num_1D + (k + 1) * num_1D * num_1D; // node_rid(i_patch, j, k+1, num_1D);
            node_lid++;

            // node_lid 2 in patch
            // index = i + (j+1)*num_1D + (k+1)*num_1D*num_1D;
            patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) = 
                        i_patch + (j + 1) * num_1D + (k + 1) * num_1D * num_1D; // node_rid(i_patch, j+1, k+1, num_1D);
            node_lid++;

            // node_lid 3 in patch
            // index = i + (j+1)*num_1D + k*num_1D*num_1D;
            patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) = 
                        i_patch + (j + 1) * num_1D + k * num_1D * num_1D; // node_rid(i_patch, j+1, k, num_1D);

        }); // end parallel for
        face_lid ++;    

        // iplus-dir patches
        i_patch = num_1D - 1;
        FOR_ALL(k, 0, num_1D-1, 
                j, 0, num_1D-1, {

            size_t node_lid = 0;
            size_t surf_patch_lid = j+k*(num_1D-1);

            // node_lid 0 in patch
            // index = i + j*num_1D + k*num_1D*num_1D;
            patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) = 
                        i_patch + j * num_1D + k * num_1D * num_1D; // node_rid(i_patch, j, k, num_1D);
            node_lid++;

            // node_lid 1 in patch
            // index = i + (j+1)*num_1D + k*num_1D*num_1D;
            patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) = 
                        i_patch + (j + 1) * num_1D + k * num_1D * num_1D; // node_rid(i_patch, j+1, k, num_1D);
            node_lid++;

            // node_lid 2 in patch
            // index = i + (j+1)*num_1D + (k+1)*num_1D*num_1D;
            patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) = 
                        i_patch + (j + 1) * num_1D + (k + 1) * num_1D * num_1D; // node_rid(i_patch, j+1, k+1, num_1D);
            node_lid++;

            // node_lid 3 in patch
            // index = i + j*num_1D + (k+1)*num_1D*num_1D;
            patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) = 
                        i_patch + j * num_1D + (k + 1) * num_1D * num_1D; // node_rid(i_patch, j, k+1, num_1D);

        }); // end parallel for
        face_lid ++;


        /*
            i,j,k layout

            k  j
            | /
            |/
            o-->i


            j=0,jmax

             (i,k+1) o--o (i+1,k+1)
                     |  |
              (i,,k) o--o (i+1,k)

        */

        j_patch = 0;
        FOR_ALL(k, 0, num_1D - 1, 
                i, 0, num_1D - 1, {

            size_t node_lid = 0;
            size_t surf_patch_lid = i+k*(num_1D-1);

            // node_lid 0 in patch
            patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) =  
                i + j_patch * num_1D + k * num_1D * num_1D; // node_rid(i, j_patch, k, num_1D);
            node_lid++;

            // node_lid 1 in patch
            patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) =  
                i + 1 + j_patch * num_1D + k * num_1D * num_1D; // node_rid(i+1, j_patch, k, num_1D);
            node_lid++;

            // node_lid 2 in patch
            patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) =  
                i + 1 + j_patch * num_1D + (k + 1) * num_1D * num_1D; // node_rid(i+1, j_patch, k+1, num_1D);
            node_lid++;

            // node_lid 3 in patch
            patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) =  
                i + j_patch * num_1D + (k + 1) * num_1D * num_1D; // node_rid(i, j_patch, k+1, num_1D);
            
        }); // end parallel for
        face_lid ++;
                

        j_patch = num_1D - 1;
        FOR_ALL(k, 0, num_1D-1, 
                i, 0, num_1D-1, {

            size_t node_lid = 0;
            size_t surf_patch_lid = i+k*(num_1D-1);

            // node_lid 0 in patch
            patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) = 
                        i + j_patch * num_1D + k * num_1D * num_1D; // node_rid(i, j_patch, k, num_1D);
            node_lid++;

            // node_lid 1 in patch
            patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) = 
                        i + j_patch * num_1D + (k + 1) * num_1D * num_1D; // node_rid(i, j_patch, k+1, num_1D);
            node_lid++;

            // node_lid 2 in patch
            patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) = 
                        i + 1 + j_patch * num_1D + (k + 1) * num_1D * num_1D; // node_rid(i+1, j_patch, k+1, num_1D);
            node_lid++;

            // node_lid 3 in patch
            patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) = 
                        i + 1 + j_patch * num_1D + k * num_1D * num_1D; // node_rid(i+1, j_patch, k, num_1D);
            
        }); // end parallel for
        face_lid ++;


        /*
        
            i,j,k layout
        
            k  j
            | /
            |/
            o-->i
        
        
            k=0,kmax
        
            (i,j+1) o--o (i+1,j+1)
                   /  /
            (i,j) o--o (i+1,j)
         
        */

        k_patch = 0;
        FOR_ALL(j, 0, num_1D-1, 
                i, 0, num_1D-1, {
                
                size_t node_lid = 0;
                size_t surf_patch_lid = i+j*(num_1D-1);
            
                // node_lid 0 in patch
                patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) = 
                            i + j * num_1D + k_patch * num_1D * num_1D; // node_rid(i, j, k_patch, num_1D);
                node_lid++;

                // node_lid 1 in patch
                patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) = 
                            i + (j + 1) * num_1D + k_patch * num_1D * num_1D; // node_rid(i, j+1, k_patch, num_1D);
                node_lid++;

                // node_lid 2 in patch
                patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) = 
                            i + 1 + (j + 1) * num_1D + k_patch * num_1D * num_1D; // node_rid(i+1, j+1, k_patch, num_1D);
                node_lid++;

                // node_lid 3 in patch
                patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) = 
                            i + 1 + j * num_1D + k_patch * num_1D * num_1D; // node_rid(i+1, j, k_patch, num_1D);

        }); // end parallel for
        face_lid ++;

        k_patch = num_1D - 1;
        FOR_ALL(j, 0, num_1D-1, 
                i, 0, num_1D-1, {
                
                size_t node_lid = 0;
                size_t surf_patch_lid = i+j*(num_1D-1);

                // node_lid 0 in patch
                patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) = 
                            i + j * num_1D + k_patch * num_1D * num_1D; // node_rid(i, j, k_patch, num_1D);
                node_lid++;

                // node_lid 1 in patch
                patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) = 
                            i + 1 + j * num_1D + k_patch * num_1D * num_1D; // node_rid(i+1, j, k_patch, num_1D);
                node_lid++;

                // node_lid 2 in patch
                patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) = 
                            i + 1 + (j + 1) * num_1D + k_patch * num_1D * num_1D; // node_rid(i+1, j+1, k_patch, num_1D);
                node_lid++;

                // node_lid 3 in patch
                patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) = 
                            i + (j + 1) * num_1D + k_patch * num_1D * num_1D; // node_rid(i, j+1, k_patch, num_1D);

        }); // end parallel for
        face_lid ++;

        if(face_lid!=6) Kokkos::abort("ERROR: wrong number of element faces in 3D when building patches.\n");

    }// end if 3D element
    else{
        // 2D arbitrary order elements

        // iminus-dir patches
        i_patch = 0;
        FOR_ALL(j, 0, num_1D-1, {

            size_t node_lid = 0;
            size_t surf_patch_lid = j;

            patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) = 
                        i_patch + j * num_1D; // node_rid(i_patch, j, num_1D;
            node_lid++;

            patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) = 
                        i_patch + (j + 1) * num_1D; // node_rid(i_patch, j+1, num_1D;
        }); // end parallel for j
        face_lid ++; 

        // i-plus-dir patches
        i_patch = num_1D - 1;
        FOR_ALL(j, 0, num_1D-1, {

            size_t node_lid = 0;
            size_t surf_patch_lid = j;

            patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) = 
                        i_patch + j * num_1D; // node_rid(i_patch, j, num_1D;
            node_lid++;

            patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) = 
                        i_patch + (j + 1) * num_1D; // node_rid(i_patch, j+1, num_1D;
        }); // end parallel for j
        face_lid ++; 

        j_patch = 0;
        FOR_ALL(i, 0, num_1D-1, {
            
            size_t node_lid = 0;
            size_t surf_patch_lid = i;

            patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) = 
                        i + j_patch * num_1D; // node_rid(i, j_patch, num_1D);
            node_lid++;

            patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) = 
                        i + 1 + j_patch * num_1D; // node_rid(i+1, j_patch, num_1D);
        }); // end parallel for i
        face_lid ++; 

        j_patch = num_1D - 1;
        FOR_ALL(i, 0, num_1D-1, {

            size_t node_lid = 0;
            size_t surf_patch_lid = i;

            patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) = 
                        i + j_patch * num_1D; // node_rid(i, j_patch, num_1D);
            node_lid++;

            patch_node_ordering_in_elem(face_lid, surf_patch_lid, node_lid) = 
                        i + 1 + j_patch * num_1D; // node_rid(i+1, j_patch, num_1D);
        }); // end parallel for i
        face_lid ++;

        if(face_lid!=4) Kokkos::abort("ERROR: wrong number of element faces in 2D when building patches.\n");

    } // end if 2D arbitrary-order element


} // end function


} // end name space

#endif