/*****************************************************************************
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
#ifndef ELEMENTS_REF_ELEMENT_H
#define ELEMENTS_REF_ELEMENT_H 
namespace elements {

    // Reference Element Informations
    class ref_element{
    private:
        
        // this assumes i,j,k structured indexing of nodes
        const int node_rlid_in_patch_in_cell_[24] = {
            // patches in xi-dir
            0,2,6,4,
            1,3,7,5,
            // patches in eta-dir
            0,4,5,1,
            2,6,7,3,
            // patches in mu-dir
            0,1,3,2,
            4,5,7,6
        };
        
        // return the local ref id for the same patch in the cell-patch neighbor
        const int patch_rlid_cell_neighbor_[6] = {
            1, // xi left side of my cell
            0, // xi right side of my cell
            3, // eta front side of my cell
            2, // eta bact side of my cell
            5, // mu bottom side of my cell
            4  // mu top side of my cell
        };
        
        // return the reference coordinate unit normal of the sides in a cell
        const real_t cell_side_unit_normals_[18] = {
            -1, 0, 0, // xi left side of my cell
            1, 0, 0, // xi right side of my cell
            0,-1, 0, // eta minus side of my cell
            0, 1, 0, // eta plus side of my cell
            0, 0,-1, // mu minus side of my cell
            0, 0, 1  // mu plus side of my cell
        };
        
        int num_dim_;
        
        int num_ref_nodes_1D_;
        int num_ref_cells_1D_;
        int num_ref_corners_1D_;
        
        int num_ref_surface_nodes_in_elem_;
        int num_ref_inside_nodes_in_elem_;
        
        // Zones
        int num_zones_1d_;
        int num_zones_in_elem_;
        
        // cells
        int num_ref_cells_in_elem_;
        
        CArray <int>cells_in_zone_list_;
        
        // nodes
        int num_ref_nodes_in_elem_;
        int num_ref_nodes_in_cell_;
        
        CArray <int> ref_nodes_in_cell_;
        CArray <int> ref_surface_nodes_in_elem_;
        CArray <int> ref_inside_nodes_in_elem_;
        
        CArray <real_t> ref_node_positions_;
        CArray <real_t> ref_node_g_weights_;
        
        CArray <int> cell_nodes_in_elem_list_;
        
        // Vertices
        int num_ref_verts_1d_;
        int num_ref_verts_in_elem_;
        
        // corners
        int num_ref_corners_in_cell_;
        int num_ref_corners_in_elem_;
        
        CArray <int> ref_corners_in_cell_;
        
        CArray <real_t> ref_corner_surf_normals_;
        CArray <real_t> ref_corner_g_weights_;
        CArray <real_t> ref_corner_surf_g_weights_;
        
        // patches and sides
        int num_patches_in_elem_;
        int num_sides_in_elem_;
        
        CArray <real_t> ref_cell_positions_;
        CArray <real_t> ref_cell_g_weights_;
        CArray <real_t> ref_patch_positions_;
        CArray <real_t> ref_patch_g_weights_;
        
        CArray <int> ref_patches_in_cell_list_;
        
        CArray <real_t> ref_patch_basis_;
        CArray <real_t> ref_patch_gradient_;  // grad basis at patch
        
        CArray <real_t> ref_cell_basis_;
        CArray <real_t> ref_cell_gradient_;   // grad basis at cell
        
        // Num basis functions
        int num_basis_;
        
        // Basis evaluation at nodes
        CArray <real_t> ref_nodal_basis_;
        
        // Pointer to HexN object
        HexN *elem_ptr;
        
    public:
        
        // Gradient of basis
        CArray <real_t> ref_nodal_gradient_;
        //real_t * ref_nodal_gradient_;
        
        // Function Declarations
        
        // Default constructor
        ref_element() : elem_ptr(nullptr) {};
        
        // Initialize reference element information
        void init(int poly_order, int num_dim, HexN &elem);
        
        int num_dim() const;
        
        int num_basis() const;
        
        int num_ref_nodes() const;
        
        int num_ref_cells_in_elem() const;
        int num_ref_corners_in_cell() const;
        
        // local i,j,k indexing
        int node_rid(int i, int j, int k) const;
        int cell_rid(int i, int j, int k) const;
        int corner_rid(int i, int j, int k) const;
        
        int ref_corners_in_cell(int cell_rid, int corner_rlid) const;
        int ref_nodes_in_cell(int cell_rid, int node_rlid) const;
        
        int ref_surface_nodes_in_elem(int node_rlid) const;
        int ref_inside_nodes_in_elem(int node_rlid) const;
        int num_ref_surface_nodes_in_elem() const;
        int num_ref_inside_nodes_in_elem() const;
        
        real_t ref_node_positions(int node_rid, int dim) const;
        
        real_t ref_corner_surface_normals(int corner_rid, int surf_rlid, int dim) const;
        
        real_t ref_corner_g_surface_weights(int corner_rid, int surf_rlid) const;
        
        real_t ref_node_g_weights(int node_rid) const;
        
        real_t ref_corner_g_weights(int corner_rid) const;
        
        real_t &ref_nodal_gradient(int node_rid, int basis_id, int dim) const;
        
        real_t &ref_nodal_basis(int node_rid, int basis_id) const;
        
        int& cell_lid_in_zone(int zone_lid, int cell_lid) const;
        
        real_t ref_cell_positions(int cell_rid, int dim) const;
        real_t ref_cell_g_weights(int cell_rid) const;
        real_t ref_cell_basis(int cell_rid, int basis_id) const;
        real_t ref_cell_gradient(int cell_rid, int basis_id, int dim) const;
        
        real_t cell_side_unit_normals(int side_rlid, int dim) const;
        
        int& cell_nodes_in_elem(int cell_lid, int node_lid) const;
        
        int node_in_patch_in_cell(int patch_lid, int node_lid) const;
        
        int patch_rlid_in_cell_neighbor(int patch_rlid) const;
        
        int ref_patches_in_cell(int cell_lid, int patch_rlid) const;
        
        real_t ref_patch_positions(int patch_rid, int dim) const;
        real_t ref_patch_g_weights(int patch_rid) const;
        real_t ref_patch_basis(int patch_rid, int basis_id) const;
        real_t ref_patch_gradient(int patch_rid, int basis_id, int dim) const;
        
        int vert_node_map(int vert_lid);
        
        // Deconstructor
        ~ref_element();
    };

}
#endif // ELEMENTS_REF_ELEMENT_H
