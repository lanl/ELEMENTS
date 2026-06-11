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
#ifndef STATE_H
#define STATE_H

#include "matar.h"
// #include "mpi_type.h"

using namespace mtr;

using REAL_T = double;

// Possible node states, used to initialize node_t
enum class node_state
{
    coords,
    coords_reference,
    coords_predicted,
    velocity,
    inv_mass,
    force_external
};


/////////////////////////////////////////////////////////////////////////////
///
/// \struct node_t
///
/// \brief Stores state information associated with a node
///
/////////////////////////////////////////////////////////////////////////////
struct Node_t
{
    // Replace with MPIDCArrayKokkos
    MPICArrayKokkos<REAL_T> coords;             ///< Nodal coordinates in the current configuration
    MPICArrayKokkos<REAL_T> coords_reference;   ///< Nodal coordinates in the reference configuration
    MPICArrayKokkos<REAL_T> coords_predicted;   ///< Predicted nodal coordinates
    MPICArrayKokkos<REAL_T> velocity;           ///< Nodal velocity
    MPICArrayKokkos<REAL_T> inv_mass;           ///< inv_mass on a node
    MPICArrayKokkos<REAL_T> force_external;     ///< Nodal force external
    
    // initialization method (num_nodes, num_dims, state to allocate)
    void initialize(size_t num_nodes, size_t num_dims, std::vector<node_state> node_states)
    {

        CommunicationPlan comm_plan;
        
        for (auto field : node_states){
            switch(field){
                case node_state::coords:
                    if (coords.size() == 0){
                        this->coords = MPICArrayKokkos<REAL_T>(num_nodes, num_dims, "node_coordinates");
                        this->coords.initialize_comm_plan(comm_plan);
                    }
                    break;
                case node_state::coords_reference:
                    if (coords_reference.size() == 0){
                        this->coords_reference = MPICArrayKokkos<REAL_T>(num_nodes, num_dims, "node_coordinates_reference");
                        this->coords_reference.initialize_comm_plan(comm_plan);
                    }
                    break;
                case node_state::coords_predicted:
                    if (coords_predicted.size() == 0){
                        this->coords_predicted = MPICArrayKokkos<REAL_T>(num_nodes, num_dims, "node_coordinates_predicted");
                        this->coords_predicted.initialize_comm_plan(comm_plan);
                    }
                    break;
                case node_state::velocity:
                    if (velocity.size() == 0) this->velocity = MPICArrayKokkos<REAL_T>(num_nodes, num_dims, "node_velocity");
                    this->velocity.initialize_comm_plan(comm_plan);
                    break;
                case node_state::inv_mass:
                    if (inv_mass.size() == 0) this->inv_mass = MPICArrayKokkos<REAL_T>(num_nodes, "node_inv_mass");
                    this->inv_mass.initialize_comm_plan(comm_plan);
                    break;
                case node_state::force_external:
                    if (force_external.size() == 0) this->force_external = MPICArrayKokkos<REAL_T>(num_nodes, num_dims, "node_force_external");
                    this->force_external.initialize_comm_plan(comm_plan);
                    break;
                default:
                    std::cout<<"Desired node state not understood in node_t initialize"<<std::endl;
                    throw std::runtime_error("**** Error in State Field Name ****");
            }
        }
    }; // end method
    
    // initialization method (num_nodes, num_dims, state to allocate)
    void initialize(size_t num_nodes, size_t num_dims, std::vector<node_state> node_states, CommunicationPlan& comm_plan)
    {
        for (auto field : node_states){
            switch(field){
                case node_state::coords:
                    if (coords.size() == 0){
                        this->coords = MPICArrayKokkos<REAL_T>(num_nodes, num_dims, "node_coordinates");
                        this->coords.initialize_comm_plan(comm_plan);
                    }
                    if (coords_reference.size() == 0){
                        this->coords_reference = MPICArrayKokkos<REAL_T>(num_nodes, num_dims, "node_coordinates_reference");
                        this->coords_reference.initialize_comm_plan(comm_plan);
                    }
                    if (coords_predicted.size() == 0){
                        this->coords_predicted = MPICArrayKokkos<REAL_T>(num_nodes, num_dims, "node_coordinates_predicted");
                        this->coords_predicted.initialize_comm_plan(comm_plan);
                    }
                    break;
                case node_state::velocity:
                    if (velocity.size() == 0) this->velocity = MPICArrayKokkos<REAL_T>(num_nodes, num_dims, "node_velocity");
                    this->velocity.initialize_comm_plan(comm_plan);
                    break;
                case node_state::inv_mass:
                    if (inv_mass.size() == 0) this->inv_mass = MPICArrayKokkos<REAL_T>(num_nodes, "node_inv_mass");
                    this->inv_mass.initialize_comm_plan(comm_plan);
                    break;
                case node_state::force_external:
                    if (force_external.size() == 0) this->force_external = MPICArrayKokkos<REAL_T>(num_nodes, num_dims, "node_force_external");
                    this->force_external.initialize_comm_plan(comm_plan);
                    break;
                default:
                    std::cout<<"Desired node state not understood in node_t initialize"<<std::endl;
                    throw std::runtime_error("**** Error in State Field Name ****");
            }
        }
    }; // end method

}; // end node_t


// Possible gauss point states, used to initialize GaussPoint_t
enum class gauss_pt_state
{
    J_inv_initial,
    jacobian_determinant,
    volume,
};

/////////////////////////////////////////////////////////////////////////////
///
/// \struct GaussPoint_t
///
/// \brief Stores state information associated with the Gauss point
///
/////////////////////////////////////////////////////////////////////////////
struct GaussPoint_t
{

    MPICArrayKokkos<REAL_T> J_inv_initial;
    MPICArrayKokkos<REAL_T> jacobian_determinant;
    MPICArrayKokkos<REAL_T> volume;

    // initialization method (num_cells, num_dims)
    void initialize(size_t num_gauss_pnts, size_t num_dims, std::vector<gauss_pt_state> gauss_pt_states, CommunicationPlan& comm_plan)
    {

        for (auto field : gauss_pt_states){
            switch(field){
                case gauss_pt_state::J_inv_initial:
                    //if (fields.size() == 0) this->fields = DCArrayKokkos<REAL_T>(num_gauss_pnts, "gauss_point_fields");
                    if (J_inv_initial.size() == 0){
                        this->J_inv_initial = MPICArrayKokkos<REAL_T>(num_gauss_pnts, num_dims, num_dims, "gauss_point_Jinv_initial");
                        this->J_inv_initial.initialize_comm_plan(comm_plan);
                    }
                    break;
                case gauss_pt_state::jacobian_determinant:
                    if (jacobian_determinant.size() == 0){
                        this->jacobian_determinant = MPICArrayKokkos<REAL_T>(num_gauss_pnts, "gauss_point_jacobian_determinant");
                        this->jacobian_determinant.initialize_comm_plan(comm_plan);
                    }
                    break;
                case gauss_pt_state::volume:
                    if (volume.size() == 0){
                        this->volume = MPICArrayKokkos<REAL_T>(num_gauss_pnts, "gauss_point_volume");
                        this->volume.initialize_comm_plan(comm_plan);
                    }
                    break;
                default:
                    std::cout<<"Desired gauss point state not understood in GaussPoint_t initialize"<<std::endl;
                    throw std::runtime_error("**** Error in State Field Name ****");
            }
        }
    }; // end method
};  // end GaussPoint_t

// Possible Element states, used to initialize Element_t
enum class element_state
{
    lambda, // lagrange multiplier

};

struct Element_t
{
    MPICArrayKokkos<REAL_T> lambda; // Lagrange multiplier

    void initialize(size_t num_elements, std::vector<element_state> element_states, CommunicationPlan& comm_plan)
    {
        for (auto field : element_states){
            switch(field){
                case element_state::lambda:
                    if (lambda.size() == 0) this->lambda = MPICArrayKokkos<REAL_T>(num_elements, "element_lambda");
                    this->lambda.initialize_comm_plan(comm_plan);
                    break;
            }
        }
    }; // end method
}; // end Element_t 

/////////////////////////////////////////////////////////////////////////////
///
/// \struct state_t
///
/// \brief Stores all state
///
/////////////////////////////////////////////////////////////////////////////
struct State_t
{
    // ---------------------------------------------------------------------
    //    state data on mesh declarations
    // ---------------------------------------------------------------------
    Node_t node;              ///< access as node.coords(node_gid,dim)
    GaussPoint_t GaussPoints; ///< access as GaussPoints.vol(gauss_pt_gid)
    Element_t Element; ///< access as Element.lambda(element_gid)
    
}; // end state_t





#endif