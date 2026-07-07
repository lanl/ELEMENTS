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
#ifndef REF_ELEM_NEW_H
#define REF_ELEM_NEW_H

#include <cmath>
#include "matar.h"
#include "ref_quadrature.h"


using namespace mtr;


namespace reference_space
{
    enum ElementType
    {
        linearElement = 0,        // single quadrature point element
        arbitraryOrderElement = 1 // fully integrated arbitrary-order element
    };

    // Location of basis DOFs
    enum BasisType
    {
        LagrangeLobatto = 0,  // GLL in Steven's solver (C0 kinematic space)
        LagrangeLegendre = 1, // GL in Steven's solver (DG thermo space)
    };

    enum QuadratureType
    {
        GaussLobatto = 0,
        GaussLegendre = 1
    };
} // end ref_elem namespace


namespace elements
{

    // Quadrature rules for surfaces and elems
    struct Quadrature_t
    {
        reference_space::QuadratureType QuadratureType;

        size_t elem_dims = 0;
        size_t num_qpts_in_elem = 0;
        size_t num_qpts_1d = 0;

        CArrayKokkos<double> qpt_positions;
        CArrayKokkos<double> qpt_weights;

        /////////////////////////////////////////////////////////////////////////////
        ///
        /// \fn initialize_quadrature
        ///
        /// \brief Set up quadrature in a volume or surface element
        ///
        /// \param TypeInp The type of quadrature (e.g., Lobatto or Legendra)
        /// \param num_qpts_1d_inp The number of quadrature in 1D, applied to each direction.
        /// \param elem_dims_in The number dimensions 
        ///
        /////////////////////////////////////////////////////////////////////////////
        void initialize_quadrature(const reference_space::QuadratureType TypeInp, 
                                   const size_t num_qpts_1d_inp, 
                                   const size_t elem_dims_in)
        {
            QuadratureType = TypeInp;

            elem_dims = elem_dims_in;
            num_qpts_1d = num_qpts_1d_inp;
            if(num_qpts_1d==0) throw std::runtime_error("ERROR: zero quadrature points specified \n");

            num_qpts_in_elem = 1;
            for(size_t dim=0; dim<elem_dims; dim++){
                num_qpts_in_elem *= num_qpts_1d;
            }
            qpt_positions = CArrayKokkos<double>(num_qpts_in_elem, elem_dims, "qpt_positions");
            qpt_weights   = CArrayKokkos<double>(num_qpts_in_elem, "qpt_weights");

            // temporary 1D variables to build 3D element
            CArrayKokkos<double> qpt_positions_1d(num_qpts_1d, "qpt_positions_1d");
            CArrayKokkos<double> qpt_weights_1d (num_qpts_1d, "qpt_weights_1d");

            if(QuadratureType == reference_space::GaussLegendre){
                RUN_CLASS({
                    get_legendre_nodes_1D(qpt_positions_1d, num_qpts_1d);
                    get_legendre_weights_1D(qpt_weights_1d, num_qpts_1d);
                });
            }
            else if(QuadratureType == reference_space::GaussLobatto){         
                RUN_CLASS({
                    get_lobatto_nodes_1D(qpt_positions_1d, num_qpts_1d);
                    get_lobatto_weights_1D(qpt_weights_1d, num_qpts_1d);
                });
            }
            else
            {
                throw std::runtime_error("ERROR: unsupported quadrature set specified \n");
            }

            // 3D volume element
            if(elem_dims==3){
                FOR_ALL_CLASS(k, 0, num_qpts_1d,
                              j, 0, num_qpts_1d,
                              i, 0, num_qpts_1d, {

                    const size_t rid = get_qpt_rid(i, j, k);

                    qpt_positions(rid, 0) = qpt_positions_1d(i);
                    qpt_positions(rid, 1) = qpt_positions_1d(j);
                    qpt_positions(rid, 2) = qpt_positions_1d(k);

                    qpt_weights(rid) = qpt_weights_1d(i) * qpt_weights_1d(j) * qpt_weights_1d(k);
                });
                Kokkos::fence();
            }
            // 2D volume or 2D surface element
            else if (elem_dims==2){
                FOR_ALL_CLASS(j, 0, num_qpts_1d,
                              i, 0, num_qpts_1d, {

                    const size_t rid = get_qpt_rid(i, j);

                    qpt_positions(rid, 0) = qpt_positions_1d(i);
                    qpt_positions(rid, 1) = qpt_positions_1d(j);

                    qpt_weights(rid) = qpt_weights_1d(i) * qpt_weights_1d(j);
                });
                Kokkos::fence();
            }
            // 1D volume, edge, or 1D surface element
            else if (elem_dims==1) {
                FOR_ALL_CLASS(i, 0, num_qpts_1d, {

                    const size_t rid = i;

                    qpt_positions(rid, 0) = qpt_positions_1d(i);

                    qpt_weights(rid) = qpt_weights_1d(i);
                });
                Kokkos::fence();
            }
            else{
            throw std::runtime_error("ERROR: unsupported quadrature elem dims \n");
            }

        } // init fcn quadrature


        /////////////////////////////////////////////////////////////////////////////
        ///
        /// \fn get_qpt_rid
        ///
        /// \brief Compute the 1D array index for a quadrature point in a 3D volume 
        ///        element.
        ///
        /// Calculates the row-major flat index corresponding to a quadrature point
        /// at the position (i, j, k) within an element. This function is typically 
        /// used for accessing basis functions, positions, and weights defined on 
        /// the tensor-product grid of quadrature points in the reference element, 
        /// which is common in high-order finite element and spectral methods.
        ///
        /// \param i Local quadrature index in the first (xi) coordinate direction.
        /// \param j Local quadrature index in the second (eta) coordinate direction.
        /// \param k Local quadrature index in the third (mu) coordinate direction.
        ///
        /// \return The row-major offset index for the quadrature point in this element.
        ///
        /////////////////////////////////////////////////////////////////////////////
        KOKKOS_INLINE_FUNCTION
        size_t get_qpt_rid(size_t i, size_t j, size_t k) const
        {
            return i + (j + k * num_qpts_1d) * num_qpts_1d;
        } // end function

        /////////////////////////////////////////////////////////////////////////////
        ///
        /// \fn get_qpt_rid
        ///
        /// \brief Compute the 1D array index for a quadrature point in an 2D volume 
        ///        or surface element.
        ///
        /// Calculates the row-major flat index corresponding to a quadrature point
        /// at the position (i, j) within a 2D volume or surface element. This 
        /// function is typically used for accessing basis functions, positions, 
        /// and weights defined on the tensor-product grid of quadrature points in 
        /// the reference element, which is common in high-order finite element and 
        /// spectral methods.
        ///
        /// \param i Local quadrature index in the first (xi) coordinate direction.
        /// \param j Local quadrature index in the second (eta) coordinate direction.
        ///
        /// \return The row-major offset index for the quadrature point in this element.
        ///
        /////////////////////////////////////////////////////////////////////////////
        KOKKOS_INLINE_FUNCTION
        size_t get_qpt_rid(size_t i, size_t j) const
        {
            return i + j * num_qpts_1d;
        } // end function

    }; // end Quadrature_t


    // reference element data structure
    struct ReferenceElement_t
    {

        reference_space::ElementType ElementType = reference_space::linearElement; ///< The type of element
        reference_space::BasisType BasisType = reference_space::LagrangeLobatto;  ///<Basis DOF location
        
        size_t elem_dims = 0;

        // Dofs
        size_t num_dofs_in_elem = 1;
        size_t num_dofs_1d = 0;

        size_t Pn = 0;
      

        // DOF positions
        CArrayKokkos<double> dof_positions;


        // Basis evaluation at quadrature points
        CArrayKokkos<double> qpt_basis;
        CArrayKokkos<double> qpt_grad_basis;


        /////////////////////////////////////////////////////////////////////////////
        ///
        /// \fn initialize_ref_elem
        ///
        /// \brief Initialize the reference element with polynomial order and dimension.
        ///        It is the companion structure for the unstructured_mesh.h
        ///
        /// Initializes internal data members for the reference element according
        /// to the given polynomial order and spatial dimension. Sets up the number of
        /// degrees of freedom (DOF), as well as associated sizes and counts for basis
        /// function evaluations, according to the dimension and polynomial order. If
        /// the polynomial order is less than element order, it is in a discontinous 
        /// space.  For that case, user must create multiple reference elements, one 
        /// that defines the position and one for the discontinous fields.
        ///
        /// The DOF's for the Lagrange basis can be at Lobatto or Legendra points, 
        /// those points are in most cases are spatially different from the quadrature
        /// points. 
        ///
        /// \param ElemTypeInp The element type (classical linear or arbitrary-order)
        /// \param BasisTypeInp The location of basis DOFs 
        /// \param Quadrature_t Quadrature set in the element
        /// \param p_order The element order e.g., constant = 0, linear = 1, quadratic = 2, ...
        ///
        /// \return void
        ///
        /////////////////////////////////////////////////////////////////////////////
        void initialize_ref_elem(const reference_space::ElementType ElemTypeInp,
                                 const reference_space::BasisType BasisTypeInp,
                                 const struct Quadrature_t Quadrature,
                                 const size_t p_order)
        {
            // set element and basis type
            ElementType = ElemTypeInp;
            BasisType = BasisTypeInp;

            elem_dims = Quadrature.elem_dims; 
            if(elem_dims==0) throw std::runtime_error("ERROR: quadrature not correctly specified \n");
            if(elem_dims>3) throw std::runtime_error("ERROR: only 1D, 2D, and 3D reference elements supported \n");

            // -----------------------------------------------------------------------
            // Step 1a: determine the number of DOFs in 3D
            // -----------------------------------------------------------------------
            Pn = p_order;
            num_dofs_1d = p_order + 1;

            num_dofs_in_elem = 1;  // Initialize to 1
            for (size_t dim = 0; dim < elem_dims; dim++) {
                num_dofs_in_elem *= num_dofs_1d;
            } // end for

            // -----------------------------------------------------------------------
            // Step 1b: get the positions in reference space for the DOFs
            // -----------------------------------------------------------------------
            dof_positions = CArrayKokkos<double>(num_dofs_in_elem, elem_dims, "dof_positions");
            CArrayKokkos<double> dof_positions_1d(num_dofs_1d, "dof_positions_1d");

            // dof positions can be at legendre or lobatto locations in elem
            if(BasisTypeInp == reference_space::LagrangeLegendre){
                RUN_CLASS({
                    get_legendre_nodes_1D(dof_positions_1d, num_dofs_1d);
                });
            }
            else if(BasisTypeInp == reference_space::LagrangeLobatto){         
                RUN_CLASS({
                    get_lobatto_nodes_1D(dof_positions_1d, num_dofs_1d);
                });
            }
            else
            {
                throw std::runtime_error("ERROR: unsupported basis DOF locations specified \n");
            }

            // 3D volume element
            if(elem_dims==3){
                FOR_ALL_CLASS(k, 0, num_dofs_1d,
                              j, 0, num_dofs_1d,
                              i, 0, num_dofs_1d, {

                        const size_t rid = get_dof_rid(i, j, k);

                        dof_positions(rid, 0) = dof_positions_1d(i);
                        dof_positions(rid, 1) = dof_positions_1d(j);
                        dof_positions(rid, 2) = dof_positions_1d(k);
                });
            } // end if 3D
            // 2D volume or 2D surface element
            else if (elem_dims==2){ 
                FOR_ALL_CLASS(j, 0, num_dofs_1d,
                              i, 0, num_dofs_1d, {

                    const size_t rid = get_dof_rid(i, j);

                    dof_positions(rid, 0) = dof_positions_1d(i);
                    dof_positions(rid, 1) = dof_positions_1d(j);
                });
            } // end if 1D
            // 1D volume, edge, or 1D surface element
            else {
                FOR_ALL_CLASS(i, 0, num_dofs_1d, {

                    const size_t rid = i;

                    dof_positions(rid, 0) = dof_positions_1d(i);
                });
            } // end if 1D
            Kokkos::fence();
            
            // -----------------------------------------------------------------------
            // Step 2: Calculate the basis values at quadrature points
            // -----------------------------------------------------------------------
            qpt_basis = CArrayKokkos<double>(Quadrature.num_qpts_in_elem, num_dofs_in_elem, "qpt_basis");           

            // temporary arrays to hold evaluations at a single point for each dof
            CArrayKokkos<double> temp_basis(num_dofs_in_elem);
            CArrayKokkos<double> temp_val_1d(num_dofs_1d);
            CArrayKokkos<double> temp_val_Nd(num_dofs_1d, elem_dims); // 2D or 3D

            CArrayKokkos<double> point(elem_dims);


            RUN_CLASS({
                for (size_t qpt_rid = 0; qpt_rid < Quadrature.num_qpts_in_elem; qpt_rid++) {
                    
                    // Get the evaluation coordinates
                    for (size_t dim = 0; dim < elem_dims; dim++) {
                        point(dim) = Quadrature.qpt_positions(qpt_rid, dim);
                    }

                    get_basis(temp_basis, dof_positions_1d, temp_val_1d, temp_val_Nd, point);

                    for (size_t basis_id = 0; basis_id < num_dofs_in_elem; basis_id++) {
                        qpt_basis(qpt_rid, basis_id) = temp_basis(basis_id);
                        temp_basis(basis_id) = 0.0;
                    }
                } // end for over qpts in elem
            });
            Kokkos::fence();

            // -----------------------------------------------------------------------
            // Step 3: Calculate the grad basis values at quadrature points
            // -----------------------------------------------------------------------
            qpt_grad_basis = CArrayKokkos<double>(Quadrature.num_qpts_in_elem, num_dofs_in_elem, elem_dims, "qpt_grad_basis");  
                

            // temporary arrays to hold evaluations at a single point for each dof
            CArrayKokkos<double> temp_partial_xi(num_dofs_in_elem);
            CArrayKokkos<double> temp_partial_eta(num_dofs_in_elem);
            CArrayKokkos<double> temp_partial_mu(num_dofs_in_elem);

            CArrayKokkos<double> temp_Dval_1d(num_dofs_1d);
            CArrayKokkos<double> temp_Dval_Nd(num_dofs_1d, elem_dims);


            RUN_CLASS({
                for (size_t qpt_rid = 0; qpt_rid < Quadrature.num_qpts_in_elem; qpt_rid++) {
                    
                    // Get the evaluation coordinates
                    for (size_t dim = 0; dim < elem_dims; dim++) {
                        point(dim) = Quadrature.qpt_positions(qpt_rid, dim);
                    }

                    partial_xi_basis(temp_partial_xi, 
                                     dof_positions_1d, 
                                     temp_val_1d, 
                                     temp_val_Nd, 
                                     temp_Dval_1d, 
                                     temp_Dval_Nd, 
                                     point);

                    if(elem_dims>1) partial_eta_basis(temp_partial_eta, 
                                                      dof_positions_1d, 
                                                      temp_val_1d, 
                                                      temp_val_Nd, 
                                                      temp_Dval_1d, 
                                                      temp_Dval_Nd, 
                                                      point);

                    if(elem_dims>2) partial_mu_basis(temp_partial_mu, 
                                                     dof_positions_1d, 
                                                     temp_val_1d, 
                                                     temp_val_Nd, 
                                                     temp_Dval_1d, 
                                                     temp_Dval_Nd, 
                                                     point);

                    for (size_t basis_id = 0; basis_id < num_dofs_in_elem; basis_id++) {
                        qpt_grad_basis(qpt_rid, basis_id, 0) = temp_partial_xi(basis_id);
                        if(elem_dims>1)qpt_grad_basis(qpt_rid, basis_id, 1) = temp_partial_eta(basis_id);
                        if(elem_dims>2)qpt_grad_basis(qpt_rid, basis_id, 2) = temp_partial_mu(basis_id);

                        temp_partial_xi(basis_id)  = 0.0;
                        if(elem_dims>1) temp_partial_eta(basis_id) = 0.0;
                        if(elem_dims>2) temp_partial_mu(basis_id)  = 0.0;
                    } // end loop over basis functions

                    
                } // end for qpts in elem
            });
            Kokkos::fence();

        } // end of member function

        
        /////////////////////////////////////////////////////////////////////////////
        ///
        /// \fn get_dof_rid
        ///
        /// \brief Compute the 1D array index for a degree of freedom (DOF) in an element.
        ///
        /// Calculates the flat row-major index corresponding to a DOF located at position (i, j, k)
        /// in the element, for continuous fields. This is used for basis functions and data fields
        /// that are continuous across element boundaries.
        ///
        /// \param i Local DOF index in the first (xi) coordinate direction.
        /// \param j Local DOF index in the second (eta) coordinate direction.
        /// \param k Local DOF index in the third (mu) coordinate direction.
        ///
        /// \return The row-major offset index for the DOF in this element.
        ///
        /////////////////////////////////////////////////////////////////////////////
        KOKKOS_INLINE_FUNCTION
        size_t get_dof_rid(size_t i, size_t j, size_t k) const
        {
            return i + (j + k * num_dofs_1d) * num_dofs_1d;
        }
        
        /////////////////////////////////////////////////////////////////////////////
        ///
        /// \fn get_dof_rid
        ///
        /// \brief Compute the 1D array index for a degree of freedom (DOF) in an element.
        ///
        /// Calculates the flat row-major index corresponding to a DOF located at position (i, j, k)
        /// in the element, for continuous fields. This is used for basis functions and data fields
        /// that are continuous across element boundaries.
        ///
        /// \param i Local DOF index in the first (xi) coordinate direction.
        /// \param j Local DOF index in the second (eta) coordinate direction.
        ///
        /// \return The row-major offset index for the DOF in this element.
        ///
        /////////////////////////////////////////////////////////////////////////////
        KOKKOS_INLINE_FUNCTION
        size_t get_dof_rid(size_t i, size_t j) const
        {
            return i + j * num_dofs_1d;
        }


        /////////////////////////////////////////////////////////////////////////////
        ///
        /// \fn get_basis
        ///
        /// \brief Computes the tensor-product nodal basis values at an arbitrary point.
        ///
        /// This function evaluates the Lagrange basis functions at a specified point within
        /// the reference element and assembles the tensor-product basis values for all degrees 
        /// of freedom (DOFs). Basis values in each coordinate direction are computed independently 
        /// using the 1D Lagrange basis, and then combined to form the full multi-dimensional basis.
        /// The results are written to the provided output array.
        ///
        /// \param basis Reference to the output CArrayKokkos<double> to hold full tensor-product basis values, sized for all DOFs in the element.
        /// \param val_1d Temporary CArrayKokkos<double> for holding 1D basis values (as workspace).
        /// \param val_3d Temporary CArrayKokkos<double> for holding basis values for each direction; shape should be (num_dofs_1d, 3).
        /// \param point Reference to CArrayKokkos<double> representing the coordinates (xi, eta, mu) at which the basis is evaluated (size 3).
        ///
        /// \return void
        ///
        /////////////////////////////////////////////////////////////////////////////
        KOKKOS_FUNCTION
        void get_basis(const CArrayKokkos<double>& basis,
            const CArrayKokkos<double>& dof_positions_1d,
            const CArrayKokkos<double>& val_1d,
            const CArrayKokkos<double>& val_Nd,
            const CArrayKokkos<double>& point) const
        {
    
            for(size_t dim=0; dim<elem_dims; dim++){

                for (size_t i = 0; i < num_dofs_1d; i++) {
                    val_1d(i) = 0.0;
                }

                // Calculate 1D basis for the this dim coordinate of the point
                lagrange_basis_1D(val_1d, dof_positions_1d, point(dim));

                // Save the basis value at the point to a temp array and zero out the temp array
                for (size_t i = 0; i < num_dofs_1d; i++) {
                    val_Nd(i, dim) = val_1d(i);
                }

            } // end loop over dims


            if(elem_dims==3){
                // Multiply the i, j, k components of the basis from each node
                // to get the tensor product basis for the node
                for (size_t k = 0; k < num_dofs_1d; k++) 
                for (size_t j = 0; j < num_dofs_1d; j++)
                for (size_t i = 0; i < num_dofs_1d; i++) {
                    const size_t rid = get_dof_rid(i, j, k);
                    basis(rid) = val_Nd(i, 0) * val_Nd(j, 1) * val_Nd(k, 2);
  
                }
            } // end if 3D
            else if(elem_dims==2){
                // Multiply the i, j components of the basis from each node
                // to get the tensor product basis for the node
                for (size_t j = 0; j < num_dofs_1d; j++) 
                for (size_t i = 0; i < num_dofs_1d; i++) {
                    const size_t rid = get_dof_rid(i, j);
                    basis(rid) = val_Nd(i, 0) * val_Nd(j, 1);  
                }
            } // end if 2D 
            else{
                for (size_t i = 0; i < num_dofs_1d; i++) {
                    const size_t rid = i;
                    basis(rid) = val_Nd(i, 0);  
                }
            } // end if 1D


            // reset values to 0.0
            for (size_t i = 0; i < num_dofs_1d; i++) {
                val_1d(i)    = 0.0;
                for(size_t dim=0; dim<elem_dims; dim++){
                    val_Nd(i, dim) = 0.0;
                }
            } // end for i

        } // end get basis function


        /////////////////////////////////////////////////////////////////////////////
        ///
        /// \fn partial_xi_basis
        ///
        /// \brief Compute the partial derivative of the basis function with respect to the xi coordinate.
        ///
        /// This function evaluates the tensor-product Lagrange basis function derivatives in the xi direction
        /// at a given point within the reference element. The computation is performed by first evaluating the
        /// 1D derivatives and basis functions for each coordinate and then combining them using the chain rule
        /// for tensor products. The result is stored in the provided array for all degrees of freedom.
        ///
        /// \param partial_xi Array to store the value of the partial derivative with respect to xi for each basis function.
        /// \param val_1d Temporary workspace array for 1D basis evaluations.
        /// \param val_Nd Temporary workspace array for 1D/2D/3D basis component evaluations.
        /// \param Dval_1d Temporary workspace array for 1D derivative evaluations.
        /// \param Dval_Nd Temporary workspace array for 1D/2D/3D derivative component evaluations.
        /// \param point Input array specifying the coordinates (xi, eta, mu) at which to evaluate the derivative.
        ///
        /// \return void
        ///
        /////////////////////////////////////////////////////////////////////////////
        KOKKOS_FUNCTION
        void partial_xi_basis(const CArrayKokkos<double>& partial_xi,
                              const CArrayKokkos<double>& dof_positions_1d,
                              const CArrayKokkos<double>& val_1d,
                              const CArrayKokkos<double>& val_Nd,
                              const CArrayKokkos<double>& Dval_1d,
                              const CArrayKokkos<double>& Dval_Nd,
                              const CArrayKokkos<double>& point) const
        {
            // get grad basis 
            for (size_t i = 0; i < num_dofs_1d; i++) {
                Dval_1d(i) = 0.0;
            }

            // Calculate 1D partial w.r.t. xi for the X coordinate of the point
            lagrange_derivative_1D(Dval_1d, dof_positions_1d, point(0));

            // Save the grad basis value at the point to a temp array and zero out the temp array
            for (size_t i = 0; i < num_dofs_1d; i++) {
                Dval_Nd(i, 0) = Dval_1d(i);
            }

            // get Y and Z basis, the latter only if elem_dims = 3
            for(size_t dim=1; dim<elem_dims; dim++){

                for (size_t i = 0; i < num_dofs_1d; i++) {
                    val_1d(i) = 0.0;
                }

                // Calculate 1D basis for the this dim coordinate of the point
                lagrange_basis_1D(val_1d, dof_positions_1d, point(dim));

                // Save the basis value at the point to a temp array and zero out the temp array
                for (size_t i = 0; i < num_dofs_1d; i++) {
                    val_Nd(i, dim) = val_1d(i);
                }

            } // end loop over dims

            if(elem_dims==3){
                // Multiply the i, j, k components of the basis and partial_xi from each node
                // to get the tensor product partial derivatives of the basis at each node
                for (size_t k = 0; k < num_dofs_1d; k++) 
                for (size_t j = 0; j < num_dofs_1d; j++) 
                for (size_t i = 0; i < num_dofs_1d; i++) {
                    const size_t rid = get_dof_rid(i, j, k);

                    // Partial w.r.t xi
                    partial_xi(rid) = Dval_Nd(i, 0) * val_Nd(j, 1) * val_Nd(k, 2);
                } // end for

                for (size_t i = 0; i < num_dofs_1d; i++) {
                    val_1d(i)     = 0.0;
                    val_Nd(i, 0)  = 0.0;
                    val_Nd(i, 1)  = 0.0;
                    val_Nd(i, 2)  = 0.0;
                    Dval_1d(i)    = 0.0;
                    Dval_Nd(i, 0) = 0.0;
                    Dval_Nd(i, 1) = 0.0;
                    Dval_Nd(i, 2) = 0.0;
                } // end for
            } //end if 3D
            else if(elem_dims==2){
                // Multiply the i, j components of the basis and partial_xi from each node
                // to get the tensor product partial derivatives of the basis at each node 
                for (size_t j = 0; j < num_dofs_1d; j++) 
                for (size_t i = 0; i < num_dofs_1d; i++) {
                    const size_t rid = get_dof_rid(i, j);

                    // Partial w.r.t xi
                    partial_xi(rid) = Dval_Nd(i, 0) * val_Nd(j, 1);
                } // end for

                for (size_t i = 0; i < num_dofs_1d; i++) {
                    val_1d(i)     = 0.0;
                    val_Nd(i, 0)  = 0.0;
                    val_Nd(i, 1)  = 0.0;
                    Dval_1d(i)    = 0.0;
                    Dval_Nd(i, 0) = 0.0;
                    Dval_Nd(i, 1) = 0.0;
                } // end for
            } // end if 2D
            else {
                // Multiply the i components of the basis and partial_xi from each node
                // to get the tensor product partial derivatives of the basis at each node  
                for (size_t i = 0; i < num_dofs_1d; i++) {
                    const size_t rid = i;

                    // Partial w.r.t xi
                    partial_xi(rid) = Dval_Nd(i, 0);
                } // end for

                for (size_t i = 0; i < num_dofs_1d; i++) {
                    val_1d(i)     = 0.0;
                    val_Nd(i, 0)  = 0.0;
                    Dval_1d(i)    = 0.0;
                    Dval_Nd(i, 0) = 0.0;
                } // end for
            } // end if 1D

        } // end partial_xi_basis

        /////////////////////////////////////////////////////////////////////////////
        ///
        /// \fn partial_eta_basis
        ///
        /// \brief Computes the partial derivative of the Lagrange basis functions with respect to eta at a given point.
        ///
        /// This function evaluates the partial derivative of the tensor-product Lagrange basis 
        /// functions with respect to the eta (second) coordinate at a given point in the reference element.
        /// It builds up the 3D basis and its derivatives using repeated calls to 1D basis and derivative 
        /// evaluations, and computes the tensor-product forms to populate the provided partial_eta array 
        /// with the values of the partials at every reference node.
        ///
        /// \param partial_eta   Output array to store the partial derivatives with respect to eta at each basis node.
        /// \param val_1d       Temporary array used for storing 1D basis values during calculation.
        /// \param val_3d       Temporary 2D array used for storing intermediate 3D basis values.
        /// \param Dval_1d      Temporary array used for storing 1D derivative values during calculation.
        /// \param Dval_3d      Temporary 2D array used for storing intermediate 3D derivative values.
        /// \param point        Input array representing the coordinates (xi, eta, zeta) at which to evaluate the derivative.
        ///
        /// \return void
        ///
        /////////////////////////////////////////////////////////////////////////////
        KOKKOS_FUNCTION
        void partial_eta_basis(const CArrayKokkos<double>& partial_eta,
                               const CArrayKokkos<double>& dof_positions_1d,
                               const CArrayKokkos<double>& val_1d,
                               const CArrayKokkos<double>& val_Nd,
                               const CArrayKokkos<double>& Dval_1d,
                               const CArrayKokkos<double>& Dval_Nd,
                               const CArrayKokkos<double>& point) const
        {


            // get X and Z basis values, the latter only if elem_dims = 3D
            for(size_t dim=0; dim<elem_dims; dim+=2){

                for (size_t i = 0; i < num_dofs_1d; i++) {
                    val_1d(i) = 0.0;
                }

                // Calculate 1D basis for the this dim coordinate of the point
                lagrange_basis_1D(val_1d, dof_positions_1d, point(dim));

                // Save the basis value at the point to a temp array and zero out the temp array
                for (size_t i = 0; i < num_dofs_1d; i++) {
                    val_Nd(i, dim) = val_1d(i);
                }

            } // end loop over dims

            // get grad basis 
            for (size_t i = 0; i < num_dofs_1d; i++) {
                Dval_1d(i) = 0.0;
            }

            // Calculate 1D partial w.r.t. eta for the Y coordinate of the point
            lagrange_derivative_1D(Dval_1d, dof_positions_1d, point(1));

            // Save the grad basis value at the point to a temp array and zero out the temp array
            for (size_t i = 0; i < num_dofs_1d; i++) {
                Dval_Nd(i, 1) = Dval_1d(i);
            }


            if(elem_dims==3){
                // Multiply the i, j, k components of the basis and partial_eta from each node
                // to get the tensor product partial derivatives of the basis at each node
                for (size_t k = 0; k < num_dofs_1d; k++) 
                for (size_t j = 0; j < num_dofs_1d; j++) 
                for (size_t i = 0; i < num_dofs_1d; i++) {
                    size_t rid = get_dof_rid(i, j, k);

                    // Partial w.r.t xi
                    partial_eta(rid) = val_Nd(i, 0) * Dval_Nd(j, 1) * val_Nd(k, 2);
                } // end for

                for (size_t i = 0; i < num_dofs_1d; i++) {
                    val_1d(i)     = 0.0;
                    val_Nd(i, 0)  = 0.0;
                    val_Nd(i, 1)  = 0.0;
                    val_Nd(i, 2)  = 0.0;
                    Dval_1d(i)    = 0.0;
                    Dval_Nd(i, 0) = 0.0;
                    Dval_Nd(i, 1) = 0.0;
                    Dval_Nd(i, 2) = 0.0;
                }
            } // end if 3D
            else if(elem_dims==2){
                // Multiply the i, j components of the basis and partial_eta from each node
                // to get the tensor product partial derivatives of the basis at each node 
                for (size_t j = 0; j < num_dofs_1d; j++) 
                for (size_t i = 0; i < num_dofs_1d; i++) {
                    size_t rid = get_dof_rid(i, j);

                    // Partial w.r.t xi
                    partial_eta(rid) = val_Nd(i, 0) * Dval_Nd(j, 1);
                } // end for

                for (size_t i = 0; i < num_dofs_1d; i++) {
                    val_1d(i)     = 0.0;
                    val_Nd(i, 0)  = 0.0;
                    val_Nd(i, 1)  = 0.0;
                    Dval_1d(i)    = 0.0;
                    Dval_Nd(i, 0) = 0.0;
                    Dval_Nd(i, 1) = 0.0;
                }
            } // end if 1D
            else {
                // Note: there is no 1D with this function, its valid for only 2D & 3D
            } // end if 1D

        } // end partial_eta_basis

        /////////////////////////////////////////////////////////////////////////////
        ///
        /// \fn partial_mu_basis
        ///
        /// \brief Compute the tensor-product partial derivative of the basis functions
        ///        with respect to the mu (z) coordinate at a given point in the reference element.
        ///
        /// This function calculates the partial derivatives of the high-order
        /// tensor-product Lagrange basis functions with respect to the mu (z) coordinate
        /// at a specified point in the reference element, using their 1D basis and derivative
        /// representations. The result is stored in the provided array. Temporary arrays
        /// for 1D and 3D basis/derivative values are passed in for workspace re-use.
        ///
        /// \param partial_mu Array to store the computed partial derivatives with respect to mu.
        /// \param val_1d Workspace array for 1D Lagrange basis evaluations.
        /// \param val_3d Workspace array for 3D tensor-product basis evaluations.
        /// \param Dval_1d Workspace array for 1D Lagrange basis derivatives.
        /// \param Dval_3d Workspace array for 3D tensor-product basis derivatives.
        /// \param point Coordinates in the reference element where the basis partials are evaluated.
        ///
        /// \return void
        ///
        /////////////////////////////////////////////////////////////////////////////
        KOKKOS_FUNCTION
        void partial_mu_basis(const CArrayKokkos<double>& partial_mu,
                              const CArrayKokkos<double>& dof_positions_1d,
                              const CArrayKokkos<double>& val_1d,
                              const CArrayKokkos<double>& val_3d,
                              const CArrayKokkos<double>& Dval_1d,
                              const CArrayKokkos<double>& Dval_3d,
                              const CArrayKokkos<double>& point) const
        {
            // this routine is only valid for 3D ref elems

            // get X and Y basis
            for(size_t dim=0; dim<elem_dims-1; dim++){

                for (size_t i = 0; i < num_dofs_1d; i++) {
                    val_1d(i) = 0.0;
                }

                // Calculate 1D basis for the this dim coordinate of the point
                lagrange_basis_1D(val_1d, dof_positions_1d, point(dim));

                // Save the basis value at the point to a temp array and zero out the temp array
                for (size_t i = 0; i < num_dofs_1d; i++) {
                    val_3d(i, dim) = val_1d(i);
                }

            } // end loop over dims


            // get grad basis
            for (size_t i = 0; i < num_dofs_1d; i++) {
                Dval_1d(i) = 0.0;
            }
            
            // Calculate 1D partial w.r.t. mu for the Z coordinate of the point
            lagrange_derivative_1D(Dval_1d, dof_positions_1d, point(2));

            // Save the basis value at the point to a temp array and zero out the temp array
            for (size_t i = 0; i < num_dofs_1d; i++) {
                Dval_3d(i, 2) = Dval_1d(i);
            }

            // Multiply the i, j, k components of the basis and partial_xi from each node
            // to get the tensor product partial derivatives of the basis at each node
            for (size_t k = 0; k < num_dofs_1d; k++) 
            for (size_t j = 0; j < num_dofs_1d; j++) 
            for (size_t i = 0; i < num_dofs_1d; i++) {
                size_t rid = get_dof_rid(i, j, k);

                // Partial w.r.t mu
                partial_mu(rid) = val_3d(i, 0) * val_3d(j, 1) * Dval_3d(k, 2);
            } // end for

            for (size_t i = 0; i < num_dofs_1d; i++) {
                val_1d(i)     = 0.0;
                val_3d(i, 0)  = 0.0;
                val_3d(i, 1)  = 0.0;
                val_3d(i, 2)  = 0.0;
                Dval_1d(i)    = 0.0;
                Dval_3d(i, 0) = 0.0;
                Dval_3d(i, 1) = 0.0;
                Dval_3d(i, 2) = 0.0;
            }
        }
        /////////////////////////////////////////////////////////////////////////////
        ///
        /// \fn lagrange_basis_1D
        ///
        /// \brief Computes the Lagrange basis functions in 1D at a given point.
        ///
        /// This function evaluates the values of the 1D Lagrange basis functions at a specified point 
        /// within the reference element, using the nodal positions. For each basis node, it computes 
        /// the interpolation value (the product over all other node positions) and stores 
        /// the result in the provided array.
        ///
        /// \param interp   Output array to store the value of each basis function at the specified point.
        /// \param dof_positions_1d the positions of the DOFs in reference coordinates
        /// \param x_point  Point at which to evaluate the basis functions.
        ///
        /// \return void
        ///
        /////////////////////////////////////////////////////////////////////////////

        KOKKOS_FUNCTION
        void lagrange_basis_1D(
            const CArrayKokkos<double>& interp,           // interpolant from each basis
            const CArrayKokkos<double>& dof_positions_1d, // location of basis DOFs in ref elem 
            const double x_point) const                   // point of interest in element
        // calculate the basis value associated with each node_i
        {
            for (size_t vert_i = 0; vert_i < num_dofs_1d; vert_i++) {
                double numerator   = 1.0;       // placeholder numerator
                double denominator = 1.0;       // placeholder denominator
                double interpolant = 1.0;       // placeholder value of numerator/denominator

                for (size_t vert_j = 0; vert_j < num_dofs_1d; vert_j++) { // looping over the verts !=vert_i
                    if (vert_j != vert_i) {
                        // Calculate the numerator
                        numerator = numerator * (x_point - dof_positions_1d(vert_j));

                        // Calculate the denominator
                        denominator = denominator * (dof_positions_1d(vert_i) - dof_positions_1d(vert_j));
                    } // end if

                    
                } // end looping over nodes != vert_i
                interpolant = numerator / denominator; // storing a single value for interpolation for node vert_i

                // writing value to vectors for later use
                interp(vert_i) = interpolant;             // Interpolant value at given point
            } // end loop over all nodes
        } // end of Lagrange_1D function

        /////////////////////////////////////////////////////////////////////////////
        ///
        /// \fn lagrange_derivative_1D
        ///
        /// \brief Computes the values of the derivatives of the 1D Lagrange basis functions at a given point.
        ///
        /// This function evaluates the first derivatives of the 1D Lagrange basis functions associated
        /// with the element's degrees of freedom at a specified point within the reference element.
        /// For each basis node, it computes the derivative of the basis function using the nodal
        /// positions and stores the results in the provided array.
        ///
        /// \param derivative Output array to store the value of each 1D basis function derivative at the given point.
        /// \param x_point    Point at which to evaluate the derivatives of the basis functions.
        ///
        /// \return void
        ///
        /////////////////////////////////////////////////////////////////////////////
        KOKKOS_INLINE_FUNCTION
        void lagrange_derivative_1D(
            const CArrayKokkos<double>& derivative,        // derivative
            const CArrayKokkos<double>& dof_positions_1d,  // location of basis DOFs in ref elem
            const double x_point) const                    // point of interest in element
        {
            for (size_t vert_i = 0; vert_i < num_dofs_1d; vert_i++) { // looping over the nodes
                double denominator  = 1.0;      // placeholder denominator
                double num_gradient = 0.0;      // placeholder for numerator of the gradient
                double gradient     = 0.0;

                for (size_t vert_j = 0; vert_j < num_dofs_1d; vert_j++) { // looping over the nodes !=vert_i
                    if (vert_j != vert_i) {
                        // Calculate the denominator that is the same for
                        // both the basis and the gradient of the basis
                        denominator = denominator * (dof_positions_1d(vert_i) - dof_positions_1d(vert_j));

                        double product_gradient = 1.0;

                        // Calculate the numerator of the gradient
                        for (size_t N = 0; N < num_dofs_1d; N++) { // looping over the nodes !=vert_i
                            if (N != vert_j && N != vert_i) {
                                product_gradient = product_gradient * (x_point - dof_positions_1d(N));
                            } // end if
                        } // end for

                        // Sum over the product of the numerator
                        // contributions from each node
                        num_gradient += product_gradient;
                    } // end if

                    
                } // end looping over nodes != vert_i
                
                gradient = (num_gradient / denominator); // storing the derivative of the interpolating function
                
                // writing value to vectors for later use
                derivative(vert_i) = gradient;     // derivative of each function
            } // end loop over all nodes
        } // end of Lagrange_1D function

    }; // end struct

} // end namespace elements

#endif