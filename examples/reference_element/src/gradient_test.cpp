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
using namespace elements;

bool Verbose = false;
size_t max_num   = 19; // max number of quadrature points to test up to, the limit is 19
size_t max_order = 19; // max polynomial order to test, limit is 19th-order with Legendre

// polynomial with terms <= p_order
KOKKOS_INLINE_FUNCTION
double polynomial(const CArrayKokkos<double> &coeff, const double x, const size_t p_order){

    double result = 0.0;
    for (int i = 0; i <= p_order; ++i) {
        result += coeff(i) * pow(x, (double)i);
    }
    return result;
} // end polynomial

KOKKOS_INLINE_FUNCTION
double polynomial(const CArrayKokkos<double> &coeff, const double x, const double y, const size_t p_order){

    double result = 0.0;
    for (int j = 0; j <= p_order; ++j) 
    for (int i = 0; i <= p_order-j; ++i) {
        result += coeff(i,j) * pow(x, (double)i) * pow(y, (double)j);
        
    } // end for
    return result;
} // end polynomial

KOKKOS_INLINE_FUNCTION
double polynomial(const CArrayKokkos<double> &coeff, const double x, const double y, const double z, const size_t p_order){

    double result = 0.0;
    for (int k = 0; k <= p_order; ++k) 
    for (int j = 0; j <= p_order-k; ++j) 
    for (int i = 0; i <= p_order-j-k; ++i) {
        result += coeff(i,j,k) * pow(x, (double)i) * pow(y, (double)j) * pow(z, (double)k); 
    } // end for

    return result;
} // end polynomial


// polynomial with terms <= p_order
KOKKOS_INLINE_FUNCTION
double grad_polynomial_x(const CArrayKokkos<double> &coeff, const double x, const size_t p_order){

    double result = 0.0;
    for (int i = 0; i <= p_order; ++i) {
        result += (double)i*coeff(i) * pow(x, (double)i-1);
    }

    return result;
} // end polynomial

KOKKOS_INLINE_FUNCTION
double grad_polynomial_x(const CArrayKokkos<double> &coeff, const double x, const double y, const size_t p_order){

    double result = 0.0;
    for (int j = 0; j <= p_order; ++j) 
    for (int i = 0; i <= p_order-j; ++i) {
        result += (double)i*coeff(i,j) * pow(x, (double)i-1) * pow(y, (double)j);
    } // end for

    return result;
} // end polynomial

KOKKOS_INLINE_FUNCTION
double grad_polynomial_y(const CArrayKokkos<double> &coeff, const double x, const double y, const size_t p_order){

    double result = 0.0;
    for (int j = 0; j <= p_order; ++j) 
    for (int i = 0; i <= p_order-j; ++i) {
        result += (double)j*coeff(i,j) * pow(x, (double)i) * pow(y, (double)j-1);
    } // end for

    return result;
} // end polynomial

KOKKOS_INLINE_FUNCTION
double grad_polynomial_x(const CArrayKokkos<double> &coeff, const double x, const double y, const double z, const size_t p_order){

    double result = 0.0;
    for (int k = 0; k <= p_order; ++k) 
    for (int j = 0; j <= p_order-k; ++j) 
    for (int i = 0; i <= p_order-j-k; ++i) {
        result += (double)i*coeff(i,j,k) * pow(x, (double)i-1.) * pow(y, (double)j) * pow(z, (double)k); 
    } // end for

    return result;
} // end polynomial

KOKKOS_INLINE_FUNCTION
double grad_polynomial_y(const CArrayKokkos<double> &coeff, const double x, const double y, const double z, const size_t p_order){

    double result = 0.0;
    for (int k = 0; k <= p_order; ++k) 
    for (int j = 0; j <= p_order-k; ++j) 
    for (int i = 0; i <= p_order-j-k; ++i) {
        result += (double)j*coeff(i,j,k) * pow(x, (double)i) * pow(y, (double)j-1.) * pow(z, (double)k); 
    } // end for

    return result;
} // end polynomial

KOKKOS_INLINE_FUNCTION
double grad_polynomial_z(const CArrayKokkos<double> &coeff, const double x, const double y, const double z, const size_t p_order){

    double result = 0.0;
    for (int k = 0; k <= p_order; ++k) 
    for (int j = 0; j <= p_order-k; ++j) 
    for (int i = 0; i <= p_order-j-k; ++i) {
        result += (double)k*coeff(i,j,k) * pow(x, (double)i) * pow(y, (double)j) * pow(z, (double)k-1.); 
    } // end for

    return result;
} // end polynomial


// Test: interpolate a polynomial that the basis can represent exactly
void test_gradient(const Quadrature_t& Quad,
                   const ReferenceElement_t& RefElem) {
    
    // For polynomial order p, Lagrange basis can represent 
    // any polynomial of degree <= p exactly
    
    const size_t p_order = RefElem.num_dofs_1d - 1;

    CArrayKokkos<double> dof_values(RefElem.num_dofs_in_elem);
    CArrayKokkos<double> coeff;


    if(RefElem.elem_dims==1){
        coeff = CArrayKokkos<double>(p_order+1);
    }
    else if (RefElem.elem_dims==2){
        coeff = CArrayKokkos<double>(p_order+1, p_order+1);
    }
    else {
        coeff = CArrayKokkos<double>(p_order+1, p_order+1, p_order+1);
    }
    coeff.set_values(0.78914567);  // a radom value

    

    // populate polynomial values at DOFs of the ref element
    FOR_ALL(dof, 0, RefElem.num_dofs_in_elem, {

        double value=0.0;
        if(RefElem.elem_dims==1){
            const double xi = RefElem.dof_positions(dof, 0);
            value = polynomial(coeff, xi, p_order);
        }
        else if (RefElem.elem_dims==2){
            const double xi  = RefElem.dof_positions(dof, 0);
            const double eta = RefElem.dof_positions(dof, 1);
            value = polynomial(coeff, xi, eta, p_order);
        }
        else {
            const double xi  = RefElem.dof_positions(dof, 0);
            const double eta = RefElem.dof_positions(dof, 1);
            const double mu  = RefElem.dof_positions(dof, 2);
            value = polynomial(coeff, xi, eta, mu, p_order);
        }

        dof_values(dof) = value; // the value at the ref elem node

    }); // end parallel


    // compare interpoloation with exact solution
    FOR_ALL(qpt, 0, Quad.num_qpts_in_elem, {
        
        double sum[3]; // gradient value
        sum[0] = 0.;
        sum[1] = 0.;
        sum[2] = 0.;

        for (size_t dof = 0; dof < RefElem.num_dofs_in_elem; dof++) {
            for(size_t dim=0; dim<RefElem.elem_dims; dim++){
                sum[dim] += RefElem.qpt_grad_basis(qpt, dof, dim)*dof_values(dof);
            }
        }

        double exact_value[3];
        exact_value[0] = 0.0;
        exact_value[1] = 0.0;
        exact_value[2] = 0.0;

        if(RefElem.elem_dims==1){
            const double xi = Quad.qpt_positions(qpt, 0);
            exact_value[0] = grad_polynomial_x(coeff, xi, p_order);
        }
        else if (RefElem.elem_dims==2){
            const double xi  = Quad.qpt_positions(qpt, 0);
            const double eta = Quad.qpt_positions(qpt, 1);
            exact_value[0] = grad_polynomial_x(coeff, xi, eta, p_order);
            exact_value[1] = grad_polynomial_y(coeff, xi, eta, p_order);
        }
        else {
            const double xi  = Quad.qpt_positions(qpt, 0);
            const double eta = Quad.qpt_positions(qpt, 1);
            const double mu  = Quad.qpt_positions(qpt, 2);
            exact_value[0] = grad_polynomial_x(coeff, xi, eta, mu, p_order);
            exact_value[1] = grad_polynomial_y(coeff, xi, eta, mu, p_order);
            exact_value[2] = grad_polynomial_z(coeff, xi, eta, mu, p_order);
        }

        // round off get bad at high p_order's due to polynomial sensativity
        for(size_t dim=0; dim<RefElem.elem_dims; dim++){
            if (fabs(sum[dim] - exact_value[dim]) > 1.e-10*(double)p_order) {
                printf("Error: interpolation failed at qpt id = %d with order = %zu \n", qpt, p_order);
                printf("interpolated = %f vs exact value = %f, error = %f \n", sum[dim], exact_value[dim], sum[dim]-exact_value[dim]);
                Kokkos::abort("Interpolation failed at quadrature point ");
            }
            if(Verbose)printf("Grad in dim = %zu: interpolated = %f vs exact value = %f \n", dim, sum[dim], exact_value[dim]);
        }

    }); // end loop over qpts
} // end function

int main(int argc, char** argv) {

MATAR_INITIALIZE(argc, argv);
{ // MATAR scope

    printf("\n--- Gradient tests ---\n");


    printf("\n--- DG element with Legendre Quadrature & Legendre DOFs ---\n");
    for(size_t num_qpts_1D = 1; num_qpts_1D<=max_num; num_qpts_1D++){

        if(Verbose)printf("num quadrature points in 1D = %zu \n", num_qpts_1D);

        Quadrature_t Quad;
        
        // elem_dims=1,2,3
        for(size_t elem_dims_test = 1; elem_dims_test<=3; elem_dims_test++){   
            Quad.initialize_quadrature(reference_space::GaussLegendre,
                                       num_qpts_1D,
                                       elem_dims_test);

            // build reference elements of varing orders, 
            for (size_t p_order = 0; p_order<max_order; p_order++){
                if(Verbose)printf("p_order = %zu: \n", p_order);
                ReferenceElement_t FERefElem;

                // p_order is the basis order for Lagrange polynomial
                FERefElem.initialize_ref_elem(reference_space::arbitraryOrderElement,
                                              reference_space::LagrangeLegendre,
                                              Quad,
                                              p_order);

                if(Verbose)printf("interpolation check: \n");
                test_gradient(Quad, FERefElem);
                if(Verbose)printf("\n");

            } // end p_order loop
        } // elem

        if(Verbose)printf("\n");
    } // end loop of num qpts 


    printf("\n--- DG element with Lobatto Quadrature & Legendre DOFs ---\n");
    for(size_t num_qpts_1D = 1; num_qpts_1D<=max_num; num_qpts_1D++){

        if(Verbose)printf("num quadrature points in 1D = %zu \n", num_qpts_1D);

        Quadrature_t Quad;
        
        // elem_dims=1,2,3
        for(size_t elem_dims_test = 1; elem_dims_test<=3; elem_dims_test++){   
            Quad.initialize_quadrature(reference_space::GaussLobatto,
                                       num_qpts_1D,
                                       elem_dims_test);

            // build reference elements of varing orders
            for (size_t p_order = 0; p_order<max_order; p_order++){
                if(Verbose)printf("p_order = %zu: \n", p_order);
                ReferenceElement_t FERefElem;

                // p_order is the basis order for Lagrange polynomial
                FERefElem.initialize_ref_elem(reference_space::arbitraryOrderElement,
                                              reference_space::LagrangeLegendre,
                                              Quad,
                                              p_order);

                if(Verbose)printf("interpolation check: \n");
                test_gradient(Quad, FERefElem);
                if(Verbose)printf("\n");

            } // end p_order loop
        } // elem

        if(Verbose)printf("\n");
    } // end loop of num qpts 



    printf("\n--- FE element with Lobatto Quadrature & Lobatto DOFs ---\n");
    for(size_t num_qpts_1D = 1; num_qpts_1D<=max_num; num_qpts_1D++){

        if(Verbose)printf("num quadrature points in 1D = %zu \n", num_qpts_1D);

        Quadrature_t Quad;
        
        // elem_dims=1,2,3
        for(size_t elem_dims_test = 1; elem_dims_test<=3; elem_dims_test++){   
            Quad.initialize_quadrature(reference_space::GaussLobatto,
                                       num_qpts_1D,
                                       elem_dims_test);

            // build reference elements of varing orders, starting at 1, its Lobatto points
            for (size_t p_order = 1; p_order<max_order; p_order++){
                if(Verbose)printf("p_order = %zu: \n", p_order);
                ReferenceElement_t FERefElem;

                // p_order is the basis order for Lagrange polynomial
                FERefElem.initialize_ref_elem(reference_space::arbitraryOrderElement,
                                              reference_space::LagrangeLobatto,
                                              Quad,
                                              p_order);

                if(Verbose)printf("interpolation check: \n");
                test_gradient(Quad, FERefElem);
                if(Verbose)printf("\n");

            } // end p_order loop
        } // elem

        if(Verbose)printf("\n");
    } // end loop of num qpts 




    printf("\n--- FE element with Legendre Quadrature & Lobatto DOFs ---\n");
   for(size_t num_qpts_1D = 1; num_qpts_1D<=max_num; num_qpts_1D++){

        if(Verbose)printf("num quadrature points in 1D = %zu \n", num_qpts_1D);

        Quadrature_t Quad;
        
        // elem_dims=1,2,3
        for(size_t elem_dims_test = 1; elem_dims_test<=3; elem_dims_test++){   
            Quad.initialize_quadrature(reference_space::GaussLegendre,
                                       num_qpts_1D,
                                       elem_dims_test);

            // build reference elements of varing orders, starting at 1, its Lobatto points
            for (size_t p_order = 1; p_order<max_order; p_order++){
                if(Verbose)printf("p_order = %zu: \n", p_order);
                ReferenceElement_t FERefElem;

                // p_order is the basis order for Lagrange polynomial
                FERefElem.initialize_ref_elem(reference_space::arbitraryOrderElement,
                                              reference_space::LagrangeLobatto,
                                              Quad,
                                              p_order);

                if(Verbose)printf("interpolation check: \n");
                test_gradient(Quad, FERefElem);
                if(Verbose)printf("\n");

            } // end p_order loop
        } // elem

        if(Verbose)printf("\n");
    } // end loop of num qpts 


    printf("\nAll gradient checks passed.\n");

}
MATAR_FINALIZE();

} // end main



