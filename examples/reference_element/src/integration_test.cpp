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

const bool Verbose = false;
const size_t max_num   = 19; // max number of quadrature points to test up to, the limit is 19
const size_t max_order = 19; // max polynomial order to test, limit is 19th-order with Legendre


// polynomial with terms <= p_order
KOKKOS_INLINE_FUNCTION
real_t polynomial(const CArrayKokkos<real_t> &coeff, const real_t x, const size_t p_order){

    real_t result = 0.0;
    for (int i = 0; i <= p_order; ++i) {
        result += coeff(i) * pow(x, (real_t)i);
    }
    return result;
} // end polynomial

// Analytical integration of polynomial over [-1, 1]
KOKKOS_INLINE_FUNCTION
real_t integrate_polynomial_1D(const CArrayKokkos<real_t> &coeff, const size_t p_order){
    real_t result = 0.0;
    for (int i = 0; i <= p_order; ++i) {
        // Integral of x^i from -1 to 1
        // = [x^(i+1)/(i+1)] from -1 to 1
        // = (1^(i+1) - (-1)^(i+1))/(i+1)
        if (i % 2 == 0) {  // even power: (-1)^i = 1
            result += coeff(i) * 2.0 / (real_t)(i + 1);
        }
        // odd powers integrate to zero over symmetric interval
    }
    return result;
} // end integrate

// 2D polynomial for QUADS
KOKKOS_INLINE_FUNCTION
real_t polynomial(const CArrayKokkos<real_t> &coeff, const real_t x, 
                  const real_t y, const size_t p_order){
    real_t result = 0.0;
    for (int j = 0; j <= p_order; ++j)      // Full range
    for (int i = 0; i <= p_order; ++i) {    // Full range - NO "-j"
        result += coeff(i,j) * pow(x, (real_t)i) * pow(y, (real_t)j);
    }
    return result;
}

// 2D analytical integration for QUADS
KOKKOS_INLINE_FUNCTION
real_t integrate_polynomial_2D_quad(const CArrayKokkos<real_t> &coeff, 
                                     const size_t p_order){
    real_t result = 0.0;
    for (int j = 0; j <= p_order; ++j)      // Full range
    for (int i = 0; i <= p_order; ++i) {    // Full range - NO "-j"
        real_t xi_integral = 0.0;
        real_t eta_integral = 0.0;
        
        if (i % 2 == 0) xi_integral = 2.0 / (real_t)(i + 1);
        if (j % 2 == 0) eta_integral = 2.0 / (real_t)(j + 1);
        
        result += coeff(i,j) * xi_integral * eta_integral;
    }
    return result;
}

// 3D polynomial for HEXES
KOKKOS_INLINE_FUNCTION
real_t polynomial(const CArrayKokkos<real_t> &coeff, const real_t x, 
                  const real_t y, const real_t z, const size_t p_order){
    real_t result = 0.0;
    for (int k = 0; k <= p_order; ++k)      // Full range
    for (int j = 0; j <= p_order; ++j)      // Full range
    for (int i = 0; i <= p_order; ++i) {    // Full range
        result += coeff(i,j,k) * pow(x, (real_t)i) * 
                                 pow(y, (real_t)j) * 
                                 pow(z, (real_t)k); 
    }
    return result;
}

// 3D analytical integration for HEXES
KOKKOS_INLINE_FUNCTION
real_t integrate_polynomial_3D_hex(const CArrayKokkos<real_t> &coeff, 
                                    const size_t p_order){
    real_t result = 0.0;
    for (int k = 0; k <= p_order; ++k)      // Full range
    for (int j = 0; j <= p_order; ++j)      // Full range  
    for (int i = 0; i <= p_order; ++i) {    // Full range
        real_t xi_integral = 0.0;
        real_t eta_integral = 0.0;
        real_t mu_integral = 0.0;
        
        if (i % 2 == 0) xi_integral = 2.0 / (real_t)(i + 1);
        if (j % 2 == 0) eta_integral = 2.0 / (real_t)(j + 1);
        if (k % 2 == 0) mu_integral = 2.0 / (real_t)(k + 1);
        
        result += coeff(i,j,k) * xi_integral * eta_integral * mu_integral;
    }
    return result;
}


// Test: interpolate a polynomial that the basis can represent exactly
void test_integration(const Quadrature_t& Quad,
                      const ReferenceElement_t& RefElem) {
    
    // For polynomial order p, Lagrange basis can represent 
    // any polynomial of degree <= p exactly
    
    const size_t p_order = RefElem.num_dofs_1d - 1;

    CArrayKokkos<real_t> coeff;

    if(RefElem.elem_dims==1){
        coeff = CArrayKokkos<real_t>(p_order+1);
    }
    else if (RefElem.elem_dims==2){
        coeff = CArrayKokkos<real_t>(p_order+1, p_order+1);
    }
    else {
        coeff = CArrayKokkos<real_t>(p_order+1, p_order+1, p_order+1);
    }
    coeff.set_values(0.78914567);  // a radom value

    // Compute numerical integral using quadratur
    real_t numerical_integral = 0.0;
    real_t sum_lcl = 0.0;

    // evaluate polynomial at quadrature point and sum 
    FOR_REDUCE_SUM(qpt, 0, Quad.num_qpts_in_elem, sum_lcl, {

        real_t value = 0.0;
        if(RefElem.elem_dims==1){
            const real_t xi = Quad.qpt_positions(qpt, 0);
            value = polynomial(coeff, xi, p_order);
        }
        else if (RefElem.elem_dims==2){
            const real_t xi  = Quad.qpt_positions(qpt, 0);
            const real_t eta = Quad.qpt_positions(qpt, 1);
            value = polynomial(coeff, xi, eta, p_order);
        }
        else {
            const real_t xi  = Quad.qpt_positions(qpt, 0);
            const real_t eta = Quad.qpt_positions(qpt, 1);
            const real_t mu  = Quad.qpt_positions(qpt, 2);
            value = polynomial(coeff, xi, eta, mu, p_order);
        }

        sum_lcl += value*Quad.qpt_weights(qpt); 

    }, numerical_integral); // end parallel


    RUN({

        // Compute analytical integral
        real_t exact_integral = 0.0;

        if(RefElem.elem_dims==1){
            exact_integral = integrate_polynomial_1D(coeff, p_order);
        }
        else if (RefElem.elem_dims==2){
            exact_integral = integrate_polynomial_2D_quad(coeff, p_order);
        }
        else {
            exact_integral = integrate_polynomial_3D_hex(coeff, p_order);
        }

        // Compare results
        const real_t error = fabs(numerical_integral - exact_integral);
        const real_t tolerance = fmax(real_tol(1.e-10), real_tol(1.e-10) * (real_t)p_order);

        if (error > tolerance) {
            printf("Error: integration failed with order = %zu \n", p_order);
            printf("tolerance = %.15e \n", tolerance);
            printf("numerical = %.15e vs exact = %.15e, error = %.15e \n", 
                numerical_integral, exact_integral, error);
            Kokkos::abort("Integration test failed");
        }
    
        if (Verbose){
            printf("Integration: numerical = %.15e vs exact = %.15e, error = %.15e \n", 
                numerical_integral, exact_integral, error);
        }

    }); // end RUN on device
} // end integration test


// Diagnostic: Check weight sums in all dimensions
void check_quadrature(const Quadrature_t& Quad) {

    real_t weight_sum = 0.0;
    real_t weight_sum_lcl = 0.0;

    FOR_REDUCE_SUM(qpt, 0, Quad.num_qpts_in_elem, weight_sum_lcl, {
        weight_sum_lcl += Quad.qpt_weights(qpt);
    }, weight_sum);
    Kokkos::fence();
    
    
    real_t expected = pow(2.0, (real_t)Quad.elem_dims);
    real_t error = fabs(weight_sum - expected);
    if (error > real_tol(1.e-12)) {
        printf("Error: quadrature weights don't correctly tally for number = %zu \n", Quad.num_qpts_in_elem);
        printf("numerical = %.15e vs exact = %.15e, error = %.15e \n", 
            weight_sum, expected, error);
        Kokkos::abort("Quadrature weights test failed");
    }

    if(Verbose){
        printf("%zuD: Sum of weights = %.15e ", Quad.elem_dims, weight_sum);
        if(Quad.elem_dims == 1) printf("(expected: 2.0)\n");
        else if(Quad.elem_dims == 2) printf("(expected: 4.0 for quad)\n");
        else if(Quad.elem_dims == 3) printf("(expected: 8.0 for hex)\n");
    }
} // end check quadrature



int main(int argc, char** argv) {

MATAR_INITIALIZE(argc, argv);
{ // MATAR scope

    printf("\n--- integration tests ---\n");


    printf("\n--- checking Lengendre quadrature weight tallies ---\n");
    for(size_t elem_dims_test = 1; elem_dims_test<=3; elem_dims_test++){   

        Quadrature_t Quad;

        for(size_t num_qpts_1D = 1; num_qpts_1D<=max_num; num_qpts_1D++){
            Quad.initialize_quadrature(reference_space::GaussLegendre,
                                       num_qpts_1D,
                                       elem_dims_test);

            check_quadrature(Quad);
        }
        
    } // end for dim


    printf("\n--- checking Lobbatto quadrature weight tallies ---\n");
    for(size_t elem_dims_test = 1; elem_dims_test<=3; elem_dims_test++){   

        Quadrature_t Quad;

        for(size_t num_qpts_1D = 2; num_qpts_1D<=max_num; num_qpts_1D++){
            Quad.initialize_quadrature(reference_space::GaussLobatto,
                                       num_qpts_1D,
                                       elem_dims_test);

            check_quadrature(Quad);
        }
        
    } // end for dim


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

            size_t p_order_ceil = 2*num_qpts_1D-1;  // Legendre
            if(max_order<p_order_ceil) p_order_ceil=max_order;
            for (size_t p_order = 0; p_order<p_order_ceil; p_order++){
                if(Verbose)printf("p_order = %zu: \n", p_order);
                ReferenceElement_t FERefElem;

                // p_order is the basis order for Lagrange polynomial
                FERefElem.initialize_ref_elem(reference_space::arbitraryOrderElement,
                                              reference_space::LagrangeLegendre,
                                              Quad,
                                              p_order);

                if(Verbose)printf("integration check: \n");
                test_integration(Quad, FERefElem);
                if(Verbose)printf("\n");

            } // end p_order loop
        } // elem

        if(Verbose)printf("\n");
    } // end loop of num qpts 


    printf("\n--- DG element with Lobatto Quadrature & Legendre DOFs ---\n");
    for(size_t num_qpts_1D = 2; num_qpts_1D<=max_num; num_qpts_1D++){

        if(Verbose)printf("num quadrature points in 1D = %zu \n", num_qpts_1D);

        Quadrature_t Quad;
        
        // elem_dims=1,2,3
        for(size_t elem_dims_test = 1; elem_dims_test<=3; elem_dims_test++){   
            Quad.initialize_quadrature(reference_space::GaussLobatto,
                                       num_qpts_1D,
                                       elem_dims_test);

            // build reference elements of varing orders

            size_t p_order_ceil = 2*num_qpts_1D-3;  // Lobatto starts at 2 points
            if(max_order<p_order_ceil) p_order_ceil=max_order;
            for (size_t p_order = 0; p_order<=p_order_ceil; p_order++){
                if(Verbose)printf("p_order = %zu: \n", p_order);
                ReferenceElement_t FERefElem;

                // p_order is the basis order for Lagrange polynomial
                FERefElem.initialize_ref_elem(reference_space::arbitraryOrderElement,
                                              reference_space::LagrangeLegendre,
                                              Quad,
                                              p_order);

                if(Verbose)printf("integration check: \n");
                test_integration(Quad, FERefElem);
                if(Verbose)printf("\n");

            } // end p_order loop
        } // elem

        if(Verbose)printf("\n");
    } // end loop of num qpts 



    printf("\n--- FE element with Lobatto Quadrature & Lobatto DOFs ---\n");
    for(size_t num_qpts_1D = 2; num_qpts_1D<=max_num; num_qpts_1D++){

        if(Verbose)printf("num quadrature points in 1D = %zu \n", num_qpts_1D);

        Quadrature_t Quad;
        
        // elem_dims=1,2,3
        for(size_t elem_dims_test = 1; elem_dims_test<=3; elem_dims_test++){   
            Quad.initialize_quadrature(reference_space::GaussLobatto,
                                       num_qpts_1D,
                                       elem_dims_test);

            // build reference elements of varing orders, starting at 1, its Lobatto points
            
            size_t p_order_ceil = 2*num_qpts_1D-3;  // Lobatto starts at 2 points
            if(max_order<p_order_ceil) p_order_ceil=max_order;
            for (size_t p_order = 1; p_order<=p_order_ceil; p_order++){
                if(Verbose)printf("p_order = %zu: \n", p_order);
                ReferenceElement_t FERefElem;

                // p_order is the basis order for Lagrange polynomial
                FERefElem.initialize_ref_elem(reference_space::arbitraryOrderElement,
                                              reference_space::LagrangeLobatto,
                                              Quad,
                                              p_order);

                if(Verbose)printf("integration check: \n");
                test_integration(Quad, FERefElem);
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
            size_t p_order_ceil = 2*num_qpts_1D-1;  // Legendre
            if(max_order<p_order_ceil) p_order_ceil=max_order;
            for (size_t p_order = 1; p_order<p_order_ceil; p_order++){
                if(Verbose)printf("p_order = %zu: \n", p_order);
                ReferenceElement_t FERefElem;

                // p_order is the basis order for Lagrange polynomial
                FERefElem.initialize_ref_elem(reference_space::arbitraryOrderElement,
                                              reference_space::LagrangeLobatto,
                                              Quad,
                                              p_order);

                if(Verbose)printf("integration check: \n");
                test_integration(Quad, FERefElem);
                if(Verbose)printf("\n");

            } // end p_order loop
        } // elem

        if(Verbose)printf("\n");
    } // end loop of num qpts 


    printf("\nAll integration checks passed.\n");

}
MATAR_FINALIZE();

} // end main

