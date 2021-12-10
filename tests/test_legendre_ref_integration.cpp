
#include "elements.h"

// random parameters to test integration routines
const real_t a0 = 1.23456;
const real_t a1 = 0.98765;
const real_t xi0 = 0.1567;



// ----------------
// a function to integrate
real_t func_term(real_t xi, real_t n){
    return a1*pow( (xi - xi0),n );
}

real_t func(real_t xi, real_t n){
    real_t polyn = a0;
    
    for(int i=1; i<=n; i++){
       polyn += func_term(xi, n);
    }
    
    return polyn;
}
// ----------------


// ----------------
// exact integral of function above
real_t integral_term(real_t n){
    return a1/(n + 1.0)*( pow((1.0 - xi0),n+1) - pow((-1.0 - xi0),n+1) );
}

real_t exact_integral(real_t n){
    real_t result = 2.0*a0;
    
    for(int i=1; i<=n; i++){
        result += integral_term(n);
    }
    
    return result;
}
// ----------------


/** Test Legendre quadrature routines */
int main() {
  
    std::cout << "testing Legendre integration" << std::endl;
    
    for (size_t num_pnts_1D = 1; num_pnts_1D<=19; num_pnts_1D++){
        
        // allocate memory for the 1D Gauss Legendre points and weights
        auto leg_points_1D = CArray <real_t> (num_pnts_1D);
        elements::legendre_nodes_1D(leg_points_1D, num_pnts_1D); // get the locations
        
        auto leg_weights_1D = CArray <real_t> (num_pnts_1D);
        elements::legendre_weights_1D(leg_weights_1D, num_pnts_1D); // get the weights
        
        real_t sum = 0;
        for (int i=0; i<num_pnts_1D; i++){
            sum += leg_weights_1D(i);
        }
    
        real_t integral = 0;
        real_t n = (real_t)(2*num_pnts_1D - 1);  // exact for 2n-1
        for (int i=0; i<num_pnts_1D; i++){
            integral += func(leg_points_1D(i), n)*leg_weights_1D(i);
        }
        
        std::cout
            << " order = " << num_pnts_1D << " , "
            << "sum of weights = " << sum << " , "
            << "relative error = " << (integral - exact_integral(n))/exact_integral(n) << " , "
            << "exact integral(fcn) = " << exact_integral(n) << "\n";
        
    }
    
    std::cout << "finished ---" << std::endl;

  return 0;
}
