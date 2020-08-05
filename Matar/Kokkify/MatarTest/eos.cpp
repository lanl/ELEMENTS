//
//  eos.h
//
//
//  Created by nmorgan on 3/18/20.
//

#include <stdio.h>
//#include <iostream>
//#include "matar.h"
#include "eos.hpp"
//#include "kokkos_alias.h"
//#include <Kokkos_Core.hpp>



//----------------------------
    KOKKOS_FUNCTION
    eos_variables::eos_variables() {};

//----------------------------
    KOKKOS_FUNCTION
    eos_models::eos_models() {};

    KOKKOS_FUNCTION
    gamma_law::gamma_law(double gamma, double cv){
        this_gamma = gamma;
        this_cv = cv;
    }
    
    KOKKOS_FUNCTION
    double gamma_law::pres_from_den_sie(double den, double sie){
        return this_gamma + this_cv;
    }

    KOKKOS_FUNCTION
    sesame::sesame() {
        this_gamma = 0.0;
        this_cv = 0.0;
    }
    
    KOKKOS_FUNCTION
    double sesame::pres_from_den_sie(double den, double sie){
        double sum = den * sie;
        return sum;
    }










