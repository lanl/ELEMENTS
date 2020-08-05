//
//  materials.h
//  
//
//  Created by nmorgan on 3/18/20.
//
// Material variables:
//   mat_vars(mat_gid).initialize_eos(num_mat_pnts)
//   mat_vars(mat_gid).eos->stress(mpnt_gid,i,j)
//   mat_vars(mat_gid).eos->p(mpnt_gid)
#ifndef MATERIAL_H
#define MATERIAL_H
#include "eos.hpp"
//#include "strength.h"


class material_variables{
    
    public:
        //eos_variables *eos_var;
        int num_mat_pnts;
        int mat_type;
        // eos variables
        double *eos_p; // pressure
        //double *eos_d; // density
        //double *eos_sie; // specific internal energy
        //double *eos_m; // mass
        //double *eos_sspd; // sound speed
        //double *eos_velgrad_matrix;  // size 9*num_mat_pts
        //double *eos_stress_matrix;  // size 9*num_mat_pts

        // hypo strength variables
        double *hypo_fake1;
        //double *hypo_fake2;

    // ...
    
        // default constructor
        KOKKOS_FUNCTION
        material_variables() {};

        // init constructor
        //material_variables(int npnts, int mtype) {
        //    num_mat_pnts = npnts;
        //    mat_type     = mtype;
        //};
    
        // deconstructor
        KOKKOS_FUNCTION
        ~material_variables(){};
    
};

class material_models{
    
    public:
        eos_models *eos;
        //hypo_strength_models *hypo_strength;
    // ...
    
    
    // deconstructor
        KOKKOS_FUNCTION
        ~material_models(){};
    
    
}; // end of material


#endif
