//
//  eos.h
//
//
//  Created by nmorgan on 3/18/20.
//
#ifndef EOS_H
#define EOS_H

#include <stdio.h>
#include <Kokkos_Core.hpp>

enum eos_type
{
    GAMMA_LAW,
    SESAME
};

class eos_variables {
    protected:

    public:
        int num_mat_pnts;
        int mat_type;
        double *p;

        KOKKOS_FUNCTION
        eos_variables();

        KOKKOS_FUNCTION
        ~eos_variables(){};
};


//----------------------------

//----------------------------
class eos_models{

    protected:
        // local parameters
        double this_gamma;
        double this_cv;
    
    public:
        KOKKOS_FUNCTION
        eos_models();

        // eos functions
        KOKKOS_FUNCTION
        virtual double pres_from_den_sie(double den, double sie) {return 0.0;};
        
        KOKKOS_FUNCTION
        virtual ~eos_models(){};
};


class gamma_law: public eos_models{
    
    public:
    
        // constructor
        KOKKOS_FUNCTION
        gamma_law(double gamma, double cv);
    
        // add member functions here to get pressure etc
        KOKKOS_FUNCTION
        double pres_from_den_sie(double den, double sie);
};

class sesame: public eos_models{

    public:
    
        // constructor
        KOKKOS_FUNCTION
        sesame();
    
        // add member functions here to get pressure etc
        KOKKOS_FUNCTION
        double pres_from_den_sie(double den, double sie);
};


#endif
