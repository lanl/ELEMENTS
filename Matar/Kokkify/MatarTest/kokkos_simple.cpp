#include "kokkos_simple.hpp"

void AllocateHost(MaterialHost1D h_mats, u_int idx, size_t size) 
{

    h_mats(idx).eos = (eos_models *) kmalloc(size);

}

void FreeHost(MaterialHost1D h_mats) 
{

    for (int mem = 0; mem < h_mats.extent(0); mem++) {
        kfree(h_mats(mem).eos);
    }
    
}

void InitEOSModels(Material1D mats, u_int idx, sesame sesame_inp) 
{

    Kokkos::parallel_for(
            "CreateObjects", 1, KOKKOS_LAMBDA(const int&) {
                //CreateEOSObjects(mats, sesame_inp, idx);
                new ((sesame *)mats(idx).eos) sesame{sesame_inp};
            });

}


void InitEOSModels(Material1D mats, u_int idx, gamma_law gamma_law_inp) 
{

    Kokkos::parallel_for(
            "CreateObjects", 1, KOKKOS_LAMBDA(const int&) {
                //CreateEOSObjects(mats, gamma_law_inp, idx);
                new ((gamma_law *)mats(idx).eos) gamma_law{gamma_law_inp};
            });

}


void ClearDeviceModels(Material1D mats)
{

    Kokkos::parallel_for(
            "DestroyObjects", 1, KOKKOS_LAMBDA(const int&) {
              mats(0).eos->~eos_models();
              mats(1).eos->~eos_models();
            });

}


// IGNORE
/*

void AllocateKokkosVarHost(MaterialVarHost1D h_mats) 
{

    h_mats(0).eos_var = (eos_variables *) kmalloc(sizeof(eos_variables));
    h_mats(0).eos_var->p = (double *) kmalloc(5*sizeof(double));

}

void FreeKokkosVarHost(MaterialVarHost1D h_mats) 
{

    //for (int mem = 0; mem < h_mats.exten; mem++) {
        kfree(h_mats(0).eos_var->p);
        kfree(h_mats(0).eos_var);
    //}
    
}

KOKKOS_FUNCTION
void CreateEOSVarObjects(MaterialVar1D mats)
{

    // new statement for each material
    // can't do a for loop because of different class types
    new ((eos_variables *)mats(0).eos_var) eos_variables();

}

void InitEOSVars(MaterialVar1D mats) 
{

    parallel_for(
            "CreateObjects", 1, KOKKOS_LAMBDA(const int&) {
                CreateEOSVarObjects(mats);
            });

}

void ClearDeviceVars(MaterialVar1D mats)
{

    parallel_for(
            "DestroyObjects", 1, KOKKOS_LAMBDA(const int&) {
              mats(0).eos_var->~eos_variables();
            });

}

#ifdef FIRST
    ProfileRegionStart("VarInit");
    int mat_pnts = 5;
    MaterialVar1D mat_vars("mat_vars", 1); 
    auto h_mat_vars = HostMirror(mat_vars);
    ProfileRegionEnd();

    ProfileRegionStart("DataInit");
    AllocateKokkosVarHost(h_mat_vars);
    ProfileRegionEnd();
    
    deep_copy(mat_vars, h_mat_vars);

    InitEOSVars(mat_vars);

    ProfileRegionStart("ValueInit");
    parallel_for(
            "InitMatVars",
            range_policy({0}, {mat_pnts}),
            KOKKOS_LAMBDA(const int idx)
            {
                mat_vars(0).eos_var->p[idx] = (real_t) idx;
            }
            );

    deep_copy(h_mat_vars, mat_vars);
    ProfileRegionEnd();

    ProfileRegionStart("VarEnd");
    for (int pr = 0; pr < mat_pnts; pr++) {
        printf("%f ", h_mat_vars(0).eos_var->p[pr]);
    }
    printf("\n");

    ClearDeviceVars(mat_vars); // Kokkos Function to call deconstructors of objects on the GPU

    FreeKokkosVarHost(h_mat_vars); // Function performed on Host to free the allocated GPU classes inside of the Host mirror
    ProfileRegionEnd();
#endif
*/
