#ifndef KOKKOS_SIMPLE_H
#define KOKKOS_SIMPLE_H

#include <stdio.h>
#include "kokkos_alias.h"

#define GAMMA_LAW_SIZE ( sizeof(gamma_law) )
#define SESAME_SIZE ( sizeof(sesame) )


void AllocateHost(MaterialHost1D h_mats, u_int idx, size_t size);

void FreeHost(MaterialHost1D h_mats); 

void InitEOSModels(Material1D mats, u_int idx, gamma_law gamma_law_inp); 
void InitEOSModels(Material1D mats, u_int idx, sesame sesame_inp); 

/*
KOKKOS_FUNCTION
void CreateEOSObjects(Material1D mats, gamma_law gamma_law_inp, u_int idx);
KOKKOS_FUNCTION
void CreateEOSObjects(Material1D mats, sesame sesame_inp, u_int idx);
*/

void ClearDeviceModels(Material1D mats); // would need modification from user by calling the correct deconstructor 

#endif
