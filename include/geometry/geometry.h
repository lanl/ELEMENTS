#ifndef GEOMETRY_H
#define GEOMETRY_H 


#include "common/utilities.h"
#include "element_types/elements.h"
#include "swage/swage.h"

#include <iostream>

using namespace utils;


void get_gauss_pt_jacobian(swage::mesh_t& mesh, elements::ref_element& ref_elem);

void get_gauss_patch_pt_jacobian(swage::mesh_t& mesh, elements::ref_element& ref_elem);

void get_gauss_cell_pt_jacobian(swage::mesh_t& mesh, elements::ref_element& ref_elem);

void get_vol_jacobi(swage::mesh_t& mesh, elements::ref_element& ref_elem);

void get_vol_hex(swage::mesh_t& mesh, elements::ref_element& ref_elem);



#endif // GEOMETRY_H
