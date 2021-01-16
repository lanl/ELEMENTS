#ifndef HEADER_H
#define HEADER_H  

#include "utilities.h"
#include "state.h"
#include "elements.h"

using namespace utils;

struct material_t {

    real_t temp;     // temperature
    real_t cond;     // Thermal conductivity
    real_t den;  // Density
    real_t cp;   // Specific heat
};


namespace region
{

// for tagging boundary faces
enum vol_tag
{
    global = 0,     // tag every cell in the mesh
    box = 1,        // tag all cells inside a box
    cylinder = 2,   // tag all cells inside a cylinder
    sphere = 3      // tag all cells inside a sphere
};

} // end of namespace


// fill instructions
struct mat_fill_t {
    
    // type
    region::vol_tag volume; // 1 is global, 2 are planes, 3 is a cylinder
    
    // material id
    int mat_id;
    
    // planes
    real_t x1;
    real_t x2;
    real_t y1;
    real_t y2;
    real_t z1;
    real_t z2;
    
    // radius
    real_t radius1;
    real_t radius2;

    real_t field1; // some field
    real_t field2; // some other field

};




namespace bdy
{

// for tagging boundary faces
enum bdy_tag
{
    x_plane  = 0,   // tag an x-plane
    y_plane  = 1,   // tag an y-plane
    z_plane  = 2,   // tag an z-plane
    cylinder = 3,   // tag an cylindrical surface
    sphere   = 4,   // tag a spherical surface
    readFile = 5    // read from a file
};



// for enforcing boundary conditions
enum bdy_thermal_conds
{
    isothermal = 0,     // constant temp
    reflected  = 1,     // reflected or wall condition
    free       = 2,     // free surface
};

} // end of bdy namespace


// tag mesh points on bdy's and set the BC type
struct boundary_t {

    // tag surface type
    bdy::bdy_tag surface;    // 0=xplane, 1=yplane, 2=zplane, 3=cylinder, 4=sphere, 5=read file
    
    // tag surface value or radius
    real_t value;
    
    // BC type
    bdy::bdy_thermal_conds thermal_bc;
    
};

#endif // end HEADER_H
