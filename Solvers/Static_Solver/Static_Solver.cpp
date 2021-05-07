/* 

Example code for smoothing some field on the mesh. 


A representative mesh is shown:

p
*---------*---------*
|         |         |
|         |         |
|    *z   |    *    |
|         |         |
|         |         |
*---------*---------*
|         |         |
|         |         |
|    *    |    *    |
|         |         |
|         |         |
*---------*---------*

The smoothing operation follows a two step process:

1. ) Loop over all the nodes (p) in a cell and 
average the field to the cell center materhial
point (z). 

2.) Loop over all of the cells (z) connected to a node (p)
and average values to the nodal field.


Each cell is within an element, and the number of cells is 
defined by the user using the p_order variable in the input

num_cells in element = (p_order*2)^3

*/


#include <iostream>
#include <stdio.h>
#include <stdlib.h> 
#include <math.h>  // fmin, fmax, abs note: fminl is long
#include <sys/stat.h>
#include <mpi.h>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Kokkos_View.hpp>
#include <Kokkos_Parallel.hpp>
#include <Kokkos_Parallel_Reduce.hpp>
#include "Tpetra_Details_makeColMap.hpp"
#include "Tpetra_Details_DefaultTypes.hpp"
#include <set>

#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"

#include "elements.h"
#include "swage.h"
#include "matar.h"
#include "utilities.h"
#include "header.h"
#include "state.h"
#include "Simulation_Parameters.h"
#include "Static_Solver.h"
#include "Teuchos_RCP.hpp"
#include "Amesos2_Version.hpp"
#include "Amesos2.hpp"

using namespace utils;

/*

Swage is a reference to a swage block used in blacksmithing.  
Its a large metal block that has multiple shaps carved into 
each surface to use for hammering metal into to form it. 

*/

Static_Solver::Static_Solver() : Solver(){
  //create parameter object
    simparam = new Simulation_Parameters();
  //create ref element object
    ref_elem = new elements::ref_element();
  //create mesh objects
    init_mesh = new swage::mesh_t();
    mesh = new swage::mesh_t();

    element_select = new elements::element_selector();
}

Static_Solver::~Static_Solver(){
   delete simparam;
   delete ref_elem;
   delete init_mesh;
   delete mesh;
   delete element_select;
}

//==============================================================================
//    Primary simulation runtime routine
//==============================================================================


void Static_Solver::run(int argc, char *argv[]){
    
    std::cout << "Running Static Solver Example" << std::endl;

    // ---- Read input file, define state and boundary conditions ---- //
    simparam->input();

    // ---- Read intial mesh, refine, and build connectivity ---- //
    read_mesh(argv[1]);
    std::cout << "Num elements = " << mesh->num_elems() << std::endl;

    // ---- Find Boundaries on mesh ---- //
    generate_bcs();
    
    //allocate and fill sparse structures needed for global solution
    init_global();

    //assemble the global solution (stiffness matrix etc. and nodal forces)
    assemble();

    //debug solve; move later after bcs applies
    int solver_exit = solve();
    if(solver_exit == EXIT_SUCCESS) return;
    std::cout << "Before free pointer  " << std::endl <<std::flush;
    //debug return to avoid printing further
    return;

    std::cout << "finished linear solver" << std::endl;

    // ---- Allocate memory for state on mesh ---- //
    allocate_state();

    initialize_state();

    std::cout << "Before boundary  " << std::endl;
    apply_boundary();


    // Calculate reference element information
    calculate_ref_elem();

    get_nodal_jacobian();

    real_t dt = simparam->dt;
    int cycle_stop = simparam->cycle_stop;
    real_t &TIME = simparam->TIME;
    real_t TFINAL = simparam->TFINAL;
    int &cycle = simparam->cycle;
    int graphics_cyc_ival = simparam->graphics_cyc_ival;
    real_t graphics_dt_ival = simparam->graphics_dt_ival;
    real_t graphics_time = simparam->graphics_dt_ival;
    real_t &percent_comp = simparam->percent_comp;


    // Solve the Laplacian


    ensight();

    for (cycle = 1; cycle <= cycle_stop; cycle++) {
        // stop calculation if flag
        if (simparam->stop_calc == 1) break;


        // Set timestep
        dt = simparam->dt;


        // smooth_cells();

        smooth_element();

        apply_boundary();

        // increment the time
        TIME += dt;



        // Runtime outputs

        if ((cycle%10) == 0) {
            percent_comp = (TIME/TFINAL)*100.;
            printf("Percent complete = %.2f \n", percent_comp); 
        }
        
        // end of calculation
        if (TIME>=TFINAL) break;

        
        if ((cycle%50) == 0) {
            printf("Step = %d   dt = %e  time = %e  \n ", cycle, dt, TIME);
        }

        
        // output a graphics dump
        if ((cycle%graphics_cyc_ival) == 0 || (TIME-graphics_time) >= 0.0) {
            
            printf("****************Graphics write*************** \n");
            ensight();
            
            // set next time to dump a graphics output
            graphics_time += graphics_dt_ival;
        }
        
        
    } // end of cycle loop

    // Data writers
    ensight();
    // vtk_writer();

    std::cout << "End of Main file" << std::endl;

}

//==============================================================================
//   Function Definitions
//==============================================================================


// Input and setup
void Static_Solver::read_mesh(char *MESH){

    FILE *in;
    char ch;
    int num_dim = simparam->num_dim;
    int p_order = simparam->p_order;

    //read the mesh    WARNING: assumes an ensight .geo file
    in = fopen(MESH,"r");  
    
    //skip 8 lines
    for (int j = 1; j <= 8; j++) {
        int i=0;
        while ((ch=(char)fgetc(in))!='\n') {
            i++;
            printf("%c",ch);
        }
        printf("\n");
    }  

    int num_nodes;

    // --- Read the number of vertices in the mesh --- //
    fscanf(in,"%d",&num_nodes);
    printf("%d\n" , num_nodes);

    // set the vertices in the mesh read in
    init_mesh->init_nodes(num_nodes); // add 1 for index starting at 1
    
    std::cout << "Num points read in = " << init_mesh->num_nodes() << std::endl;


    // read the initial mesh coordinates
    // x-coords
    for (int node_gid = 0; node_gid < init_mesh->num_nodes(); node_gid++) {
        fscanf(in,"%le", &init_mesh->node_coords(node_gid, 0));
        //std::cout<<" "<< init_mesh->node_coords(node_gid, 0)<<std::endl;
    }
    
    if(num_dim>1)
    // y-coords
    for (int node_gid = 0; node_gid < init_mesh->num_nodes(); node_gid++) {
        fscanf(in,"%le", &init_mesh->node_coords(node_gid, 1));
        //std::cout<<" "<< init_mesh->node_coords(node_gid, 1)<<std::endl;
    }  

    // z-coords
    if(num_dim>2)
    for (int node_gid = 0; node_gid < init_mesh->num_nodes(); node_gid++) {
        fscanf(in,"%le", &init_mesh->node_coords(node_gid, 2));
        //std::cout<<" "<< init_mesh->node_coords(node_gid, 2)<<std::endl;
    }

    /*
    ch = (char)fgetc(in);
    printf("%c",ch);

    //skip 1 line
    for (int j=1; j<=1; j++) {
        int i=0;
        while ((ch=(char)fgetc(in))!='\n') {
            i++;
            printf("%c",ch);
        }
        printf("\n");
    }
    */
    char skip_string[20];
    fscanf(in,"%s",skip_string);
    std::cout << skip_string <<std::endl;
    int num_elem = 0;
    
    // --- read the number of cells in the mesh ---
    fscanf(in,"%d",&num_elem);
    printf("Num elements read in %d\n" , num_elem);

    std::cout<<"before initial mesh initialization"<<std::endl;

    init_mesh->init_element(0, 3, num_elem);
    init_mesh->init_cells(num_elem);

    Element_Types = CArray<size_t>(num_elem); 



    // for each cell read the list of associated nodes
    for (int cell_gid = 0; cell_gid < num_elem; cell_gid++) {
        for (int node_lid = 0; node_lid < 8; node_lid++){
            
            fscanf(in,"%d",&init_mesh->nodes_in_cell(cell_gid, node_lid));
            //std::cout<<" "<< init_mesh->nodes_in_cell(cell_gid, node_lid);
            // shift to start vertex index space at 0
            init_mesh->nodes_in_cell(cell_gid,node_lid) = init_mesh->nodes_in_cell(cell_gid, node_lid) - 1;
        }
        //std::cout<<" "<<std::endl;
        Element_Types(cell_gid) = 4; //Hex8 default
    }

    //element type selection (subject to change)
    // ---- Set Element Type ---- //
    // allocate element type memory
    elements::elem_type_t* elem_choice;

    int NE = 1; // number of element types in problem
    elem_choice = new elements::elem_type_t[NE];
    
    //current default
    elem_choice->type = elements::elem_types::elem_type::Hex8;
    
    //set base type pointer to one of the existing derived type object references
    if(simparam->num_dim==2)
    element_select->choose_2Delem_type(elem_choice[0], elem2D);
    else if(simparam->num_dim==3)
    element_select->choose_3Delem_type(elem_choice[0], elem);

    // Convert ijk index system to the finite element numbering convention
    // for vertices in cell
    int convert_ensight_to_ijk[8];
    convert_ensight_to_ijk[0] = 0;
    convert_ensight_to_ijk[1] = 1;
    convert_ensight_to_ijk[2] = 3;
    convert_ensight_to_ijk[3] = 2;
    convert_ensight_to_ijk[4] = 4;
    convert_ensight_to_ijk[5] = 5;
    convert_ensight_to_ijk[6] = 7;
    convert_ensight_to_ijk[7] = 6;
    
    int tmp_ijk_indx[8];

    for (int cell_gid = 0; cell_gid < num_elem; cell_gid++) {
        for (int node_lid = 0; node_lid < 8; node_lid++){
    
            tmp_ijk_indx[node_lid] = init_mesh->nodes_in_cell(cell_gid, convert_ensight_to_ijk[node_lid]);
        }   
        
        for (int node_lid = 0; node_lid < 8; node_lid++){

            init_mesh->nodes_in_cell(cell_gid, node_lid) = tmp_ijk_indx[node_lid];
        }
    }

    // Build all connectivity in initial mesh
    std::cout<<"Before initial mesh connectivity"<<std::endl;

    if(num_elem < 0) std::cout << "ERROR, NO ELEMENTS IN MESH" << std::endl;
    if(num_elem > 1) {

        // -- NODE TO CELL CONNECTIVITY -- //
        init_mesh->build_node_cell_connectivity(); 

        // -- CORNER CONNECTIVITY -- //
        init_mesh->build_corner_connectivity(); 

        // -- CELL TO CELL CONNECTIVITY -- //
        init_mesh->build_cell_cell_connectivity(); 

        // -- PATCHES -- //
        init_mesh->build_patch_connectivity(); 
    }

    std::cout<<"refine mesh"<<std::endl;

    // refine subcell mesh to desired order
    swage::refine_mesh(*init_mesh, *mesh, 0, num_dim);
    std::cout<<"done refining"<< simparam->num_dim <<std::endl;
    // Close mesh input file
    fclose(in);
    
    // Create reference element
    ref_elem->init(p_order, num_dim, elem->num_basis());
    std::cout<<"done with ref elem"<<std::endl;
     
    std::cout << "number of patches = " << mesh->num_patches() << std::endl;
    std::cout << "End of setup " << std::endl;     
} // end read_mesh


void Static_Solver::generate_bcs(){
    
    // build boundary mesh patches
    mesh->build_bdy_patches();
    std::cout << "number of boundary patches = " << mesh->num_bdy_patches() << std::endl;
    std::cout << "building boundary sets " << std::endl;
    // set the number of boundary sets
    
    int num_bdy_sets = simparam->NB;
    int num_surface_force_sets = simparam->NBSF;
    int num_surface_disp_sets = simparam->NBD;
    int num_dim = simparam->num_dim;
    int num_nodes = mesh->num_nodes();
    int current_bdy_id = 0;
    int bdy_set_id;
    int surf_force_set_id = 0;
    int surf_disp_set_id = 0;

    mesh->init_bdy_sets(num_bdy_sets);
    Boundary_Condition_Type_List = CArray<int>(num_bdy_sets); 
    Boundary_Surface_Force_Densities = CArray<real_t>(num_surface_force_sets,3);
    Boundary_Surface_Displacements = CArray<real_t>(num_surface_disp_sets,3);
    //initialize
    for(int ibdy=0; ibdy < num_bdy_sets; ibdy++) Boundary_Condition_Type_List(ibdy) = NONE;
    
    // tag the x=0 plane,  (Direction, value, bdy_set)
    std::cout << "tagging x = 0 " << std::endl;
    int bc_tag = 0;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
    real_t value = 0.0;
    bdy_set_id = current_bdy_id++;
    mesh->tag_bdys(bc_tag, value, bdy_set_id);
    Boundary_Condition_Type_List(bdy_set_id) = DISPLACEMENT_CONDITION;
    Boundary_Surface_Displacements(surf_disp_set_id,0) = 0;
    Boundary_Surface_Displacements(surf_disp_set_id,1) = 0;
    Boundary_Surface_Displacements(surf_disp_set_id,2) = 0;
    surf_disp_set_id++;
    
    std::cout << "tagged a set " << std::endl;
    std::cout << "number of bdy patches in this set = " << mesh->num_bdy_patches_in_set(bdy_set_id) << std::endl;
    std::cout << std::endl;
    /*
    //This part should be changed so it interfaces with simparam to handle multiple input cases
    // tag the y=0 plane,  (Direction, value, bdy_set)
    std::cout << "tagging y = 0 " << std::endl;
    bc_tag = 1;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
    value = 0.0;
    bdy_set_id = 1;
    mesh->tag_bdys(bc_tag, value, bdy_set_id);
    Boundary_Condition_Type_List(bdy_set_id) = DISPLACEMENT_CONDITION;
    std::cout << "tagged a set " << std::endl;
    std::cout << "number of bdy patches in this set = " << mesh->num_bdy_patches_in_set(bdy_set_id) << std::endl;
    std::cout << std::endl;
    

    // tag the z=0 plane,  (Direction, value, bdy_set)
    std::cout << "tagging z = 0 " << std::endl;
    bc_tag = 2;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
    value = 0.0;
    bdy_set_id = 2;
    mesh->tag_bdys(bc_tag, value, bdy_set_id);
    Boundary_Condition_Type_List(bdy_set_id) = DISPLACEMENT_CONDITION;
    std::cout << "tagged a set " << std::endl;
    std::cout << "number of bdy patches in this set = " << mesh->num_bdy_patches_in_set(bdy_set_id) << std::endl;
    std::cout << std::endl;
    */

    std::cout << "tagging x = 2 " << std::endl;
    bc_tag = 0;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
    value = 1.0;
    bdy_set_id = current_bdy_id++;
    mesh->tag_bdys(bc_tag, value, bdy_set_id);
    Boundary_Condition_Type_List(bdy_set_id) = LOADING_CONDITION;
    Boundary_Surface_Force_Densities(surf_force_set_id,0) = 0;
    Boundary_Surface_Force_Densities(surf_force_set_id,1) = 2;
    Boundary_Surface_Force_Densities(surf_force_set_id,2) = 0;
    surf_force_set_id++;
    std::cout << "tagged a set " << std::endl;
    std::cout << "number of bdy patches in this set = " << mesh->num_bdy_patches_in_set(bdy_set_id) << std::endl;
    std::cout << std::endl;

    /*
    std::cout << "tagging y = 2 " << std::endl;
    bc_tag = 1;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
    value = 2.0;
    bdy_set_id = 4;
    mesh->tag_bdys(bc_tag, value, bdy_set_id);
    Boundary_Condition_Type_List(bdy_set_id) = LOADING_CONDITION;
    std::cout << "tagged a set " << std::endl;
    std::cout << "number of bdy patches in this set = " << mesh->num_bdy_patches_in_set(bdy_set_id) << std::endl;
    std::cout << std::endl;

    std::cout << "tagging z = 2 " << std::endl;
    bc_tag = 2;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
    value = 2.0;
    bdy_set_id = 5;
    mesh->tag_bdys(bc_tag, value, bdy_set_id);
    Boundary_Condition_Type_List(bdy_set_id) = LOADING_CONDITION;
    std::cout << "tagged a set " << std::endl;
    std::cout << "number of bdy patches in this set = " << mesh->num_bdy_patches_in_set(bdy_set_id) << std::endl;
    std::cout << std::endl;
    */

    //allocate nodal data
    Node_DOF_Boundary_Condition_Type = CArray<int>(num_nodes *num_dim);
    Node_DOF_Displacement_Boundary_Conditions = CArray<real_t>(num_nodes *num_dim);
    Node_DOF_Force_Boundary_Conditions = CArray<real_t>(num_nodes *num_dim);

    //initialize
    for(int init=0; init < num_nodes*num_dim; init++)
      Node_DOF_Boundary_Condition_Type(init) = NONE;
} // end generate_bcs


void Static_Solver::allocate_state(){

    int rk_storage = simparam->rk_storage;
    int num_dim = simparam->num_dim;
    node_t *node = simparam->node;
    mat_pt_t *mat_pt = simparam->mat_pt;

    std::cout << "Allocate and Initialize"  << std::endl;
    std::cout << "RK num stages = "<< rk_storage  << std::endl;
    // --- allocate and initialize the defaults for the problem ---

    // ---- Node initialization ---- //
    node->init_node_state(num_dim, *mesh, rk_storage);
    std::cout << "Node state allocated and initialized to zero"  << std::endl;
    std::cout << std::endl;

    // ---- Material point initialization ---- //
    mat_pt->init_mat_pt_state(num_dim, *mesh, rk_storage);
    std::cout << "Material point state allocated and initialized to zero"  << std::endl;
    std::cout << std::endl;
} // end allocate_state


void Static_Solver::initialize_state(){
    int NF = simparam->NF;
    mat_pt_t *mat_pt = simparam->mat_pt;
    mat_fill_t *mat_fill = simparam->mat_fill;
    int rk_stage = simparam->rk_stage;

    std::cout << "Before fill instructions"  << std::endl;
    //--- apply the fill instructions ---//
    for (int f_id = 0; f_id < NF; f_id++){
        
        for (int cell_gid = 0; cell_gid < mesh->num_cells(); cell_gid++) {
            
            // calculate the coordinates and radius of the cell
            real_t cell_coords_x = 0.0;
            real_t cell_coords_y = 0.0;
            real_t cell_coords_z = 0.0;
            
            for (int node_lid = 0; node_lid < 8; node_lid++){
                
                // increment the number of cells attached to this vertex
                int vert_gid = mesh->nodes_in_cell(cell_gid, node_lid); // get the global_id
                
                cell_coords_x += mesh->node_coords(vert_gid, 0);
                cell_coords_y += mesh->node_coords(vert_gid, 1);
                cell_coords_z += mesh->node_coords(vert_gid, 2);
                
            }// end for loop over node_lid
            
            cell_coords_x = cell_coords_x/8.0;
            cell_coords_y = cell_coords_y/8.0;
            cell_coords_z = cell_coords_z/8.0;

            // Material points at cell center

            mat_pt->coords(rk_stage, cell_gid, 0) = cell_coords_x;
            mat_pt->coords(rk_stage, cell_gid, 1) = cell_coords_y;
            mat_pt->coords(rk_stage, cell_gid, 2) = cell_coords_z;

            
            // spherical radius
            real_t radius = sqrt( cell_coords_x*cell_coords_x +
                                  cell_coords_y*cell_coords_y +
                                  cell_coords_z*cell_coords_z );
            
            // cylinderical radius
            real_t radius_cyl = sqrt( cell_coords_x*cell_coords_x +
                                      cell_coords_y*cell_coords_y);
            
            
            // default is not to fill the cell
            int fill_this = 0;
            
            // check to see if this cell should be filled
            switch(mat_fill[f_id].volume)
            {
                case region::global:
                {
                    fill_this = 1;
                }
                case region::box:
                {
                    if ( cell_coords_x >= mat_fill[f_id].x1
                        && cell_coords_x <= mat_fill[f_id].x2
                        && cell_coords_y >= mat_fill[f_id].y1
                        && cell_coords_y <= mat_fill[f_id].y2
                        && cell_coords_z >= mat_fill[f_id].z1
                        && cell_coords_z <= mat_fill[f_id].z2 ) fill_this = 1;
                }
                case region::cylinder:
                {
                    if ( radius_cyl >= mat_fill[f_id].radius1
                      && radius_cyl <= mat_fill[f_id].radius2 ) fill_this = 1;
                }
                case region::sphere:
                {
                    if ( radius >= mat_fill[f_id].radius1
                      && radius <= mat_fill[f_id].radius2 ) fill_this = 1;
                }
            } // end of switch

            if (fill_this == 1){    
                
                mat_pt->field(cell_gid) = mat_fill[f_id].field2;


            } // end if fill volume
        } // end for cell loop
    } // end for fills

    std::cout << "After fill instructions"  << std::endl;
} // end initialize_state


void Static_Solver::calculate_ref_elem(){
    int num_dim = simparam->num_dim;

    real_t partial_xia[elem->num_basis()];
    auto partial_xi = ViewCArray <real_t> (partial_xia, elem->num_basis());

    real_t partial_etaa[elem->num_basis()];
    auto partial_eta = ViewCArray <real_t> (partial_etaa, elem->num_basis());

    real_t partial_mua[elem->num_basis()];
    auto partial_mu = ViewCArray <real_t> (partial_mua, elem->num_basis());


    std::cout << "::::  Getting partials of basis  ::::" << std::endl;
    std::cout << "Num Basis =  "<< elem->num_basis() << std::endl;

    std::cout << "Num_ref_nodes  "<< ref_elem->num_ref_nodes() << std::endl;

    for(int node_rid = 0; node_rid < ref_elem->num_ref_nodes(); node_rid++){

        // make temp array of ref node positions
        real_t ref_node_loc_a[mesh->num_dim()];
        auto ref_node_loc = ViewCArray<real_t> (ref_node_loc_a, mesh->num_dim());

        for(int dim = 0; dim < mesh->num_dim(); dim++){

            ref_node_loc(dim) = ref_elem->ref_node_positions(node_rid, dim);
        }
    
        std::cout << "Local Node =  "<< node_rid << std::endl;

        // Calculate array of partials of each basis at the point ref_node
        elem->partial_xi_basis(partial_xi, ref_node_loc);
        elem->partial_eta_basis(partial_eta, ref_node_loc);
        elem->partial_mu_basis(partial_mu, ref_node_loc);
        
        // Save the partials of each basis to the reference node
        
        for(int basis_id = 0; basis_id < elem->num_basis(); basis_id++){
            
            ref_elem->ref_nodal_gradient(node_rid, basis_id, 0) = partial_xi(basis_id);
            ref_elem->ref_nodal_gradient(node_rid, basis_id, 1) = partial_eta(basis_id);
            ref_elem->ref_nodal_gradient(node_rid, basis_id, 2) = partial_mu(basis_id);

            // std::cout << "Partial Xi for basis "<< basis_id <<" = " << partial_xi(basis_id) << std::endl;
            // std::cout << "Partial Eta for basis "<< basis_id <<" = " << partial_eta(basis_id) << std::endl;
            // std::cout << "Partial Mu for basis "<< basis_id <<" = " << partial_eta(basis_id) << std::endl;
            // std::cout << "Basis " << basis_id  << std::endl;
            // std::cout << "XX Partial Xi = "  << ref_elem.ref_nodal_gradient(node_rid, basis_id, 0) << std::endl;
            // std::cout << "XX Partial Eta = " << ref_elem.ref_nodal_gradient(node_rid, basis_id, 1) << std::endl;
            // std::cout << "XX Partial Mu = "  << ref_elem.ref_nodal_gradient(node_rid, basis_id, 2) << std::endl;
            



            partial_xi(basis_id) = 0.0;
            partial_eta(basis_id)= 0.0;
            partial_mu(basis_id) = 0.0;

        }


        std::cout << std::endl;

    } // end for node_rid


    // for(int node_rid = 0; node_rid < ref_elem.num_ref_nodes(); node_rid++){

    //     std::cout << "Local Reference Node =  "<< node_rid << std::endl;

    //     for(int basis_id = 0; basis_id < elem->num_basis(); basis_id++){
            
    //         std::cout << "Basis " << basis_id  << std::endl;
    //         std::cout << "Partial Xi = "  << ref_elem.ref_nodal_gradient(node_rid, basis_id, 0) << std::endl;
    //         std::cout << "Partial Eta = " << ref_elem.ref_nodal_gradient(node_rid, basis_id, 1) << std::endl;
    //         std::cout << "Partial Mu = "  << ref_elem.ref_nodal_gradient(node_rid, basis_id, 2) << std::endl;
    //     }

    //     std::cout << std::endl;
    // } // end vert lid

    std::cout << "Finished gradient vector" << std::endl;

}


// Runtime Functions

void Static_Solver::apply_boundary(){

    node_t *node = simparam->node;

    for (int bdy_patch_gid = 0; bdy_patch_gid < mesh->num_bdy_patches_in_set(3); bdy_patch_gid++){
                
        // get the global id for this boundary patch
        int patch_gid = mesh->bdy_patches_in_set(3, bdy_patch_gid);

        // apply boundary condition at nodes on boundary
        for(int node_lid = 0; node_lid < 4; node_lid++){
            
            int node_gid = mesh->node_in_patch(patch_gid, node_lid);

            // Set nodal temp to zero
            node->field(node_gid) = 0.0;

        }
    }

    // std::cout << "Apply temp" << std::endl;
    // Apply constant temp of 1 to x=0 plane of mesh
    for (int bdy_patch_gid = 0; bdy_patch_gid < mesh->num_bdy_patches_in_set(0); bdy_patch_gid++){
                
        // get the global id for this boundary patch
        int patch_gid = mesh->bdy_patches_in_set(0, bdy_patch_gid);

        // apply boundary condition at nodes on boundary
        for(int node_lid = 0; node_lid < 4; node_lid++){
            
            int node_gid = mesh->node_in_patch(patch_gid, node_lid);

            // Set nodal temp to zero
            node->field(node_gid) = 40.0;

        }
    }
}


void Static_Solver::smooth_cells(){
    mat_pt_t *mat_pt = simparam->mat_pt;
    mat_fill_t *mat_fill = simparam->mat_fill;
    node_t *node = simparam->node;

    // Walk over cells 
    for(int cell_gid = 0; cell_gid < mesh->num_cells(); cell_gid++){

        // Temporary holder variable
        real_t temp_avg = 0.0;

        // Walk over nodes in cell and calculate average
        for(int node_lid = 0; node_lid < mesh->num_nodes_in_cell(); node_lid++){

            // Get global index of the node
            int node_gid = mesh->nodes_in_cell(cell_gid, node_lid);

            temp_avg += node->field(node_gid)/8.0;

        }

        // Save average to material point at cell center
        mat_pt->field(cell_gid) = temp_avg;

    } // end of loop over cells

    // Walk over all the nodes
    for(int node_gid = 0; node_gid < mesh->num_nodes(); node_gid++){

        real_t temp_avg = 0.0;
        
        // Walk over all the cells connected to the node and average values
        for(int cell_lid = 0; cell_lid < mesh->num_cells_in_node(node_gid); cell_lid++){

            // Get global index of the cell
            int cell_gid = mesh->cells_in_node(node_gid, cell_lid);

            temp_avg += mat_pt->field(cell_gid)/ (real_t)mesh->num_cells_in_node(node_gid);
        
        }

        // Save average to the node
        node->field(node_gid) = temp_avg;

    } // end of loop over nodes
}


void Static_Solver::smooth_element(){
    mat_pt_t *mat_pt = simparam->mat_pt;
    mat_fill_t *mat_fill = simparam->mat_fill;
    node_t *node = simparam->node;

    // Walk over each element in the mesh
    for(int elem_gid = 0; elem_gid < mesh->num_elems(); elem_gid++){

        // Walk over each cell in the element
        for(int cell_lid = 0; cell_lid < mesh->num_cells_in_elem(); cell_lid++){

            // Get the global ID of the cell
            int cell_gid = mesh->cells_in_elem(elem_gid, cell_lid);

            real_t temp_avg = 0.0;

            // Loop over nodes in the cell
            for(int node_lid = 0; node_lid < mesh->num_nodes_in_cell(); node_lid++){

                // Get global ID for this node
                int node_gid = mesh->nodes_in_cell(cell_gid, node_lid);

                temp_avg += node->field(node_gid)/mesh->num_nodes_in_cell();

            }

            // Save averaged values to cell centered material point
            mat_pt->field(cell_gid) = temp_avg;
        }// end loop over nodes
    }// end loop over elements


    // Walk over each element in the mesh
    for(int elem_gid = 0; elem_gid < mesh->num_elems(); elem_gid++){
       
        // Walk over each node in the element
        for(int node_lid = 0; node_lid < mesh->num_nodes_in_elem(); node_lid++){

            // Get global ID of the node
            int node_gid = mesh->nodes_in_elem(elem_gid, node_lid);

            real_t temp_avg = 0.0;

            // Walk over all cell connected to the node
            for(int cell_lid = 0; cell_lid < mesh->num_cells_in_node(node_gid); cell_lid++){

                // Get globa ID for the cell
                int cell_gid = mesh->cells_in_node(node_gid, cell_lid);

                temp_avg += mat_pt->field(cell_gid)/ (real_t)mesh->num_cells_in_node(node_gid);
            }

            // Save averaged field to node
            node->field(node_gid) = temp_avg;

        }// end loop over nodes
    }// end loop over elements
}

void Static_Solver::get_nodal_jacobian(){
  mat_pt_t *mat_pt = simparam->mat_pt;
  mat_fill_t *mat_fill = simparam->mat_fill;
  int num_dim = simparam->num_dim;

  // loop over the mesh

    for(int elem_gid = 0; elem_gid < mesh->num_elems(); elem_gid++){
        
        for(int cell_lid = 0; cell_lid < mesh->num_cells_in_elem(); cell_lid++){ // 1 for linear elements

            int cell_gid = mesh->cells_in_elem(elem_gid, cell_lid);

            // Loop over nodes in cell and initialize jacobian to zero
            for(int node_lid = 0; node_lid < mesh->num_nodes_in_cell(); node_lid++){

                int node_gid = mesh->nodes_in_cell(cell_gid, node_lid);
                
                for(int dim_i = 0; dim_i < mesh->num_dim(); dim_i++){
                    for(int dim_j = 0; dim_j < mesh->num_dim(); dim_j++){

                        mesh->node_jacobian(node_gid, dim_i, dim_j) = 0.0;
                    }
                }
            }

            // Calculate the actual jacobian for that node
            for(int node_lid = 0; node_lid < mesh->num_nodes_in_cell(); node_lid++){
                
                int node_gid = mesh->nodes_in_cell(cell_gid, node_lid);



                for(int dim_i = 0; dim_i < mesh->num_dim(); dim_i++){
                    for(int dim_j = 0; dim_j < mesh->num_dim(); dim_j++){

                        // Sum over the basis functions and nodes where they are defined
                        for(int basis_id = 0; basis_id < mesh->num_nodes_in_cell(); basis_id++){

                            int ref_node_gid = mesh->nodes_in_cell(cell_gid, basis_id);

                            mesh->node_jacobian(node_gid, dim_i, dim_j) += 
                                mesh->node_coords(ref_node_gid , dim_i) * ref_elem->ref_nodal_gradient(node_lid, basis_id, dim_j);
                        }
                    }
                }
            }

        }
    }

// NOTE: Save only J^inverse and det_J
#pragma omp simd //Modified by Daniel
    // loop over the nodes of the mesh
    for (int node_gid = 0; node_gid < mesh->num_nodes(); node_gid++) {

        mesh->node_det_j(node_gid) = 
            mesh->node_jacobian(node_gid, 0, 0) 
          * ( (mesh->node_jacobian(node_gid, 1, 1)*mesh->node_jacobian(node_gid, 2, 2)) - (mesh->node_jacobian(node_gid, 2, 1)*mesh->node_jacobian(node_gid, 1, 2)) )  //group 1
          - mesh->node_jacobian(node_gid, 0, 1) 
          * ( (mesh->node_jacobian(node_gid, 1, 0)*mesh->node_jacobian(node_gid, 2, 2)) - (mesh->node_jacobian(node_gid, 2, 0)*mesh->node_jacobian(node_gid, 1, 2)) ) // group 2
          + mesh->node_jacobian(node_gid, 0, 2) 
          * ( (mesh->node_jacobian(node_gid, 1, 0)*mesh->node_jacobian(node_gid, 2, 1)) - (mesh->node_jacobian(node_gid, 2, 0)*mesh->node_jacobian(node_gid, 1, 1)) ); // group 3

    }

    for(int cell_gid = 0; cell_gid < mesh->num_cells(); cell_gid++){

        mat_pt->volume(cell_gid) = 0;

        for(int node_lid = 0; node_lid < mesh->num_nodes_in_cell(); node_lid++){

            int node_gid = mesh->nodes_in_cell(cell_gid, node_lid);

            mat_pt->volume(cell_gid) += mesh->node_det_j(node_gid);
        }

        std::cout<< "Volume for cell  "<< cell_gid << " = "<< mat_pt->volume(cell_gid) << std::endl;
    }

    // NOTE: Invert and save J^inverse here!



    for (int node_gid = 0; node_gid < mesh->num_nodes(); node_gid++) {
        
        auto jacobian_view = ViewCArray <real_t> (&mesh->node_jacobian(node_gid, 0, 0), mesh->num_dim(), mesh->num_dim());
        auto jacobian_inverse_view = ViewCArray <real_t> (&mesh->node_jacobian_inv(node_gid, 0, 0), mesh->num_dim(), mesh->num_dim());

        elements::jacobian_inverse_3d(jacobian_inverse_view, jacobian_view);
    }


    // Check J^-1 * J = I


    for (int node_gid = 0; node_gid < mesh->num_nodes(); node_gid++) {
        
        real_t test_array[3][3];

        for(int dim_k = 0; dim_k < mesh->num_dim(); dim_k++){
            for(int dim_i = 0; dim_i < mesh->num_dim(); dim_i++){

                test_array[dim_i][dim_k] = 0.0;
            }
        }

        for(int dim_k = 0; dim_k < mesh->num_dim(); dim_k++){
            for(int dim_j = 0; dim_j < mesh->num_dim(); dim_j++){
                for(int dim_i = 0; dim_i < mesh->num_dim(); dim_i++){
                    test_array[dim_i][dim_k] += mesh->node_jacobian(node_gid, dim_i, dim_j) * mesh->node_jacobian_inv(node_gid, dim_j, dim_k);
                }
            }
        }
        std::cout << std::fixed;
        std::cout<< "J * J^-1 for node  "<< node_gid << " = " << std::endl;
        std::cout  << test_array[0][0] << "     " << test_array[0][1] << "     " << test_array[0][2] << std::endl;
        std::cout  << test_array[1][0] << "     " << test_array[1][1] << "     " << test_array[1][2] << std::endl;
        std::cout  << test_array[2][0] << "     " << test_array[2][1] << "     " << test_array[2][2] << std::endl;
    }


    
} // end subroutine






// Output writers

void Static_Solver::vtk_writer(){
    int graphics_id = simparam->graphics_id;
    int num_dim = simparam->num_dim;

    const int num_scalar_vars = 2;
    const int num_vec_vars = 1;
    const int num_cell_vars = 1;

    const char name[10] = {"Testing"};
    const char scalar_var_names[num_scalar_vars][15] = {
        "fake1", "elem"
    };
    const char vec_var_names[num_vec_vars][15] = {
        "velocity"
    };

    const char cell_var_names[num_vec_vars][15] = {
        "cell_gid"
    };
    //  ---------------------------------------------------------------------------
    //  Setup of file and directoring for exporting
    //  ---------------------------------------------------------------------------

    
    FILE *out[20];   // the output files that are written to
    char filename[128];
    
    struct stat st;
    
    if(stat("vtk",&st) != 0)
        system("mkdir vtk");
    
    if(stat("vtk/data",&st) != 0)
        system("mkdir vtk/data");


    //  ---------------------------------------------------------------------------
    //  Write the Geometry file
    //  ---------------------------------------------------------------------------
    

    
    sprintf(filename, "vtk/data/%s_%05d_%i.vtu", name, graphics_id, 0);
    // filename has the full string
    
    out[0] = fopen(filename, "w");
    
    int num_nodes = mesh->num_nodes();
    int num_cells = mesh->num_cells();


    fprintf(out[0],"<?xml version=\"1.0\"?> \n");
    fprintf(out[0],"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(out[0],"<UnstructuredGrid>\n");
    fprintf(out[0],"<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", num_nodes, num_cells);

    

    //  ---------------------------------------------------------------------------
    //  Write point data
    //  ---------------------------------------------------------------------------


    fprintf(out[0],"<PointData> \n");


    fprintf(out[0],"</PointData> \n");
    fprintf(out[0],"\n");



    //  ---------------------------------------------------------------------------
    //  Write cell data
    //  ---------------------------------------------------------------------------

    fprintf(out[0],"<CellData> \n");

    for(int cell_var = 0; cell_var < num_cell_vars; cell_var++){
        
        fprintf(out[0],"<DataArray type=\"Float64\" Name=\"%s\" Format=\"ascii\">\n", cell_var_names[cell_var]);
        
        for (int cell_gid = 0; cell_gid < num_cells; cell_gid++){
            
            fprintf(out[0],"%f\n",(float) cell_gid);

        } // end for k over cells
    }
    fprintf(out[0],"</DataArray> \n");
    fprintf(out[0],"</CellData> \n");
    fprintf(out[0],"\n");


    //  ---------------------------------------------------------------------------
    //  Write node positions
    //  ---------------------------------------------------------------------------

    real_t min_coord = 0;
    real_t max_coord = 2.0;
    fprintf(out[0],"<Points> \n");

        
    fprintf(out[0],"<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"%i\" format=\"ascii\">\n", num_dim);
    
    for (int node_gid = 0; node_gid < num_nodes; node_gid++){
        
        fprintf(out[0],"%f   %f   %f   \n",mesh->node_coords(node_gid, 0),
                                           mesh->node_coords(node_gid, 1),
                                           mesh->node_coords(node_gid, 2));

    } 

    fprintf(out[0],"</DataArray> \n");
    fprintf(out[0],"</Points> \n");
    fprintf(out[0],"\n");


    //  ---------------------------------------------------------------------------
    //  Write cell type definitions
    //  ---------------------------------------------------------------------------

    fprintf(out[0],"<Cells> \n");
    fprintf(out[0],"\n");


    // Cell connectivity

    fprintf(out[0],"<DataArray type=\"Int64\" Name=\"connectivity\">\n");

    // write nodes in a cell
    for (int cell_gid = 0; cell_gid < num_cells; cell_gid++){
        
        for(int node = 0; node < 8; node++){
            fprintf(out[0],"%i  ", mesh->nodes_in_cell(cell_gid, node));
        }

        fprintf(out[0],"\n");

    } // end for k over cells
    fprintf(out[0],"</DataArray> \n");
    fprintf(out[0],"\n");


    fprintf(out[0],"<DataArray type=\"Int64\" Name=\"offsets\">\n");
    for (int cell_gid = 0; cell_gid < num_cells; cell_gid++){
        fprintf(out[0],"%i  \n", 8*(cell_gid+1));
    } // end for k over cells
    fprintf(out[0],"</DataArray> \n");
    fprintf(out[0],"\n");


    fprintf(out[0],"<DataArray type=\"UInt64\" Name=\"types\">\n");
    for (int cell_gid = 0; cell_gid < num_cells; cell_gid++){
        fprintf(out[0],"%i  \n", 42);
    } // end for k over cells
    fprintf(out[0],"</DataArray> \n");
    fprintf(out[0],"\n");


    fprintf(out[0],"<DataArray type=\"Int64\" Name=\"faces\">\n");
    
    for (int cell_gid = 0; cell_gid < num_cells; cell_gid++){
        fprintf(out[0],"%i  \n", 6);

        for(int patch_lid = 0; patch_lid < 6; patch_lid++){

            fprintf(out[0],"4  ");
            for(int node_lid = 0; node_lid < 4; node_lid++){
                fprintf(out[0],"%i  ", mesh->node_in_patch_in_cell(cell_gid, patch_lid, node_lid));
            }

            fprintf(out[0],"\n");
        }


    } // end for k over cells
    fprintf(out[0],"</DataArray> \n");
    fprintf(out[0],"\n");


    int faceoffsets = 31;
    fprintf(out[0],"<DataArray type=\"Int64\" Name=\"faceoffsets\">\n");
    for (int cell_gid = 0; cell_gid < num_cells; cell_gid++){
        fprintf(out[0],"%i  \n", faceoffsets*(cell_gid+1));
    } // end for k over cells
    
    fprintf(out[0],"</DataArray> \n");
    fprintf(out[0],"\n");


    fprintf(out[0],"</Cells> \n");
    fprintf(out[0],"\n");


    fprintf(out[0],"\n");
    fprintf(out[0],"</Piece> \n");
    fprintf(out[0],"</UnstructuredGrid> \n");
    fprintf(out[0],"</VTKFile> \n");

    fclose(out[0]);
} // end vtk_writer


void Static_Solver::ensight(){
    mat_pt_t *mat_pt = simparam->mat_pt;
    int &graphics_id = simparam->graphics_id;
    real_t *graphics_times = simparam->graphics_times;
    real_t &TIME = simparam->TIME;

    auto convert_vert_list_ord_Ensight = CArray <int> (8);
    convert_vert_list_ord_Ensight(0) = 1;
    convert_vert_list_ord_Ensight(1) = 0;
    convert_vert_list_ord_Ensight(2) = 2;
    convert_vert_list_ord_Ensight(3) = 3;
    convert_vert_list_ord_Ensight(4) = 5;
    convert_vert_list_ord_Ensight(5) = 4;
    convert_vert_list_ord_Ensight(6) = 6;
    convert_vert_list_ord_Ensight(7) = 7;


    const int num_scalar_vars = 4;
    const int num_vec_vars = 1;

    const char name[10] = {"Testing"};
    const char scalar_var_names[num_scalar_vars][15] = {
       "cell_field1", "elem", "elem_id", "cell_field2"
    };
    const char vec_var_names[num_vec_vars][15] = {
        "position"
    };


    int num_nodes = mesh->num_nodes();
    int num_cells = mesh->num_cells();

    // save the cell state to an array for exporting to graphics files
    auto cell_fields = CArray <real_t> (num_cells, num_scalar_vars);

    int cell_cnt = 0;
    int c_in_e = mesh->num_cells_in_elem();
    int elem_val = 1;


    for (int cell_gid=0; cell_gid<num_cells; cell_gid++){
        cell_fields(cell_gid, 0) = mat_pt->field(cell_gid);    
    } // end for k over cells

    int num_elem = mesh->num_elems();

    int num_sub_1d;

    if(mesh->elem_order() == 0){
        num_sub_1d = 1;
    }
    
    else{
        num_sub_1d = mesh->elem_order()*2;
    }



    for (int elem_gid = 0; elem_gid < num_elem; elem_gid++){

        for(int k = 0; k < num_sub_1d; k++){
            for(int j = 0; j < num_sub_1d; j++){
                for(int i = 0; i < num_sub_1d; i++){

                    int cell_index = i + j*num_sub_1d + k*num_sub_1d*num_sub_1d;
                    int cell_mesh_index = cell_index + num_sub_1d*num_sub_1d*num_sub_1d*(elem_gid);
                    cell_fields(cell_mesh_index, 1) = elem_val;

                }

            }
        }
        elem_val *= -1;
    }

    for (int elem_gid = 0; elem_gid < num_elem; elem_gid++){

        for(int k = 0; k < num_sub_1d; k++){
            for(int j = 0; j < num_sub_1d; j++){
                for(int i = 0; i < num_sub_1d; i++){

                    int cell_index = i + j*num_sub_1d + k*num_sub_1d*num_sub_1d;
                    int cell_mesh_index = cell_index + num_sub_1d*num_sub_1d*num_sub_1d*(elem_gid);
                    cell_fields(cell_mesh_index, 2) = elem_gid;

                }

            }
        }
    }

    // Use average temp from each node as cell temp
    for (int cell_gid = 0; cell_gid < num_cells; cell_gid++){

        cell_fields(cell_gid, 3) = mat_pt->field(cell_gid);    

    } // end for k over cells


    // save the vertex vector fields to an array for exporting to graphics files
    auto vec_fields = CArray <real_t> (num_nodes, num_vec_vars, 3);

    for (int node_gid = 0; node_gid < num_nodes; node_gid++){
        
        vec_fields(node_gid, 0, 0) = mesh->node_coords(node_gid, 0); 
        vec_fields(node_gid, 0, 1) = mesh->node_coords(node_gid, 1);
        vec_fields(node_gid, 0, 2) = mesh->node_coords(node_gid, 2);

    } // end for loop over vertices

    //  ---------------------------------------------------------------------------
    //  Setup of file and directoring for exporting
    //  ---------------------------------------------------------------------------

    
    FILE *out[20];   // the output files that are written to
    char filename[128];
    
    struct stat st;
    
    if(stat("ensight",&st) != 0)
        system("mkdir ensight");
    
    if(stat("ensight/data",&st) != 0)
        system("mkdir ensight/data");


    //  ---------------------------------------------------------------------------
    //  Write the Geometry file
    //  ---------------------------------------------------------------------------
    
    
    sprintf(filename, "ensight/data/%s.%05d.geo", name, graphics_id);
    // filename has the full string
    
    out[0] = fopen(filename, "w");
    
    
    fprintf(out[0],"A graphics dump by Cercion \n");
    fprintf(out[0],"%s","EnSight Gold geometry\n");
    fprintf(out[0],"%s","node id assign\n");
    fprintf(out[0],"%s","element id assign\n");
    
    fprintf(out[0],"part\n");
    fprintf(out[0],"%10d\n",1);
    fprintf(out[0],"Mesh\n");
    
    
    // --- vertices ---
    fprintf(out[0],"coordinates\n");
    fprintf(out[0],"%10d\n",mesh->num_nodes());
    
    // write all components of the point coordinates
    for (int node_gid=0; node_gid<num_nodes; node_gid++){
        fprintf(out[0],"%12.5e\n",mesh->node_coords(node_gid, 0));
    }
    
    for (int node_gid=0; node_gid<num_nodes; node_gid++){
        fprintf(out[0],"%12.5e\n",mesh->node_coords(node_gid, 1));
    }
    
    for (int node_gid=0; node_gid<num_nodes; node_gid++){
        fprintf(out[0],"%12.5e\n",mesh->node_coords(node_gid, 2));
    }
    
    // convert_vert_list_ord_Ensight
    // --- cells ---
    fprintf(out[0],"hexa8\n");
    fprintf(out[0],"%10d\n",num_cells);
    int this_index;
    
    // write all global point numbers for this cell
    for (int cell_gid = 0; cell_gid<num_cells; cell_gid++) {

        for (int j=0; j<8; j++){
            this_index = convert_vert_list_ord_Ensight(j);
            fprintf(out[0],"%10d\t",mesh->nodes_in_cell(cell_gid, this_index)+1); // note node_gid starts at 1

        }
        fprintf(out[0],"\n");
    }
  
    fclose(out[0]);
    
 
    // ---------------------------------------------------------------------------
    // Write the Scalar variable files
    // ---------------------------------------------------------------------------
    // ensight_vars = (den, pres,...)
    for (int var = 0; var < num_scalar_vars; var++){
        
        // write a scalar value
        sprintf(filename,"ensight/data/%s.%05d.%s", name, graphics_id, scalar_var_names[var]);

        out[0]=fopen(filename,"w");
        
        fprintf(out[0],"Per_elem scalar values\n");
        fprintf(out[0],"part\n");
        fprintf(out[0],"%10d\n",1);
        
        fprintf(out[0],"hexa8\n");  // e.g., hexa8
        
        for (int cell_gid=0; cell_gid<num_cells; cell_gid++) {
            fprintf(out[0],"%12.5e\n", cell_fields(cell_gid, var));
        }
        
        fclose(out[0]);
        
    } // end for var

    //  ---------------------------------------------------------------------------
    //  Write the Vector variable files
    //  ---------------------------------------------------------------------------
    
    // ensight vector vars = (position, velocity, force)
    for (int var=0; var<num_vec_vars; var++){
        
        sprintf(filename,"ensight/data/%s.%05d.%s", name, graphics_id, vec_var_names[var]);
        
        out[0]=fopen(filename,"w");
        fprintf(out[0],"Per_node vector values\n");
        fprintf(out[0],"part\n");
        fprintf(out[0],"%10d\n",1);
        fprintf(out[0],"coordinates\n");
        
        for (int node_gid=0; node_gid<num_nodes; node_gid++){
            fprintf(out[0],"%12.5e\n",vec_fields(node_gid, var, 0));
        }
        
        for (int node_gid=0; node_gid<num_nodes; node_gid++){
            fprintf(out[0],"%12.5e\n",vec_fields(node_gid, var, 1));
        }
        
        for (int node_gid=0; node_gid<num_nodes; node_gid++){
            fprintf(out[0],"%12.5e\n",vec_fields(node_gid, var, 2));
        }
        
        fclose(out[0]);
    } // end for var

    // ---------------------------------------------------------------------------
    // Write the case file
    // ---------------------------------------------------------------------------
    
    sprintf(filename,"ensight/%s.case",name);
    out[0]=fopen(filename,"w");
    
    fprintf(out[0],"FORMAT\n");
    fprintf(out[0],"type: ensight gold\n");
    fprintf(out[0],"GEOMETRY\n");
    
    sprintf(filename,"model: data/%s.*****.geo\n",name);
    fprintf(out[0],"%s",filename);
    fprintf(out[0],"VARIABLE\n");
    
    for (int var=0; var<num_scalar_vars; var++){
        sprintf(filename,
                "scalar per element: %s data/%s.*****.%s\n",
                scalar_var_names[var], name, scalar_var_names[var]);
        fprintf(out[0],"%s",filename);
    }
    
    for (int var=0; var<num_vec_vars; var++){
        
        sprintf(filename,
                "vector per node: %s data/%s.*****.%s\n",
                vec_var_names[var], name, vec_var_names[var]);
        fprintf(out[0],"%s",filename);
    }
    
    fprintf(out[0],"TIME\n");
    fprintf(out[0],"time set: 1\n");
    fprintf(out[0],"number of steps: %4d\n",graphics_id+1);
    fprintf(out[0],"filename start number: 0\n");
    fprintf(out[0],"filename increment: 1\n");
    fprintf(out[0],"time values: \n");

    
    graphics_times[graphics_id]=TIME;
    
    for (int i=0;i<=graphics_id;i++) {
        fprintf(out[0],"%12.5e\n",graphics_times[i]);
    }
    fclose(out[0]);
    
    
    // ---------------------------------------------------------------------------
    // Done writing the graphics dump
    // ---------------------------------------------------------------------------

    // increment graphics id counter
    graphics_id++;
} // end ensight

/* ----------------------------------------------------------------------
   Initialize global vectors and array maps needed for matrix assembly
------------------------------------------------------------------------- */

void Static_Solver::init_global(){

  int num_dim = simparam->num_dim;
  int num_elems = mesh->num_elems();
  int num_nodes = mesh->num_nodes();
  int nodes_per_elem = mesh->num_nodes_in_elem();
  Stiffness_Matrix_strides = CArray <size_t> (num_nodes*num_dim);
  CArray <int> Graph_Fill = CArray <int> (num_nodes);
  //CArray <int> nodes_in_cell_list_ = mesh->nodes_in_cell_list_;
  CArray <size_t> current_row_nodes_scanned;
  int current_row_n_nodes_scanned;
  int local_node_index, current_column_index;
  int max_stride = 0;

  //allocate right hand side vector of nodal forces
  Nodal_Forces = CArray <real_t> (num_nodes*num_dim);
  
  //allocate stride arrays
  CArray <size_t> Graph_Matrix_strides_initial = CArray <size_t> (num_nodes);
  Graph_Matrix_strides = CArray <size_t> (num_nodes);

  //allocate storage for the sparse stiffness matrix map used in the assembly process
  Global_Stiffness_Matrix_Assembly_Map = CArray <size_t> (num_elems,nodes_per_elem,nodes_per_elem);

  //allocate array used to determine global node repeats in the sparse graph later
  CArray <int> global_nodes_used = CArray <int> (num_nodes);

  /*allocate array that stores which column the global node index occured on for the current row
    when removing repeats*/
  CArray <size_t> column_index = CArray <size_t> (num_nodes);
  
  //initialize stride arrays
  for(int inode = 0; inode < num_nodes; inode++){
      Graph_Matrix_strides_initial(inode) = 0;
      Graph_Matrix_strides(inode) = 0;
      global_nodes_used(inode) = 0;
      column_index(inode) = 0;
      Graph_Fill(inode) = 0;
  }
  
  //count strides for Sparse Pattern Graph with global repeats
  for (int ielem = 0; ielem < num_elems; ielem++)
    for (int lnode = 0; lnode < nodes_per_elem; lnode++){
      local_node_index = mesh->nodes_in_cell_list_(ielem, lnode);
      Graph_Matrix_strides_initial(local_node_index) += nodes_per_elem;
    }
  
  //equate strides for later
  for(int inode = 0; inode < num_nodes; inode++)
    Graph_Matrix_strides(inode) = Graph_Matrix_strides_initial(inode);
  
  //compute maximum stride
  for(int inode = 0; inode < num_nodes; inode++)
    if(Graph_Matrix_strides_initial(inode) > max_stride) max_stride = Graph_Matrix_strides_initial(inode);
  
  //allocate array used in the repeat removal process
  current_row_nodes_scanned = CArray <size_t> (max_stride);

  //allocate sparse graph with node repeats
  RaggedRightArray <size_t> Repeat_Graph_Matrix = RaggedRightArray <size_t> (Graph_Matrix_strides_initial);
  RaggedRightArrayofVectors <size_t> Element_local_indices = RaggedRightArrayofVectors <size_t> (Graph_Matrix_strides_initial,3);
  
  //Fill the initial Graph with repeats
  
  for (int ielem = 0; ielem < num_elems; ielem++)
    for (int lnode = 0; lnode < nodes_per_elem; lnode++){
      local_node_index = mesh->nodes_in_cell_list_(ielem, lnode);
      for (int jnode = 0; jnode < nodes_per_elem; jnode++){
        current_column_index = Graph_Fill(local_node_index)+jnode;
        Repeat_Graph_Matrix(local_node_index, current_column_index) = mesh->nodes_in_cell_list_(ielem,jnode);

        //fill inverse map
        Element_local_indices(local_node_index,current_column_index,0) = ielem;
        Element_local_indices(local_node_index,current_column_index,1) = lnode;
        Element_local_indices(local_node_index,current_column_index,2) = jnode;

        //fill forward map
        Global_Stiffness_Matrix_Assembly_Map(ielem,lnode,jnode) = current_column_index;
      }
      Graph_Fill(local_node_index) += nodes_per_elem;
    }
  
  //debug statement
  //std::cout << "started run" << std::endl;

  //remove repeats from the inital graph setup
  int current_node, current_element_index, element_row_index, element_column_index, current_stride;
  for (int inode = 0; inode < num_nodes; inode++){
    current_row_n_nodes_scanned = 0;
    for (int istride = 0; istride < Graph_Matrix_strides(inode); istride++){
      current_node = Repeat_Graph_Matrix(inode,istride);
      if(global_nodes_used(current_node)){
        //set forward map index to the index where this global node was first found
        current_element_index = Element_local_indices(inode,istride,0);
        element_row_index = Element_local_indices(inode,istride,1);
        element_column_index = Element_local_indices(inode,istride,2);

        Global_Stiffness_Matrix_Assembly_Map(current_element_index,element_row_index, element_column_index) 
            = column_index(current_node);   

        
        //swap current node with the end of the current row and shorten the stride of the row
        //first swap information about the inverse and forward maps

        current_stride = Graph_Matrix_strides(inode);
        if(istride!=current_stride-1){
        Element_local_indices(inode,istride,0) = Element_local_indices(inode,current_stride-1,0);
        Element_local_indices(inode,istride,1) = Element_local_indices(inode,current_stride-1,1);
        Element_local_indices(inode,istride,2) = Element_local_indices(inode,current_stride-1,2);
        current_element_index = Element_local_indices(inode,istride,0);
        element_row_index = Element_local_indices(inode,istride,1);
        element_column_index = Element_local_indices(inode,istride,2);

        Global_Stiffness_Matrix_Assembly_Map(current_element_index,element_row_index, element_column_index) 
            = istride;

        //now that the element map information has been copied, copy the global node index and delete the last index

        Repeat_Graph_Matrix(inode,istride) = Repeat_Graph_Matrix(inode,current_stride-1);
        }
        istride--;
        Graph_Matrix_strides(inode)--;
      }
      else{
        /*this node hasn't shown up in the row before; add it to the list of nodes
          that have been scanned uniquely. Use this list to reset the flag array
          afterwards without having to loop over all the nodes in the system*/
        global_nodes_used(current_node) = 1;
        column_index(current_node) = istride;
        current_row_nodes_scanned(current_row_n_nodes_scanned) = current_node;
        current_row_n_nodes_scanned++;
      }
    }
    //reset nodes used list for the next row of the sparse list
    for(int node_reset = 0; node_reset < current_row_n_nodes_scanned; node_reset++)
      global_nodes_used(current_row_nodes_scanned(node_reset)) = 0;

  }
  
  //copy reduced content to non_repeat storage
  Graph_Matrix = RaggedRightArray <size_t> (Graph_Matrix_strides);
  for(int inode = 0; inode < num_nodes; inode++)
    for(int istride = 0; istride < Graph_Matrix_strides(inode); istride++)
      Graph_Matrix(inode,istride) = Repeat_Graph_Matrix(inode,istride);

  //deallocate repeat matrix
  
  /*At this stage the sparse graph should have unique global indices on each row.
    The constructed Assembly map (to the global sparse matrix)
    is used to loop over each element's local stiffness matrix in the assembly process.*/
  
  //expand strides for stiffness matrix by multipling by dim
  for(int inode = 0; inode < num_nodes; inode++){
    for (int idim = 0; idim < num_dim; idim++)
    Stiffness_Matrix_strides(num_dim*inode + idim) = num_dim*Graph_Matrix_strides(inode);
  }

  Stiffness_Matrix = RaggedRightArray <real_t> (Stiffness_Matrix_strides);
  DOF_Graph_Matrix = RaggedRightArray <size_t> (Stiffness_Matrix_strides);

  //initialize stiffness Matrix entries to 0
  for (int idof = 0; idof < num_dim*num_nodes; idof++)
    for (int istride = 0; istride < Stiffness_Matrix_strides(idof); istride++){
      Stiffness_Matrix(idof,istride) = 0;
      DOF_Graph_Matrix(idof,istride) = Graph_Matrix(idof/num_dim,istride/num_dim)*num_dim + istride%num_dim;
    }
  
  /*
  //debug print nodal positions and indices
  std::cout << " ------------NODAL POSITIONS--------------"<<std::endl;
  for (int inode = 0; inode < num_nodes; inode++){
      std::cout << "node: " << inode + 1 << " { ";
    for (int istride = 0; istride < num_dim; istride++){
        std::cout << mesh->node_coords(inode,istride) << " , ";
    }
    std::cout << " }"<< std::endl;
  }
  //debug print element edof
  
  std::cout << " ------------ELEMENT EDOF--------------"<<std::endl;

  for (int ielem = 0; ielem < num_elems; ielem++){
    std::cout << "elem:  " << ielem+1 << std::endl;
    for (int lnode = 0; lnode < nodes_per_elem; lnode++){
        std::cout << "{ ";
          std::cout << lnode+1 << " = " << mesh->nodes_in_cell_list_(ielem,lnode) + 1 << " ";
        
        std::cout << " }"<< std::endl;
    }
    std::cout << std::endl;
  }
  

  //debug section; print stiffness matrix graph and per element map
  std::cout << " ------------SPARSE GRAPH MATRIX--------------"<<std::endl;
  for (int inode = 0; inode < num_nodes; inode++){
      std::cout << "row: " << inode + 1 << " { ";
    for (int istride = 0; istride < Graph_Matrix_strides(inode); istride++){
        std::cout << istride + 1 << " = " << Repeat_Graph_Matrix(inode,istride) + 1 << " , " ;
    }
    std::cout << " }"<< std::endl;
  }

  std::cout << " ------------ELEMENT ASSEMBLY MAP--------------"<<std::endl;

  for (int ielem = 0; ielem < num_elems; ielem++){
    std::cout << "elem:  " << ielem+1 << std::endl;
    for (int lnode = 0; lnode < nodes_per_elem; lnode++){
        std::cout << "{ "<< std::endl;
        for (int jnode = 0; jnode < nodes_per_elem; jnode++){
          std::cout <<"(" << lnode+1 << "," << jnode+1 << "," << mesh->nodes_in_cell(ielem,lnode)+1 << ")"<< " = " << Global_Stiffness_Matrix_Assembly_Map(ielem,lnode, jnode) + 1 << " ";
        }
        std::cout << " }"<< std::endl;
    }
    std::cout << std::endl;
  }
  */
  
}

/* ----------------------------------------------------------------------
   Assemble the Sparse Stiffness Matrix and force vector
------------------------------------------------------------------------- */

void Static_Solver::assemble(){
  int num_dim = simparam->num_dim;
  int num_elems = mesh->num_elems();
  int num_nodes = mesh->num_nodes();
  int nodes_per_elem = mesh->num_nodes_in_elem();
  int current_row_n_nodes_scanned;
  int local_node_index, current_row, current_column;
  int max_stride = 0;
  CArray <real_t> Local_Stiffness_Matrix = CArray <real_t> (num_dim*nodes_per_elem,num_dim*nodes_per_elem);

  //assemble the global stiffness matrix
  for (int ielem = 0; ielem < num_elems; ielem++){
    //construct local stiffness matrix for this element
    local_matrix_multiply(ielem, Local_Stiffness_Matrix);
    //assign entries of this local matrix to the sparse global matrix storage;
    for (int inode = 0; inode < nodes_per_elem; inode++){
      current_row = num_dim*mesh->nodes_in_cell_list_(ielem,inode);
      for(int jnode = 0; jnode < nodes_per_elem; jnode++){
        
        current_column = num_dim*Global_Stiffness_Matrix_Assembly_Map(ielem,inode,jnode);
        for (int idim = 0; idim < num_dim; idim++){
          for (int jdim = 0; jdim < num_dim; jdim++){

            //debug print
            //if(current_row + idim==15&&current_column + jdim==4)
            //std::cout << " Local stiffness matrix contribution for row " << current_row + idim +1 << " and column " << current_column + jdim + 1 << " : " <<
            //Local_Stiffness_Matrix(num_dim*inode + idim,num_dim*jnode + jdim) << " from " << ielem +1 << " i: " << num_dim*inode+idim+1 << " j: " << num_dim*jnode + jdim +1 << std::endl << std::endl;
            //end debug

            Stiffness_Matrix(current_row + idim, current_column + jdim) += Local_Stiffness_Matrix(num_dim*inode + idim,num_dim*jnode + jdim);
          }
        }
      }
    }
  }

  //debug print of stiffness matrix
  //debug section; print stiffness matrix graph and per element map
  /*
  std::cout << " ------------SPARSE STIFFNESS MATRIX--------------"<<std::endl;
  for (int idof = 0; idof < num_nodes*num_dim; idof++){
      std::cout << "row: " << idof + 1 << " { ";
    for (int istride = 0; istride < Stiffness_Matrix_strides(idof); istride++){
        std::cout << istride + 1 << " = " << Stiffness_Matrix(idof,istride) << " , " ;
    }
    std::cout << " }"<< std::endl;
  }
  */
  //force vector construction
  
  //initialize
  for(int i=0; i < num_dim*num_nodes; i++)
    Nodal_Forces(i) = 0;
  
  //Tag nodes for Boundary conditions such as displacements
  Displacement_Boundary_Conditions();

  //Construct applied nodal force vector with quadrature
  Force_Vector_Construct();

}

/* ----------------------------------------------------------------------
   Retrieve material properties associated with a finite element
------------------------------------------------------------------------- */

void Static_Solver::Element_Material_Properties(size_t ielem, real_t &Element_Modulus, real_t &Poisson_Ratio){

  Element_Modulus = simparam->Elastic_Modulus;
  Poisson_Ratio = simparam->Poisson_Ratio;

}

/* ----------------------------------------------------------------------
   Construct the local stiffness matrix
------------------------------------------------------------------------- */

void Static_Solver::local_matrix(int ielem, CArray <real_t> &Local_Matrix){
  int num_dim = simparam->num_dim;
  int num_elems = mesh->num_elems();
  int num_nodes = mesh->num_nodes();
  int nodes_per_elem = mesh->num_nodes_in_elem();
  int num_gauss_points = simparam->num_gauss_points;
  int z_quad,y_quad,x_quad, direct_product_count;

  direct_product_count = std::pow(num_gauss_points,num_dim);
  real_t Elastic_Constant, Shear_Term, Pressure_Term, matrix_term;
  real_t matrix_subterm1, matrix_subterm2, matrix_subterm3, Jacobian;
  real_t Element_Modulus, Poisson_Ratio;
  CArray<real_t> legendre_nodes_1D(num_gauss_points);
  CArray<real_t> legendre_weights_1D(num_gauss_points);
  real_t pointer_quad_coordinate[num_dim];
  real_t pointer_quad_coordinate_weight[num_dim];
  real_t pointer_interpolated_point[num_dim];
  real_t pointer_JT_row1[num_dim];
  real_t pointer_JT_row2[num_dim];
  real_t pointer_JT_row3[num_dim];
  ViewCArray<real_t> quad_coordinate(pointer_quad_coordinate,num_dim);
  ViewCArray<real_t> quad_coordinate_weight(pointer_quad_coordinate_weight,num_dim);
  ViewCArray<real_t> interpolated_point(pointer_interpolated_point,num_dim);
  ViewCArray<real_t> JT_row1(pointer_JT_row1,num_dim);
  ViewCArray<real_t> JT_row2(pointer_JT_row2,num_dim);
  ViewCArray<real_t> JT_row3(pointer_JT_row3,num_dim);

  real_t pointer_basis_values[elem->num_basis()];
  real_t pointer_basis_derivative_s1[elem->num_basis()];
  real_t pointer_basis_derivative_s2[elem->num_basis()];
  real_t pointer_basis_derivative_s3[elem->num_basis()];
  ViewCArray<real_t> basis_values(pointer_basis_values,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s1(pointer_basis_derivative_s1,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s2(pointer_basis_derivative_s2,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s3(pointer_basis_derivative_s3,elem->num_basis());
  CArray<real_t> nodal_positions(elem->num_basis(),num_dim);

  //initialize weights
  elements::legendre_nodes_1D(legendre_nodes_1D,num_gauss_points);
  elements::legendre_weights_1D(legendre_weights_1D,num_gauss_points);

  //acquire set of nodes for this local element
  for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
    nodal_positions(node_loop,0) = mesh->node_coords((mesh->nodes_in_cell_list_(ielem, node_loop)),0);
    nodal_positions(node_loop,1) = mesh->node_coords((mesh->nodes_in_cell_list_(ielem, node_loop)),1);
    nodal_positions(node_loop,2) = mesh->node_coords((mesh->nodes_in_cell_list_(ielem, node_loop)),2);
  }

  //look up element material properties
  Element_Material_Properties((size_t) ielem,Element_Modulus,Poisson_Ratio);
  Elastic_Constant = Element_Modulus/(1+Poisson_Ratio)/(1-2*Poisson_Ratio);
  Shear_Term = 0.5-Poisson_Ratio;
  Pressure_Term = 1 - Poisson_Ratio;

  //initialize local stiffness matrix storage
  for(int ifill=0; ifill < num_dim*nodes_per_elem; ifill++)
      for(int jfill=0; jfill < num_dim*nodes_per_elem; jfill++)
      Local_Matrix(ifill,jfill) = 0;

  //loop over quadrature points
  for(int iquad=0; iquad < direct_product_count; iquad++){

    //set current quadrature point
    if(num_dim==3) z_quad = iquad/(num_gauss_points*num_gauss_points);
    y_quad = (iquad % (num_gauss_points*num_gauss_points))/num_gauss_points;
    x_quad = iquad % num_gauss_points;
    quad_coordinate(0) = legendre_nodes_1D(x_quad);
    quad_coordinate(1) = legendre_nodes_1D(y_quad);
    if(num_dim==3)
    quad_coordinate(2) = legendre_nodes_1D(z_quad);

    //set current quadrature weight
    quad_coordinate_weight(0) = legendre_weights_1D(x_quad);
    quad_coordinate_weight(1) = legendre_weights_1D(y_quad);
    if(num_dim==3)
    quad_coordinate_weight(2) = legendre_weights_1D(z_quad);

    //compute shape functions at this point for the element type
    elem->basis(basis_values,quad_coordinate);

    //compute all the necessary coordinates and derivatives at this point

    //compute shape function derivatives
    elem->partial_xi_basis(basis_derivative_s1,quad_coordinate);
    elem->partial_eta_basis(basis_derivative_s2,quad_coordinate);
    elem->partial_mu_basis(basis_derivative_s3,quad_coordinate);

    //compute derivatives of x,y,z w.r.t the s,t,w isoparametric space needed by JT (Transpose of the Jacobian)
    //derivative of x,y,z w.r.t s
    JT_row1(0) = 0;
    JT_row1(1) = 0;
    JT_row1(2) = 0;
    for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
      JT_row1(0) += nodal_positions(node_loop,0)*basis_derivative_s1(node_loop);
      JT_row1(1) += nodal_positions(node_loop,1)*basis_derivative_s1(node_loop);
      JT_row1(2) += nodal_positions(node_loop,2)*basis_derivative_s1(node_loop);
    }

    //derivative of x,y,z w.r.t t
    JT_row2(0) = 0;
    JT_row2(1) = 0;
    JT_row2(2) = 0;
    for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
      JT_row2(0) += nodal_positions(node_loop,0)*basis_derivative_s2(node_loop);
      JT_row2(1) += nodal_positions(node_loop,1)*basis_derivative_s2(node_loop);
      JT_row2(2) += nodal_positions(node_loop,2)*basis_derivative_s2(node_loop);
    }

    //derivative of x,y,z w.r.t w
    JT_row3(0) = 0;
    JT_row3(1) = 0;
    JT_row3(2) = 0;
    for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
      JT_row3(0) += nodal_positions(node_loop,0)*basis_derivative_s3(node_loop);
      JT_row3(1) += nodal_positions(node_loop,1)*basis_derivative_s3(node_loop);
      JT_row3(2) += nodal_positions(node_loop,2)*basis_derivative_s3(node_loop);
    }

    //compute the determinant of the Jacobian
    Jacobian = JT_row1(0)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
               JT_row1(1)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
               JT_row1(2)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1));
    if(Jacobian<0) Jacobian = -Jacobian;

    //compute the contributions of this quadrature point to all the matrix elements
    int index_x,index_y,basis_index_x,basis_index_y,swap1,swap2;
    for(int ifill=0; ifill < num_dim*nodes_per_elem; ifill++)
      for(int jfill=0; jfill < num_dim*nodes_per_elem; jfill++){
        index_x = ifill%num_dim;
        index_y = jfill%num_dim;
        basis_index_x = ifill/num_dim;
        basis_index_y = jfill/num_dim;

        //compute stiffness matrix terms involving derivatives of the shape function and cofactor determinants from cramers rule
        if(index_x==0&&index_y==0){
          matrix_subterm1 = Pressure_Term*(basis_derivative_s1(basis_index_x)*
          (JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(basis_index_x)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(basis_index_x)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)))*
          (basis_derivative_s1(basis_index_y)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(basis_index_y)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(basis_index_y)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));

          matrix_subterm2 = Shear_Term*(-basis_derivative_s1(basis_index_x)*
          (JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(basis_index_x)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(basis_index_x)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)))*
          (-basis_derivative_s1(basis_index_y)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(basis_index_y)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(basis_index_y)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));

          matrix_subterm3 = Shear_Term*(basis_derivative_s1(basis_index_x)*
          (JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(basis_index_x)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(basis_index_x)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)))*
          (basis_derivative_s1(basis_index_y)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(basis_index_y)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(basis_index_y)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));

          matrix_term = matrix_subterm1 + matrix_subterm2 + matrix_subterm3;

           //debug print block
           /*
          if(ielem==0&&((jfill==3&&ifill==3))){
           std::cout << " ------------quadrature point "<< iquad + 1 <<"--------------"<<std::endl;
           std::cout << "row: " << ifill + 1 << " { ";
           std::cout << jfill + 1 << " = " << matrix_subterm1 << " , " << " bi:JT11 "<<JT_row1(0) << " bi:JT12 " <<  JT_row1(1) << " bi:JT13 " << JT_row1(2)
           <<  " bi:JT21 "<<JT_row2(0) << " bi:JT22 " <<  JT_row2(1) << " bi:JT23 " << JT_row2(2) <<  " bi:JT31 "<<JT_row3(0) << " bi:JT32 " <<  JT_row3(1) << " bi:JT33 " << JT_row3(2);
           std::cout << " }"<< std::endl;
          }
          
          if(ielem==0&&((jfill==3&&ifill==3))){
           std::cout << " ------------quadrature point "<< iquad + 1 <<"--------------"<<std::endl;
           std::cout << "row: " << ifill + 1 << " { ";
           std::cout << jfill + 1 << " = " << matrix_term*Elastic_Constant/Jacobian << " , " << " basis index x s1 "<< basis_derivative_s1(basis_index_x) << " quad x " <<  quad_coordinate(0)
           <<  " quad y "<< quad_coordinate(1) << " quad z " <<  quad_coordinate(2) << "Force Vector " << Local_Matrix(3,3);
           std::cout << " }"<< std::endl;
          }
          */
          
        }

        if(index_x==1&&index_y==1){
          
          matrix_subterm1 = Pressure_Term*(-basis_derivative_s1(basis_index_x)*
          (JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(basis_index_x)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(basis_index_x)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)))*
          (-basis_derivative_s1(basis_index_y)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(basis_index_y)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(basis_index_y)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));

          matrix_subterm2 = Shear_Term*(basis_derivative_s1(basis_index_x)*
          (JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(basis_index_x)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(basis_index_x)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)))*
          (basis_derivative_s1(basis_index_y)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(basis_index_y)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(basis_index_y)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));

          matrix_subterm3 = Shear_Term*(basis_derivative_s1(basis_index_x)*
          (JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(basis_index_x)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(basis_index_x)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)))*
          (basis_derivative_s1(basis_index_y)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(basis_index_y)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(basis_index_y)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));

          matrix_term = matrix_subterm1 + matrix_subterm2 + matrix_subterm3;
        }

        if(index_x==2&&index_y==2){
          
          matrix_subterm1 = Pressure_Term*(basis_derivative_s1(basis_index_x)*
          (JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(basis_index_x)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(basis_index_x)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)))*
          (basis_derivative_s1(basis_index_y)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(basis_index_y)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(basis_index_y)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));

          matrix_subterm2 = Shear_Term*(basis_derivative_s1(basis_index_x)*
          (JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(basis_index_x)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(basis_index_x)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)))*
          (basis_derivative_s1(basis_index_y)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(basis_index_y)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(basis_index_y)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));

          matrix_subterm3 = Shear_Term*(-basis_derivative_s1(basis_index_x)*
          (JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(basis_index_x)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(basis_index_x)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)))*
          (-basis_derivative_s1(basis_index_y)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(basis_index_y)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(basis_index_y)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));

          matrix_term = matrix_subterm1 + matrix_subterm2 + matrix_subterm3;
        }

        if((index_x==0&&index_y==1)||(index_x==1&&index_y==0)){
          if(index_x==1&&index_y==0){
            swap1 = basis_index_x;
            swap2 = basis_index_y;
          }
          else{
            swap1 = basis_index_y;
            swap2 = basis_index_x;
          }
          
          matrix_subterm1 = Poisson_Ratio*(-basis_derivative_s1(swap1)*
          (JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(swap1)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(swap1)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)))*
          (basis_derivative_s1(swap2)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(swap2)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(swap2)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
          
          
          matrix_subterm2 = Shear_Term*(basis_derivative_s1(swap1)*
          (JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(swap1)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(swap1)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)))*
          (-basis_derivative_s1(swap2)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(swap2)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(swap2)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));

          matrix_term = matrix_subterm1 + matrix_subterm2;

          /* debug print block
          if(iquad==0&&((jfill==4&&ifill==0)||(jfill==0&&ifill==4))){
           std::cout << " ------------LOCAL STIFFNESS MATRIX "<< ielem + 1 <<"--------------"<<std::endl;
           std::cout << "row: " << ifill + 1 << " { ";
           std::cout << jfill + 1 << " = " << matrix_subterm2 << " , " << " bi:JT11 "<<JT_row1(0) << " bi:JT12 " <<  JT_row1(1) << " bi:JT13 " << JT_row1(2)
           <<  " bi:JT21 "<<JT_row2(0) << " bi:JT22 " <<  JT_row2(1) << " bi:JT23 " << JT_row2(2) <<  " bi:JT31 "<<JT_row3(0) << " bi:JT32 " <<  JT_row3(1) << " bi:JT33 " << JT_row3(2);
           std::cout << " }"<< std::endl;
          }
          */
        }

        if((index_x==0&&index_y==2)||(index_x==2&&index_y==0)){
          if(index_x==2&index_y==0){
            swap1 = basis_index_x;
            swap2 = basis_index_y;
          }
          else{
            swap1 = basis_index_y;
            swap2 = basis_index_x;
          }
          matrix_subterm1 = Poisson_Ratio*(basis_derivative_s1(swap1)*
          (JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(swap1)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(swap1)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)))*
          (basis_derivative_s1(swap2)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(swap2)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(swap2)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));

          matrix_subterm2 = Shear_Term*(basis_derivative_s1(swap1)*
          (JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(swap1)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(swap1)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)))*
          (basis_derivative_s1(swap2)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(swap2)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(swap2)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));

          matrix_term = matrix_subterm1 + matrix_subterm2;
        }

        if((index_x==1&&index_y==2)||(index_x==2&&index_y==1)){
          if(index_x==2&&index_y==1){
            swap1 = basis_index_x;
            swap2 = basis_index_y;
          }
          else{
            swap1 = basis_index_y;
            swap2 = basis_index_x;
          }
          matrix_subterm1 = Poisson_Ratio*(basis_derivative_s1(swap1)*
          (JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(swap1)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(swap1)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)))*
          (-basis_derivative_s1(swap2)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(swap2)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(swap2)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));

          matrix_subterm2 = Shear_Term*(-basis_derivative_s1(swap1)*
          (JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(swap1)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(swap1)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)))*
          (basis_derivative_s1(swap2)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(swap2)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(swap2)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));

          matrix_term = matrix_subterm1 + matrix_subterm2;
        }
        
        Local_Matrix(ifill,jfill) += Elastic_Constant*quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2)*matrix_term/Jacobian;
      }
      
    }

    //debug print of local stiffness matrix
      /*
      std::cout << " ------------LOCAL STIFFNESS MATRIX "<< ielem + 1 <<"--------------"<<std::endl;
      for (int idof = 0; idof < num_dim*nodes_per_elem; idof++){
        std::cout << "row: " << idof + 1 << " { ";
        for (int istride = 0; istride < num_dim*nodes_per_elem; istride++){
          std::cout << istride + 1 << " = " << Local_Matrix(idof,istride) << " , " ;
        }
        std::cout << " }"<< std::endl;
        }
      */
}

/* ----------------------------------------------------------------------
   Construct the local stiffness matrix
------------------------------------------------------------------------- */

void Static_Solver::local_matrix_multiply(int ielem, CArray <real_t> &Local_Matrix){
  int num_dim = simparam->num_dim;
  int num_elems = mesh->num_elems();
  int num_nodes = mesh->num_nodes();
  int nodes_per_elem = elem->num_basis();
  int num_gauss_points = simparam->num_gauss_points;
  int z_quad,y_quad,x_quad, direct_product_count;

  direct_product_count = std::pow(num_gauss_points,num_dim);
  real_t Elastic_Constant, Shear_Term, Pressure_Term, matrix_term;
  real_t matrix_subterm1, matrix_subterm2, matrix_subterm3, Jacobian;
  real_t Element_Modulus, Poisson_Ratio;
  CArray<real_t> legendre_nodes_1D(num_gauss_points);
  CArray<real_t> legendre_weights_1D(num_gauss_points);
  real_t pointer_quad_coordinate[num_dim];
  real_t pointer_quad_coordinate_weight[num_dim];
  real_t pointer_interpolated_point[num_dim];
  real_t pointer_JT_row1[num_dim];
  real_t pointer_JT_row2[num_dim];
  real_t pointer_JT_row3[num_dim];
  ViewCArray<real_t> quad_coordinate(pointer_quad_coordinate,num_dim);
  ViewCArray<real_t> quad_coordinate_weight(pointer_quad_coordinate_weight,num_dim);
  ViewCArray<real_t> interpolated_point(pointer_interpolated_point,num_dim);
  ViewCArray<real_t> JT_row1(pointer_JT_row1,num_dim);
  ViewCArray<real_t> JT_row2(pointer_JT_row2,num_dim);
  ViewCArray<real_t> JT_row3(pointer_JT_row3,num_dim);

  real_t pointer_basis_values[elem->num_basis()];
  real_t pointer_basis_derivative_s1[elem->num_basis()];
  real_t pointer_basis_derivative_s2[elem->num_basis()];
  real_t pointer_basis_derivative_s3[elem->num_basis()];
  ViewCArray<real_t> basis_values(pointer_basis_values,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s1(pointer_basis_derivative_s1,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s2(pointer_basis_derivative_s2,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s3(pointer_basis_derivative_s3,elem->num_basis());
  CArray<real_t> nodal_positions(elem->num_basis(),num_dim);
  size_t Brows;
  if(num_dim==2) Brows = 3;
  if(num_dim==3) Brows = 6;
  CArray<real_t> B_matrix_contribution(Brows,num_dim*elem->num_basis());
  CArray<real_t> B_matrix(Brows,num_dim*elem->num_basis());
  CArray<real_t> CB_matrix_contribution(Brows,num_dim*elem->num_basis());
  CArray<real_t> CB_matrix(Brows,num_dim*elem->num_basis());
  CArray<real_t> C_matrix(Brows,Brows);

  //initialize weights
  elements::legendre_nodes_1D(legendre_nodes_1D,num_gauss_points);
  elements::legendre_weights_1D(legendre_weights_1D,num_gauss_points);

  //acquire set of nodes for this local element
  for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
    nodal_positions(node_loop,0) = mesh->node_coords((mesh->nodes_in_cell_list_(ielem, node_loop)),0);
    nodal_positions(node_loop,1) = mesh->node_coords((mesh->nodes_in_cell_list_(ielem, node_loop)),1);
    nodal_positions(node_loop,2) = mesh->node_coords((mesh->nodes_in_cell_list_(ielem, node_loop)),2);
  }

  //look up element material properties
  Element_Material_Properties((size_t) ielem,Element_Modulus,Poisson_Ratio);
  Elastic_Constant = Element_Modulus/(1+Poisson_Ratio)/(1-2*Poisson_Ratio);
  Shear_Term = 0.5-Poisson_Ratio;
  Pressure_Term = 1 - Poisson_Ratio;

  //initialize C matrix
  for(int irow = 0; irow < Brows; irow++)
    for(int icol = 0; icol < Brows; icol++)
      C_matrix(irow,icol) = 0;

  //compute Elastic (C) matrix
  if(num_dim==2){
  C_matrix(0,0) = Pressure_Term;
  C_matrix(1,1) = Pressure_Term;
  C_matrix(0,1) = Poisson_Ratio;
  C_matrix(1,0) = Poisson_Ratio;
  C_matrix(2,2) = Shear_Term;
  }
  if(num_dim==3){
  C_matrix(0,0) = Pressure_Term;
  C_matrix(1,1) = Pressure_Term;
  C_matrix(2,2) = Pressure_Term;
  C_matrix(0,1) = Poisson_Ratio;
  C_matrix(0,2) = Poisson_Ratio;
  C_matrix(1,0) = Poisson_Ratio;
  C_matrix(1,2) = Poisson_Ratio;
  C_matrix(2,0) = Poisson_Ratio;
  C_matrix(2,1) = Poisson_Ratio;
  C_matrix(3,3) = Shear_Term;
  C_matrix(4,4) = Shear_Term;
  C_matrix(5,5) = Shear_Term;
  }
  
  /*
  //debug print of elasticity matrix
  std::cout << " ------------ELASTICITY MATRIX "<< ielem + 1 <<"--------------"<<std::endl;
  for (int idof = 0; idof < Brows; idof++){
    std::cout << "row: " << idof + 1 << " { ";
    for (int istride = 0; istride < Brows; istride++){
      std::cout << istride + 1 << " = " << C_matrix(idof,istride) << " , " ;
    }
    std::cout << " }"<< std::endl;
  }
  //end debug block
  */

  //initialize local stiffness matrix storage
  for(int ifill=0; ifill < num_dim*nodes_per_elem; ifill++)
      for(int jfill=0; jfill < num_dim*nodes_per_elem; jfill++)
      Local_Matrix(ifill,jfill) = 0;

  //B matrix initialization
  for(int irow=0; irow < Brows; irow++)
      for(int icol=0; icol < num_dim*nodes_per_elem; icol++){
        CB_matrix(irow,icol) = B_matrix(irow,icol) = 0;
      }

  //loop over quadrature points
  for(int iquad=0; iquad < direct_product_count; iquad++){

    //set current quadrature point
    if(num_dim==3) z_quad = iquad/(num_gauss_points*num_gauss_points);
    y_quad = (iquad % (num_gauss_points*num_gauss_points))/num_gauss_points;
    x_quad = iquad % num_gauss_points;
    quad_coordinate(0) = legendre_nodes_1D(x_quad);
    quad_coordinate(1) = legendre_nodes_1D(y_quad);
    if(num_dim==3)
    quad_coordinate(2) = legendre_nodes_1D(z_quad);

    //set current quadrature weight
    quad_coordinate_weight(0) = legendre_weights_1D(x_quad);
    quad_coordinate_weight(1) = legendre_weights_1D(y_quad);
    if(num_dim==3)
    quad_coordinate_weight(2) = legendre_weights_1D(z_quad);

    //compute shape functions at this point for the element type
    elem->basis(basis_values,quad_coordinate);

    //compute all the necessary coordinates and derivatives at this point

    //compute shape function derivatives
    elem->partial_xi_basis(basis_derivative_s1,quad_coordinate);
    elem->partial_eta_basis(basis_derivative_s2,quad_coordinate);
    elem->partial_mu_basis(basis_derivative_s3,quad_coordinate);

    //compute derivatives of x,y,z w.r.t the s,t,w isoparametric space needed by JT (Transpose of the Jacobian)
    //derivative of x,y,z w.r.t s
    JT_row1(0) = 0;
    JT_row1(1) = 0;
    JT_row1(2) = 0;
    for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
      JT_row1(0) += nodal_positions(node_loop,0)*basis_derivative_s1(node_loop);
      JT_row1(1) += nodal_positions(node_loop,1)*basis_derivative_s1(node_loop);
      JT_row1(2) += nodal_positions(node_loop,2)*basis_derivative_s1(node_loop);
    }

    //derivative of x,y,z w.r.t t
    JT_row2(0) = 0;
    JT_row2(1) = 0;
    JT_row2(2) = 0;
    for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
      JT_row2(0) += nodal_positions(node_loop,0)*basis_derivative_s2(node_loop);
      JT_row2(1) += nodal_positions(node_loop,1)*basis_derivative_s2(node_loop);
      JT_row2(2) += nodal_positions(node_loop,2)*basis_derivative_s2(node_loop);
    }

    //derivative of x,y,z w.r.t w
    JT_row3(0) = 0;
    JT_row3(1) = 0;
    JT_row3(2) = 0;
    for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
      JT_row3(0) += nodal_positions(node_loop,0)*basis_derivative_s3(node_loop);
      JT_row3(1) += nodal_positions(node_loop,1)*basis_derivative_s3(node_loop);
      JT_row3(2) += nodal_positions(node_loop,2)*basis_derivative_s3(node_loop);
    }

    //compute the determinant of the Jacobian
    Jacobian = JT_row1(0)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
               JT_row1(1)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
               JT_row1(2)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1));
    if(Jacobian<0) Jacobian = -Jacobian;

    //compute the contributions of this quadrature point to the B matrix
    if(num_dim==2)
    for(int ishape=0; ishape < nodes_per_elem; ishape++){
      B_matrix_contribution(0,ishape*num_dim) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(1,ishape*num_dim) = 0;
      B_matrix_contribution(2,ishape*num_dim) = 0;
      B_matrix_contribution(3,ishape*num_dim) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
      B_matrix_contribution(4,ishape*num_dim) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(5,ishape*num_dim) = 0;
      B_matrix_contribution(0,ishape*num_dim+1) = 0;
      B_matrix_contribution(1,ishape*num_dim+1) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
      B_matrix_contribution(2,ishape*num_dim+1) = 0;
      B_matrix_contribution(3,ishape*num_dim+1) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(4,ishape*num_dim+1) = 0;
      B_matrix_contribution(5,ishape*num_dim+1) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(0,ishape*num_dim+2) = 0;
      B_matrix_contribution(1,ishape*num_dim+2) = 0;
      B_matrix_contribution(2,ishape*num_dim+2) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(3,ishape*num_dim+2) = 0;
      B_matrix_contribution(4,ishape*num_dim+2) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(5,ishape*num_dim+2) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
    }
    if(num_dim==3)
    for(int ishape=0; ishape < nodes_per_elem; ishape++){
      B_matrix_contribution(0,ishape*num_dim) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(1,ishape*num_dim) = 0;
      B_matrix_contribution(2,ishape*num_dim) = 0;
      B_matrix_contribution(3,ishape*num_dim) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
      B_matrix_contribution(4,ishape*num_dim) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(5,ishape*num_dim) = 0;
      B_matrix_contribution(0,ishape*num_dim+1) = 0;
      B_matrix_contribution(1,ishape*num_dim+1) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
      B_matrix_contribution(2,ishape*num_dim+1) = 0;
      B_matrix_contribution(3,ishape*num_dim+1) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(4,ishape*num_dim+1) = 0;
      B_matrix_contribution(5,ishape*num_dim+1) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(0,ishape*num_dim+2) = 0;
      B_matrix_contribution(1,ishape*num_dim+2) = 0;
      B_matrix_contribution(2,ishape*num_dim+2) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(3,ishape*num_dim+2) = 0;
      B_matrix_contribution(4,ishape*num_dim+2) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(5,ishape*num_dim+2) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
    }
    /*
    //debug print of B matrix per quadrature point
    std::cout << " ------------B MATRIX QUADRATURE CONTRIBUTION"<< ielem + 1 <<"--------------"<<std::endl;
    for (int idof = 0; idof < Brows; idof++){
      std::cout << "row: " << idof + 1 << " { ";
      for (int istride = 0; istride < nodes_per_elem*num_dim; istride++){
        std::cout << istride + 1 << " = " << B_matrix_contribution(idof,istride) << " , " ;
      }
      std::cout << " }"<< std::endl;
    }
    //end debug block
    */
    //accumulate B matrix
    for(int irow=0; irow < Brows; irow++)
      for(int icol=0; icol < num_dim*nodes_per_elem; icol++)
      B_matrix(irow,icol) += B_matrix_contribution(irow,icol);

    //compute the previous multiplied by the Elastic (C) Matrix
    for(int irow=0; irow < Brows; irow++){
      for(int icol=0; icol < num_dim*nodes_per_elem; icol++){
        CB_matrix_contribution(irow,icol) = 0;
        for(int span=0; span < Brows; span++){
          CB_matrix_contribution(irow,icol) += C_matrix(irow,span)*B_matrix_contribution(span,icol);
        }
      }
    }

    //accumulate CB matrix
    for(int irow=0; irow < Brows; irow++)
      for(int icol=0; icol < num_dim*nodes_per_elem; icol++)
      CB_matrix(irow,icol) += CB_matrix_contribution(irow,icol);

    //compute the contributions of this quadrature point to all the local stiffness matrix elements
    for(int ifill=0; ifill < num_dim*nodes_per_elem; ifill++)
      for(int jfill=0; jfill < num_dim*nodes_per_elem; jfill++){
        matrix_term = 0;
        for(int span = 0; span < Brows; span++){
          matrix_term += B_matrix_contribution(span,ifill)*CB_matrix_contribution(span,jfill);
        }
        Local_Matrix(ifill,jfill) += Elastic_Constant*quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2)*matrix_term/Jacobian;
      }
      
    }
    
    /*
    //debug print of B matrix per quadrature point
    std::cout << " ------------B MATRIX"<< ielem + 1 <<"--------------"<<std::endl;
    for (int idof = 0; idof < Brows; idof++){
      std::cout << "row: " << idof + 1 << " { ";
      for (int istride = 0; istride < nodes_per_elem*num_dim; istride++){
        std::cout << istride + 1 << " = " << B_matrix(idof,istride) << " , " ;
      }
      std::cout << " }"<< std::endl;
    }
    //end debug block

    //debug print of B matrix per quadrature point
    std::cout << " ------------CB MATRIX"<< ielem + 1 <<"--------------"<<std::endl;
    for (int idof = 0; idof < Brows; idof++){
      std::cout << "row: " << idof + 1 << " { ";
      for (int istride = 0; istride < nodes_per_elem*num_dim; istride++){
        std::cout << istride + 1 << " = " << CB_matrix(idof,istride) << " , " ;
      }
      std::cout << " }"<< std::endl;
    }
    //end debug block
    */
}

/* ----------------------------------------------------------------------
   Loop through applied boundary conditions and tag node ids to redecule 
   necessary rows and columns from the assembled linear system
------------------------------------------------------------------------- */

void Static_Solver::Displacement_Boundary_Conditions(){
  int num_bdy_patches_in_set;
  int warning_flag = 0;
  int current_node_index, patch_gid, current_node_id;
  int num_nodes = mesh->num_nodes();
  int num_bdy_sets = mesh->num_bdy_sets();
  int surface_disp_set_id = 0;
  int num_dim = simparam->num_dim;
  int bc_option, bc_dim_set[3];
  CArray<real_t> displacement(num_dim);
  CArray<size_t> Displacement_Conditions(num_dim);
  CArray<size_t> first_condition_per_node(num_nodes*num_dim);
  Number_DOF_BCS = 0;
  Displacement_Conditions(0) = X_DISPLACEMENT_CONDITION;
  Displacement_Conditions(1) = Y_DISPLACEMENT_CONDITION;
  Displacement_Conditions(2) = Z_DISPLACEMENT_CONDITION;

  //initialize to -1 (DO NOT make -1 an index for bdy sets)
  for(int inode = 0 ; inode < num_nodes*num_dim; inode++)
    first_condition_per_node(inode) = -1;
  
  //scan for surface method of setting fixed nodal displacements
  for(int iboundary = 0; iboundary < num_bdy_sets; iboundary++){
    
    if(Boundary_Condition_Type_List(iboundary)==DISPLACEMENT_CONDITION){bc_option=3;}
    else if(Boundary_Condition_Type_List(iboundary)==X_DISPLACEMENT_CONDITION){bc_option=0;}
    else if(Boundary_Condition_Type_List(iboundary)==Y_DISPLACEMENT_CONDITION){bc_option=1;}
    else if(Boundary_Condition_Type_List(iboundary)==Z_DISPLACEMENT_CONDITION){bc_option=2;}
    else{
      continue;
    }
      
      //debug print of surface conditions
      //std::cout << " Surface BC types " << Boundary_Condition_Type_List(iboundary) <<std::endl;

      num_bdy_patches_in_set = mesh->num_bdy_patches_in_set(iboundary);
      if(bc_option==0) {
        bc_dim_set[0]=1;
        displacement(0) = Boundary_Surface_Displacements(surface_disp_set_id,0);
      }
      else if(bc_option==1) {
        bc_dim_set[1]=1;
        displacement(1) = Boundary_Surface_Displacements(surface_disp_set_id,1);
      }
      else if(bc_option==2) {
        bc_dim_set[2]=1;
        displacement(2) = Boundary_Surface_Displacements(surface_disp_set_id,2);
      }
      else if(bc_option==3) {
        bc_dim_set[0]=1;
        bc_dim_set[1]=1;
        bc_dim_set[2]=1;
        displacement(0) = Boundary_Surface_Displacements(surface_disp_set_id,0);
        displacement(1) = Boundary_Surface_Displacements(surface_disp_set_id,1);
        displacement(2) = Boundary_Surface_Displacements(surface_disp_set_id,2);
      }
      surface_disp_set_id++;
      
      //loop over boundary set patches for their respective nodes
      for (int bdy_patch_gid = 0; bdy_patch_gid < num_bdy_patches_in_set; bdy_patch_gid++){
        // get the global id for this boundary patch
        patch_gid = mesh->bdy_patches_in_set(iboundary, bdy_patch_gid);
        for(int inode = 0; inode < 4; inode++){
          current_node_id = mesh->node_in_patch(patch_gid, inode);
          
          /*
          //debug print of nodal bc settings
          std::cout << " Node BC types " << Node_DOF_Boundary_Condition_Type(current_node_id);
          std::cout << " node: " << inode << " { ";
          for (int istride = 0; istride < num_dim; istride++){
            std::cout << mesh->node_coords(current_node_id,istride) << " , ";
          }
          std::cout << " }"<< std::endl;
          */

          for(int idim=0; idim < num_dim; idim++){
          //warning for reapplied a displacement boundary condition (For now there is an output flag being set that triggers output later)
          if(Node_DOF_Boundary_Condition_Type(current_node_id*num_dim + idim)==DISPLACEMENT_CONDITION||
          Node_DOF_Boundary_Condition_Type(current_node_id*num_dim + idim)==Displacement_Conditions(idim)){
            //if overlap is just due to the loop over patches, a warning is not needed
            if(first_condition_per_node(current_node_id*num_dim + idim)!=iboundary) warning_flag = 1;
          }
          else{
            if(bc_dim_set[idim]){
              first_condition_per_node(current_node_id*num_dim + idim) = iboundary;
              Node_DOF_Boundary_Condition_Type(current_node_id*num_dim+idim) = Boundary_Condition_Type_List(iboundary);
              Node_DOF_Displacement_Boundary_Conditions(current_node_id*num_dim+idim) = displacement(idim);
              Number_DOF_BCS++;
            }
          }
          }
        }
      }
  }

  //scan for direct setting of nodal displacements from input
  //indices for nodal BC settings referred to here start at num_bdy_sets

  //debug print of nodal bc settings
  /*
  std::cout << " ------------NODE BC SETTINGS--------------"<<std::endl;
  for(int inode=0; inode < num_nodes*num_dim; inode++)
  std::cout << " Node BC types " << Node_DOF_Boundary_Condition_Type(inode) <<std::endl;
  //end debug block
  */

  //print warning for overlapping boundary conditions
  if(warning_flag)
  std::cout << std::endl << "One or more displacement boundary conditions overlap on a subset of nodes; please revise input" << std::endl << std::endl;

}

/* ----------------------------------------------------------------------
   Construct the global applied force vector
------------------------------------------------------------------------- */

void Static_Solver::Force_Vector_Construct(){

  int num_dim = simparam->num_dim;
  int num_elems = mesh->num_elems();
  int num_nodes = mesh->num_nodes();
  int nodes_per_elem = mesh->num_nodes_in_elem();
  int num_gauss_points = simparam->num_gauss_points;
  int z_quad,y_quad,x_quad, direct_product_count;
  int current_element_index, local_surface_id, surf_dim1, surf_dim2;
  CArray<real_t> legendre_nodes_1D(num_gauss_points);
  CArray<real_t> legendre_weights_1D(num_gauss_points);
  real_t pointer_quad_coordinate[num_dim];
  real_t pointer_quad_coordinate_weight[num_dim];
  real_t pointer_interpolated_point[num_dim];
  ViewCArray<real_t> quad_coordinate(pointer_quad_coordinate,num_dim);
  ViewCArray<real_t> quad_coordinate_weight(pointer_quad_coordinate_weight,num_dim);
  ViewCArray<real_t> interpolated_point(pointer_interpolated_point,num_dim);
  real_t force_density[3], wedge_product;
  
  CArray<real_t> JT_row1(num_dim);
  CArray<real_t> JT_row2(num_dim);
  CArray<real_t> JT_row3(num_dim);

  int num_bdy_patches_in_set;
  int current_node_index;
  int num_bdy_sets = mesh->num_bdy_sets();
  int surface_force_set_id = 0;

  /*Loop through boundary sets and check if they apply surface forces.
  These sets can have overlapping nodes since applied loading conditions
  are assumed to be additive*/
  for(int iboundary = 0; iboundary < num_bdy_sets; iboundary++){
    if(Boundary_Condition_Type_List(iboundary)!=LOADING_CONDITION) continue;
    //std::cout << "I REACHED THE LOADING BOUNDARY CONDITION" <<std::endl;
    num_bdy_patches_in_set = mesh->num_bdy_patches_in_set(iboundary);
    
    force_density[0] = Boundary_Surface_Force_Densities(surface_force_set_id,0);
    force_density[1] = Boundary_Surface_Force_Densities(surface_force_set_id,1);
    force_density[2] = Boundary_Surface_Force_Densities(surface_force_set_id,2);
    surface_force_set_id++;

    real_t pointer_basis_values[elem->num_basis()];
    ViewCArray<real_t> basis_values(pointer_basis_values,elem->num_basis());

    //initialize weights
    elements::legendre_nodes_1D(legendre_nodes_1D,num_gauss_points);
    elements::legendre_weights_1D(legendre_weights_1D,num_gauss_points);

    direct_product_count = std::pow(num_gauss_points,num_dim-1);
    //loop over boundary sets for their applied forces; use quadrature for distributed forces
    for (int bdy_patch_gid = 0; bdy_patch_gid < num_bdy_patches_in_set; bdy_patch_gid++){
                
    // get the global id for this boundary patch
    int patch_gid = mesh->bdy_patches_in_set(iboundary, bdy_patch_gid);
    current_element_index = mesh->cells_in_patch(patch_gid,0);
    if(current_element_index==-1) current_element_index = mesh->cells_in_patch(patch_gid,1);
    local_surface_id = mesh->cells_in_patch_local_id(patch_gid,0);
    if(local_surface_id==-1) local_surface_id = mesh->cells_in_patch_local_id(patch_gid,1);

    //loop over quadrature points if this is a distributed force
    for(int iquad=0; iquad < direct_product_count; iquad++){
      
      if(Element_Types(current_element_index)==4){

      real_t pointer_basis_values[8];
      real_t pointer_basis_derivative_s1[8];
      real_t pointer_basis_derivative_s2[8];
      ViewCArray<real_t> basis_values(pointer_basis_values,8);
      ViewCArray<real_t> basis_derivative_s1(pointer_basis_derivative_s1,8);
      ViewCArray<real_t> basis_derivative_s2(pointer_basis_derivative_s2,8);
      CArray<real_t> nodal_positions(4,num_dim);
      CArray<real_t> surf_basis_derivative_s1(4,num_dim);
      CArray<real_t> surf_basis_derivative_s2(4,num_dim);
      CArray<real_t> surf_basis_values(4,num_dim);
      int local_nodes[4];
      //set current quadrature point
      y_quad = iquad / num_gauss_points;
      x_quad = iquad % num_gauss_points;
      
      if(local_surface_id<2){
      surf_dim1 = 0;
      surf_dim2 = 1;
      quad_coordinate(0) = legendre_nodes_1D(x_quad);
      quad_coordinate(1) = legendre_nodes_1D(y_quad);
      //set to -1 or 1 for an isoparametric space
        if(local_surface_id%2==0)
        quad_coordinate(2) = -1;
        else
        quad_coordinate(2) = 1;
      }
      else if(local_surface_id<4){
      surf_dim1 = 0;
      surf_dim2 = 2;
      quad_coordinate(0) = legendre_nodes_1D(x_quad);
      quad_coordinate(2) = legendre_nodes_1D(y_quad);
      //set to -1 or 1 for an isoparametric space
        if(local_surface_id%2==0)
        quad_coordinate(1) = -1;
        else
        quad_coordinate(1) = 1;
      }
      else if(local_surface_id<6){
      surf_dim1 = 1;
      surf_dim2 = 2;
      quad_coordinate(1) = legendre_nodes_1D(x_quad);
      quad_coordinate(2) = legendre_nodes_1D(y_quad);
      //set to -1 or 1 for an isoparametric space
        if(local_surface_id%2==0)
        quad_coordinate(0) = -1;
        else
        quad_coordinate(0) = 1;
      }
      
      //set current quadrature weight
      quad_coordinate_weight(0) = legendre_weights_1D(x_quad);
      quad_coordinate_weight(1) = legendre_weights_1D(y_quad);

      //find local dof set for this surface
      local_nodes[0] = elem->surface_to_dof_lid(local_surface_id,0);
      local_nodes[1] = elem->surface_to_dof_lid(local_surface_id,1);
      local_nodes[2] = elem->surface_to_dof_lid(local_surface_id,2);
      local_nodes[3] = elem->surface_to_dof_lid(local_surface_id,3);

      //acquire set of nodes for this face
      for(int node_loop=0; node_loop < 4; node_loop++){
        current_node_index = mesh->nodes_in_cell(current_element_index, local_nodes[node_loop]);
        nodal_positions(node_loop,0) = mesh->node_coords(current_node_index,0);
        nodal_positions(node_loop,1) = mesh->node_coords(current_node_index,1);
        nodal_positions(node_loop,2) = mesh->node_coords(current_node_index,2);
      }

      if(local_surface_id<2){
        //compute shape function derivatives
        elem->partial_xi_basis(basis_derivative_s1,quad_coordinate);
        elem->partial_eta_basis(basis_derivative_s2,quad_coordinate);
      }
      else if(local_surface_id<4){
        //compute shape function derivatives
        elem->partial_xi_basis(basis_derivative_s1,quad_coordinate);
        elem->partial_mu_basis(basis_derivative_s2,quad_coordinate);
      }
      else if(local_surface_id<6){
        //compute shape function derivatives
        elem->partial_eta_basis(basis_derivative_s1,quad_coordinate);
        elem->partial_mu_basis(basis_derivative_s2,quad_coordinate);
      }

      //set values relevant to this surface
      for(int node_loop=0; node_loop < 4; node_loop++){
        surf_basis_derivative_s1(node_loop) = basis_derivative_s1(local_nodes[node_loop]);
        surf_basis_derivative_s2(node_loop) = basis_derivative_s2(local_nodes[node_loop]);
      }

      //compute derivatives of x,y,z w.r.t the s,t coordinates of this surface; needed to compute dA in surface integral
      //derivative of x,y,z w.r.t s
      JT_row1(0) = 0;
      JT_row1(1) = 0;
      JT_row1(2) = 0;
      for(int node_loop=0; node_loop < 4; node_loop++){
        JT_row1(0) += nodal_positions(node_loop,0)*surf_basis_derivative_s1(node_loop);
        JT_row1(1) += nodal_positions(node_loop,1)*surf_basis_derivative_s1(node_loop);
        JT_row1(2) += nodal_positions(node_loop,2)*surf_basis_derivative_s1(node_loop);
      }

      //derivative of x,y,z w.r.t t
      JT_row2(0) = 0;
      JT_row2(1) = 0;
      JT_row2(2) = 0;
      for(int node_loop=0; node_loop < 4; node_loop++){
        JT_row2(0) += nodal_positions(node_loop,0)*surf_basis_derivative_s2(node_loop);
        JT_row2(1) += nodal_positions(node_loop,1)*surf_basis_derivative_s2(node_loop);
        JT_row2(2) += nodal_positions(node_loop,2)*surf_basis_derivative_s2(node_loop);
      }
      

      //compute jacobian for this surface
      //compute the determinant of the Jacobian
      wedge_product = sqrt(pow(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2),2)+
               pow(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2),2)+
               pow(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1),2));

      //compute shape functions at this point for the element type
      elem->basis(basis_values,quad_coordinate);

      // loop over nodes of this face and 
      for(int node_count = 0; node_count < 4; node_count++){
            
        int node_gid = mesh->nodes_in_cell(current_element_index, local_nodes[node_count]);
        
        /*
        //debug print block
        std::cout << " ------------Element "<< current_element_index + 1 <<"--------------"<<std::endl;
        std::cout <<  " = , " << " Wedge Product: " << wedge_product << " local node " << local_nodes[node_count] << " node " << node_gid + 1<< " : s " 
        << quad_coordinate(0) << " t " << quad_coordinate(1) << " w " << quad_coordinate(2) << " basis value  "<< basis_values(local_nodes[node_count])
        << " Nodal Force value"<< Nodal_Forces(num_dim*node_gid)+wedge_product*quad_coordinate_weight(0)*quad_coordinate_weight(1)*force_density[0]*basis_values(local_nodes[node_count]);
        
        std::cout << " }"<< std::endl;
        //end debug print block
        */

        // Accumulate force vector contribution from this quadrature point
        for(int idim = 0; idim < num_dim; idim++){
          if(force_density[idim]!=0)
          //Nodal_Forces(num_dim*node_gid + idim) += wedge_product*quad_coordinate_weight(0)*quad_coordinate_weight(1)*force_density[idim]*basis_values(local_nodes[node_count]);
          Nodal_Forces(num_dim*node_gid + idim) += wedge_product*quad_coordinate_weight(0)*quad_coordinate_weight(1)*force_density[idim]*basis_values(local_nodes[node_count]);
        }
      }
      }
    }
  }
  }

    //apply line distribution of forces

    //apply point forces

    //apply contribution from non-zero displacement boundary conditions

    //debug print of force vector
    /*
    std::cout << "---------FORCE VECTOR-------------" << std::endl;
    for(int iforce=0; iforce < num_nodes*num_dim; iforce++)
      std::cout << " DOF: "<< iforce+1 << ", "<< Nodal_Forces(iforce) << std::endl;
    */

}

/* ----------------------------------------------------------------------
   Solve the FEA linear system
------------------------------------------------------------------------- */

int Static_Solver::solve(){

  typedef double Scalar;
  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;

  typedef Tpetra::CrsMatrix<Scalar,LO,GO> MAT;
  typedef Tpetra::MultiVector<Scalar,LO,GO> MV;

  typedef Kokkos::ViewTraits<LO*, Kokkos::LayoutLeft, void, void>::size_type SizeType;
  typedef Tpetra::Details::DefaultTypes::node_type node_type;
  using traits = Kokkos::ViewTraits<LO*, Kokkos::LayoutLeft, void, void>;
  
  using array_layout    = typename traits::array_layout;
  using execution_space = typename traits::execution_space;
  using device_type     = typename traits::device_type;
  using memory_traits   = typename traits::memory_traits;

  using Tpetra::global_size_t;
  using Teuchos::tuple;
  using Teuchos::RCP;
  using Teuchos::rcp;

  int num_dim = simparam->num_dim;
  int num_elems = mesh->num_elems();
  int num_nodes = mesh->num_nodes();
  int nodes_per_elem = mesh->num_nodes_in_elem();
  int current_row_n_nodes_scanned;
  int local_node_index, current_row, current_column;
  int max_stride = 0;
  global_size_t global_index, reduced_row_count;

  Teuchos::RCP<const Teuchos::Comm<int> > comm =
  Tpetra::getDefaultComm();
  
  typedef Kokkos::View<Scalar*, Kokkos::LayoutRight, device_type, memory_traits> values_array;
  typedef Kokkos::View<GO*, array_layout, device_type, memory_traits> global_indices_array;
  typedef Kokkos::View<LO*, array_layout, device_type, memory_traits> indices_array;
  typedef Kokkos::View<SizeType*, array_layout, device_type, memory_traits> row_pointers;
 
  using vec_map_type = Tpetra::Map<LO, GO>;
  using vec_device_type = typename vec_map_type::device_type;
  typedef Kokkos::DualView<Scalar**, Kokkos::LayoutLeft, device_type>::t_dev vec_array;

  // Before we do anything, check that KLU2 is enabled
  if( !Amesos2::query("KLU2") ){
    std::cerr << "KLU2 not enabled in this run.  Exiting..." << std::endl;
    return EXIT_SUCCESS;        // Otherwise CTest will pick it up as
                                // failure, which it isn't really
  }

  size_t myRank = comm->getRank();

  std::ostream &out = std::cout;

  out << Amesos2::version() << std::endl << std::endl;

  const size_t numVectors = 1;
  
  //construct global data for example; normally this would be from file input and distributed
  //according to the row map at that point
  global_size_t nrows = num_nodes*num_dim;
 
  //check for number of boundary conditions on this mpi rank
  global_size_t nboundaries = Number_DOF_BCS;
  global_size_t nrows_reduced = nrows - nboundaries;

  //Rebalance distribution of the global stiffness matrix rows here later since
  //rows and columns are being removed.

  // create a Map for the global stiffness matrix (move this up in the code later)
  RCP<Tpetra::Map<LO,GO,node_type> > map
    = rcp( new Tpetra::Map<LO,GO,node_type>(nrows,0,comm) );
  
  global_size_t local_nrows = map->getNodeNumElements();
  
  //populate local row offset data from global data
  global_size_t min_gid = map->getMinGlobalIndex();
  global_size_t max_gid = map->getMaxGlobalIndex();
  global_size_t index_base = map->getIndexBase();

  global_size_t *entries_per_row = new global_size_t[local_nrows];
  row_pointers row_offsets = row_pointers("row_offsets", local_nrows+1);

  //init row_offsets
  for(int i=0; i < local_nrows+1; i++){
    row_offsets(i) = 0;
  }

  //allocate upper bound storage for global indices this task owns for the reduced stiffness matrix
  CArrayKokkos<GO, array_layout, device_type, memory_traits> my_reduced_global_indices_buffer(local_nrows);
  
  //find out which rows to remove from the global stiffness matrix due to boundary conditions
  size_t index_counter = 0;
  for(LO i=0; i < local_nrows; i++){
    global_index = map->getGlobalElement(i);
    //change to local index access when rowmap moves up i.e. Node_DOF_Boundary_Condition_Type stores local node data
    if((Node_DOF_Boundary_Condition_Type(global_index)!=DISPLACEMENT_CONDITION)&&
       (Node_DOF_Boundary_Condition_Type(global_index)!=X_DISPLACEMENT_CONDITION)&&
       (Node_DOF_Boundary_Condition_Type(global_index)!=Y_DISPLACEMENT_CONDITION)&&
       (Node_DOF_Boundary_Condition_Type(global_index)!=Z_DISPLACEMENT_CONDITION))
    my_reduced_global_indices_buffer(index_counter++) = global_index;
  }
  
  //allocate storage with no extra spaces for global indices this task owns for the reduced stiffness matrix
  CArrayKokkos<GO, array_layout, device_type, memory_traits> my_reduced_global_indices(index_counter);

  //copy values from buffer to compressed array
  for(int ifill=0; ifill < index_counter; ifill++)
    my_reduced_global_indices(ifill) = my_reduced_global_indices_buffer(ifill);
  
  //debug print of nodal bc settings
  //for(int inode=0; inode < index_counter; inode++)
  //std::cout << " my_reduced_global_indices " << my_reduced_global_indices(inode) <<std::endl;
  
  // create a Map for the reduced global stiffness matrix (BC rows removed)
  RCP<Tpetra::Map<LO,GO,node_type> > reduced_map
    = rcp( new Tpetra::Map<LO,GO,node_type>(nrows_reduced,my_reduced_global_indices.get_kokkos_view(),0,comm) );
  
  //local rows for the reduced stiffness matrix
  global_size_t local_nrows_reduced = reduced_map->getNodeNumElements();

  global_size_t *reduced_entries_per_row = new global_size_t[local_nrows_reduced];
  row_pointers reduced_row_offsets = row_pointers("reduced_row_offsets", local_nrows_reduced+1);

  //init row_offsets
  for(int i=0; i < local_nrows_reduced+1; i++){
    reduced_row_offsets(i) = 0;
  }
  
  //count entries per row (doesn't count nodal indices that correspond to displacement boundary conditions)
  for(LO i=0; i < local_nrows_reduced; i++){
    global_index = reduced_map->getGlobalElement(i);
    reduced_row_count = 0;
    reduced_entries_per_row[i] = Stiffness_Matrix_strides(global_index);
    for(LO j=0; j < reduced_entries_per_row[i]; j++){
      if((Node_DOF_Boundary_Condition_Type(DOF_Graph_Matrix(global_index,j))!=DISPLACEMENT_CONDITION)&&
         (Node_DOF_Boundary_Condition_Type(DOF_Graph_Matrix(global_index,j))!=X_DISPLACEMENT_CONDITION)&&
         (Node_DOF_Boundary_Condition_Type(DOF_Graph_Matrix(global_index,j))!=Y_DISPLACEMENT_CONDITION)&&
         (Node_DOF_Boundary_Condition_Type(DOF_Graph_Matrix(global_index,j))!=Z_DISPLACEMENT_CONDITION)){
        reduced_row_count++;
      }
    }
    reduced_row_offsets(i+1) = reduced_row_count + reduced_row_offsets(i);
  }
  
  //remove number of boundary conditions
  global_size_t nnz = reduced_row_offsets(local_nrows_reduced);

  //indices_array all_indices = indices_array("indices_array",nnz);
  //values_array all_values = values_array("values_array",nnz);
  CArrayKokkos<GO, array_layout, device_type, memory_traits> all_global_indices(nnz);
  CArrayKokkos<LO, array_layout, device_type, memory_traits> all_indices(nnz);
  CArrayKokkos<Scalar, Kokkos::LayoutRight, device_type, memory_traits> all_values(nnz);
  CArrayKokkos<Scalar, Kokkos::LayoutLeft, device_type, memory_traits> Bview(local_nrows_reduced);
  CArrayKokkos<Scalar, Kokkos::LayoutLeft, device_type, memory_traits> Xview(local_nrows_reduced);

  //set Kokkos view data
  LO entrycount = 0;
  for(LO i=0; i < local_nrows_reduced; i++){
    global_index = reduced_map->getGlobalElement(i);
    for(LO j=0; j < reduced_entries_per_row[i]; j++){
      if((Node_DOF_Boundary_Condition_Type(DOF_Graph_Matrix(global_index,j))!=DISPLACEMENT_CONDITION)&&
         (Node_DOF_Boundary_Condition_Type(DOF_Graph_Matrix(global_index,j))!=X_DISPLACEMENT_CONDITION)&&
         (Node_DOF_Boundary_Condition_Type(DOF_Graph_Matrix(global_index,j))!=Y_DISPLACEMENT_CONDITION)&&
         (Node_DOF_Boundary_Condition_Type(DOF_Graph_Matrix(global_index,j))!=Z_DISPLACEMENT_CONDITION)){
        all_global_indices(entrycount) = DOF_Graph_Matrix(reduced_map->getGlobalElement(i),j);
        all_values(entrycount) = Stiffness_Matrix(reduced_map->getGlobalElement(i),j);
        entrycount++;
      }
    }
  }
  
  //build column map
  RCP<const Tpetra::Map<LO,GO,node_type> > colmap;
  const RCP<const Tpetra::Map<LO,GO,node_type> > dommap = reduced_map;
  
  Tpetra::Details::makeColMap<LO,GO,node_type>(colmap,dommap,all_global_indices.get_kokkos_view(), nullptr);

  //debug print of global indices
  std::cout << " DEBUG PRINT FOR GLOBAL INDICES" << std::endl;
  for(int debug = 0; debug < entrycount; debug++)
  std::cout <<  all_global_indices.get_kokkos_view()(debug) << " ";
  std::cout << std::endl;
  //end debug print

  //debug print of reduced row offsets
  std::cout << " DEBUG PRINT FOR ROW OFFSETS" << std::endl;
  for(int debug = 0; debug < local_nrows_reduced+1; debug++)
  std::cout <<  reduced_row_offsets(debug) << " ";
  std::cout << std::endl;
  //end debug print
  
  //convert global indices to local indices using column map
  entrycount = 0;
  for(LO i=0; i < local_nrows_reduced; i++){
    for(LO j=0; j < reduced_row_offsets(i+1) - reduced_row_offsets(i); j++){
    all_indices(entrycount) = colmap->getLocalElement(all_global_indices(entrycount));
    entrycount++;
    }
  }

  //debug print of local indices
  std::cout << " DEBUG PRINT FOR LOCAL INDICES" << std::endl;
  for(int debug = 0; debug < entrycount; debug++)
  std::cout <<  all_indices(debug) << " ";
  std::cout << std::endl;
  //end debug print
  
  RCP<MAT> A = rcp( new MAT(reduced_map, colmap, reduced_row_offsets, all_indices.get_kokkos_view(), all_values.get_kokkos_view()) );
  
   RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
  A->fillComplete();
  //This completes the setup for A matrix of the linear system

  //debug print of A matrix
  *fos << "Matrix :" << std::endl;
  A->describe(*fos,Teuchos::VERB_EXTREME);
  *fos << std::endl;
  
  
  /*
  //Matlab compatible debug print (checks solver result)
  std::cout << "Matlab A matrix" << std::endl;
  CArray<real_t> row_set(local_nrows_reduced);
  int column_index;
  
  for(int irow = 0 ; irow < local_nrows_reduced; irow++){
    global_index = reduced_map->getGlobalElement(irow);
    //initialize
    for(int icol = 0 ; icol < local_nrows_reduced; icol++) row_set(icol) = 0;
    //set non-zero elements
    for(int istride=0; istride < reduced_entries_per_row[irow]; istride++){
      if((Node_DOF_Boundary_Condition_Type(DOF_Graph_Matrix(global_index,istride))!=DISPLACEMENT_CONDITION)&&
         (Node_DOF_Boundary_Condition_Type(DOF_Graph_Matrix(global_index,istride))!=X_DISPLACEMENT_CONDITION)&&
         (Node_DOF_Boundary_Condition_Type(DOF_Graph_Matrix(global_index,istride))!=Y_DISPLACEMENT_CONDITION)&&
         (Node_DOF_Boundary_Condition_Type(DOF_Graph_Matrix(global_index,istride))!=Z_DISPLACEMENT_CONDITION)){
      column_index = colmap->getLocalElement(DOF_Graph_Matrix(reduced_map->getGlobalElement(irow),istride));
      //std::cout << "COLUMN INDEX " << column_index << std::endl;
      row_set(column_index) = istride+1;
      }
    }

    for(int icol = 0 ; icol < local_nrows_reduced; icol++){
      if(row_set(icol))
      std::cout << Stiffness_Matrix(reduced_map->getGlobalElement(irow),row_set(icol)-1)<<" " ;
      else
      std::cout << " 0" ;

      if(icol<local_nrows_reduced-1) std::cout << " ,";
    }
    std::cout << " ; " <<std::endl;
  }

  //end debug print
  */
  
  // Create random X vector
  vec_array Xview_pass = vec_array("Xview_pass", local_nrows_reduced,1);
  Xview_pass.assign_data(Xview.pointer());
  RCP<MV> X = rcp(new MV(reduced_map, Xview_pass));
  X->randomize();
  
  //print allocation of the solution vector to check distribution
  *fos << "ALLOCATED SOLUTION VECTOR :" << std::endl;
  X->describe(*fos,Teuchos::VERB_EXTREME);
  *fos << std::endl;
  
  // Create Kokkos view of RHS B vector (Force Vector)  
  vec_array Bview_pass = vec_array("Bview_pass", local_nrows_reduced,1);

  //set bview to force vector data
  for(LO i=0; i < local_nrows_reduced; i++)
    Bview(i) = Nodal_Forces(reduced_map->getGlobalElement(i));
  
  Bview_pass.assign_data(Bview.pointer());
  
  RCP<MV> B = rcp(new MV(reduced_map, Bview_pass));

  *fos << "RHS :" << std::endl;
  B->describe(*fos,Teuchos::VERB_EXTREME);
  *fos << std::endl;
  
  //debug block to print nodal positions corresponding to nodes of reduced system
  CArrayKokkos<Scalar, Kokkos::LayoutLeft, device_type, memory_traits> Node_Positions_View(local_nrows_reduced);

  vec_array Node_Positions_pass = vec_array("Node_Positions", local_nrows_reduced,1);

  //print the nodal positions
  for(LO i=0; i < local_nrows_reduced; i++){
    Node_Positions_View(i) = mesh->node_coords(reduced_map->getGlobalElement(i)/num_dim,reduced_map->getGlobalElement(i)%num_dim);
  }

  Node_Positions_pass.assign_data(Node_Positions_View.pointer());

  RCP<MV> Node_Positions = rcp(new MV(reduced_map, Node_Positions_pass));

  *fos << "Node_Positions :" << std::endl;
  Node_Positions->describe(*fos,Teuchos::VERB_EXTREME);
  *fos << std::endl;
  //end of debug block
  
  //std::fflush(stdout);
  //std::cout << " PRINT BEFORE ERROR " << std::endl;
  // Create solver interface to KLU2 with Amesos2 factory method
  RCP<Amesos2::Solver<MAT,MV> > solver = Amesos2::create<MAT,MV>("KLU2", A, X, B);
  
  //declare non-contiguous map
  // Create a Teuchos::ParameterList to hold solver parameters
  Teuchos::ParameterList amesos2_params("Amesos2");
  amesos2_params.sublist("KLU2").set("IsContiguous", false, "Are GIDs Contiguous");
  solver->setParameters( Teuchos::rcpFromRef(amesos2_params) );
  //Solve the system
  solver->symbolicFactorization().numericFactorization().solve();
  
  //Print solution vector
  *fos << "Solution :" << std::endl;
  X->describe(*fos,Teuchos::VERB_EXTREME);
  *fos << std::endl;
  
  return 0;
}