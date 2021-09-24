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
#include <fstream>
#include <string>
#include <sstream>
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
#include "Tpetra_Details_FixedHashTable.hpp"
#include "Tpetra_Import.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include <set>

#include "elements.h"
#include "swage.h"
#include "matar.h"
#include "utilities.h"
#include "header.h"
#include "state.h"
#include "node_combination.h"
#include "Simulation_Parameters.h"
#include "Amesos2_Version.hpp"
#include "Amesos2.hpp"
#include "Static_Solver_Parallel.h"

//debug and performance includes
#include <sys/time.h>
#include <sys/resource.h>
#include <ctime>

#define BUFFER_LINES 1000
#define MAX_WORD 30
#define MAX_ELEM_NODES 8

using namespace utils;

/*

Swage is a reference to a swage block used in blacksmithing.  
Its a large metal block that has multiple shaps carved into 
each surface to use for hammering metal into to form it. 

*/

Static_Solver_Parallel::Static_Solver_Parallel() : Solver(){
  //create parameter object
    simparam = new Simulation_Parameters();
  // ---- Read input file, define state and boundary conditions ---- //
    simparam->input();
  //create ref element object
    ref_elem = new elements::ref_element();
  //create mesh objects
    init_mesh = new swage::mesh_t(simparam);
    mesh = new swage::mesh_t(simparam);

    element_select = new elements::element_selector();
    num_nodes = 0;
}

Static_Solver_Parallel::~Static_Solver_Parallel(){
   delete simparam;
   delete ref_elem;
   delete mesh;
   delete element_select;
   if(myrank==0)
   delete in;
}

//==============================================================================
//    Primary simulation runtime routine
//==============================================================================


void Static_Solver_Parallel::run(int argc, char *argv[]){
    
    std::cout << "Running Static Solver Example" << std::endl;
    //MPI info
    world = MPI_COMM_WORLD; //used for convenience to represent all the ranks in the job
    MPI_Comm_rank(world,&myrank);
    MPI_Comm_size(world,&nranks);

    //initialize Trilinos communicator class
    comm = Tpetra::getDefaultComm();

    //error handle for file input name
    //if(argc < 2)

    // ---- Read intial mesh, refine, and build connectivity ---- //
    read_mesh(argv[1]);
    
    std::cout << "Num elements = " << mesh->num_elems() << std::endl;
    
    //initialize timing
    if(simparam->report_runtime_flag)
    init_clock();
    
    // ---- Find Boundaries on mesh ---- //
    generate_bcs();
    
    //allocate and fill sparse structures needed for global solution
    init_global();
    
    //assemble the global solution (stiffness matrix etc. and nodal forces)
    assemble();
    
    int solver_exit = solve();
    if(solver_exit == EXIT_SUCCESS){
      std::cout << "Before free pointer" << std::endl <<std::flush;
      return;
    }
    //CPU time
    double current_cpu = CPU_Time();
    std::cout << " RUNTIME OF CODE ON TASK " << myrank << " is "<< current_cpu-initial_CPU_time <<std::endl;
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
void Static_Solver_Parallel::read_mesh(char *MESH){

  char ch;
  int num_dim = simparam->num_dim;
  int p_order = simparam->p_order;
  std::string skip_line, read_line, substring;
  std::stringstream line_parse;
  CArray <char> read_buffer;
  //Nodes_Per_Element_Type =  elements::elem_types::Nodes_Per_Element_Type;

  //read the mesh
  //PLACEHOLDER: ensight_format(MESH);
  // abaqus_format(MESH);
  // vtk_format(MESH)

  //task 0 reads file
  if(myrank==0){
  in = new std::ifstream();
  in->open(MESH);  
    
  //skip 8 lines
  for (int j = 1; j <= 8; j++) {
    getline(*in, skip_line);
    std::cout << skip_line << std::endl;
  }
  }

  // --- Read the number of nodes in the mesh --- //
  if(myrank==0){
    getline(*in, read_line);
    line_parse.str(read_line);
    line_parse >> num_nodes;
    std::cout << "declared node count: " << num_nodes << std::endl;
  }
  
  //broadcast number of nodes
  MPI_Bcast(&num_nodes,1,MPI_INT,0,world);
  
  //construct contiguous parallel row map now that we know the number of nodes
  map = Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(num_nodes,0,comm));

  // set the vertices in the mesh read in
  global_size_t local_nrows = map->getNodeNumElements();
  nlocal_nodes = local_nrows;
  //populate local row offset data from global data
  global_size_t min_gid = map->getMinGlobalIndex();
  global_size_t max_gid = map->getMaxGlobalIndex();
  global_size_t index_base = map->getIndexBase();
  //debug print
  //std::cout << "local node count on task: " << " " << nlocal_nodes << std::endl;
  //construct dof map that follows from the node map (used for distributed matrix and vector objects later)
  CArrayKokkos<GO, array_layout, device_type, memory_traits> local_dof_indices(nlocal_nodes*num_dim, "local_dof_indices");
  for(int i = 0; i < nlocal_nodes; i++){
    for(int j = 0; j < num_dim; j++)
    local_dof_indices(i*num_dim + j) = map->getGlobalElement(i)*num_dim + j;
  }
  
  local_dof_map = Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(num_nodes*num_dim,local_dof_indices.get_kokkos_view(),0,comm) );

  //allocate node storage with dual view
  dual_node_data = dual_vec_array("dual_node_data", nlocal_nodes,num_dim);
  dual_nodal_forces = dual_vec_array("dual_nodal_forces", nlocal_nodes*num_dim,1);

  //local variable for host view in the dual view
  host_vec_array node_data = dual_node_data.view_host();
  //notify that the host view is going to be modified in the file readin
  dual_node_data.modify_host();

  //old swage method
  //mesh->init_nodes(local_nrows); // add 1 for index starting at 1
    
  std::cout << "Num nodes assigned to task " << myrank << " = " << nlocal_nodes << std::endl;

  // read the initial mesh coordinates
  // x-coords
  /*only task 0 reads in nodes and elements from the input file
  stores node data in a buffer and communicates once the buffer cap is reached
  or the data ends*/

  words_per_line = simparam->words_per_line;
  elem_words_per_line = simparam->elem_words_per_line;

  //allocate read buffer
  read_buffer = CArray<char>(BUFFER_LINES,words_per_line,MAX_WORD);

  //depends on file input format; this is for ensight
  int dof_limit = num_nodes;
  int buffer_loop, buffer_iteration, scan_loop;
  size_t read_index_start, node_gid, node_rid, elem_gid;
  real_t dof_value;
  int buffer_iterations = dof_limit/BUFFER_LINES;
  if(dof_limit%BUFFER_LINES!=0) buffer_iterations++;
  
  //x-coords
  read_index_start = 0;
  for(buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++){
    //pack buffer on rank 0
    if(myrank==0&&buffer_iteration<buffer_iterations-1){
      for (buffer_loop = 0; buffer_loop < BUFFER_LINES; buffer_loop++) {
        getline(*in,read_line);
        line_parse.clear();
        line_parse.str(read_line);
        
        for(int iword = 0; iword < words_per_line; iword++){
        //read portions of the line into the substring variable
        line_parse >> substring;
        //debug print
        //std::cout<<" "<< substring <<std::endl;
        //assign the substring variable as a word of the read buffer
        strcpy(&read_buffer(buffer_loop,iword,0),substring.c_str());
        }
      }
    }
    else if(myrank==0){
      buffer_loop=0;
      while(buffer_iteration*BUFFER_LINES+buffer_loop < num_nodes) {
        getline(*in,read_line);
        line_parse.clear();
        line_parse.str(read_line);
        for(int iword = 0; iword < words_per_line; iword++){
        //read portions of the line into the substring variable
        line_parse >> substring;
        //assign the substring variable as a word of the read buffer
        strcpy(&read_buffer(buffer_loop,iword,0),substring.c_str());
        }
        buffer_loop++;
      }
      
    }

    //broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
    MPI_Bcast(read_buffer.pointer(),BUFFER_LINES*words_per_line*MAX_WORD,MPI_CHAR,0,world);
    //broadcast how many nodes were read into this buffer iteration
    MPI_Bcast(&buffer_loop,1,MPI_INT,0,world);

    //debug_print
    //std::cout << "NODE BUFFER LOOP IS: " << buffer_loop << std::endl;
    //for(int iprint=0; iprint < buffer_loop; iprint++)
      //std::cout<<"buffer packing: " << std::string(&read_buffer(iprint,0,0)) << std::endl;
    //return;

    //determine which data to store in the swage mesh members (the local node data)
    //loop through read buffer
    for(scan_loop = 0; scan_loop < buffer_loop; scan_loop++){
      //set global node id (ensight specific order)
      node_gid = read_index_start + scan_loop;
      //let map decide if this node id belongs locally; if yes store data
      if(map->isNodeGlobalElement(node_gid)){
        //set local node index in this mpi rank
        node_rid = map->getLocalElement(node_gid);
        //extract nodal position from the read buffer
        //for ensight format this is just one coordinate per line
        dof_value = atof(&read_buffer(scan_loop,0,0));
        node_data(node_rid, 0) = dof_value;
      }
    }
    read_index_start+=BUFFER_LINES;
  }

  
  // y-coords
  read_index_start = 0;
  for(buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++){
    
    //pack buffer on rank 0
    if(myrank==0&&buffer_iteration<buffer_iterations-1){
      for (buffer_loop = 0; buffer_loop < BUFFER_LINES; buffer_loop++) {
        getline(*in,read_line);
        line_parse.clear();
        line_parse.str(read_line);
        for(int iword = 0; iword < words_per_line; iword++){
        //read portions of the line into the substring variable
        line_parse >> substring;
        //assign the substring variable as a word of the read buffer
        strcpy(&read_buffer(buffer_loop,iword,0),substring.c_str());
        }
      }
    }
    else if(myrank==0){
      buffer_loop=0;
      while(buffer_iteration*BUFFER_LINES+buffer_loop < num_nodes) {
        getline(*in,read_line);
        line_parse.clear();
        line_parse.str(read_line);
        for(int iword = 0; iword < words_per_line; iword++){
        //read portions of the line into the substring variable
        line_parse >> substring;
        //assign the substring variable as a word of the read buffer
        strcpy(&read_buffer(buffer_loop,iword,0),substring.c_str());
        }
        buffer_loop++;
        //std::cout<<" "<< node_data(node_gid, 0)<<std::endl;
      }
    }

    //broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
    MPI_Bcast(read_buffer.pointer(),BUFFER_LINES*words_per_line*MAX_WORD,MPI_CHAR,0,world);
    //broadcast how many nodes were read into this buffer iteration
    MPI_Bcast(&buffer_loop,1,MPI_INT,0,world);

    //determine which data to store in the swage mesh members (the local node data)
    //loop through read buffer
    for(scan_loop = 0; scan_loop < buffer_loop; scan_loop++){
      //set global node id (ensight specific order)
      node_gid = read_index_start + scan_loop;
      //let map decide if this node id belongs locally; if yes store data
      if(map->isNodeGlobalElement(node_gid)){
        //set local node index in this mpi rank
        node_rid = map->getLocalElement(node_gid);
        //extract nodal position from the read buffer
        //for ensight format this is just one coordinate per line
        dof_value = atof(&read_buffer(scan_loop,0,0));
        node_data(node_rid, 1) = dof_value;
      }
    }
    read_index_start+=BUFFER_LINES;
  }

  // z-coords
  read_index_start = 0;
  if(num_dim==3)
  for(buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++){
    
    //pack buffer on rank 0
    if(myrank==0&&buffer_iteration<buffer_iterations-1){
      for (buffer_loop = 0; buffer_loop < BUFFER_LINES; buffer_loop++) {
        getline(*in,read_line);
        line_parse.clear();
        line_parse.str(read_line);
        for(int iword = 0; iword < words_per_line; iword++){
        //read portions of the line into the substring variable
        line_parse >> substring;
        //assign the substring variable as a word of the read buffer
        strcpy(&read_buffer(buffer_loop,iword,0),substring.c_str());
        }
      }
    }
    else if(myrank==0){
      buffer_loop=0;
      while(buffer_iteration*BUFFER_LINES+buffer_loop < num_nodes) {
        getline(*in,read_line);
        line_parse.clear();
        line_parse.str(read_line);
        for(int iword = 0; iword < words_per_line; iword++){
        //read portions of the line into the substring variable
        line_parse >> substring;
        //assign the substring variable as a word of the read buffer
        strcpy(&read_buffer(buffer_loop,iword,0),substring.c_str());
        }
        buffer_loop++;
        //std::cout<<" "<< node_data(node_gid, 0)<<std::endl;
      }
    }

    //broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
    MPI_Bcast(read_buffer.pointer(),BUFFER_LINES*words_per_line*MAX_WORD,MPI_CHAR,0,world);
    //broadcast how many nodes were read into this buffer iteration
    MPI_Bcast(&buffer_loop,1,MPI_INT,0,world);

    //determine which data to store in the swage mesh members (the local node data)
    //loop through read buffer
    for(scan_loop = 0; scan_loop < buffer_loop; scan_loop++){
      //set global node id (ensight specific order)
      node_gid = read_index_start + scan_loop;
      //let map decide if this node id belongs locally; if yes store data
      if(map->isNodeGlobalElement(node_gid)){
        //set local node index in this mpi rank
        node_rid = map->getLocalElement(node_gid);
        //extract nodal position from the read buffer
        //for ensight format this is just one coordinate per line
        dof_value = atof(&read_buffer(scan_loop,0,0));
        node_data(node_rid, 2) = dof_value;
      }
    }
    read_index_start+=BUFFER_LINES;
  }
  
  
  //debug print of nodal data
  
  //debug print nodal positions and indices
  /*
  std::cout << " ------------NODAL POSITIONS ON TASK " << myrank << " --------------"<<std::endl;
  for (int inode = 0; inode < local_nrows; inode++){
      std::cout << "node: " << map->getGlobalElement(inode) + 1 << " { ";
    for (int istride = 0; istride < num_dim; istride++){
        std::cout << node_data(inode,istride) << " , ";
    }
    std::cout << " }"<< std::endl;
  }
  */

  //check that local assignments match global total

  
  //read in element info (ensight file format is organized in element type sections)
  //loop over this later for several element type sections

  num_elem = 0;
  rnum_elem = 0;
  CArray<int> node_store(elem_words_per_line);

  if(myrank==0){
  //skip element type name line
    getline(*in, skip_line);
    std::cout << skip_line << std::endl;
  }
    
  // --- read the number of cells in the mesh ---
  // --- Read the number of vertices in the mesh --- //
  if(myrank==0){
    getline(*in, read_line);
    line_parse.clear();
    line_parse.str(read_line);
    line_parse >> num_elem;
    std::cout << "declared element count: " << num_elem << std::endl;
    if(num_elem <= 0) std::cout << "ERROR, NO ELEMENTS IN MESH" << std::endl;
  }
  
  //broadcast number of elements
  MPI_Bcast(&num_elem,1,MPI_INT,0,world);

  std::cout<<"before initial mesh initialization"<<std::endl;
  
  //read in element connectivity
  //we're gonna reallocate for the words per line expected for the element connectivity
  read_buffer = CArray<char>(BUFFER_LINES,elem_words_per_line,MAX_WORD);

  //calculate buffer iterations to read number of lines
  buffer_iterations = num_elem/BUFFER_LINES;
  int assign_flag;
  //dynamic buffer used to store elements before we know how many this rank needs
  std::vector<size_t> element_temp(BUFFER_LINES*elem_words_per_line);
  if(num_elem%BUFFER_LINES!=0) buffer_iterations++;

  read_index_start = 0;
  //std::cout << "ELEMENT BUFFER ITERATIONS: " << buffer_iterations << std::endl;
  rnum_elem = 0;
  for(buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++){
    
    //pack buffer on rank 0
    if(myrank==0&&buffer_iteration<buffer_iterations-1){
      for (buffer_loop = 0; buffer_loop < BUFFER_LINES; buffer_loop++) {
        getline(*in,read_line);
        line_parse.clear();
        line_parse.str(read_line);
        for(int iword = 0; iword < elem_words_per_line; iword++){
        //read portions of the line into the substring variable
        line_parse >> substring;
        //assign the substring variable as a word of the read buffer
        strcpy(&read_buffer(buffer_loop,iword,0),substring.c_str());
        }
      }
    }
    else if(myrank==0){
      buffer_loop=0;
      while(buffer_iteration*BUFFER_LINES+buffer_loop < num_elem) {
        getline(*in,read_line);
        line_parse.clear();
        line_parse.str(read_line);
        for(int iword = 0; iword < elem_words_per_line; iword++){
        //read portions of the line into the substring variable
        line_parse >> substring;
        //assign the substring variable as a word of the read buffer
        strcpy(&read_buffer(buffer_loop,iword,0),substring.c_str());
        }
        buffer_loop++;
        //std::cout<<" "<< node_data(node_gid, 0)<<std::endl;
      }
    }

    //broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
    MPI_Bcast(read_buffer.pointer(),BUFFER_LINES*elem_words_per_line*MAX_WORD,MPI_CHAR,0,world);
    //broadcast how many nodes were read into this buffer iteration
    MPI_Bcast(&buffer_loop,1,MPI_INT,0,world);
    
    //determine which data to store in the swage mesh members (the local element data)
    //loop through read buffer

    //std::cout << "ELEMENT BUFFER LOOP IS: " << buffer_loop << std::endl;
    
    for(scan_loop = 0; scan_loop < buffer_loop; scan_loop++){
      //set global node id (ensight specific order)
      elem_gid = read_index_start + scan_loop;
      //add this element to the local list if any of its nodes belong to this rank according to the map
      //get list of nodes for each element line and check if they belong to the map
      assign_flag = 0;
      for(int inode = 0; inode < elem_words_per_line; inode++){
        //as we loop through the nodes belonging to this element we store them
        //if any of these nodes belongs to this rank this list is used to store the element locally
        node_gid = atoi(&read_buffer(scan_loop,inode,0));
        node_store(inode) = node_gid - 1; //subtract 1 since file index start is 1 but code expects 0
        //first we add the elements to a dynamically allocated list
        if(map->isNodeGlobalElement(node_gid-1)&&!assign_flag){
          assign_flag = 1;
          rnum_elem++;
        }
      }

      if(assign_flag)
      for(int inode = 0; inode < elem_words_per_line; inode++){
        if((rnum_elem-1)*elem_words_per_line + inode>=BUFFER_LINES*elem_words_per_line) 
          element_temp.resize((rnum_elem-1)*elem_words_per_line + inode + BUFFER_LINES*elem_words_per_line);
        element_temp[(rnum_elem-1)*elem_words_per_line + inode] = node_store(inode); 
        //std::cout << "VECTOR STORAGE FOR ELEM " << rnum_elem << " ON TASK " << myrank << " NODE " << inode+1 << " IS " << node_store(inode) + 1 << std::endl;
      }
    }
    read_index_start+=BUFFER_LINES;
  }
  
  std::cout << "RNUM ELEMENTS IS: " << rnum_elem << std::endl;
  //copy temporary element storage to swage arrays
  mesh->init_element(0, 3, rnum_elem);
  mesh->init_cells(rnum_elem);
  Element_Types = CArray<elements::elem_types::elem_type>(rnum_elem);
  for(int ielem = 0; ielem < rnum_elem; ielem++)
    for(int inode = 0; inode < elem_words_per_line; inode++)
      mesh->nodes_in_cell(ielem, inode) = element_temp[ielem*elem_words_per_line + inode];

  //delete temporary element connectivity
  //element_temp.~vector();
  std::vector<size_t>().swap(element_temp);

  //debug print element edof
  /*
  std::cout << " ------------ELEMENT EDOF ON TASK " << myrank << " --------------"<<std::endl;

  for (int ielem = 0; ielem < rnum_elem; ielem++){
    std::cout << "elem:  " << ielem+1 << std::endl;
    for (int lnode = 0; lnode < 8; lnode++){
        std::cout << "{ ";
          std::cout << lnode+1 << " = " << mesh->nodes_in_cell_list_(ielem,lnode) + 1 << " ";
        
        std::cout << " }"<< std::endl;
    }
    std::cout << std::endl;
  }
  */
 
  //simplified for now
  for(int ielem = 0; ielem < rnum_elem; ielem++)
    Element_Types(ielem) = elements::elem_types::Hex8;

  //element type selection (subject to change)
  // ---- Set Element Type ---- //
  // allocate element type memory
  //elements::elem_type_t* elem_choice;

  int NE = 1; // number of element types in problem
    
  //set base type pointer to one of the existing derived type object references
  if(simparam->num_dim==2)
  element_select->choose_2Delem_type(Element_Types(0), elem2D);
  else if(simparam->num_dim==3)
  element_select->choose_3Delem_type(Element_Types(0), elem);

  // Convert ijk index system to the finite element numbering convention
  // for vertices in cell
  max_nodes_per_element = MAX_ELEM_NODES;
  CArray<size_t> convert_ensight_to_ijk(max_nodes_per_element);
  CArray<size_t> tmp_ijk_indx(max_nodes_per_element);
  convert_ensight_to_ijk(0) = 0;
  convert_ensight_to_ijk(1) = 1;
  convert_ensight_to_ijk(2) = 3;
  convert_ensight_to_ijk(3) = 2;
  convert_ensight_to_ijk(4) = 4;
  convert_ensight_to_ijk(5) = 5;
  convert_ensight_to_ijk(6) = 7;
  convert_ensight_to_ijk(7) = 6;
    
  int nodes_per_element;
  
  if(num_dim==2)
  for (int cell_rid = 0; cell_rid < rnum_elem; cell_rid++) {
    //set nodes per element
    element_select->choose_2Delem_type(Element_Types(cell_rid), elem2D);
    nodes_per_element = elem2D->num_nodes();
    for (int node_lid = 0; node_lid < nodes_per_element; node_lid++){
      tmp_ijk_indx(node_lid) = mesh->nodes_in_cell(cell_rid, convert_ensight_to_ijk(node_lid));
    }   
        
    for (int node_lid = 0; node_lid < nodes_per_element; node_lid++){
      mesh->nodes_in_cell(cell_rid, node_lid) = tmp_ijk_indx(node_lid);
    }
  }

  if(num_dim==3)
  for (int cell_rid = 0; cell_rid < rnum_elem; cell_rid++) {
    //set nodes per element
    element_select->choose_3Delem_type(Element_Types(cell_rid), elem);
    nodes_per_element = elem->num_nodes();
    for (int node_lid = 0; node_lid < nodes_per_element; node_lid++){
      tmp_ijk_indx(node_lid) = mesh->nodes_in_cell(cell_rid, convert_ensight_to_ijk(node_lid));
    }   
        
    for (int node_lid = 0; node_lid < nodes_per_element; node_lid++){
      mesh->nodes_in_cell(cell_rid, node_lid) = tmp_ijk_indx(node_lid);
    }
  }

  // Build all connectivity in initial mesh
  std::cout<<"Before initial mesh connectivity"<<std::endl;
  
  if(rnum_elem >= 1) {

    // -- NODE TO CELL CONNECTIVITY -- //
    //mesh->build_node_cell_connectivity(); 

    // -- CORNER CONNECTIVITY -- //
    //mesh->build_corner_connectivity(); 

    // -- CELL TO CELL CONNECTIVITY -- //
    //mesh->build_cell_cell_connectivity(); 

    // -- PATCHES -- //
    //mesh->build_patch_connectivity();

    //Construct set of ghost nodes; start with a buffer with upper limit
    size_t buffer_limit = 0;
    if(num_dim==2)
    for(int ielem = 0; ielem < rnum_elem; ielem++){
      element_select->choose_2Delem_type(Element_Types(ielem), elem2D);
      buffer_limit += elem2D->num_nodes();
    }

    if(num_dim==3)
    for(int ielem = 0; ielem < rnum_elem; ielem++){
      element_select->choose_3Delem_type(Element_Types(ielem), elem);
      buffer_limit += elem->num_nodes();
    }

    CArray<size_t> ghost_node_buffer(buffer_limit);
    std::set<GO> ghost_node_set;

    //search through local elements for global node indices not owned by this MPI rank
    if(num_dim==2)
    for (int cell_rid = 0; cell_rid < rnum_elem; cell_rid++) {
      //set nodes per element
      element_select->choose_2Delem_type(Element_Types(cell_rid), elem2D);
      nodes_per_element = elem2D->num_nodes();  
      for (int node_lid = 0; node_lid < nodes_per_element; node_lid++){
        node_gid = mesh->nodes_in_cell(cell_rid, node_lid);
        if(!map->isNodeGlobalElement(node_gid)) ghost_node_set.insert(node_gid);
      }
    }

    if(num_dim==3)
    for (int cell_rid = 0; cell_rid < rnum_elem; cell_rid++) {
      //set nodes per element
      element_select->choose_3Delem_type(Element_Types(cell_rid), elem);
      nodes_per_element = elem->num_nodes();  
      for (int node_lid = 0; node_lid < nodes_per_element; node_lid++){
        node_gid = mesh->nodes_in_cell(cell_rid, node_lid);
        if(!map->isNodeGlobalElement(node_gid)) ghost_node_set.insert(node_gid);
      }
    }

    //by now the set contains, with no repeats, all the global node indices that are ghosts for this rank
    //now pass the contents of the set over to a CArray, then create a map to find local ghost indices from global ghost indices
    nghost_nodes = ghost_node_set.size();
    ghost_nodes = CArrayKokkos<GO, Kokkos::LayoutLeft, node_type::device_type>(nghost_nodes, "ghost_nodes");
    ghost_node_ranks = CArrayKokkos<int, array_layout, device_type, memory_traits>(nghost_nodes, "ghost_node_ranks");
    int ighost = 0;
    auto it = ghost_node_set.begin();
    while(it!=ghost_node_set.end()){
      ghost_nodes(ighost++) = *it;
      it++;
    }

    //debug print of ghost nodes
    //std::cout << " GHOST NODE SET ON TASK " << myrank << std::endl;
    //for(int i = 0; i < nghost_nodes; i++)
      //std::cout << "{" << i + 1 << "," << ghost_nodes(i) + 1 << "}" << std::endl;

    //create map object; note that kokkos view members in our map won't lose contents when local object is destroyed
    //since kokkos views use a reference counting mechanism to deallocate contents
    global2local_map = global_to_local_table_host_type(ghost_nodes.get_kokkos_view());

    //find which mpi rank each ghost node belongs to and store the information in a CArray
    //allocate Teuchos Views since they are the only input available at the moment in the map definitions
    ghost_nodes_pass = Teuchos::ArrayView<GO>(ghost_nodes.get_kokkos_view().data(), nghost_nodes);
    ghost_node_ranks_pass = Teuchos::ArrayView<int>(ghost_node_ranks.get_kokkos_view().data(), nghost_nodes);
    map->getRemoteIndexList(ghost_nodes_pass, ghost_node_ranks_pass);
    
    //debug print of ghost nodes
    //std::cout << " GHOST NODE MAP ON TASK " << myrank << std::endl;
    //for(int i = 0; i < nghost_nodes; i++)
      //std::cout << "{" << i + 1 << "," << global2local_map.get(ghost_nodes(i)) + 1 << "}" << std::endl;

  }

  // Close mesh input file
  if(myrank==0)
  in->close();

  //synchronize device data
  dual_node_data.sync_device();
    
  // Create reference element
  ref_elem->init(p_order, num_dim, elem->num_basis());
  std::cout<<"done with ref elem"<<std::endl;

  //communicate ghost node positions; construct multivector distributed object using local node data

  //construct array for all indices (ghost + local)
  nall_nodes = nlocal_nodes + nghost_nodes;
  CArrayKokkos<GO, array_layout, device_type, memory_traits> all_node_indices(nall_nodes, "all_node_indices");
  for(int i = 0; i < nall_nodes; i++){
    if(i<nlocal_nodes) all_node_indices(i) = map->getGlobalElement(i);
    else all_node_indices(i) = ghost_nodes(i-nlocal_nodes);
  }
  
  //debug print of node indices
  //for(int inode=0; inode < index_counter; inode++)
  //std::cout << " my_reduced_global_indices " << my_reduced_global_indices(inode) <<std::endl;
  
  // create a Map for all the node indices (ghost + local)
  all_node_map = Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(Teuchos::OrdinalTraits<GO>::invalid(),all_node_indices.get_kokkos_view(),0,comm));

  //construct dof map that follows from the all_node map (used for distributed matrix and vector objects later)
  CArrayKokkos<GO, array_layout, device_type, memory_traits> all_dof_indices(nall_nodes*num_dim, "all_dof_indices");
  for(int i = 0; i < nall_nodes; i++){
    for(int j = 0; j < num_dim; j++)
    all_dof_indices(i*num_dim + j) = all_node_map->getGlobalElement(i)*num_dim + j;
  }
  
  //pass invalid global count so the map reduces the global count automatically
  all_dof_map = Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(Teuchos::OrdinalTraits<GO>::invalid(),all_dof_indices.get_kokkos_view(),0,comm) );

  //debug print of map
  //debug print
  std::ostream &out = std::cout;
  Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
  //if(myrank==0)
  //*fos << "Ghost Node Map :" << std::endl;
  //all_node_map->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;
  //std::fflush(stdout);

  //create local dof map for multivector of local node data

  //create distributed multivector of the local node data and all (local + ghost) node storage
  node_data_distributed = Teuchos::rcp(new MV(map, dual_node_data));
  all_node_data_distributed = Teuchos::rcp(new MV(all_node_map, num_dim));
  
  //debug print
  //std::ostream &out = std::cout;
  //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
  //if(myrank==0)
  //*fos << "Node Data :" << std::endl;
  //node_data_distributed->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;
  //std::fflush(stdout);

  //create import object using local node indices map and all indices map
  Tpetra::Import<LO, GO> importer(map, all_node_map);

  //comms to get ghosts
  all_node_data_distributed->doImport(*node_data_distributed, importer, Tpetra::INSERT);

  //debug print
  //std::ostream &out = std::cout;
  //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
  //if(myrank==0)
  //*fos << "Node Data with Ghosts :" << std::endl;
  //all_node_data_distributed->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;
  //std::fflush(stdout);

  //get view of all node data (local + ghost) on the device (multivector function forces sync of dual view)

  vec_array all_node_data_device = all_node_data_distributed->getLocalView<device_type> (Tpetra::Access::ReadWrite);
  host_vec_array all_node_data_host = all_node_data_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);

  //allocate node storage with dual view
  dual_all_node_data = dual_vec_array(all_node_data_device, all_node_data_host);

  //debug print of views node indices
  //std::cout << "Local View of All Nodes on Task " << myrank <<std::endl;
  //for(int inode=0; inode < all_node_map->getNodeNumElements(); inode++){
    //std::cout << "node "<<all_node_map->getGlobalElement(inode) << " } " ;
    //std::cout << dual_all_node_data.view_host()(inode,0) << " " << dual_all_node_data.view_host()(inode,1) << " " << dual_all_node_data.view_host()(inode,2) << " " << std::endl;
  //}
     
  //std::cout << "number of patches = " << mesh->num_patches() << std::endl;
  std::cout << "End of setup " << std::endl;
} // end read_mesh

/* ----------------------------------------------------------------------
   Find boundary surface segments that belong to this MPI rank
------------------------------------------------------------------------- */

void Static_Solver_Parallel::Get_Boundary_Patches(){
  size_t npatches_repeat, npatches, element_npatches, num_nodes_in_patch, node_gid;
  int local_node_id;
  int num_dim = simparam->num_dim;
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Surface_Nodes;
  
  std::set<Node_Combination> my_patches;
  //inititializes type for the pair variable (finding the iterator type is annoying)
  std::pair<std::set<Node_Combination>::iterator, bool> current_combination;
  std::set<Node_Combination>::iterator it;
  
  //compute the number of patches in this MPI rank with repeats for adjacent cells
  npatches_repeat = 0;

  if(num_dim==2)
  for(int ielem = 0; ielem < rnum_elem; ielem++){
    element_select->choose_2Delem_type(Element_Types(ielem), elem2D);
    element_npatches = elem2D->nsurfaces;
    npatches_repeat += element_npatches;
  }
  
  else if(num_dim==3)
  for(int ielem = 0; ielem < rnum_elem; ielem++){
    element_select->choose_3Delem_type(Element_Types(ielem), elem);
    element_npatches = elem->nsurfaces;
    npatches_repeat += element_npatches;
  }
  
  //information for all patches on this rank
  CArrayKokkos<Node_Combination,array_layout, device_type, memory_traits> Patch_Nodes(npatches_repeat, "Patch_Nodes");
  CArrayKokkos<size_t,array_layout, device_type, memory_traits> Patch_Boundary_Flags(npatches_repeat, "Patch_Boundary_Flags");
  
  //initialize boundary patch flags
  for(int init = 0; init < npatches_repeat; init++)
    Patch_Boundary_Flags(init) = 1;

  //use set of nodal combinations to find boundary set
  //boundary patches will not try to add nodal combinations twice
  //loop through elements in this rank to find boundary patches
  npatches_repeat = 0;
  if(num_dim==2)
  for(int ielem = 0; ielem < rnum_elem; ielem++){
    element_select->choose_2Delem_type(Element_Types(ielem), elem2D);
    element_npatches = elem2D->nsurfaces;
    //loop through local surfaces
    for(int isurface = 0; isurface < element_npatches; isurface++){
      num_nodes_in_patch = elem2D->surface_to_dof_lid.stride(isurface);
      Surface_Nodes = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(num_nodes_in_patch, "Surface_Nodes");
      for(int inode = 0; inode < num_nodes_in_patch; inode++){
        local_node_id = elem2D->surface_to_dof_lid(isurface,inode);
        Surface_Nodes(inode) = mesh->nodes_in_cell(ielem, local_node_id);
      }
      Node_Combination temp(Surface_Nodes);
      //construct Node Combination object for this surface
      Patch_Nodes(npatches_repeat) = temp;
      Patch_Nodes(npatches_repeat).patch_id = npatches_repeat;
      Patch_Nodes(npatches_repeat).element_id = ielem;
      Patch_Nodes(npatches_repeat).local_patch_id = isurface;
      //test if this patch has already been added; if yes set boundary flags to 0
      current_combination = my_patches.insert(Patch_Nodes(npatches_repeat));
      //if the set determines this is a duplicate access the original element's patch id and set flag to 0
      if(current_combination.second==false){
      //set original element flag to 0
      Patch_Boundary_Flags((*current_combination.first).patch_id) = 0;
      //set this current flag to 0 for the duplicate as well
      Patch_Boundary_Flags(npatches_repeat) = 0;
      npatches_repeat++;
      }

    }
  }

  if(num_dim==3)
  for(int ielem = 0; ielem < rnum_elem; ielem++){
    element_select->choose_3Delem_type(Element_Types(ielem), elem);
    element_npatches = elem->nsurfaces;
    //loop through local surfaces
    for(int isurface = 0; isurface < element_npatches; isurface++){
      num_nodes_in_patch = elem->surface_to_dof_lid.stride(isurface);
      //debug print
      //std::cout << "NUMBER OF PATCH NODES FOR ELEMENT " << ielem+1 << " ON LOCAL SURFACE " << isurface+1 << " IS " << elem->surface_to_dof_lid.start_index_[isurface+1] << std::endl;
      Surface_Nodes = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(num_nodes_in_patch, "Surface_Nodes");
      for(int inode = 0; inode < num_nodes_in_patch; inode++){
        local_node_id = elem->surface_to_dof_lid(isurface,inode);
        Surface_Nodes(inode) = mesh->nodes_in_cell(ielem, local_node_id);
      }
      Node_Combination temp(Surface_Nodes);
      //construct Node Combination object for this surface
      Patch_Nodes(npatches_repeat) = temp;
      Patch_Nodes(npatches_repeat).patch_id = npatches_repeat;
      Patch_Nodes(npatches_repeat).element_id = ielem;
      Patch_Nodes(npatches_repeat).local_patch_id = isurface;
      //test if this patch has already been added; if yes set boundary flags to 0
      current_combination = my_patches.insert(Patch_Nodes(npatches_repeat));
      //if the set determines this is a duplicate access the original element's patch id and set flag to 0
      if(current_combination.second==false){
        //set original element flag to 0
        Patch_Boundary_Flags((*current_combination.first).patch_id) = 0;
        //set this current flag to 0 for the duplicate as well
        Patch_Boundary_Flags(npatches_repeat) = 0;
      }
      npatches_repeat++;
    }
  }

  //debug print of all patches
  /*
  std::cout << " ALL PATCHES " << npatches_repeat <<std::endl;
  for(int iprint = 0; iprint < npatches_repeat; iprint++){
    std::cout << "Patch " << iprint + 1 << " ";
    for(int j = 0; j < Patch_Nodes(iprint).node_set.size(); j++)
      std::cout << Patch_Nodes(iprint).node_set(j) << " ";
    std::cout << std::endl;
  }
  */

  //loop through patch boundary flags to isolate boundary patches
  nboundary_patches = 0;
  for(int iflags = 0 ; iflags < npatches_repeat; iflags++){
    if(Patch_Boundary_Flags(iflags)) nboundary_patches++;
  }
  //upper bound that is not much larger
  Boundary_Patches = CArrayKokkos<Node_Combination, array_layout, device_type, memory_traits>(nboundary_patches, "Boundary_Patches");
  nboundary_patches = 0;
  bool my_rank_flag;
  size_t remote_count;
  for(int ipatch = 0 ; ipatch < npatches_repeat; ipatch++){
    if(Patch_Boundary_Flags(ipatch)){
      /*check if Nodes in the combination for this patch belong to this MPI rank.
        If all are local then this is a boundary patch belonging to this rank.
        If all nodes are remote then another rank must decide if that patch is a boundary.
        If only a subset of the nodes are local it must be a boundary patch; this
        case assigns the patch to the lowest mpi rank index the nodes in this patch belong to */
      num_nodes_in_patch = Patch_Nodes(ipatch).node_set.size();
      my_rank_flag = true;
      remote_count = 0;

      //assign as a local boundary patch if any of the nodes on the patch are local
      //only the local nodes on the patch will contribute to the equation assembly on this rank
      for(int inode = 0; inode < num_nodes_in_patch; inode++){
        node_gid = Patch_Nodes(ipatch).node_set(inode);
        if(!map->isNodeGlobalElement(node_gid)){
          //if(ghost_node_ranks(global2local_map.get(node_gid))<myrank)
          //my_rank_flag = false;
          remote_count++;
          //test
        } 
      }
      //all nodes were remote
      //if(remote_count == num_nodes_in_patch) my_rank_flag = false;

      //if all nodes were not local
      if(my_rank_flag)
        Boundary_Patches(nboundary_patches++) = Patch_Nodes(ipatch);
    }
  }

  //debug print of boundary patches
  /*std::cout << " BOUNDARY PATCHES ON TASK " << myrank << " = " << nboundary_patches <<std::endl;
  for(int iprint = 0; iprint < nboundary_patches; iprint++){
    std::cout << "Patch " << iprint + 1 << " ";
    for(int j = 0; j < Boundary_Patches(iprint).node_set.size(); j++)
      std::cout << Boundary_Patches(iprint).node_set(j) << " ";
    std::cout << std::endl;
  }
  std::fflush(stdout);
  */
}

/* ----------------------------------------------------------------------
   Assign sets of element boundary surfaces corresponding to user BCs
------------------------------------------------------------------------- */

void Static_Solver_Parallel::generate_bcs(){
    
    // build boundary mesh patches
    //mesh->build_bdy_patches();
    Get_Boundary_Patches();

    std::cout << "number of boundary patches on task " << myrank << " = " << nboundary_patches << std::endl;
    std::cout << "building boundary sets " << std::endl;
    // set the number of boundary sets
    
    int num_boundary_sets = simparam->NB;
    int num_surface_force_sets = simparam->NBSF;
    int num_surface_disp_sets = simparam->NBD;
    int num_dim = simparam->num_dim;
    int current_bdy_id = 0;
    int bdy_set_id;
    int surf_force_set_id = 0;
    int surf_disp_set_id = 0;

    init_boundary_sets(num_boundary_sets);
    Boundary_Condition_Type_List = CArray<int>(num_boundary_sets); 
    Boundary_Surface_Force_Densities = CArray<real_t>(num_surface_force_sets,3);
    Boundary_Surface_Displacements = CArray<real_t>(num_surface_disp_sets,3);
    //initialize
    for(int ibdy=0; ibdy < num_boundary_sets; ibdy++) Boundary_Condition_Type_List(ibdy) = NONE;
    
    // tag the x=0 plane,  (Direction, value, bdy_set)
    std::cout << "tagging x = 0 " << std::endl;
    int bc_tag = 0;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
    real_t value = 0.0;
    bdy_set_id = current_bdy_id++;
    tag_boundaries(bc_tag, value, bdy_set_id);
    Boundary_Condition_Type_List(bdy_set_id) = DISPLACEMENT_CONDITION;
    Boundary_Surface_Displacements(surf_disp_set_id,0) = 0;
    Boundary_Surface_Displacements(surf_disp_set_id,1) = 0;
    Boundary_Surface_Displacements(surf_disp_set_id,2) = 0;
    surf_disp_set_id++;
    
    std::cout << "tagged a set " << std::endl;
    std::cout << "number of bdy patches in this set = " << NBoundary_Condition_Patches(bdy_set_id) << std::endl;
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

    std::cout << "tagging x = 1.2 " << std::endl;
    bc_tag = 0;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
    value = 1.2;
    //value = 2;
    bdy_set_id = current_bdy_id++;
    //find boundary patches this BC corresponds to
    tag_boundaries(bc_tag, value, bdy_set_id);
    Boundary_Condition_Type_List(bdy_set_id) = LOADING_CONDITION;
    Boundary_Surface_Force_Densities(surf_force_set_id,0) = 2;
    Boundary_Surface_Force_Densities(surf_force_set_id,1) = 0;
    Boundary_Surface_Force_Densities(surf_force_set_id,2) = 0;
    surf_force_set_id++;
    std::cout << "tagged a set " << std::endl;
    std::cout << "number of bdy patches in this set = " << NBoundary_Condition_Patches(bdy_set_id) << std::endl;
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
    Node_DOF_Boundary_Condition_Type = CArrayKokkos<int, array_layout, device_type, memory_traits>(nall_nodes*num_dim, "Node_DOF_Boundary_Condition_Type");
    Node_DOF_Displacement_Boundary_Conditions = CArray<real_t>(nall_nodes*num_dim);
    Node_DOF_Force_Boundary_Conditions = CArray<real_t>(nall_nodes*num_dim);

    //initialize
    for(int init=0; init < nall_nodes*num_dim; init++)
      Node_DOF_Boundary_Condition_Type(init) = NONE;
} // end generate_bcs

/* ----------------------------------------------------------------------
   initialize storage for element boundary surfaces corresponding to user BCs
------------------------------------------------------------------------- */

void Static_Solver_Parallel::init_boundary_sets (int num_sets){
    
  num_boundary_conditions = num_sets;
  if(num_sets == 0){
    std::cout << " Warning: number of boundary conditions = 0";
    return;
  }
  Boundary_Condition_Patches_strides = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(num_sets, "Boundary_Condition_Patches_strides");
  NBoundary_Condition_Patches = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(num_sets, "NBoundary_Condition_Patches");
  Boundary_Condition_Patches = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(num_sets, nboundary_patches, "Boundary_Condition_Patches");

  //initialize data
  for(int iset = 0; iset < num_sets; iset++) NBoundary_Condition_Patches(iset) = 0;
}

/* ----------------------------------------------------------------------
   find which boundary patches correspond to the given BC.
   bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
   val = plane value, cylinder radius, shell radius
------------------------------------------------------------------------- */

void Static_Solver_Parallel::tag_boundaries(int bc_tag, real_t val, int bdy_set){
  
  int num_boundary_sets = simparam->NB;
  int is_on_set;
  /*
  if (bdy_set == num_bdy_sets_){
    std::cout << " ERROR: number of boundary sets must be increased by "
      << bdy_set-num_bdy_sets_+1 << std::endl;
    exit(0);
  }
  */
    
  // save the boundary vertices to this set that are on the plane
  int counter = 0;
  for (int iboundary_patch = 0; iboundary_patch < nboundary_patches; iboundary_patch++) {

    // check to see if this patch is on the specified plane
    is_on_set = check_boundary(Boundary_Patches(iboundary_patch), bc_tag, val); // no=0, yes=1
        
    if (is_on_set == 1){
      Boundary_Condition_Patches(bdy_set,counter) = iboundary_patch;
      counter ++;
    }
  } // end for bdy_patch
    
  // save the number of bdy patches in the set
  NBoundary_Condition_Patches(bdy_set) = counter;
    
  std::cout << " tagged boundary patches " << std::endl;
}

/* ----------------------------------------------------------------------
   routine for checking to see if a patch is on a boundary set
   bc_tag = 0 xplane, 1 yplane, 3 zplane, 4 cylinder, 5 is shell
   val = plane value, radius, radius
------------------------------------------------------------------------- */


int Static_Solver_Parallel::check_boundary(Node_Combination &Patch_Nodes, int bc_tag, real_t val){
  
  int is_on_set = 1;
  vec_array all_node_data = dual_all_node_data.view_device();

  //Nodes on the Patch
  auto node_list = Patch_Nodes.node_set;
  size_t nnodes = node_list.size();
  size_t node_rid;
  real_t node_coord;
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> node_on_flags(nnodes, "node_on_flags");

  //initialize
  for(int inode = 0; inode < nnodes; inode++) node_on_flags(inode) = 0;
    
  //test for planes
  if(bc_tag < 3)
  for(int inode = 0; inode < nnodes; inode++){

    node_rid = all_node_map->getLocalElement(node_list(inode));
    node_coord = all_node_data(node_rid,bc_tag);
    if ( fabs(node_coord - val) <= 1.0e-8 ) node_on_flags(inode) = 1;
    //debug print of node id and node coord
    //std::cout << "node coords on task " << myrank << " for node " << node_rid << std::endl;
    //std::cout << "coord " <<node_coord << " flag " << node_on_flags(inode) << " bc_tag " << bc_tag << std::endl;
  }
    
    /*
    // cylinderical shell where radius = sqrt(x^2 + y^2)
    else if (this_bc_tag == 3){
        
        real_t R = sqrt(these_patch_coords[0]*these_patch_coords[0] +
                        these_patch_coords[1]*these_patch_coords[1]);
        
        if ( fabs(R - val) <= 1.0e-8 ) is_on_bdy = 1;

        
    }// end if on type
    
    // spherical shell where radius = sqrt(x^2 + y^2 + z^2)
    else if (this_bc_tag == 4){
        
        real_t R = sqrt(these_patch_coords[0]*these_patch_coords[0] +
                        these_patch_coords[1]*these_patch_coords[1] +
                        these_patch_coords[2]*these_patch_coords[2]);
        
        if ( fabs(R - val) <= 1.0e-8 ) is_on_bdy = 1;
        
    } // end if on type
    */
    //check if all nodes lie on the boundary set
  for(int inode = 0; inode < nnodes; inode++)
    if(!node_on_flags(inode)) is_on_set = 0;
  
  //debug print of return flag
  //std::cout << "patch flag on task " << myrank << " is " << is_on_set << std::endl;
  return is_on_set;
    
} // end method to check bdy

void Static_Solver_Parallel::allocate_state(){
    //local variable for host view in the dual view
    host_vec_array node_data = dual_node_data.view_host();
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


void Static_Solver_Parallel::initialize_state(){
  //local variable for host view in the dual view
  host_vec_array node_data = dual_node_data.view_host();
    int NF = simparam->NF;
    mat_pt_t *mat_pt = simparam->mat_pt;
    mat_fill_t *mat_fill = simparam->mat_fill;
    int rk_stage = simparam->rk_stage;

    std::cout << "Before fill instructions"  << std::endl;
    //--- apply the fill instructions ---//
    for (int f_id = 0; f_id < NF; f_id++){
        
        for (int cell_rid = 0; cell_rid < mesh->num_cells(); cell_rid++) {
            
            // calculate the coordinates and radius of the cell
            real_t cell_coords_x = 0.0;
            real_t cell_coords_y = 0.0;
            real_t cell_coords_z = 0.0;
            
            for (int node_lid = 0; node_lid < 8; node_lid++){
                
                // increment the number of cells attached to this vertex
                int vert_gid = mesh->nodes_in_cell(cell_rid, node_lid); // get the global_id
                
                cell_coords_x += node_data(vert_gid, 0);
                cell_coords_y += node_data(vert_gid, 1);
                cell_coords_z += node_data(vert_gid, 2);
                
            }// end for loop over node_lid
            
            cell_coords_x = cell_coords_x/8.0;
            cell_coords_y = cell_coords_y/8.0;
            cell_coords_z = cell_coords_z/8.0;

            // Material points at cell center

            mat_pt->coords(rk_stage, cell_rid, 0) = cell_coords_x;
            mat_pt->coords(rk_stage, cell_rid, 1) = cell_coords_y;
            mat_pt->coords(rk_stage, cell_rid, 2) = cell_coords_z;

            
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
                
                mat_pt->field(cell_rid) = mat_fill[f_id].field2;


            } // end if fill volume
        } // end for cell loop
    } // end for fills

    std::cout << "After fill instructions"  << std::endl;
} // end initialize_state


void Static_Solver_Parallel::calculate_ref_elem(){
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

void Static_Solver_Parallel::apply_boundary(){

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


void Static_Solver_Parallel::smooth_cells(){
    mat_pt_t *mat_pt = simparam->mat_pt;
    mat_fill_t *mat_fill = simparam->mat_fill;
    node_t *node = simparam->node;

    // Walk over cells 
    for(int cell_rid = 0; cell_rid < mesh->num_cells(); cell_rid++){

        // Temporary holder variable
        real_t temp_avg = 0.0;

        // Walk over nodes in cell and calculate average
        for(int node_lid = 0; node_lid < mesh->num_nodes_in_cell(); node_lid++){

            // Get global index of the node
            int node_gid = mesh->nodes_in_cell(cell_rid, node_lid);

            temp_avg += node->field(node_gid)/8.0;

        }

        // Save average to material point at cell center
        mat_pt->field(cell_rid) = temp_avg;

    } // end of loop over cells

    // Walk over all the nodes
    for(int node_gid = 0; node_gid < mesh->num_nodes(); node_gid++){

        real_t temp_avg = 0.0;
        
        // Walk over all the cells connected to the node and average values
        for(int cell_lid = 0; cell_lid < mesh->num_cells_in_node(node_gid); cell_lid++){

            // Get global index of the cell
            int cell_rid = mesh->cells_in_node(node_gid, cell_lid);

            temp_avg += mat_pt->field(cell_rid)/ (real_t)mesh->num_cells_in_node(node_gid);
        
        }

        // Save average to the node
        node->field(node_gid) = temp_avg;

    } // end of loop over nodes
}


void Static_Solver_Parallel::smooth_element(){
    mat_pt_t *mat_pt = simparam->mat_pt;
    mat_fill_t *mat_fill = simparam->mat_fill;
    node_t *node = simparam->node;

    // Walk over each element in the mesh
    for(int elem_gid = 0; elem_gid < mesh->num_elems(); elem_gid++){

        // Walk over each cell in the element
        for(int cell_lid = 0; cell_lid < mesh->num_cells_in_elem(); cell_lid++){

            // Get the global ID of the cell
            int cell_rid = mesh->cells_in_elem(elem_gid, cell_lid);

            real_t temp_avg = 0.0;

            // Loop over nodes in the cell
            for(int node_lid = 0; node_lid < mesh->num_nodes_in_cell(); node_lid++){

                // Get global ID for this node
                int node_gid = mesh->nodes_in_cell(cell_rid, node_lid);

                temp_avg += node->field(node_gid)/mesh->num_nodes_in_cell();

            }

            // Save averaged values to cell centered material point
            mat_pt->field(cell_rid) = temp_avg;
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
                int cell_rid = mesh->cells_in_node(node_gid, cell_lid);

                temp_avg += mat_pt->field(cell_rid)/ (real_t)mesh->num_cells_in_node(node_gid);
            }

            // Save averaged field to node
            node->field(node_gid) = temp_avg;

        }// end loop over nodes
    }// end loop over elements
}

void Static_Solver_Parallel::get_nodal_jacobian(){
  //local variable for host view in the dual view
  host_vec_array node_data = dual_node_data.view_host();
  mat_pt_t *mat_pt = simparam->mat_pt;
  mat_fill_t *mat_fill = simparam->mat_fill;
  int num_dim = simparam->num_dim;

  // loop over the mesh

    for(int elem_gid = 0; elem_gid < mesh->num_elems(); elem_gid++){
        
        for(int cell_lid = 0; cell_lid < mesh->num_cells_in_elem(); cell_lid++){ // 1 for linear elements

            int cell_rid = mesh->cells_in_elem(elem_gid, cell_lid);

            // Loop over nodes in cell and initialize jacobian to zero
            for(int node_lid = 0; node_lid < mesh->num_nodes_in_cell(); node_lid++){

                int node_gid = mesh->nodes_in_cell(cell_rid, node_lid);
                
                for(int dim_i = 0; dim_i < mesh->num_dim(); dim_i++){
                    for(int dim_j = 0; dim_j < mesh->num_dim(); dim_j++){

                        mesh->node_jacobian(node_gid, dim_i, dim_j) = 0.0;
                    }
                }
            }

            // Calculate the actual jacobian for that node
            for(int node_lid = 0; node_lid < mesh->num_nodes_in_cell(); node_lid++){
                
                int node_gid = mesh->nodes_in_cell(cell_rid, node_lid);



                for(int dim_i = 0; dim_i < mesh->num_dim(); dim_i++){
                    for(int dim_j = 0; dim_j < mesh->num_dim(); dim_j++){

                        // Sum over the basis functions and nodes where they are defined
                        for(int basis_id = 0; basis_id < mesh->num_nodes_in_cell(); basis_id++){

                            int ref_node_gid = mesh->nodes_in_cell(cell_rid, basis_id);

                            mesh->node_jacobian(node_gid, dim_i, dim_j) += 
                                node_data(ref_node_gid , dim_i) * ref_elem->ref_nodal_gradient(node_lid, basis_id, dim_j);
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

    for(int cell_rid = 0; cell_rid < mesh->num_cells(); cell_rid++){

        mat_pt->volume(cell_rid) = 0;

        for(int node_lid = 0; node_lid < mesh->num_nodes_in_cell(); node_lid++){

            int node_gid = mesh->nodes_in_cell(cell_rid, node_lid);

            mat_pt->volume(cell_rid) += mesh->node_det_j(node_gid);
        }

        std::cout<< "Volume for cell  "<< cell_rid << " = "<< mat_pt->volume(cell_rid) << std::endl;
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

void Static_Solver_Parallel::vtk_writer(){
    //local variable for host view in the dual view
    host_vec_array node_data = dual_node_data.view_host();
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
        "cell_rid"
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
        
        for (int cell_rid = 0; cell_rid < num_cells; cell_rid++){
            
            fprintf(out[0],"%f\n",(float) cell_rid);

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
        
        fprintf(out[0],"%f   %f   %f   \n",node_data(node_gid, 0),
                                           node_data(node_gid, 1),
                                           node_data(node_gid, 2));

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
    for (int cell_rid = 0; cell_rid < num_cells; cell_rid++){
        
        for(int node = 0; node < 8; node++){
            fprintf(out[0],"%i  ", mesh->nodes_in_cell(cell_rid, node));
        }

        fprintf(out[0],"\n");

    } // end for k over cells
    fprintf(out[0],"</DataArray> \n");
    fprintf(out[0],"\n");


    fprintf(out[0],"<DataArray type=\"Int64\" Name=\"offsets\">\n");
    for (int cell_rid = 0; cell_rid < num_cells; cell_rid++){
        fprintf(out[0],"%i  \n", 8*(cell_rid+1));
    } // end for k over cells
    fprintf(out[0],"</DataArray> \n");
    fprintf(out[0],"\n");


    fprintf(out[0],"<DataArray type=\"UInt64\" Name=\"types\">\n");
    for (int cell_rid = 0; cell_rid < num_cells; cell_rid++){
        fprintf(out[0],"%i  \n", 42);
    } // end for k over cells
    fprintf(out[0],"</DataArray> \n");
    fprintf(out[0],"\n");


    fprintf(out[0],"<DataArray type=\"Int64\" Name=\"faces\">\n");
    
    for (int cell_rid = 0; cell_rid < num_cells; cell_rid++){
        fprintf(out[0],"%i  \n", 6);

        for(int patch_lid = 0; patch_lid < 6; patch_lid++){

            fprintf(out[0],"4  ");
            for(int node_lid = 0; node_lid < 4; node_lid++){
                fprintf(out[0],"%i  ", mesh->node_in_patch_in_cell(cell_rid, patch_lid, node_lid));
            }

            fprintf(out[0],"\n");
        }


    } // end for k over cells
    fprintf(out[0],"</DataArray> \n");
    fprintf(out[0],"\n");


    int faceoffsets = 31;
    fprintf(out[0],"<DataArray type=\"Int64\" Name=\"faceoffsets\">\n");
    for (int cell_rid = 0; cell_rid < num_cells; cell_rid++){
        fprintf(out[0],"%i  \n", faceoffsets*(cell_rid+1));
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


void Static_Solver_Parallel::ensight(){
    //local variable for host view in the dual view
    host_vec_array node_data = dual_node_data.view_host();
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


    for (int cell_rid=0; cell_rid<num_cells; cell_rid++){
        cell_fields(cell_rid, 0) = mat_pt->field(cell_rid);    
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
    for (int cell_rid = 0; cell_rid < num_cells; cell_rid++){

        cell_fields(cell_rid, 3) = mat_pt->field(cell_rid);    

    } // end for k over cells


    // save the vertex vector fields to an array for exporting to graphics files
    auto vec_fields = CArray <real_t> (num_nodes, num_vec_vars, 3);

    for (int node_gid = 0; node_gid < num_nodes; node_gid++){
        
        vec_fields(node_gid, 0, 0) = node_data(node_gid, 0); 
        vec_fields(node_gid, 0, 1) = node_data(node_gid, 1);
        vec_fields(node_gid, 0, 2) = node_data(node_gid, 2);

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
        fprintf(out[0],"%12.5e\n",node_data(node_gid, 0));
    }
    
    for (int node_gid=0; node_gid<num_nodes; node_gid++){
        fprintf(out[0],"%12.5e\n",node_data(node_gid, 1));
    }
    
    for (int node_gid=0; node_gid<num_nodes; node_gid++){
        fprintf(out[0],"%12.5e\n",node_data(node_gid, 2));
    }
    
    // convert_vert_list_ord_Ensight
    // --- cells ---
    fprintf(out[0],"hexa8\n");
    fprintf(out[0],"%10d\n",num_cells);
    int this_index;
    
    // write all global point numbers for this cell
    for (int cell_rid = 0; cell_rid<num_cells; cell_rid++) {

        for (int j=0; j<8; j++){
            this_index = convert_vert_list_ord_Ensight(j);
            fprintf(out[0],"%10d\t",mesh->nodes_in_cell(cell_rid, this_index)+1); // note node_gid starts at 1

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
        
        for (int cell_rid=0; cell_rid<num_cells; cell_rid++) {
            fprintf(out[0],"%12.5e\n", cell_fields(cell_rid, var));
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

void Static_Solver_Parallel::init_global(){

  int num_dim = simparam->num_dim;
  int num_elems = mesh->num_elems();
  int num_nodes = mesh->num_nodes();
  Stiffness_Matrix_Strides = CArrayKokkos<size_t, array_layout, device_type, memory_traits> (nlocal_nodes*num_dim, "Stiffness_Matrix_Strides");
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Graph_Fill(nall_nodes, "nall_nodes");
  //CArray <int> nodes_in_cell_list_ = mesh->nodes_in_cell_list_;
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> current_row_nodes_scanned;
  int current_row_n_nodes_scanned;
  int local_node_index, global_node_index, current_column_index;
  int max_stride = 0;
  size_t nodes_per_element;
  
  //allocate stride arrays
  CArrayKokkos <size_t, array_layout, device_type, memory_traits> Graph_Matrix_Strides_initial(nlocal_nodes, "Graph_Matrix_Strides_initial");
  Graph_Matrix_Strides = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(nlocal_nodes, "Graph_Matrix_Strides");

  //allocate storage for the sparse stiffness matrix map used in the assembly process
  Global_Stiffness_Matrix_Assembly_Map = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(rnum_elem,
                                         max_nodes_per_element,max_nodes_per_element, "Global_Stiffness_Matrix_Assembly_Map");

  //allocate array used to determine global node repeats in the sparse graph later
  CArrayKokkos <int, array_layout, device_type, memory_traits> node_indices_used(nall_nodes, "node_indices_used");

  /*allocate array that stores which column the node index occured on for the current row
    when removing repeats*/
  CArrayKokkos <size_t, array_layout, device_type, memory_traits> column_index(nall_nodes, "column_index");
  
  //initialize nlocal arrays
  for(int inode = 0; inode < nlocal_nodes; inode++){
    Graph_Matrix_Strides_initial(inode) = 0;
    Graph_Matrix_Strides(inode) = 0;
    Graph_Fill(inode) = 0;
  }

  //initialize nall arrays
  for(int inode = 0; inode < nall_nodes; inode++){
    node_indices_used(inode) = 0;
    column_index(inode) = 0;
  }
  
  //count upper bound of strides for Sparse Pattern Graph by allowing repeats due to connectivity
  if(num_dim == 2)
  for (int ielem = 0; ielem < rnum_elem; ielem++){
    element_select->choose_2Delem_type(Element_Types(ielem), elem2D);
    nodes_per_element = elem2D->num_nodes();
    for (int lnode = 0; lnode < nodes_per_element; lnode++){
      global_node_index = mesh->nodes_in_cell_list_(ielem, lnode);
      if(map->isNodeGlobalElement(global_node_index)){
        local_node_index = map->getLocalElement(global_node_index);
        Graph_Matrix_Strides_initial(local_node_index) += nodes_per_element;
      }
    }
  }

  if(num_dim == 3)
  for (int ielem = 0; ielem < rnum_elem; ielem++){
    element_select->choose_3Delem_type(Element_Types(ielem), elem);
    nodes_per_element = elem->num_nodes();
    for (int lnode = 0; lnode < nodes_per_element; lnode++){
      global_node_index = mesh->nodes_in_cell_list_(ielem, lnode);
      if(map->isNodeGlobalElement(global_node_index)){
        local_node_index = map->getLocalElement(global_node_index);
        Graph_Matrix_Strides_initial(local_node_index) += nodes_per_element;
      }
    }
  }
  
  //equate strides for later
  for(int inode = 0; inode < nlocal_nodes; inode++)
    Graph_Matrix_Strides(inode) = Graph_Matrix_Strides_initial(inode);
  
  //compute maximum stride
  for(int inode = 0; inode < nlocal_nodes; inode++)
    if(Graph_Matrix_Strides_initial(inode) > max_stride) max_stride = Graph_Matrix_Strides_initial(inode);
  
  //allocate array used in the repeat removal process
  current_row_nodes_scanned = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(max_stride, "current_row_nodes_scanned");

  //allocate sparse graph with node repeats
  RaggedRightArrayKokkos<size_t, array_layout, device_type, memory_traits> Repeat_Graph_Matrix(Graph_Matrix_Strides_initial);
  RaggedRightArrayofVectorsKokkos<size_t, array_layout, device_type, memory_traits> Element_local_indices(Graph_Matrix_Strides_initial,num_dim);
  
  //Fill the initial Graph with repeats
  if(num_dim == 2)
  for (int ielem = 0; ielem < rnum_elem; ielem++){
    element_select->choose_2Delem_type(Element_Types(ielem), elem2D);
    nodes_per_element = elem2D->num_nodes();
    for (int lnode = 0; lnode < nodes_per_element; lnode++){
      global_node_index = mesh->nodes_in_cell_list_(ielem, lnode);
      if(map->isNodeGlobalElement(global_node_index)){
        local_node_index = map->getLocalElement(global_node_index);
        for (int jnode = 0; jnode < nodes_per_element; jnode++){
          current_column_index = Graph_Fill(local_node_index)+jnode;
          Repeat_Graph_Matrix(local_node_index, current_column_index) = mesh->nodes_in_cell_list_(ielem,jnode);

          //fill inverse map
          Element_local_indices(local_node_index,current_column_index,0) = ielem;
          Element_local_indices(local_node_index,current_column_index,1) = lnode;
          Element_local_indices(local_node_index,current_column_index,2) = jnode;

          //fill forward map
          Global_Stiffness_Matrix_Assembly_Map(ielem,lnode,jnode) = current_column_index;
        }
        Graph_Fill(local_node_index) += nodes_per_element;
      }
    }
  }
  
  if(num_dim == 3)
  for (int ielem = 0; ielem < rnum_elem; ielem++){
    element_select->choose_3Delem_type(Element_Types(ielem), elem);
    nodes_per_element = elem->num_nodes();
    for (int lnode = 0; lnode < nodes_per_element; lnode++){
      global_node_index = mesh->nodes_in_cell_list_(ielem, lnode);
      if(map->isNodeGlobalElement(global_node_index)){
        local_node_index = map->getLocalElement(global_node_index);
        for (int jnode = 0; jnode < nodes_per_element; jnode++){
          current_column_index = Graph_Fill(local_node_index)+jnode;
          Repeat_Graph_Matrix(local_node_index, current_column_index) = mesh->nodes_in_cell_list_(ielem,jnode);

          //fill inverse map
          Element_local_indices(local_node_index,current_column_index,0) = ielem;
          Element_local_indices(local_node_index,current_column_index,1) = lnode;
          Element_local_indices(local_node_index,current_column_index,2) = jnode;

          //fill forward map
          Global_Stiffness_Matrix_Assembly_Map(ielem,lnode,jnode) = current_column_index;
        }
        Graph_Fill(local_node_index) += nodes_per_element;
      }
    }
  }
  
  //debug statement
  //std::cout << "started run" << std::endl;
  //std::cout << "Graph Matrix Strides Repeat on task " << myrank << std::endl;
  //for (int inode = 0; inode < nlocal_nodes; inode++)
    //std::cout << Graph_Matrix_Strides(inode) << std::endl;
  
  //remove repeats from the inital graph setup
  int current_node, current_element_index, element_row_index, element_column_index, current_stride;
  for (int inode = 0; inode < nlocal_nodes; inode++){
    current_row_n_nodes_scanned = 0;
    for (int istride = 0; istride < Graph_Matrix_Strides(inode); istride++){
      //convert global index in graph to its local index for the flagging array
      current_node = all_node_map->getLocalElement(Repeat_Graph_Matrix(inode,istride));
      //debug
      //if(current_node==-1)
      //std::cout << "Graph Matrix node access on task " << myrank << std::endl;
      //std::cout << Repeat_Graph_Matrix(inode,istride) << std::endl;
      if(node_indices_used(current_node)){
        //set global assembly map index to the location in the graph matrix where this global node was first found
        current_element_index = Element_local_indices(inode,istride,0);
        element_row_index = Element_local_indices(inode,istride,1);
        element_column_index = Element_local_indices(inode,istride,2);
        Global_Stiffness_Matrix_Assembly_Map(current_element_index,element_row_index, element_column_index) 
            = column_index(current_node);   

        
        //swap current node with the end of the current row and shorten the stride of the row
        //first swap information about the inverse and forward maps

        current_stride = Graph_Matrix_Strides(inode);
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
        Graph_Matrix_Strides(inode)--;
      }
      else{
        /*this node hasn't shown up in the row before; add it to the list of nodes
          that have been scanned uniquely. Use this list to reset the flag array
          afterwards without having to loop over all the nodes in the system*/
        node_indices_used(current_node) = 1;
        column_index(current_node) = istride;
        current_row_nodes_scanned(current_row_n_nodes_scanned) = current_node;
        current_row_n_nodes_scanned++;
      }
    }
    //reset nodes used list for the next row of the sparse list
    for(int node_reset = 0; node_reset < current_row_n_nodes_scanned; node_reset++)
      node_indices_used(current_row_nodes_scanned(node_reset)) = 0;

  }

  //copy reduced content to non_repeat storage
  Graph_Matrix = RaggedRightArrayKokkos<size_t, array_layout, device_type, memory_traits>(Graph_Matrix_Strides);
  for(int inode = 0; inode < nlocal_nodes; inode++)
    for(int istride = 0; istride < Graph_Matrix_Strides(inode); istride++)
      Graph_Matrix(inode,istride) = Repeat_Graph_Matrix(inode,istride);

  //deallocate repeat matrix
  
  /*At this stage the sparse graph should have unique global indices on each row.
    The constructed Assembly map (to the global sparse matrix)
    is used to loop over each element's local stiffness matrix in the assembly process.*/
  
  //expand strides for stiffness matrix by multipling by dim
  for(int inode = 0; inode < nlocal_nodes; inode++){
    for (int idim = 0; idim < num_dim; idim++)
    Stiffness_Matrix_Strides(num_dim*inode + idim) = num_dim*Graph_Matrix_Strides(inode);
  }

  Stiffness_Matrix = RaggedRightArrayKokkos<real_t, Kokkos::LayoutRight, device_type, memory_traits, array_layout>(Stiffness_Matrix_Strides);
  DOF_Graph_Matrix = RaggedRightArrayKokkos<GO, array_layout, device_type, memory_traits> (Stiffness_Matrix_Strides);

  //initialize stiffness Matrix entries to 0
  //debug print
    //std::cout << "DOF GRAPH MATRIX ENTRIES ON TASK " << myrank << std::endl;
  for (int idof = 0; idof < num_dim*nlocal_nodes; idof++){
    for (int istride = 0; istride < Stiffness_Matrix_Strides(idof); istride++){
      Stiffness_Matrix(idof,istride) = 0;
      DOF_Graph_Matrix(idof,istride) = Graph_Matrix(idof/num_dim,istride/num_dim)*num_dim + istride%num_dim;
      //debug print
      //std::cout << "{" <<istride + 1 << "," << DOF_Graph_Matrix(idof,istride) << "} ";
    }
    //debug print
    //std::cout << std::endl;
  }
  
  /*
  //debug print nodal positions and indices
  std::cout << " ------------NODAL POSITIONS--------------"<<std::endl;
  for (int inode = 0; inode < num_nodes; inode++){
      std::cout << "node: " << inode + 1 << " { ";
    for (int istride = 0; istride < num_dim; istride++){
        std::cout << node_data(inode,istride) << " , ";
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
    for (int istride = 0; istride < Graph_Matrix_Strides(inode); istride++){
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

void Static_Solver_Parallel::assemble(){
  int num_dim = simparam->num_dim;
  int num_elems = mesh->num_elems();
  int nodes_per_element;
  int current_row_n_nodes_scanned;
  int local_dof_index, global_node_index, current_row, current_column;
  int max_stride = 0;
  //local variable for host view in the dual view
  host_vec_array Nodal_Forces = dual_nodal_forces.view_host();
  //signal modification of the host view of nodal forces
  dual_nodal_forces.modify_host();
  CArray <real_t> Local_Stiffness_Matrix = CArray <real_t> (num_dim*max_nodes_per_element,num_dim*max_nodes_per_element);

  //assemble the global stiffness matrix
  if(num_dim==2)
  for (int ielem = 0; ielem < num_elems; ielem++){
    element_select->choose_2Delem_type(Element_Types(ielem), elem2D);
    nodes_per_element = elem2D->num_nodes();
    //construct local stiffness matrix for this element
    local_matrix_multiply(ielem, Local_Stiffness_Matrix);
    //assign entries of this local matrix to the sparse global matrix storage;
    for (int inode = 0; inode < nodes_per_element; inode++){
      //see if this node is local
      global_node_index = mesh->nodes_in_cell_list_(ielem,inode);
      if(!map->isNodeGlobalElement(global_node_index)) continue;
      //set dof row start index
      current_row = num_dim*map->getLocalElement(global_node_index);
      for(int jnode = 0; jnode < nodes_per_element; jnode++){
        
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

  if(num_dim==3)
  for (int ielem = 0; ielem < num_elems; ielem++){
    element_select->choose_3Delem_type(Element_Types(ielem), elem);
    nodes_per_element = elem->num_nodes();
    //construct local stiffness matrix for this element
    local_matrix_multiply(ielem, Local_Stiffness_Matrix);
    //assign entries of this local matrix to the sparse global matrix storage;
    for (int inode = 0; inode < nodes_per_element; inode++){
      //see if this node is local
      global_node_index = mesh->nodes_in_cell_list_(ielem,inode);
      if(!map->isNodeGlobalElement(global_node_index)) continue;
      //set dof row start index
      current_row = num_dim*map->getLocalElement(global_node_index);
      for(int jnode = 0; jnode < nodes_per_element; jnode++){
        
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

  //force vector initialization
  
  //initialize
  for(int i=0; i < num_dim*nlocal_nodes; i++)
    Nodal_Forces(i,0) = 0;
  
  //Tag nodes for Boundary conditions such as displacements
  Displacement_Boundary_Conditions();

  //Construct applied nodal force vector with quadrature
  Force_Vector_Construct();

  //sync device data for the nodal forces dual view
  dual_nodal_forces.modify_host();

  //construct distributed stiffness matrix and force vector from local kokkos data
  
  //build column map for the global stiffness matrix
  Teuchos::RCP<const Tpetra::Map<LO,GO,node_type> > colmap;
  const Teuchos::RCP<const Tpetra::Map<LO,GO,node_type> > dommap = local_dof_map;
  
  //debug print
  /*
    std::cout << "DOF GRAPH MATRIX ENTRIES ON TASK " << myrank << std::endl;
  for (int idof = 0; idof < num_dim*nlocal_nodes; idof++){
    for (int istride = 0; istride < Stiffness_Matrix_Strides(idof); istride++){
      //debug print
      std::cout << "{" <<istride + 1 << "," << DOF_Graph_Matrix(idof,istride) << "} ";
    }
    //debug print
    std::cout << std::endl;
  } */

  //debug print of stiffness matrix
  /*
  std::cout << " ------------SPARSE STIFFNESS MATRIX ON TASK"<< myrank << std::endl;
  for (int idof = 0; idof < num_dim*nlocal_nodes; idof++){
      std::cout << "row: " << idof + 1 << " { ";
    for (int istride = 0; istride < Stiffness_Matrix_Strides(idof); istride++){
        std::cout << istride + 1 << " = " << Stiffness_Matrix(idof,istride) << " , " ;
    }
    std::cout << " }"<< std::endl;
  }
  */
  Tpetra::Details::makeColMap<LO,GO,node_type>(colmap,dommap,DOF_Graph_Matrix.get_kokkos_view(), nullptr);

  size_t nnz = DOF_Graph_Matrix.size();

  //debug print
  //std::cout << "DOF GRAPH SIZE ON RANK " << myrank << " IS " << nnz << std::endl;
  
  //local indices in the graph using the constructed column map
  CArrayKokkos<LO, array_layout, device_type, memory_traits> stiffness_local_indices(nnz, "stiffness_local_indices");
  
  //row offsets with compatible template arguments
    row_pointers row_offsets = DOF_Graph_Matrix.start_index_;
    row_pointers row_offsets_pass("row_offsets", nlocal_nodes*num_dim+1);
    for(int ipass = 0; ipass < nlocal_nodes*num_dim + 1; ipass++){
      row_offsets_pass(ipass) = row_offsets(ipass);
    }

  size_t entrycount = 0;
  for(int irow = 0; irow < nlocal_nodes*num_dim; irow++){
    for(int istride = 0; istride < Stiffness_Matrix_Strides(irow); istride++){
      stiffness_local_indices(entrycount) = colmap->getLocalElement(DOF_Graph_Matrix(irow,istride));
      entrycount++;
    }
  }
  
  Global_Stiffness_Matrix = Teuchos::rcp(new MAT(local_dof_map, colmap, row_offsets_pass, stiffness_local_indices.get_kokkos_view(), Stiffness_Matrix.get_kokkos_view()));

  Global_Stiffness_Matrix->fillComplete();
  //This completes the setup for A matrix of the linear system
  
  //file to debug print
  std::ostream &out = std::cout;
  Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));

  //debug print of A matrix
  //*fos << "Global Stiffness Matrix :" << std::endl;
  //Global_Stiffness_Matrix->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;
  
  //allocate global nodal force vector
  //Global_Nodal_Forces = Teuchos::rcp(new MV(local_dof_map, dual_nodal_forces));

  //Print solution vector
  //*fos << "Global Nodal Forces :" << std::endl;
  //Global_Nodal_Forces->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;
}

/* ----------------------------------------------------------------------
   Retrieve material properties associated with a finite element
------------------------------------------------------------------------- */

void Static_Solver_Parallel::Element_Material_Properties(size_t ielem, real_t &Element_Modulus, real_t &Poisson_Ratio){

  Element_Modulus = simparam->Elastic_Modulus;
  Poisson_Ratio = simparam->Poisson_Ratio;

}

/* ----------------------------------------------------------------------
   Construct the local stiffness matrix
------------------------------------------------------------------------- */

void Static_Solver_Parallel::local_matrix(int ielem, CArray <real_t> &Local_Matrix){
  //local variable for host view in the dual view
  host_vec_array all_node_data = dual_all_node_data.view_device();
  int num_dim = simparam->num_dim;
  int num_elems = mesh->num_elems();
  int nodes_per_elem = mesh->num_nodes_in_elem();
  int num_gauss_points = simparam->num_gauss_points;
  int z_quad,y_quad,x_quad, direct_product_count;
  size_t local_node_id;

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
    local_node_id = all_node_map->getLocalElement(mesh->nodes_in_cell_list_(ielem, node_loop));
    nodal_positions(node_loop,0) = all_node_data(local_node_id,0);
    nodal_positions(node_loop,1) = all_node_data(local_node_id,1);
    nodal_positions(node_loop,2) = all_node_data(local_node_id,2);
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

void Static_Solver_Parallel::local_matrix_multiply(int ielem, CArray <real_t> &Local_Matrix){
  //local variable for host view in the dual view
  host_vec_array all_node_data = dual_all_node_data.view_device();
  int num_dim = simparam->num_dim;
  int nodes_per_elem = elem->num_basis();
  int num_gauss_points = simparam->num_gauss_points;
  int z_quad,y_quad,x_quad, direct_product_count;
  size_t local_node_id;

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
    local_node_id = all_node_map->getLocalElement(mesh->nodes_in_cell_list_(ielem, node_loop));
    nodal_positions(node_loop,0) = all_node_data(local_node_id,0);
    nodal_positions(node_loop,1) = all_node_data(local_node_id,1);
    nodal_positions(node_loop,2) = all_node_data(local_node_id,2);
    /*
    if(myrank==1&&nodal_positions(node_loop,2)>10000000){
      std::cout << " LOCAL MATRIX DEBUG ON TASK " << myrank << std::endl;
      std::cout << node_loop+1 <<" " << local_node_id <<" "<< mesh->nodes_in_cell_list_(ielem, node_loop) << " "<< nodal_positions(node_loop,2) << std::endl;
      std::fflush(stdout);
    }
    */
    //std::cout << local_node_id << " " << mesh->nodes_in_cell_list_(ielem, node_loop) << " " << nodal_positions(node_loop,0) << " " << nodal_positions(node_loop,1) << " "<< nodal_positions(node_loop,2) <<std::endl;
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
      //debug print
    /*if(myrank==1&&nodal_positions(node_loop,2)*basis_derivative_s3(node_loop)<-10000000){
      std::cout << " LOCAL MATRIX DEBUG ON TASK " << myrank << std::endl;
      std::cout << node_loop+1 << " " << JT_row3(2) << " "<< nodal_positions(node_loop,2) <<" "<< basis_derivative_s3(node_loop) << std::endl;
      std::fflush(stdout);
    }*/
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

void Static_Solver_Parallel::Displacement_Boundary_Conditions(){
  int num_bdy_patches_in_set, patch_id;
  int warning_flag = 0;
  int local_flag;
  int current_node_index, current_node_id;
  int num_boundary_sets = num_boundary_conditions;
  int surface_disp_set_id = 0;
  int num_dim = simparam->num_dim;
  int bc_option, bc_dim_set[3];
  CArray<real_t> displacement(num_dim);
  CArray<size_t> Displacement_Conditions(num_dim);
  CArray<size_t> first_condition_per_node(nall_nodes*num_dim);
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Surface_Nodes;
  Number_DOF_BCS = 0;
  Displacement_Conditions(0) = X_DISPLACEMENT_CONDITION;
  Displacement_Conditions(1) = Y_DISPLACEMENT_CONDITION;
  Displacement_Conditions(2) = Z_DISPLACEMENT_CONDITION;

  //initialize to -1 (DO NOT make -1 an index for bdy sets)
  for(int inode = 0 ; inode < nlocal_nodes*num_dim; inode++)
    first_condition_per_node(inode) = -1;
  
  //scan for surface method of setting fixed nodal displacements
  for(int iboundary = 0; iboundary < num_boundary_sets; iboundary++){
    
    if(Boundary_Condition_Type_List(iboundary)==DISPLACEMENT_CONDITION){bc_option=3;}
    else if(Boundary_Condition_Type_List(iboundary)==X_DISPLACEMENT_CONDITION){bc_option=0;}
    else if(Boundary_Condition_Type_List(iboundary)==Y_DISPLACEMENT_CONDITION){bc_option=1;}
    else if(Boundary_Condition_Type_List(iboundary)==Z_DISPLACEMENT_CONDITION){bc_option=2;}
    else{
      continue;
    }
      
      //debug print of surface conditions
      //std::cout << " Surface BC types " << Boundary_Condition_Type_List(iboundary) <<std::endl;

      num_bdy_patches_in_set = NBoundary_Condition_Patches(iboundary);
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
        //get number
        patch_id = Boundary_Condition_Patches(iboundary, bdy_patch_gid);
        Surface_Nodes = Boundary_Patches(patch_id).node_set;
        for(int inode = 0; inode < Surface_Nodes.size(); inode++){
          local_flag = 0;
          current_node_id = Surface_Nodes(inode);
          if(map->isNodeGlobalElement(current_node_id)) local_flag = 1;
          current_node_id = all_node_map->getLocalElement(current_node_id);
          
          /*
          //debug print of nodal bc settings
          std::cout << " Node BC types " << Node_DOF_Boundary_Condition_Type(current_node_id);
          std::cout << " node: " << inode << " { ";
          for (int istride = 0; istride < num_dim; istride++){
            std::cout << node_data(current_node_id,istride) << " , ";
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
              //counts local DOF being constrained
              if(local_flag)
              Number_DOF_BCS++;
            }
          }
          }
        }
      }
  }

  //scan for direct setting of nodal displacements from input
  //indices for nodal BC settings referred to here start at num_boundary_sets

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

void Static_Solver_Parallel::Force_Vector_Construct(){
  //local variable for host view in the dual view
  host_vec_array all_node_data = dual_all_node_data.view_host();
  //local variable for host view in the dual view
  host_vec_array Nodal_Forces = dual_nodal_forces.view_host();
  int num_bdy_patches_in_set;
  size_t current_node_index, node_id, patch_id;
  int num_boundary_sets = num_boundary_conditions;
  int surface_force_set_id = 0;
  int num_dim = simparam->num_dim;
  int num_elems = mesh->num_elems();
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
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Surface_Nodes;
  
  CArray<real_t> JT_row1(num_dim);
  CArray<real_t> JT_row2(num_dim);
  CArray<real_t> JT_row3(num_dim);

  /*Loop through boundary sets and check if they apply surface forces.
  These sets can have overlapping nodes since applied loading conditions
  are assumed to be additive*/
  for(int iboundary = 0; iboundary < num_boundary_sets; iboundary++){
    if(Boundary_Condition_Type_List(iboundary)!=LOADING_CONDITION) continue;
    //std::cout << "I REACHED THE LOADING BOUNDARY CONDITION" <<std::endl;
    num_bdy_patches_in_set = NBoundary_Condition_Patches(iboundary);
    
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
    patch_id = Boundary_Condition_Patches(iboundary, bdy_patch_gid);
    Surface_Nodes = Boundary_Patches(patch_id).node_set;
    //find element index this boundary patch is on
    current_element_index = Boundary_Patches(patch_id).element_id;
    local_surface_id = Boundary_Patches(patch_id).local_patch_id;
    //debug print of local surface ids
    //std::cout << " LOCAL SURFACE IDS " << std::endl;
    //std::cout << local_surface_id << std::endl;

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
        current_node_index = Surface_Nodes(node_loop);
        current_node_index = all_node_map->getLocalElement(current_node_index);
        nodal_positions(node_loop,0) = all_node_data(current_node_index,0);
        nodal_positions(node_loop,1) = all_node_data(current_node_index,1);
        nodal_positions(node_loop,2) = all_node_data(current_node_index,2);
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
            
        node_id = mesh->nodes_in_cell(current_element_index, local_nodes[node_count]);
        //check if node is local to alter Nodal Forces vector
        if(!map->isNodeGlobalElement(node_id)) continue;
        node_id = map->getLocalElement(node_id);
        
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
          Nodal_Forces(num_dim*node_id + idim,0) += wedge_product*quad_coordinate_weight(0)*quad_coordinate_weight(1)*force_density[idim]*basis_values(local_nodes[node_count]);
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

int Static_Solver_Parallel::solve(){
  //local variable for host view in the dual view
  host_vec_array node_data = dual_node_data.view_host();
  //local variable for host view in the dual view
  host_vec_array Nodal_Forces = dual_nodal_forces.view_host();
  int num_dim = simparam->num_dim;
  int num_elems = mesh->num_elems();
  int nodes_per_elem = mesh->num_nodes_in_elem();
  int local_node_index, current_row, current_column;
  size_t reduced_index;
  int max_stride = 0;
  global_size_t global_index, reduced_row_count;
  
  // Before we do anything, check that KLU2 is enabled
  if( !Amesos2::query("KLU2") ){
    std::cerr << "SuperLUDist not enabled in this run.  Exiting..." << std::endl;
    return EXIT_SUCCESS;        // Otherwise CTest will pick it up as
                                // failure, which it isn't really
  }
  
  std::ostream &out = std::cout;
  Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));

  *fos << Amesos2::version() << std::endl << std::endl;

  bool printTiming   = true;
  bool verbose       = false;
  std::string filename("arc130.mtx");
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("print-timing","no-print-timing",&printTiming,"Print solver timing statistics");

  const size_t numVectors = 1;
  
  //construct global data for example; normally this would be from file input and distributed
  //according to the row map at that point
  global_size_t nrows = num_nodes*num_dim;
  
 
  //number of boundary conditions on this mpi rank
  global_size_t local_nboundaries = Number_DOF_BCS;
  size_t local_nrows_reduced = nlocal_nodes*num_dim - local_nboundaries;

  //obtain total number of boundary conditions on all ranks
  global_size_t global_nboundaries = 0;
  MPI_Allreduce(&local_nboundaries,&global_nboundaries,1,MPI_INT,MPI_SUM,world);
  global_size_t nrows_reduced = nrows - global_nboundaries;

  //Rebalance distribution of the global stiffness matrix rows here later since
  //rows and columns are being removed.

  //global_size_t *entries_per_row = new global_size_t[local_nrows];
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Reduced_Stiffness_Matrix_Strides(local_nrows_reduced,"Reduced_Stiffness_Matrix_Strides");
  row_pointers reduced_row_offsets_pass = row_pointers("reduced_row_offsets_pass", local_nrows_reduced+1);

  //init row_offsets
  for(int i=0; i < local_nrows_reduced+1; i++){
    reduced_row_offsets_pass(i) = 0;
  }
  
  //stores global indices belonging to this MPI rank from the non-reduced map corresponding to the reduced system
  CArrayKokkos<GO, array_layout, device_type, memory_traits> Free_Indices(local_nrows_reduced,"Free_Indices");
  reduced_index = 0;
  for(LO i=0; i < nlocal_nodes*num_dim; i++)
    if((Node_DOF_Boundary_Condition_Type(i)!=DISPLACEMENT_CONDITION)&&
       (Node_DOF_Boundary_Condition_Type(i)!=X_DISPLACEMENT_CONDITION)&&
       (Node_DOF_Boundary_Condition_Type(i)!=Y_DISPLACEMENT_CONDITION)&&
       (Node_DOF_Boundary_Condition_Type(i)!=Z_DISPLACEMENT_CONDITION)){
        Free_Indices(reduced_index) = local_dof_map->getGlobalElement(i);
        reduced_index++;
      }
    
  
  //compute reduced matrix strides
  size_t access_index, row_access_index, row_counter;
  GO global_dof_index, reduced_local_dof_index;
  for(LO i=0; i < local_nrows_reduced; i++){
    access_index = local_dof_map->getLocalElement(Free_Indices(i));
    reduced_row_count = 0;
    for(LO j=0; j < Stiffness_Matrix_Strides(access_index); j++){
      global_dof_index = DOF_Graph_Matrix(access_index,j);
      row_access_index = all_dof_map->getLocalElement(global_dof_index);
      if((Node_DOF_Boundary_Condition_Type(row_access_index)!=DISPLACEMENT_CONDITION)&&
         (Node_DOF_Boundary_Condition_Type(row_access_index)!=X_DISPLACEMENT_CONDITION)&&
         (Node_DOF_Boundary_Condition_Type(row_access_index)!=Y_DISPLACEMENT_CONDITION)&&
         (Node_DOF_Boundary_Condition_Type(row_access_index)!=Z_DISPLACEMENT_CONDITION)){
        reduced_row_count++;
      }
    }
    Reduced_Stiffness_Matrix_Strides(i) = reduced_row_count;
  }
  
  RaggedRightArrayKokkos<GO, array_layout, device_type, memory_traits> Reduced_DOF_Graph_Matrix(Reduced_Stiffness_Matrix_Strides); //stores global dof indices
  RaggedRightArrayKokkos<LO, array_layout, device_type, memory_traits> Reduced_Local_DOF_Graph_Matrix(Reduced_Stiffness_Matrix_Strides); //stores local dof indices
  RaggedRightArrayKokkos<real_t, Kokkos::LayoutRight, device_type, memory_traits, array_layout> Reduced_Stiffness_Matrix(Reduced_Stiffness_Matrix_Strides);

  //template compatible row offsets (may change to shallow copy later if it works on device types etc.)
  row_pointers reduced_row_offsets = Reduced_DOF_Graph_Matrix.start_index_;
    for(int ipass = 0; ipass < local_nrows_reduced+1; ipass++){
      reduced_row_offsets_pass(ipass) = reduced_row_offsets(ipass);
    }
  
  //construct maps to define set of global indices for the reduced node set
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > local_reduced_dof_original_map =
    Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(nrows_reduced,Free_Indices.get_kokkos_view(),0,comm) );
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > local_reduced_dof_map = 
    Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(nrows_reduced,local_nrows_reduced,0,comm));
  
  //dual view of the local global index to reduced global index map
  dual_vec_array dual_reduced_dof_original("dual_reduced_dof_original",local_nrows_reduced,1);

  //local variable for host view in the dual view
  host_vec_array reduced_dof_original = dual_reduced_dof_original.view_host();
  //notify that the host view is going to be modified
  dual_reduced_dof_original.modify_host();

  //set contents
  for(LO i=0; i < local_nrows_reduced; i++){
    reduced_dof_original(i,0) = local_reduced_dof_map->getGlobalElement(i);
  }
  
  //create a multivector where each remaining local index entry stores the new reduced global index associated with each old global index
  Teuchos::RCP<MV> local_reduced_global_indices = Teuchos::rcp(new MV(local_reduced_dof_original_map, dual_reduced_dof_original));
  
  //construct map of all indices including ghosts for the reduced system
  //stores global indices belonging to this MPI rank and ghosts from the non-reduced map corresponding to the reduced system
  size_t all_nrows_reduced = local_nrows_reduced + nghost_nodes*num_dim;
  for(LO i=nlocal_nodes*num_dim; i < nall_nodes*num_dim; i++){
      if((Node_DOF_Boundary_Condition_Type(i)==DISPLACEMENT_CONDITION)||
         (Node_DOF_Boundary_Condition_Type(i)==X_DISPLACEMENT_CONDITION)||
         (Node_DOF_Boundary_Condition_Type(i)==Y_DISPLACEMENT_CONDITION)||
         (Node_DOF_Boundary_Condition_Type(i)==Z_DISPLACEMENT_CONDITION))
      all_nrows_reduced--;
  }
  
  CArrayKokkos<GO, array_layout, device_type, memory_traits> All_Free_Indices(all_nrows_reduced,"All_Free_Indices");

  //debug print
  /*
  if(myrank==0||myrank==4){
  std::cout << "DOF flags global :" << std::endl;
  std::cout << "Reduced DOF Graph Matrix on Rank " << myrank << std::endl;
  for(LO i=0; i < nall_nodes*num_dim; i++){
    std::cout << all_dof_map->getGlobalElement(i) << " " << Node_DOF_Boundary_Condition_Type(i) <<" ";   
    std::cout << std::endl;
  }
  std::fflush(stdout);
  }
  */
  
  reduced_index = 0;
  for(LO i=0; i < nall_nodes*num_dim; i++)
    if((Node_DOF_Boundary_Condition_Type(i)!=DISPLACEMENT_CONDITION)&&
       (Node_DOF_Boundary_Condition_Type(i)!=X_DISPLACEMENT_CONDITION)&&
       (Node_DOF_Boundary_Condition_Type(i)!=Y_DISPLACEMENT_CONDITION)&&
       (Node_DOF_Boundary_Condition_Type(i)!=Z_DISPLACEMENT_CONDITION)){
        All_Free_Indices(reduced_index) = all_dof_map->getGlobalElement(i);
        reduced_index++;
      }
  
  //construct map to define set of global indices for the reduced node set including ghosts
  //passing invalid forces the map to count the global elements
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > all_reduced_dof_original_map =
    Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(Teuchos::OrdinalTraits<GO>::invalid(),All_Free_Indices.get_kokkos_view(),0,comm) );

  //debug print
  /*
  if(myrank==0)
  *fos << "All reduced dof original indices :" << std::endl;
  local_reduced_dof_original_map->describe(*fos,Teuchos::VERB_EXTREME);
  *fos << std::endl;
  std::fflush(stdout);

  //debug print
  if(myrank==0)
  *fos << "All reduced dof original indices :" << std::endl;
  all_reduced_dof_original_map->describe(*fos,Teuchos::VERB_EXTREME);
  *fos << std::endl;
  std::fflush(stdout);
  */

  //communicate the new reduced global indices for ghost dof indices using the local information on other ranks through import
  //create import object using local node indices map and all indices map
  Tpetra::Import<LO, GO> importer(local_reduced_dof_original_map, all_reduced_dof_original_map);

  Teuchos::RCP<MV> all_reduced_global_indices = Teuchos::rcp(new MV(all_reduced_dof_original_map, 1));

  //comms to get ghosts
  all_reduced_global_indices->doImport(*local_reduced_global_indices, importer, Tpetra::INSERT);

  //get views from multivector object

  vec_array all_reduced_global_indices_device = all_reduced_global_indices->getLocalView<device_type> (Tpetra::Access::ReadWrite);
  host_vec_array all_reduced_global_indices_host = all_reduced_global_indices->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);

  //allocate node storage with dual view
  dual_vec_array dual_all_reduced_global_indices = dual_vec_array(all_reduced_global_indices_device, all_reduced_global_indices_host);

  //debug print
  /*
  if(myrank==0)
  *fos << "All reduced global indices :" << std::endl;
  all_reduced_global_indices->describe(*fos,Teuchos::VERB_EXTREME);
  *fos << std::endl;
  std::fflush(stdout);
  */
  
  //store the new global indices for the reduced matrix graph
  for(LO i=0; i < local_nrows_reduced; i++){
    access_index = local_dof_map->getLocalElement(Free_Indices(i));
    row_counter = 0;
    for(LO j=0; j < Stiffness_Matrix_Strides(access_index); j++){
      global_dof_index = DOF_Graph_Matrix(access_index,j);
      row_access_index = all_dof_map->getLocalElement(global_dof_index);
      if((Node_DOF_Boundary_Condition_Type(row_access_index)!=DISPLACEMENT_CONDITION)&&
         (Node_DOF_Boundary_Condition_Type(row_access_index)!=X_DISPLACEMENT_CONDITION)&&
         (Node_DOF_Boundary_Condition_Type(row_access_index)!=Y_DISPLACEMENT_CONDITION)&&
         (Node_DOF_Boundary_Condition_Type(row_access_index)!=Z_DISPLACEMENT_CONDITION)){
        reduced_local_dof_index = all_reduced_dof_original_map->getLocalElement(global_dof_index);
        //std::cout << "REDUCED LOCAL INDEX ON TASK " << myrank << " is " << Reduced_Stiffness_Matrix_Strides(i) << Reduced_DOF_Graph_Matrix(i,row_counter++) << std::endl;
        Reduced_DOF_Graph_Matrix(i,row_counter++) = all_reduced_global_indices_host(reduced_local_dof_index,0);
      }
    }
  }
  
  //store reduced stiffness matrix values
  for(LO i=0; i < local_nrows_reduced; i++){
    access_index = local_dof_map->getLocalElement(Free_Indices(i));
    row_counter = 0;
    for(LO j=0; j < Stiffness_Matrix_Strides(access_index); j++){
      global_dof_index = DOF_Graph_Matrix(access_index,j);
      row_access_index = all_dof_map->getLocalElement(global_dof_index);
      if((Node_DOF_Boundary_Condition_Type(row_access_index)!=DISPLACEMENT_CONDITION)&&
         (Node_DOF_Boundary_Condition_Type(row_access_index)!=X_DISPLACEMENT_CONDITION)&&
         (Node_DOF_Boundary_Condition_Type(row_access_index)!=Y_DISPLACEMENT_CONDITION)&&
         (Node_DOF_Boundary_Condition_Type(row_access_index)!=Z_DISPLACEMENT_CONDITION)){
        Reduced_Stiffness_Matrix(i,row_counter++) = Stiffness_Matrix(access_index,j);
      }
    }
  }
  
  //indices_array all_indices = indices_array("indices_array",nnz);
  //values_array all_values = values_array("values_array",nnz);
  CArrayKokkos<real_t, Kokkos::LayoutLeft, device_type, memory_traits> Bview(local_nrows_reduced);
  CArrayKokkos<real_t, Kokkos::LayoutLeft, device_type, memory_traits> Xview(local_nrows_reduced);
  
  // create a Map for the reduced global stiffness matrix that is evenly distributed amongst mpi ranks
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > local_balanced_reduced_dof_map = 
    Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(nrows_reduced,0,comm));

  //build column map
  Teuchos::RCP<const Tpetra::Map<LO,GO,node_type> > colmap;
  const Teuchos::RCP<const Tpetra::Map<LO,GO,node_type> > dommap = local_reduced_dof_map;

  //debug print
  /*
  if(myrank==4){
  std::cout << "Reduced DOF Graph Matrix on Rank " << myrank << std::endl;
  for(LO i=0; i < local_nrows_reduced; i++){
    for(LO j=0; j < Reduced_Stiffness_Matrix_Strides(i); j++){
      std::cout << Reduced_DOF_Graph_Matrix(i,j) <<" ";
    }
    std::cout << std::endl;
  }
  }
  */

  Tpetra::Details::makeColMap<LO,GO,node_type>(colmap,dommap,Reduced_DOF_Graph_Matrix.get_kokkos_view(), nullptr);

  /*//debug print of reduced row offsets
  std::cout << " DEBUG PRINT FOR ROW OFFSETS" << std::endl;
  for(int debug = 0; debug < local_nrows_reduced+1; debug++)
  std::cout <<  reduced_row_offsets_pass(debug) << " ";
  std::cout << std::endl;
  //end debug print
  */

  //convert global indices to local indices using column map
  for(LO i=0; i < local_nrows_reduced; i++)
    for(LO j=0; j < Reduced_Stiffness_Matrix_Strides(i); j++)
      Reduced_Local_DOF_Graph_Matrix(i,j) = colmap->getLocalElement(Reduced_DOF_Graph_Matrix(i,j));
  
  Teuchos::RCP<MAT> unbalanced_A = Teuchos::rcp(new MAT(local_reduced_dof_map, colmap, reduced_row_offsets_pass,
                   Reduced_Local_DOF_Graph_Matrix.get_kokkos_view(), Reduced_Stiffness_Matrix.get_kokkos_view()));
  unbalanced_A->fillComplete();
  Teuchos::RCP<const_MAT> const_unbalanced_A(new const_MAT(*unbalanced_A));
  
  //This completes the setup for A matrix of the linear system

  //debug print of A matrix before balancing
  //*fos << "Reduced Stiffness Matrix :" << std::endl;
  //const_unbalanced_A->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;

  //communicate reduced stiffness matrix entries for better load balancing
  //create import object using the unbalanced map and the balanced map
  Tpetra::Import<LO, GO> matrix_importer(local_reduced_dof_map, local_balanced_reduced_dof_map);
  Teuchos::RCP<MAT> balanced_A = Tpetra::importAndFillCompleteCrsMatrix(const_unbalanced_A, matrix_importer, local_balanced_reduced_dof_map, local_balanced_reduced_dof_map);

  //debug print of map
  //if(myrank==0)
  //*fos << "Reduced DOF Map :" << std::endl;
  //local_balanced_reduced_dof_map->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;
  //std::fflush(stdout);

  //debug print of A matrix after balancing
  //if(myrank==0)
  //*fos << "Reduced Stiffness Matrix :" << std::endl;
  //balanced_A->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;
  //std::fflush(stdout);
  
  // Create random X vector
  size_t balanced_local_nrows = local_balanced_reduced_dof_map->getNodeNumElements();
  vec_array Xview_pass = vec_array("Xview_pass", balanced_local_nrows, 1);
  Xview_pass.assign_data(Xview.pointer());
  X = Teuchos::rcp(new MV(local_balanced_reduced_dof_map, Xview_pass));
  //return !EXIT_SUCCESS;
  //X->randomize();
  
  
  //print allocation of the solution vector to check distribution
  /*
  if(myrank==0)
  *fos << "ALLOCATED SOLUTION VECTOR:" << std::endl;
  X->describe(*fos,Teuchos::VERB_EXTREME);
  *fos << std::endl;
  */

  // Create Kokkos view of RHS B vector (Force Vector)  
  vec_array Bview_pass = vec_array("Bview_pass", local_nrows_reduced,1);

  //set bview to force vector data
  for(LO i=0; i < local_nrows_reduced; i++){
    access_index = local_dof_map->getLocalElement(Free_Indices(i));
    Bview(i) = Nodal_Forces(access_index,0);
  }
  
  Bview_pass.assign_data(Bview.pointer());
  
  Teuchos::RCP<MV> unbalanced_B = Teuchos::rcp(new MV(local_reduced_dof_map, Bview_pass));
  
  //import object to rebalance force vector
  Tpetra::Import<LO, GO> Bvec_importer(local_reduced_dof_map, local_balanced_reduced_dof_map);

  Teuchos::RCP<MV> balanced_B = Teuchos::rcp(new MV(local_balanced_reduced_dof_map, 1));
  
  //comms to rebalance force vector
  balanced_B->doImport(*unbalanced_B, Bvec_importer, Tpetra::INSERT);
  
  //if(myrank==0)
  //*fos << "RHS :" << std::endl;
  //balanced_B->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;

  //debug print
  //Tpetra::MatrixMarket::Writer<MAT> market_writer();
  //Tpetra::MatrixMarket::Writer<MAT>::writeSparseFile("A_matrix.txt", *balanced_A, "A_matrix", "Stores stiffness matrix values");
  //return !EXIT_SUCCESS;
  // Create solver interface to KLU2 with Amesos2 factory method
  
  Teuchos::RCP<Amesos2::Solver<MAT,MV>> solver = Amesos2::create<MAT,MV>("KLU2", balanced_A, X, balanced_B);
  
  //declare non-contiguous map
  //Create a Teuchos::ParameterList to hold solver parameters
  //Teuchos::ParameterList amesos2_params("Amesos2");
  //amesos2_params.sublist("KLU2").set("IsContiguous", false, "Are GIDs Contiguous");
  //solver->setParameters( Teuchos::rcpFromRef(amesos2_params) );

  //Solve the system
  solver->symbolicFactorization().numericFactorization().solve();

  //timing statistics for LU solver
  //solver->printTiming(*fos);
  
  //Print solution vector
  //if(myrank==0)
  //*fos << "Solution :" << std::endl;
  //X->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;
  
  return !EXIT_SUCCESS;
}

/* ----------------------------------------------------------------------
   Return the CPU time for the current process in seconds very
   much in the same way as MPI_Wtime() returns the wall time.
------------------------------------------------------------------------- */

double Static_Solver_Parallel::CPU_Time()
{
  double rv = 0.0;
/*
#ifdef _WIN32

  // from MSD docs.
  FILETIME ct,et,kt,ut;
  union { FILETIME ft; uint64_t ui; } cpu;
  if (GetProcessTimes(GetCurrentProcess(),&ct,&et,&kt,&ut)) {
    cpu.ft = ut;
    rv = cpu.ui * 0.0000001;
  }
*/

  struct rusage ru;
  if (getrusage(RUSAGE_SELF, &ru) == 0) {
    rv = (double) ru.ru_utime.tv_sec;
    rv += (double) ru.ru_utime.tv_usec * 0.000001;
  }


  return rv;
}

/* ----------------------------------------------------------------------
   Clock variable initialization
------------------------------------------------------------------------- */

void Static_Solver_Parallel::init_clock(){
  double current_cpu = 0;
  initial_CPU_time = CPU_Time();
}
