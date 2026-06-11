#ifndef MESH_IO_H
#define MESH_IO_H

#include "ELEMENTS.h"
#include "state.h"

using namespace mtr;

#include <map>
#include <memory>
#include <cstring>
#include <sys/stat.h>
#include <iostream>
#include <regex>    // for string pattern recoginition
#include <fstream>
#include <sstream>
#include <vector>
#include <string>   
#include <mpi.h>
#include <string>



/////////////////////////////////////////////////////////////////////////////
///
/// \fn split
///
/// \brief Splits a string by a given delimiter
///
/// \param Input string
/// \param delimiter
///
/// \return Vector of split string values
///
/////////////////////////////////////////////////////////////////////////////
inline std::vector<std::string> split(std::string s, std::string delimiter)
{
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> res;

    while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
        token     = s.substr(pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back(token);
    }

    res.push_back(s.substr(pos_start));
    return res;
} // end of split

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_id
///
/// \brief This gives the index value of the point or the elem
///
/// Assumes that the grid has an i,j,k structure
/// the elem = i + (j)*(num_points_i-1) + (k)*(num_points_i-1)*(num_points_j-1)
/// the point = i + (j)*num_points_i + (k)*num_points_i*num_points_j
///
/// \param i index
/// \param j index
/// \param k index
/// \param Number of i indices
/// \param Number of j indices
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_INLINE_FUNCTION
size_t get_id(int i, int j, int k, int num_i, int num_j)
{
    return i + j * num_i + k * num_i * num_j;
} // end get_id

/////////////////////////////////////////////////////////////////////////////
///
/// \fn PointIndexFromIJK
///
/// \brief Given (i,j,k) coordinates within the Lagrange hex, return an 
/// offset into the local connectivity (PointIds) array. The order parameter
/// must point to an array of 3 integers specifying the order along each 
/// axis of the hexahedron.
///
/////////////////////////////////////////////////////////////////////////////
int PointIndexFromIJK(int i, int j, int k, const int* order)
{
    bool ibdy = (i == 0 || i == order[0]);
    bool jbdy = (j == 0 || j == order[1]);
    bool kbdy = (k == 0 || k == order[2]);
    // How many boundaries do we lie on at once?
    int nbdy = (ibdy ? 1 : 0) + (jbdy ? 1 : 0) + (kbdy ? 1 : 0);
  
    if (nbdy == 3) // Vertex DOF
      { // ijk is a corner node. Return the proper index (somewhere in [0,7]):
      return (i ? (j ? 2 : 1) : (j ? 3 : 0)) + (k ? 4 : 0);
      }
  
    int offset = 8;
    if (nbdy == 2) // Edge DOF
      {
      if (!ibdy)
        { // On i axis
        return (i - 1) +
          (j ? order[0] - 1 + order[1] - 1 : 0) +
          (k ? 2 * (order[0] - 1 + order[1] - 1) : 0) +
          offset;
        }
      if (!jbdy)
        { // On j axis
        return (j - 1) +
          (i ? order[0] - 1 : 2 * (order[0] - 1) + order[1] - 1) +
          (k ? 2 * (order[0] - 1 + order[1] - 1) : 0) +
          offset;
        }
      // !kbdy, On k axis
      offset += 4 * (order[0] - 1) + 4 * (order[1] - 1);
      return (k - 1) + (order[2] - 1) * (i ? (j ? 3 : 1) : (j ? 2 : 0)) + offset;
      }
  
    offset += 4 * (order[0] - 1 + order[1] - 1 + order[2] - 1);
    if (nbdy == 1) // Face DOF
      {
      if (ibdy) // On i-normal face
        {
        return (j - 1) + ((order[1] - 1) * (k - 1)) + (i ? (order[1] - 1) * (order[2] - 1) : 0) + offset;
        }
      offset += 2 * (order[1] - 1) * (order[2] - 1);
      if (jbdy) // On j-normal face
        {
        return (i - 1) + ((order[0] - 1) * (k - 1)) + (j ? (order[2] - 1) * (order[0] - 1) : 0) + offset;
        }
      offset += 2 * (order[2] - 1) * (order[0] - 1);
      // kbdy, On k-normal face
      return (i - 1) + ((order[0] - 1) * (j - 1)) + (k ? (order[0] - 1) * (order[1] - 1) : 0) + offset;
      }
  
    // nbdy == 0: Body DOF
    offset += 2 * (
      (order[1] - 1) * (order[2] - 1) +
      (order[2] - 1) * (order[0] - 1) +
      (order[0] - 1) * (order[1] - 1));
    return offset +
      (i - 1) + (order[0] - 1) * (
        (j - 1) + (order[1] - 1) * (
          (k - 1)));
  }
  
  




/////////////////////////////////////////////////////////////////////////////
///
/// \fn build_3d_box
///
/// \brief Builds an unstructured 3D rectilinear mesh and initializes node coordinates
///
/// \param mesh         Reference to the mesh object to be constructed
/// \param node_coords  Reference to the array where node coordinates will be stored
/// \param origin       Array specifying the origin of the mesh (x,y,z)
/// \param length       Array specifying the length of the mesh in each direction (x,y,z)
/// \param num_elems_dim Array specifying the number of elements in each direction (i,j,k)
///
/////////////////////////////////////////////////////////////////////////////
void build_3d_box(
    swage::Mesh& mesh,
    MPICArrayKokkos<double>& node_coords,
    double origin[3],
    double length[3],
    int num_elems_dim[3],
    int Pn_order)
{
    printf("Creating a 3D box mesh \n");

    const int num_dim = 3;

    mesh.Pn = Pn_order;
    mesh.num_dims = num_dim;

    // Note: In fierro, these come from the simulation parameters
    const double lx = length[0];
    const double ly = length[1];
    const double lz = length[2];

    // Note: In fierro, these come from the simulation parameters
    const int num_elems_i = num_elems_dim[0];
    const int num_elems_j = num_elems_dim[1];
    const int num_elems_k = num_elems_dim[2];

    int num_points_in_elem = std::pow(Pn_order + 1, num_dim);

    // Calculate the number of unique node points along the X axis,
    // including shared nodes at element boundaries
    int num_points_i = num_elems_i * Pn_order + 1;
    int num_points_j = num_elems_j * Pn_order + 1;
    int num_points_k = num_elems_k * Pn_order + 1;

    const int num_nodes = num_points_i * num_points_j * num_points_k;

    // Note: this should be modified to account for quadrature spacing, but not yet. 
    const double dx = lx / ((double)num_points_i - 1);  // len/(num_points_i)
    const double dy = ly / ((double)num_points_j - 1);  // len/(num_points_j)
    const double dz = lz / ((double)num_points_k - 1);  // len/(num_points_k)

    const int num_elems = num_elems_i * num_elems_j * num_elems_k;

    // --- 3D parameters ---
    // const int num_faces_in_elem  = 6;  // number of faces in elem
    // const int num_points_in_elem = 8;  // number of points in elem
    // const int num_points_in_face = 4;  // number of points in a face
    // const int num_edges_in_elem  = 12; // number of edges in a elem

    // initialize mesh node variables
    mesh.initialize_nodes(num_nodes);


    // initialize node coordinates
    node_coords = MPICArrayKokkos<double>(num_nodes, num_dim, "node_coordinates");

    // --- Build nodes ---

    CArrayDual<double> origin_mtr(3, "origin_mtr");
    origin_mtr.host(0) = origin[0];
    origin_mtr.host(1) = origin[1];
    origin_mtr.host(2) = origin[2];
    origin_mtr.update_device();

    // populate the point data structures
    FOR_ALL(k, 0, num_points_k,
            j, 0, num_points_j,
            i, 0, num_points_i,{

        // global id for the point
        size_t node_gid = get_id(i, j, k, num_points_i, num_points_j);

        // store the point coordinates
        node_coords(node_gid, 0) = origin_mtr(0) + (double)i * dx;
        node_coords(node_gid, 1) = origin_mtr(1) + (double)j * dy;
        node_coords(node_gid, 2) = origin_mtr(2) + (double)k * dz;
    });
    // Update the host side
    node_coords.update_host();

    // initialize elem variables
    if (Pn_order == 0){
        mesh.initialize_elems(num_elems, num_dim); 
        Pn_order = 1;
    } else {
        mesh.initialize_elems_Pn(num_elems, num_dim, Pn_order);
    }

    // populate the point data structures
    FOR_ALL(k, 0, num_elems_k,
            j, 0, num_elems_j,
            i, 0, num_elems_i,{

        // global id for the elem
        size_t elem_gid = get_id(i, j, k, num_elems_i, num_elems_j);

        // size_t k_local_max = num_dim < 3 ? 0 : Pn_order;


        // store the point IDs for this elem where the range is
        // (i:i+1, j:j+1, k:k+1) for a linear hexahedron
        int this_point = 0;

        size_t i_offset = i * Pn_order;
        size_t j_offset = j * Pn_order;
        size_t k_offset = k * Pn_order;


        for (int kcount = 0; kcount <= Pn_order; kcount++) {
            for (int jcount = 0; jcount <= Pn_order; jcount++) {
                for (int icount = 0; icount <= Pn_order; icount++) {
                    // global id for the points

                    size_t node_gid = get_id(
                        (icount + i_offset) % num_points_i, 
                        (jcount + j_offset) % num_points_j, 
                        (kcount + k_offset) % num_points_k,
                        num_points_i, num_points_j
                    );

                    // convert this_point index to the FE index convention
                    int this_index = this_point; //convert_point_number_in_Hex(this_point);

                    // store the points in this elem according the the finite
                    // element numbering convention
                    mesh.nodes_in_elem(elem_gid, this_index) = node_gid;

                    // increment the point counting index
                    this_point++;
                } // end for icount
            } // end for jcount
        }  // end for kcount
    }); // end parallel for


    // Update the host side
    mesh.nodes_in_elem.update_host();

    Kokkos::fence();

    mesh.initialize_corners(num_elems * mesh.num_nodes_in_elem);

    // Build connectivity
    mesh.build_connectivity();
} // end build_3d_box


 /////////////////////////////////////////////////////////////////////////////
///
/// \fn build_2d_polar
///
/// \brief Builds an unstructured 2D polar mesh
///
/// \param Simulation mesh that is built
/// \param Element state data
/// \param Node state data
/// \param Corner state data
/// \param Simulation parameters
///
/////////////////////////////////////////////////////////////////////////////
void build_2d_polar(
    swage::Mesh& mesh,
    MPICArrayKokkos<double>& node_coords,
    double& inner_radius,
    double& outer_radius,
    double& start_angle,
    double& end_angle,
    int& num_elems_i,
    int& num_elems_j)
{
    printf("Creating a 2D polar mesh \n");

    int num_dim     = 2;

    const int num_points_i = num_elems_i + 1; // num points in x
    const int num_points_j = num_elems_j + 1; // num points in y

    const int num_nodes = num_points_i * num_points_j;

    const double dx = (outer_radius - inner_radius) / ((double)num_elems_i);  // len/(elems)



    // Convert degrees to radians
    start_angle = start_angle * M_PI / 180.0;
    end_angle = end_angle * M_PI / 180.0;
    const double dy = (end_angle - start_angle) / ((double)num_elems_j);  // len/(elems)

    const int num_elems = num_elems_i * num_elems_j;

    std::vector<double> origin(num_dim);

    for (int i = 0; i < num_dim; i++) { origin[i] = 0.0; }

    // --- 2D parameters ---
    // const int num_faces_in_elem  = 4;  // number of faces in elem
    // const int num_points_in_elem = 4;  // number of points in elem
    // const int num_points_in_face = 2;  // number of points in a face
    // const int num_edges_in_elem  = 4;  // number of edges in a elem

    // --- mesh node ordering ---
    // Convert ijk index system to the finite element numbering convention
    // for vertices in elem
    auto convert_point_number_in_quad = CArray<int>(4);
    convert_point_number_in_quad(0) = 0;
    convert_point_number_in_quad(1) = 1;
    convert_point_number_in_quad(2) = 3;
    convert_point_number_in_quad(3) = 2;

    // intialize node variables
    mesh.initialize_nodes(num_nodes);

    node_coords = MPICArrayKokkos<double>(num_nodes, num_dim, "node_coordinates");

    // populate the point data structures
    for (int j = 0; j < num_points_j; j++) {
        for (int i = 0; i < num_points_i; i++) {
            // global id for the point
            int node_gid = get_id(i, j, 0, num_points_i, num_points_j);

            double r_i     = inner_radius + (double)i * dx;
            double theta_j = start_angle + (double)j * dy;

            // store the point coordinates
            node_coords.host(node_gid, 0) = origin[0] + r_i * cos(theta_j);
            node_coords.host(node_gid, 1) = origin[1] + r_i * sin(theta_j);

            if(node_coords.host(node_gid, 0) < 0.0){
                throw std::runtime_error("**** NODE RADIUS FOR RZ MESH MUST BE POSITIVE ****");
            }

        } // end for i
    } // end for j


    node_coords.update_device();

    // initialize elem variables
    mesh.initialize_elems(num_elems, num_dim);

    // populate the elem center data structures
    for (int j = 0; j < num_elems_j; j++) {
        for (int i = 0; i < num_elems_i; i++) {
            // global id for the elem
            int elem_gid = get_id(i, j, 0, num_elems_i, num_elems_j);

            // store the point IDs for this elem where the range is
            // (i:i+1, j:j+1 for a linear quad
            int this_point = 0;

            for (int jcount = j; jcount <= j + 1; jcount++) {
                for (int icount = i; icount <= i + 1; icount++) {
                    // global id for the points
                    int node_gid = get_id(icount, jcount, 0, num_points_i, num_points_j);

                    // convert this_point index to the FE index convention
                    int this_index = convert_point_number_in_quad(this_point);

                    // store the points in this elem according the the finite
                    // element numbering convention
                    mesh.nodes_in_elem.host(elem_gid, this_index) = node_gid;

                    // increment the point counting index
                    this_point = this_point + 1;
                } // end for icount
            } // end for jcount
        } // end for i
    } // end for j

    // update device side
    mesh.nodes_in_elem.update_device();

    // intialize corner variables
    int num_corners = num_elems * mesh.num_nodes_in_elem;
    mesh.initialize_corners(num_corners);
    // corner.initialize(num_corners, num_dim);

    // Build connectivity
    mesh.build_connectivity();
} // end build_2d_box


/////////////////////////////////////////////////////////////////////////////
///
/// \fn write_vtu
///
/// \brief Writes a VTU (XML VTK) output file per MPI rank and a PVTU file
///        for parallel visualization in ParaView
///
/// Node PointData (when allocated): coords, coords_reference, coords_predicted,
/// velocity, force_external, inv_mass
///
/// Gauss point CellData (element-averaged when allocated): jacobian_determinant,
/// volume
///
/// Element CellData (when allocated): lambda
///
/// \param mesh mesh
/// \param node node state data
/// \param gauss_point gauss point state data
/// \param element element state data
/// \param num_quadrature_points_per_elem number of gauss points per element
/// \param rank MPI rank
/// \param comm MPI communicator
/// \param graphics_id output frame index
///
/////////////////////////////////////////////////////////////////////////////
void write_vtu(swage::Mesh& mesh,
               Node_t& node,
               GaussPoint_t& gauss_point,
               Element_t& element,
               size_t num_quadrature_points_per_elem,
               int rank,
               MPI_Comm comm,
               int graphics_id = 0)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);

    // ---- Update host data ----
    node.coords.update_host();
    if (node.coords_reference.size() > 0) {
        node.coords_reference.update_host();
    }
    if (node.coords_predicted.size() > 0) {
        node.coords_predicted.update_host();
    }
    if (node.velocity.size() > 0) {
        node.velocity.update_host();
    }
    if (node.force_external.size() > 0) {
        node.force_external.update_host();
    }
    if (node.inv_mass.size() > 0) {
        node.inv_mass.update_host();
    }
    if (gauss_point.jacobian_determinant.size() > 0) {
        gauss_point.jacobian_determinant.update_host();
    }
    if (gauss_point.volume.size() > 0) {
        gauss_point.volume.update_host();
    }
    if (element.lambda.size() > 0) {
        element.lambda.update_host();
    }
    Kokkos::fence();

    const bool has_coords_reference = node.coords_reference.size() > 0;
    const bool has_coords_predicted = node.coords_predicted.size() > 0;
    const bool has_velocity = node.velocity.size() > 0;
    const bool has_force_external = node.force_external.size() > 0;
    const bool has_inv_mass = node.inv_mass.size() > 0;
    const bool has_jacobian_determinant = gauss_point.jacobian_determinant.size() > 0;
    const bool has_volume = gauss_point.volume.size() > 0;
    const bool has_lambda = element.lambda.size() > 0;

    // short hand (write owned only)
    const size_t num_nodes = mesh.num_owned_nodes;
    const size_t num_elems = mesh.num_owned_elems;
    const size_t num_dims  = mesh.num_dims;

    size_t num_qp = num_quadrature_points_per_elem;
    if (num_qp == 0 && mesh.num_elems > 0) {
        if (has_jacobian_determinant) {
            num_qp = gauss_point.jacobian_determinant.size() / mesh.num_elems;
        } else if (has_volume) {
            num_qp = gauss_point.volume.size() / mesh.num_elems;
        }
    }

    CArray<double> avg_jacobian_determinant;
    if (has_jacobian_determinant && num_qp > 0) {
        avg_jacobian_determinant = CArray<double>(num_elems);
        for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
            double sum = 0.0;
            for (size_t quadrature_point_lid = 0; quadrature_point_lid < num_qp; ++quadrature_point_lid) {
                const size_t gauss_gid = elem_gid * num_qp + quadrature_point_lid;
                sum += gauss_point.jacobian_determinant.host(gauss_gid);
            }
            avg_jacobian_determinant(elem_gid) = sum / static_cast<double>(num_qp);
        }
    }

    CArray<double> avg_volume;
    if (has_volume && num_qp > 0) {
        avg_volume = CArray<double>(num_elems);
        for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
            double sum = 0.0;
            for (size_t quadrature_point_lid = 0; quadrature_point_lid < num_qp; ++quadrature_point_lid) {
                const size_t gauss_gid = elem_gid * num_qp + quadrature_point_lid;
                sum += gauss_point.volume.host(gauss_gid);
            }
            avg_volume(elem_gid) = sum / static_cast<double>(num_qp);
        }
    }

    CArray<double> elem_lambda;
    if (has_lambda) {
        elem_lambda = CArray<double>(num_elems);
        for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
            elem_lambda(elem_gid) = element.lambda.host(elem_gid);
        }
    }

    // File management
    char filename[200];
    int max_len = sizeof filename;
    int str_output_len;

    struct stat st;
    if (stat("vtk", &st) != 0) {
        system("mkdir vtk");
    }

    // Create VTU filename for this rank
    str_output_len = snprintf(filename, max_len, "vtk/Fierro.%05d_rank%d.vtu", graphics_id, rank);
    if (str_output_len >= max_len) { fputs("Filename length exceeded; string truncated", stderr); }

    FILE* vtu_file = fopen(filename, "w");
    if (!vtu_file) {
        std::cerr << "[rank " << rank << "] Failed to open VTU file: " << filename << std::endl;
        return;
    }

    // Write VTU XML header
    fprintf(vtu_file, "<?xml version=\"1.0\"?>\n");
    fprintf(vtu_file, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(vtu_file, "  <UnstructuredGrid>\n");
    fprintf(vtu_file, "    <Piece NumberOfPoints=\"%zu\" NumberOfCells=\"%zu\">\n", num_nodes, num_elems);


    std::cout<<"Writing points to VTU file for rank "<<rank<<std::endl;
    // Write Points (coordinates) — VTK expects 3 components; pad z for 2D
    fprintf(vtu_file, "      <Points>\n");
    fprintf(vtu_file, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for (size_t node_gid = 0; node_gid < num_nodes; node_gid++) {
        double x = node.coords.host(node_gid, 0);
        double y = (num_dims > 1) ? node.coords.host(node_gid, 1) : 0.0;
        double z = (num_dims > 2) ? node.coords.host(node_gid, 2) : 0.0;
        fprintf(vtu_file, "          %f %f %f\n", x, y, z);
    }
    fprintf(vtu_file, "        </DataArray>\n");
    fprintf(vtu_file, "      </Points>\n");

    // Write Cells (connectivity)
    std::cout<<"Writing cells to VTU file for rank "<<rank<<std::endl;
    fprintf(vtu_file, "      <Cells>\n");
    
    // Connectivity array - all node indices for all cells, space-separated
    fprintf(vtu_file, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
    int Pn_order = mesh.Pn;
    // WARNING: look into high-order Pn 2D elements with paraview
    int Pn_order_z = 0;
    if (num_dims == 3){
        Pn_order_z = Pn_order;
    }
    int order[3] = {Pn_order, Pn_order, Pn_order_z};

    CArray<int> convert_fierro_to_vtk((Pn_order+1)*(Pn_order+1)*(Pn_order+1));
    int this_node_lid = 0;
    for (int k = 0; k <= Pn_order_z; k++) {
        for (int j = 0; j <= Pn_order; j++) {
            for (int i = 0; i <= Pn_order; i++) {
                
                // convert this_point index to the FE index convention
                size_t vtk_index = PointIndexFromIJK(i, j, k, order);
                
                // store the points in this elem according the the finite
                // element numbering convention
                convert_fierro_to_vtk(vtk_index) = this_node_lid;
                // increment the point counting index
                this_node_lid ++;
                
            } // end for icount
        } // end for jcount
    }  // end for kcount


    
    // Write connectivity: all node IDs for all elements, space-separated
    std::cout<<"Writing connectivity to VTU file for rank "<<rank<<std::endl;
    for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
        fprintf(vtu_file, "          ");  // adding indentation before printing nodes in element
        
        int this_node_lid = 0;
        for (int k = 0; k <= Pn_order_z; k++) {
            for (int j = 0; j <= Pn_order; j++) {
                for (int i = 0; i <= Pn_order; i++) {
                    //  size_t node_lid = PointIndexFromIJK(i, j, k, order);
                    int node_lid = convert_fierro_to_vtk(this_node_lid);
                    if (num_dims == 2) node_lid = this_node_lid;
                    fprintf(vtu_file, " %zu", static_cast<unsigned long>(mesh.nodes_in_elem.host(elem_gid, node_lid)));
                    this_node_lid ++;
                }
            }
        }
    }


    fprintf(vtu_file, "\n");
    fprintf(vtu_file, "        </DataArray>\n");

    // Offsets array - cumulative index where each cell's connectivity ends
    fprintf(vtu_file, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
    int offset = 0;
    for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
        offset += static_cast<int>(mesh.num_nodes_in_elem);
        fprintf(vtu_file, " %d", offset);
    }
    fprintf(vtu_file, "\n");
    fprintf(vtu_file, "        </DataArray>\n");

    // Types array (72 = VTK_LAGRANGE_HEXAHEDRON)
    fprintf(vtu_file, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
        if(mesh.num_dims == 3) {
            fprintf(vtu_file, " 72");
        } else {
            fprintf(vtu_file, " 9");
        }
    }
    fprintf(vtu_file, "\n");
    fprintf(vtu_file, "        </DataArray>\n");
    fprintf(vtu_file, "      </Cells>\n");

    // Write PointData (node fields)
    fprintf(vtu_file, "      <PointData>\n");

    fprintf(vtu_file, "        <DataArray type=\"Float32\" Name=\"coords\" NumberOfComponents=\"%d\" format=\"ascii\">\n",
            static_cast<int>(num_dims));
    for (size_t node_gid = 0; node_gid < num_nodes; node_gid++) {
        for (int dim = 0; dim < static_cast<int>(num_dims); dim++) {
            fprintf(vtu_file, "          %f\n", node.coords.host(node_gid, dim));
        }
    }
    fprintf(vtu_file, "        </DataArray>\n");

    if (has_coords_reference) {
        fprintf(vtu_file, "        <DataArray type=\"Float32\" Name=\"coords_reference\" NumberOfComponents=\"%d\" format=\"ascii\">\n",
                static_cast<int>(num_dims));
        for (size_t node_gid = 0; node_gid < num_nodes; node_gid++) {
            for (int dim = 0; dim < static_cast<int>(num_dims); dim++) {
                fprintf(vtu_file, "          %f\n", node.coords_reference.host(node_gid, dim));
            }
        }
        fprintf(vtu_file, "        </DataArray>\n");
    }

    if (has_coords_predicted) {
        fprintf(vtu_file, "        <DataArray type=\"Float32\" Name=\"coords_predicted\" NumberOfComponents=\"%d\" format=\"ascii\">\n",
                static_cast<int>(num_dims));
        for (size_t node_gid = 0; node_gid < num_nodes; node_gid++) {
            for (int dim = 0; dim < static_cast<int>(num_dims); dim++) {
                fprintf(vtu_file, "          %f\n", node.coords_predicted.host(node_gid, dim));
            }
        }
        fprintf(vtu_file, "        </DataArray>\n");
    }

    if (has_velocity) {
        fprintf(vtu_file, "        <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"%d\" format=\"ascii\">\n",
                static_cast<int>(num_dims));
        for (size_t node_gid = 0; node_gid < num_nodes; node_gid++) {
            for (int dim = 0; dim < static_cast<int>(num_dims); dim++) {
                fprintf(vtu_file, "          %f\n", node.velocity.host(node_gid, dim));
            }
        }
        fprintf(vtu_file, "        </DataArray>\n");
    }

    if (has_force_external) {
        fprintf(vtu_file, "        <DataArray type=\"Float32\" Name=\"force_external\" NumberOfComponents=\"%d\" format=\"ascii\">\n",
                static_cast<int>(num_dims));
        for (size_t node_gid = 0; node_gid < num_nodes; node_gid++) {
            for (int dim = 0; dim < static_cast<int>(num_dims); dim++) {
                fprintf(vtu_file, "          %f\n", node.force_external.host(node_gid, dim));
            }
        }
        fprintf(vtu_file, "        </DataArray>\n");
    }

    if (has_inv_mass) {
        fprintf(vtu_file, "        <DataArray type=\"Float32\" Name=\"inv_mass\" format=\"ascii\">\n");
        for (size_t node_gid = 0; node_gid < num_nodes; node_gid++) {
            fprintf(vtu_file, "          %f\n", node.inv_mass.host(node_gid));
        }
        fprintf(vtu_file, "        </DataArray>\n");
    }
    fprintf(vtu_file, "      </PointData>\n");

    // Write CellData (element and averaged gauss point fields)
    fprintf(vtu_file, "      <CellData>\n");

    if (has_lambda) {
        fprintf(vtu_file, "        <DataArray type=\"Float32\" Name=\"lambda\" format=\"ascii\">\n");
        for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
            fprintf(vtu_file, "          %f\n", elem_lambda(elem_gid));
        }
        fprintf(vtu_file, "        </DataArray>\n");
    }

    if (has_jacobian_determinant && num_qp > 0) {
        fprintf(vtu_file, "        <DataArray type=\"Float32\" Name=\"jacobian_determinant\" format=\"ascii\">\n");
        for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
            fprintf(vtu_file, "          %f\n", avg_jacobian_determinant(elem_gid));
        }
        fprintf(vtu_file, "        </DataArray>\n");
    }

    if (has_volume && num_qp > 0) {
        fprintf(vtu_file, "        <DataArray type=\"Float32\" Name=\"volume\" format=\"ascii\">\n");
        for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
            fprintf(vtu_file, "          %f\n", avg_volume(elem_gid));
        }
        fprintf(vtu_file, "        </DataArray>\n");
    }

    fprintf(vtu_file, "      </CellData>\n");

    // Close VTU file
    fprintf(vtu_file, "    </Piece>\n");
    fprintf(vtu_file, "  </UnstructuredGrid>\n");
    fprintf(vtu_file, "</VTKFile>\n");
    fclose(vtu_file);

    // Write PVTU file (only rank 0, after all ranks have written their VTU files)
    MPI_Barrier(comm);
    
    if (rank == 0) {
        str_output_len = snprintf(filename, max_len, "vtk/Fierro.%05d.pvtu", graphics_id);
        if (str_output_len >= max_len) { fputs("Filename length exceeded; string truncated", stderr); }

        FILE* pvtu_file = fopen(filename, "w");
        if (!pvtu_file) {
            std::cerr << "[rank 0] Failed to open PVTU file: " << filename << std::endl;
            return;
        }

        // Write PVTU XML header
        fprintf(pvtu_file, "<?xml version=\"1.0\"?>\n");
        fprintf(pvtu_file, "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
        fprintf(pvtu_file, "  <PUnstructuredGrid GhostLevel=\"0\">\n");
        
        // Write PPoints
        fprintf(pvtu_file, "    <PPoints>\n");
        fprintf(pvtu_file, "      <PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n");
        fprintf(pvtu_file, "    </PPoints>\n");

        // Write PCells
        fprintf(pvtu_file, "    <PCells>\n");
        fprintf(pvtu_file, "      <PDataArray type=\"Int32\" Name=\"connectivity\"/>\n");
        fprintf(pvtu_file, "      <PDataArray type=\"Int32\" Name=\"offsets\"/>\n");
        fprintf(pvtu_file, "      <PDataArray type=\"UInt8\" Name=\"types\"/>\n");
        fprintf(pvtu_file, "    </PCells>\n");

        // Write PPointData
        fprintf(pvtu_file, "    <PPointData>\n");
        fprintf(pvtu_file, "      <PDataArray type=\"Float32\" Name=\"coords\" NumberOfComponents=\"%d\"/>\n",
                static_cast<int>(num_dims));
        if (has_coords_reference) {
            fprintf(pvtu_file, "      <PDataArray type=\"Float32\" Name=\"coords_reference\" NumberOfComponents=\"%d\"/>\n",
                    static_cast<int>(num_dims));
        }
        if (has_coords_predicted) {
            fprintf(pvtu_file, "      <PDataArray type=\"Float32\" Name=\"coords_predicted\" NumberOfComponents=\"%d\"/>\n",
                    static_cast<int>(num_dims));
        }
        if (has_velocity) {
            fprintf(pvtu_file, "      <PDataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"%d\"/>\n",
                    static_cast<int>(num_dims));
        }
        if (has_force_external) {
            fprintf(pvtu_file, "      <PDataArray type=\"Float32\" Name=\"force_external\" NumberOfComponents=\"%d\"/>\n",
                    static_cast<int>(num_dims));
        }
        if (has_inv_mass) {
            fprintf(pvtu_file, "      <PDataArray type=\"Float32\" Name=\"inv_mass\"/>\n");
        }
        fprintf(pvtu_file, "    </PPointData>\n");

        // Write PCellData
        fprintf(pvtu_file, "    <PCellData>\n");
        if (has_lambda) {
            fprintf(pvtu_file, "      <PDataArray type=\"Float32\" Name=\"lambda\"/>\n");
        }
        if (has_jacobian_determinant && num_qp > 0) {
            fprintf(pvtu_file, "      <PDataArray type=\"Float32\" Name=\"jacobian_determinant\"/>\n");
        }
        if (has_volume && num_qp > 0) {
            fprintf(pvtu_file, "      <PDataArray type=\"Float32\" Name=\"volume\"/>\n");
        }
        fprintf(pvtu_file, "    </PCellData>\n");

        // Write Piece references for each rank
        for (int r = 0; r < world_size; r++) {
            fprintf(pvtu_file, "    <Piece Source=\"Fierro.%05d_rank%d.vtu\"/>\n", graphics_id, r);
        }

        // Close PVTU file
        fprintf(pvtu_file, "  </PUnstructuredGrid>\n");
        fprintf(pvtu_file, "</VTKFile>\n");
        fclose(pvtu_file);
    }

} // end write_vtu




#endif