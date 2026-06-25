

#include <cmath> // for sin


#include "ELEMENTS.h"
#include "state.h"
#include "mesh_io.h"


KOKKOS_INLINE_FUNCTION
void get_shape_grad_at_center(int node_idx, double dN[3]) {
    // node_idx is 0-7, following the i + 2j + 4k mapping
    
    // Axis 0 (i / xi)
    // If bit 0 is set (1, 3, 5, 7), node is at i=1 (positive gradient)
    dN[0] = (node_idx & 1) ? 0.25 : -0.25;

    // Axis 1 (j / eta)
    // If bit 1 is set (2, 3, 6, 7), node is at j=1
    dN[1] = (node_idx & 2) ? 0.25 : -0.25;

    // Axis 2 (k / zeta)
    // If bit 2 is set (4, 5, 6, 7), node is at k=1
    dN[2] = (node_idx & 4) ? 0.25 : -0.25;
}


int main(int argc, char** argv) {

    MATAR_INITIALIZE(argc, argv);
    { // MATAR scope



    int num_dims = 3;
    int Pn_order = 0;


    // Mesh size for 3D box
    double origin[3] = {0.0, 0.0, 0.0};
    double length[3] = {1.0, 1.0, 1.0};
    int num_elems_dim[3] = {1, 1, 1};


    // --- XPBD Parameters ---
    double dt = 0.01;           
    const int num_iterations = 30;
    int repeats = 2;
    const int num_steps = 50*repeats;

    // Stiffness controls: Higher compliance = Softer
    // Note, these are starting value, they will be updates inthe relaxation loop
    float boundary_compliance = 0.0000001f; 
    float volume_compliance = 0.001f; 
    float edge_compliance = 0.00001f;
   
    

    // Initial mesh built on rank zero
    swage::Mesh mesh;
    node_t node;
    DCArrayKokkos<double> node_coords;

    GaussPoint_t gauss_point;

// ********************************************************  
//              Build the initial mesh
// ********************************************************  


    std::cout<<"Initializing mesh"<<std::endl;
    std::cout<<"Pn order: "<<Pn_order<<std::endl;
    build_3d_box(mesh,  node_coords, origin, length, num_elems_dim, Pn_order);
    
    // Read the mesh from a file
   //  read_vtk_mesh(mesh, node_coords, 3, "/home/jacobmoore/Desktop/repos/ELEMENTS/examples/mesh_form/stl_files/lattice3.vtk");

    std::cout << "Initial mesh has " << mesh.num_elems << " elements and " << mesh.num_nodes << " nodes" << std::endl;
    std::cout <<" Num_nodes_in_elem: " << mesh.num_nodes_in_elem << std::endl;
    std::cout <<" Num_dims: " << mesh.num_dims << std::endl;
    std::cout <<" Num_elems: " << mesh.num_elems << std::endl;
    std::cout <<" Num_nodes: " << mesh.num_nodes << std::endl;
    std::cout <<" Num_nodes_in_elem: " << mesh.num_nodes_in_elem << std::endl;
    std::cout <<" Num boundary patches: " << mesh.num_bdy_patches << std::endl;
    std::cout <<" Num boundary nodes: " << mesh.num_bdy_nodes << std::endl;

    
    std::vector<node_state> node_states = {node_state::coords, node_state::scalar_field, node_state::normal, node_state::closest_facet_point};
    node.initialize(mesh.num_nodes, mesh.num_dims, node_states);
    
    node.coords = node_coords;
    node.coords.update_device();

    std::vector<gauss_pt_state> gauss_pt_states = {gauss_pt_state::fields, gauss_pt_state::position};
    gauss_point.initialize(mesh.num_elems, mesh.num_dims, gauss_pt_states);



    // Read the STL file
    stl_data stl_data;
    binary_stl_reader("/home/jacobmoore/Desktop/repos/ELEMENTS/examples/mesh_form/stl_files/Sphere_medium.stl", stl_data);

    std::cout << "Number of facets: " << stl_data.num_facets << std::endl;
    
    stl_data.buildAABBTree();
    stl_data.verifyAABBTree();

    // Compute the gauss point position for each element
    std::cout<<"Computing gauss point position for each element"<<std::endl;
    FOR_ALL(elem_gid, 0, mesh.num_elems, {
        for (size_t gauss_pt_lid = 0; gauss_pt_lid < 1; gauss_pt_lid++) {
            size_t gauss_pt_gid = mesh.gauss_in_elem(elem_gid, gauss_pt_lid);

            // verify that the gauss gid matches the element git
            if(gauss_pt_gid != elem_gid) {
                std::cout<<"Error: gauss gid "<<gauss_pt_gid<<" does not match element gid "<<elem_gid<<std::endl;
                throw std::runtime_error("**** Error in Gauss Point ID ****");
            }
            
            double position[mesh.num_dims];

            position[0] = 0.0;
            position[1] = 0.0;
            position[2] = 0.0;
            
            // Average the node coordinates to get the gauss point position
            for( int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
                size_t node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
                
                
                for( int dim = 0; dim < mesh.num_dims; dim++) {
                    position[dim] += node.coords(node_gid, dim);
                }
            }

            for( int dim = 0; dim < mesh.num_dims; dim++) {
                gauss_point.position(gauss_pt_gid, dim) = position[dim] / mesh.num_nodes_in_elem;
            }
        }
    });
    gauss_point.position.update_host();

    node.normal.set_values(0.0);

    // Compute an outward normal for each surface node using the gauss point position
    std::cout<<"Computing outward normal for each surface node"<<std::endl;
    FOR_ALL(node_lid, 0, mesh.num_bdy_nodes, {
        size_t node_gid = mesh.bdy_nodes(node_lid);
        double normal[mesh.num_dims] = {0.0};

        for(int elem_lid = 0; elem_lid < mesh.num_corners_in_node(node_gid); elem_lid++) {
            size_t elem_gid = mesh.elems_in_node(node_gid, elem_lid);
            
            double vector[mesh.num_dims] = {0.0};
            for(int gauss_pt_lid = 0; gauss_pt_lid < mesh.num_gauss_in_elem; gauss_pt_lid++) {
                size_t gauss_pt_gid = mesh.gauss_in_elem(elem_gid, gauss_pt_lid);
                
                
                // Compute a vector pointing from the gauss point to the node
                for(int dim = 0; dim < mesh.num_dims; dim++) {
                    vector[dim] += node.coords(node_gid, dim) - gauss_point.position(gauss_pt_gid, dim);
                }
            }
            for(int dim = 0; dim < mesh.num_dims; dim++) {
                normal[dim] += vector[dim];
            }
        }
        for(int dim = 0; dim < mesh.num_dims; dim++) {
            normal[dim] /= (double)mesh.num_corners_in_node(node_gid);
        }

        // Normalize the normal vector
        double normal_magnitude = 0.0;
        for(int dim = 0; dim < mesh.num_dims; dim++) {
            normal_magnitude += normal[dim] * normal[dim];
        }
        normal_magnitude = std::sqrt(normal_magnitude);
        for(int dim = 0; dim < mesh.num_dims; dim++) {
            normal[dim] /= normal_magnitude;
        }

        for(int dim = 0; dim < mesh.num_dims; dim++) {
            node.normal(node_gid, dim) = normal[dim];
        }
    });


    // Compute the shortest distance between any node in the mesh
    std::cout<<"Computing the shortest distance between any node in the mesh"<<std::endl;
    double distance_lcl = 1.0e32;
    double min_distance_calc = 1.0e32;
    FOR_REDUCE_MIN(elem_gid, 0, mesh.num_elems, distance_lcl, {

        double coords0[24];  // element coords
        ViewCArrayKokkos<double> coords(coords0, 8, 3);

        double distance0[28];  // array for holding distances between each node
        ViewCArrayKokkos<double> dist(distance0, 28);

        // Getting the coordinates of the element
        for (size_t node_lid = 0; node_lid < 8; node_lid++) {
            for (size_t dim = 0; dim < mesh.num_dims; dim++) {
                coords(node_lid, dim) = node_coords(mesh.nodes_in_elem(elem_gid, node_lid), dim);
            } // end for dim
        } // end for loop over node_lid

        // loop conditions needed for distance calculation
        size_t countA = 0;
        size_t countB = 1;
        size_t a;
        size_t b;
        size_t loop = 0;

        // Only works for 3D
        // Solving for the magnitude of distance between each node
        for (size_t i = 0; i < 28; i++) {
            a = countA;
            b = countB;

            // returns magnitude of distance between each node, 28 total options
            dist(i) = fabs(sqrt((pow((coords(b, 0) - coords(a, 0)), 2.0)
            + pow((coords(b, 1) - coords(a, 1)), 2.0)
            + pow((coords(b, 2) - coords(a, 2)), 2.0))));

            countB++;
            countA++;

            // tricky indexing
            if (countB > 7) {
                loop++;
                countB = 1 + loop;
                countA = 0;
            }
        } // end for i

        double dist_min = dist(0);

        for (int i = 0; i < 28; ++i) {
            dist_min = fmin(dist(i), dist_min);
        }

        if (dist_min < distance_lcl) {
            distance_lcl = dist_min;
        }
    }, min_distance_calc);  // end parallel reduction
    Kokkos::fence();

    std::cout<<"Minimum distance: "<<min_distance_calc<<std::endl;
    
    DCArrayKokkos<double> nearest_facet_point(mesh.num_bdy_nodes, 3, "nearest_facet_point");

    // Walk over all of the nodes and compute the distance to the nearest faces
    FOR_ALL(node_lid, 0, mesh.num_bdy_nodes, {
        size_t node_gid = mesh.bdy_nodes(node_lid);
        float p[3];
        for(int i=0; i<3; i++) {
            p[i] = static_cast<float>(node.coords(node_gid, i));
        }
        
        float closest[3];
        for(int i=0; i<3; i++) closest[i] = 0.0f;
        float dist_sq = 0.0f;
        int facet_idx = stl_data.query_nearest_facet(p, dist_sq, closest);


        nearest_facet_point(node_lid, 0) = closest[0];
        nearest_facet_point(node_lid, 1) = closest[1];
        nearest_facet_point(node_lid, 2) = closest[2];
        // std::cout<<"Nearest facet: "<<facet_idx<<std::endl;
        // std::cout<<"Distance: "<<dist_sq<<std::endl;
        // std::cout<<"Closest point: "<<closest[0]<<", "<<closest[1]<<", "<<closest[2]<<std::endl;

        node.scalar_field(node_gid) = sqrt(dist_sq);

        node.closest_facet_point(node_gid, 0) = closest[0] - node.coords(node_gid, 0);
        node.closest_facet_point(node_gid, 1) = closest[1] - node.coords(node_gid, 1);
        node.closest_facet_point(node_gid, 2) = closest[2] - node.coords(node_gid, 2);
    });

    // Compute inverse reference matrix for each element

    DCArrayKokkos<float> inverse_reference_matrix(mesh.num_elems, 3, 3, "inverse_reference_matrix");
    DCArrayKokkos<float> reference_volume(mesh.num_elems, "reference_volume");


    DCArrayKokkos<float> corner_delta(mesh.num_corners, 3, "corner_delta");
    std::cout<<"NUmber of corners: "<<mesh.num_corners<<std::endl;

    // Index maps for the 8-node hex trilinear gradients at center
    // These represent the nodes on the "positive" and "negative" faces for each axis
    static const int pos_nodes[3][4] = {
        {1, 3, 5, 7}, // i = 1 face
        {2, 3, 6, 7}, // j = 1 face
        {4, 5, 6, 7}  // k = 1 face
    };
    static const int neg_nodes[3][4] = {
        {0, 2, 4, 6}, // i = 0 face
        {0, 1, 4, 5}, // j = 0 face
        {0, 1, 2, 3}  // k = 0 face
    };

    FOR_ALL(elem_gid, 0, mesh.num_elems, {
        
        // Compute the B matrix
        double B[3][3];
        for(int i=0; i<3; i++) {
            for(int j=0; j<3; j++) {
                B[i][j] = 0.0;
            }
        }

        for (int axis = 0; axis < 3; axis++) {
            for (int i = 0; i < 4; i++) {
                int p_idx = pos_nodes[axis][i];
                int p_node_gid = mesh.nodes_in_elem(elem_gid, p_idx);
                int n_idx = neg_nodes[axis][i];
                int n_node_gid = mesh.nodes_in_elem(elem_gid, n_idx);
                for (int d = 0; d < 3; d++) {
                    // The 0.25 weight averages the 4 edges along the specific axis
                    B[d][axis] += 0.25 * (node.coords(p_node_gid, d) - node.coords(n_node_gid, d));
                }
            }
        }

        // 2. Compute the Determinant of B
        // In the reference configuration, this represents the "local" volume scaling.
        double detB = B[0][0] * (B[1][1] * B[2][2] - B[1][2] * B[2][1]) -
        B[0][1] * (B[1][0] * B[2][2] - B[1][2] * B[2][0]) +
        B[0][2] * (B[1][0] * B[2][1] - B[1][1] * B[2][0]);

        // 3. The Sanity Check
        if (detB <= 1e-12) {
            // If this triggers, your pos_nodes/neg_nodes do not match build_3d_box
            // or the initial element is inverted/flat.
            printf("CRITICAL ERROR: Element %ld is invalid!\n", (long)elem_gid);
            printf("  detB: %e (Should be positive)\n", detB);

            // Output B matrix columns to see which axis is flipped
            for(int ax=0; ax<3; ax++) {
                printf("  Axis %d: [%f, %f, %f]\n", ax, B[0][ax], B[1][ax], B[2][ax]);
            }
        }

        // Invert the B matrix (3x3) and store in inverse_reference_matrix
        // Use Gauss-Jordan elimination for 3x3

        double inv[3][3];
        double temp[3][3];

        // Copy B to temp
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                temp[i][j] = B[i][j];

        // Initialize inv as identity
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                inv[i][j] = (i == j) ? 1.0 : 0.0;

        // Gauss-Jordan elimination
        for (int i = 0; i < 3; i++) {
            // Pivot
            double pivot = temp[i][i];
            if (fabs(pivot) < 1e-12) {
                // Handle singular matrix: fill with zeros
                for (int r = 0; r < 3; r++)
                    for (int c = 0; c < 3; c++)
                        inv[r][c] = 0.0;
                break;
            }
            for (int j = 0; j < 3; j++) {
                temp[i][j] /= pivot;
                inv[i][j] /= pivot;
            }
            // Eliminate other rows
            for (int row = 0; row < 3; row++) {
                if (row == i) continue;
                double factor = temp[row][i];
                for (int col = 0; col < 3; col++) {
                    temp[row][col] -= factor * temp[i][col];
                    inv[row][col]  -= factor * inv[i][col];
                }
            }
        }

        // Store result into inverse_reference_matrix
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                inverse_reference_matrix(elem_gid, i, j) = static_cast<float>(inv[i][j]);
            }
        }

        // Compute the volume of the element
        const size_t num_nodes = 8;

        double x_array[8];
        double y_array[8];
        double z_array[8];

        // x, y, z coordinates of elem vertices
        auto x = ViewCArrayKokkos<double>(x_array, num_nodes);
        auto y = ViewCArrayKokkos<double>(y_array, num_nodes);
        auto z = ViewCArrayKokkos<double>(z_array, num_nodes);

        // get the coordinates of the nodes(rk,elem,node) in this element
        for (int node_lid = 0; node_lid < num_nodes; node_lid++) {
            size_t node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
            x(node_lid) = node_coords(node_gid, 0);
            y(node_lid) = node_coords(node_gid, 1);
            z(node_lid) = node_coords(node_gid, 2);
        }     // end for

        double twelth = 1. / 12.;

        // element volume
        double volume =
            (x(1) * (y(2) * (-z(0) + z(3)) + y(4) * (z(0) - z(5)) + y(0) * (z(3) + z(2) - z(4) - z(5)) + y(7) * (-z(3) + z(5)) + y(5) * (z(0) - z(3) + z(4) - z(7)) + y(3) * (-z(0) - z(2) + z(5) + z(7))) +
            x(6) * (y(0) * (-z(2) + z(4)) + y(7) * (z(3) + z(2) - z(4) - z(5)) + y(3) * (z(2) - z(7)) + y(2) * (z(0) - z(3) + z(4) - z(7)) + y(5) * (-z(4) + z(7)) + y(4) * (-z(0) - z(2) + z(5) + z(7))) +
            x(2) * (y(1) * (z(0) - z(3)) + y(6) * (-z(0) + z(3) - z(4) + z(7)) + y(7) * (z(3) - z(6)) + y(3) * (z(0) + z(1) - z(7) - z(6)) + y(4) * (-z(0) + z(6)) + y(0) * (-z(1) - z(3) + z(4) + z(6))) +
            x(5) * (y(0) * (z(1) - z(4)) + y(6) * (z(4) - z(7)) + y(3) * (-z(1) + z(7)) + y(1) * (-z(0) + z(3) - z(4) + z(7)) + y(4) * (z(0) + z(1) - z(7) - z(6)) + y(7) * (-z(1) - z(3) + z(4) + z(6))) +
            x(7) * (y(1) * (z(3) - z(5)) + y(6) * (-z(3) - z(2) + z(4) + z(5)) + y(5) * (z(1) + z(3) - z(4) - z(6)) + y(4) * (z(5) - z(6)) + y(2) * (-z(3) + z(6)) + y(3) * (-z(1) + z(2) - z(5) + z(6))) +
            x(0) * (y(3) * (z(1) - z(2)) + y(6) * (z(2) - z(4)) + y(5) * (-z(1) + z(4)) + y(1) * (-z(3) - z(2) + z(4) + z(5)) + y(2) * (z(1) + z(3) - z(4) - z(6)) + y(4) * (-z(1) + z(2) - z(5) + z(6))) +
            x(3) * (y(0) * (-z(1) + z(2)) + y(5) * (z(1) - z(7)) + y(1) * (z(0) + z(2) - z(5) - z(7)) + y(6) * (-z(2) + z(7)) + y(7) * (z(1) - z(2) + z(5) - z(6)) + y(2) * (-z(0) - z(1) + z(7) + z(6))) +
            x(4) *
            (y(1) * (-z(0) + z(5)) + y(6) * (z(0) + z(2) - z(5) - z(7)) + y(2) * (z(0) - z(6)) + y(0) * (z(1) - z(2) + z(5) - z(6)) + y(7) * (-z(5) + z(6)) + y(5) * (-z(0) - z(1) + z(7) + z(6)))) *
            twelth;

        reference_volume(elem_gid) = static_cast<float>(volume);
    });


    
    

    DCArrayKokkos<double> volume_lambda(mesh.num_elems, "volume_lambda");
    DCArrayKokkos<double> boundary_lambda(mesh.num_bdy_nodes, "boundary_lambda");
    DCArrayKokkos<double> element_lambda_iso(mesh.num_elems, "element_lambda_iso");
    DCArrayKokkos<double> edge_lambda(mesh.num_nodes, "edge_lambda");
    DCArrayKokkos<double> element_lambda_diag(mesh.num_elems, "element_lambda_diag");
    
    
    // Reset lambdas at the START of the relaxation
    volume_lambda.set_values(0.0);
    boundary_lambda.set_values(0.0);
    element_lambda_iso.set_values(0.0);
    edge_lambda.set_values(0.0);
    element_lambda_diag.set_values(0.0);

    std::cout << "Starting XPBD Relaxation..." << std::endl;



    int step = 0;
    std::cout<<"Writing VTU file "<<std::endl;
    int rank = 0;
    write_vtu(mesh, node, gauss_point, rank, MPI_COMM_WORLD, step, 0.0);


    return 0;

    for(step = 1; step < num_steps/repeats; step++) {

        float cutoff = 0.8;


        // Linearly drop boundary_compliance to zero after cutoff% of steps
        if(step > cutoff * num_steps/repeats) {
            // fraction = (step - cutoff*num_steps/repeats) / (0.2*num_steps/repeats), in [0,1] after cutoff% progress
            float t_frac = float(step - cutoff * num_steps/repeats) / float(0.2 * num_steps/repeats);
            if (t_frac > 1.0f) t_frac = 1.0f;
            if (t_frac < 0.0f) t_frac = 0.0f;
            // Linear decay to zero
            boundary_compliance = boundary_compliance * (1.0f - t_frac);
            // Linearly drop volume_compliance to zero after cutoff% of steps
            volume_compliance = 0.0001f + volume_compliance * (1.0f - t_frac);
            // Linearly drop edge_compliance to zero after cutoff% of steps
            edge_compliance = 0.0001f + edge_compliance * (1.0f - t_frac);
        }

        double sub_step = dt/num_steps;
        
        double alpha_tilde_bdy = boundary_compliance / (sub_step * sub_step);
        double alpha_tilde_vol =  volume_compliance / (sub_step * sub_step);
        double alpha_tilde_edge =  edge_compliance / (sub_step * sub_step);
        
        volume_lambda.set_values(0.0);
        boundary_lambda.set_values(0.0);
        element_lambda_iso.set_values(0.0);
        element_lambda_diag.set_values(0.0);
        
        for(int iter = 0; iter < num_iterations; iter++) {
            
            // --- 1. Boundary Snap Constraint ---
            FOR_ALL(node_lid, 0, mesh.num_bdy_nodes, {
                size_t node_gid = mesh.bdy_nodes(node_lid);
                
                float p[3];
                for(int i=0; i<3; i++) p[i] = static_cast<float>(node.coords(node_gid, i));
                
                float closest[3];
                float dist_sq;
                // int facet_idx = stl_data.query_nearest_facet(p, dist_sq, closest);
                
                // float dist = sqrtf(dist_sq);

                // Calculate the distance using the previously saved nearest facet point
                dist_sq = 0.0f;
                for(int d = 0; d < 3; d++) {
                    dist_sq += (node.coords(node_gid, d) - nearest_facet_point(node_lid, d))*(node.coords(node_gid, d) - nearest_facet_point(node_lid, d));
                }

                float dist = sqrt(dist_sq);

                if (dist < 1e-9) return; 

                node.scalar_field(node_gid) = sqrt(dist_sq);

                // node.closest_facet_point(node_gid, 0) = closest[0] - node.coords(node_gid, 0);
                // node.closest_facet_point(node_gid, 1) = closest[1] - node.coords(node_gid, 1);
                // node.closest_facet_point(node_gid, 2) = closest[2] - node.coords(node_gid, 2);

                closest[0] = nearest_facet_point(node_lid, 0);
                closest[1] = nearest_facet_point(node_lid, 1);
                closest[2] = nearest_facet_point(node_lid, 2);

                node.closest_facet_point(node_gid, 0) = nearest_facet_point(node_lid, 0) - node.coords(node_gid, 0);
                node.closest_facet_point(node_gid, 1) = nearest_facet_point(node_lid, 1) - node.coords(node_gid, 1);
                node.closest_facet_point(node_gid, 2) = nearest_facet_point(node_lid, 2) - node.coords(node_gid, 2);


                if (dist < 1e-9) return; 
            
                float grad[3];
                for(int i=0; i<3; i++) grad[i] = (p[i] - closest[i]) / dist;
            
                // XPBD update for boundary
                double old_lambda = boundary_lambda(node_lid);
                double delta_lambda = (-dist - alpha_tilde_bdy * old_lambda) / (1.0 + alpha_tilde_bdy);
                boundary_lambda(node_lid) += delta_lambda;
                
                for(int d=0; d<3; d++) {
                    node.coords(node_gid, d) += delta_lambda * grad[d];
                }
            });


            // --- 4. Edge Length Constraint (Relaxation) ---
            FOR_ALL(node_gid, 0, mesh.num_nodes, {

                size_t node_a = node_gid;
                for(int node_lid = 0; node_lid < mesh.num_nodes_in_node(node_gid); node_lid++){


                    size_t node_b = mesh.nodes_in_node(node_gid, node_lid);
                    double L0 = min_distance_calc;

                    double dx[3];
                    double dist_sq = 0.0;
                    for(int d=0; d<3; d++) {
                        dx[d] = node.coords(node_a, d) - node.coords(node_b, d);
                        dist_sq += dx[d] * dx[d];
                    }
                    double dist = sqrt(dist_sq);
                    
                    if (dist < 1e-12) return; // Prevent division by zero

                    double C = dist - L0;
                    
                    // XPBD Denominator: sum(w*grad^2) + alpha_tilde
                    // w_a = 1.0, w_b = 1.0 -> denominator = 2.0 + alpha_tilde
                    double denominator = 2.0 + alpha_tilde_edge;

                    double old_lambda = edge_lambda(node_gid);
                    double delta_lambda = (-C - alpha_tilde_edge * old_lambda) / denominator;
                    edge_lambda(node_gid) += delta_lambda/mesh.num_nodes_in_node(node_gid);

                    // Unit gradient vector n
                    double n[3];
                    for(int d=0; d<3; d++) n[d] = dx[d] / dist;

                    // Update positions
                    for(int d=0; d<3; d++) {
                        double corr = delta_lambda * n[d];
                        Kokkos::atomic_add(&node.coords(node_a, d), corr);
                        Kokkos::atomic_add(&node.coords(node_b, d), -corr);
                    }
                }
            });

            // Element internal Diagonal Constraints
            FOR_ALL(elem_gid, 0, mesh.num_elems, {

                int num_sets = 4;

                int set_[num_sets*2];
                auto set = ViewCArrayKokkos<int>(set_, num_sets, 2);


                // Daigonals within the elelent
                set(0,0) = 0;
                set(0,1) = 7;
                set(1,0) = 1;
                set(1,1) = 6;
                set(2,0) = 2;
                set(2,1) = 5;
                set(3,0) = 3;
                set(3,1) = 4;

                double target = min_distance_calc*min_distance_calc + min_distance_calc*min_distance_calc + min_distance_calc*min_distance_calc;
                target = sqrt(target);

                for(int i = 0; i < 4; i++){

                    size_t node_a = mesh.nodes_in_elem(elem_gid, set(i,0));
                    size_t node_b = mesh.nodes_in_elem(elem_gid, set(i,1));

                    
                    double L0 = target;

                    double dx[3];
                    double dist_sq = 0.0;
                    for(int d=0; d<3; d++) {
                        dx[d] = node.coords(node_a, d) - node.coords(node_b, d);
                        dist_sq += dx[d] * dx[d];
                    }
                    double dist = sqrt(dist_sq);
                    
                    if (dist < 1e-12) return; // Prevent division by zero

                    double C = dist - L0;
                    
                    // XPBD Denominator: sum(w*grad^2) + alpha_tilde
                    // w_a = 1.0, w_b = 1.0 -> denominator = 2.0 + alpha_tilde
                    double denominator = 2.0 + alpha_tilde_edge;

                    double old_lambda = element_lambda_diag(elem_gid);
                    double delta_lambda = (-C - alpha_tilde_edge * old_lambda) / denominator;
                    element_lambda_diag(elem_gid) += delta_lambda / 4;

                    // Unit gradient vector n
                    double n[3];
                    for(int d=0; d<3; d++) n[d] = dx[d] / dist;

                    // Update positions
                    for(int d=0; d<3; d++) {
                        double corr = delta_lambda * n[d];
                        Kokkos::atomic_add(&node.coords(node_a, d), corr);
                        Kokkos::atomic_add(&node.coords(node_b, d), -corr);
                    }
                }
            
            });

            // Element face Diagonal Constraints
            FOR_ALL(elem_gid, 0, mesh.num_elems, {

                int num_sets = 8;

                int set_[num_sets*2];
                auto set = ViewCArrayKokkos<int>(set_, num_sets, 2);

                // Diagonals on the element face
                // -i face face
                set(0,0) = 0;
                set(0,1) = 6;
                set(1,0) = 3;
                set(1,1) = 4;
                // +i face face
                set(2,0) = 1;
                set(2,1) = 7;
                set(3,0) = 5;
                set(3,1) = 3;
                // -j face face
                set(4,0) = 0;
                set(4,1) = 5;
                set(5,0) = 1;
                set(5,1) = 4;
                // +j face face
                set(6,0) = 2;
                set(6,1) = 7;
                set(7,0) = 3;
                set(7,1) = 6;


                double target = min_distance_calc*min_distance_calc + min_distance_calc*min_distance_calc;
                target = sqrt(target);

                for(int i = 0; i < num_sets; i++){

                    size_t node_a = mesh.nodes_in_elem(elem_gid, set(i,0));
                    size_t node_b = mesh.nodes_in_elem(elem_gid, set(i,1));

                    
                    double L0 = target;

                    double dx[3];
                    double dist_sq = 0.0;
                    for(int d=0; d<3; d++) {
                        dx[d] = node.coords(node_a, d) - node.coords(node_b, d);
                        dist_sq += dx[d] * dx[d];
                    }
                    double dist = sqrt(dist_sq);
                    
                    if (dist < 1e-12) return; // Prevent division by zero

                    double C = dist - L0;
                    
                    // XPBD Denominator: sum(w*grad^2) + alpha_tilde
                    // w_a = 1.0, w_b = 1.0 -> denominator = 2.0 + alpha_tilde
                    double denominator = 2.0 + alpha_tilde_edge;

                    double old_lambda = element_lambda_diag(elem_gid);
                    double delta_lambda = (-C - alpha_tilde_edge * old_lambda) / denominator;
                    element_lambda_diag(elem_gid) += delta_lambda / 4;

                    // Unit gradient vector n
                    double n[3];
                    for(int d=0; d<3; d++) n[d] = dx[d] / dist;

                    // Update positions
                    for(int d=0; d<3; d++) {
                        double corr = delta_lambda * n[d];
                        Kokkos::atomic_add(&node.coords(node_a, d), corr);
                        Kokkos::atomic_add(&node.coords(node_b, d), -corr);
                    }
                }
            
            });

            // --- 2. Element Volumetric Constraint ---
            // FOR_ALL(elem_gid, 0, mesh.num_elems, {
            // //for(int elem_gid=0; elem_gid<mesh.num_elems; elem_gid++) {
            //     double x1[8][3];
            //     for (int n = 0; n < 8; n++) {
            //         size_t node_gid = mesh.nodes_in_elem(elem_gid, n);
            //         for (int d = 0; d < 3; d++) x1[n][d] = node.coords(node_gid, d);
            //     }
            
            //     // Shape Matrix A (using i,j,k mapping)
            //     double A_mem[9];
            //     ViewCArrayKokkos<double> A(A_mem, 3, 3);
            //     A.set_values(0.0);

            //     for (int axis = 0; axis < 3; axis++) {
            //         for (int i = 0; i < 4; i++) {
            //             for (int d = 0; d < 3; d++) {
            //                 A(d, axis) += 0.25 * (x1[pos_nodes[axis][i]][d] - x1[neg_nodes[axis][i]][d]);
            //             }
            //         }
            //     }
            
            //     // Deformation Gradient F = A * invB
            //     double F_mem[9];
            //     ViewCArrayKokkos<double> F(F_mem, 3, 3);
            //     F.set_values(0.0);
                
            //     for(int i=0; i<3; i++)
            //         for(int j=0; j<3; j++)
            //             for(int k=0; k<3; k++)
            //                 F(i, j) += A(i, k) * inverse_reference_matrix(elem_gid, k, j);
            
            //     // C = det(F) - 1
            //     double detF = F(0,0)*(F(1,1)*F(2,2) - F(1,2)*F(2,1))
            //                 - F(0,1)*(F(1,0)*F(2,2) - F(1,2)*F(2,0))
            //                 + F(0,2)*(F(1,0)*F(2,1) - F(1,1)*F(2,0));
                
            //     double C = detF - 1.0;
            
            //     // Gradient of C (Cofactor matrix H)
            //     double H_mem[9];
            //     ViewCArrayKokkos<double> H(H_mem, 3, 3);
            //     H(0, 0) = F(1, 1)*F(2, 2) - F(1, 2)*F(2, 1); H(0, 1) = F(1, 2)*F(2, 0) - F(1, 0)*F(2, 2); H(0, 2) = F(1, 0)*F(2, 1) - F(1, 1)*F(2, 0);
            //     H(1, 0) = F(0, 2)*F(2, 1) - F(0, 1)*F(2, 2); H(1, 1) = F(0, 0)*F(2, 2) - F(0, 2)*F(2, 0); H(1, 2) = F(0, 1)*F(2, 0) - F(0, 0)*F(2, 1);
            //     H(2, 0) = F(0, 1)*F(1, 2) - F(0, 2)*F(1, 1); H(2, 1) = F(0, 2)*F(1, 0) - F(0, 0)*F(1, 2); H(2, 2) = F(0, 0)*F(1, 1) - F(0, 1)*F(1, 0);

            //     double grad_mem[24];
            //     ViewCArrayKokkos<double> grad(grad_mem, 8, 3);
            //     grad.set_values(0.0);
            //     double sum_grad_sq = 0;

            //     for (int n = 0; n < 8; n++) {
            //         double dN_dxi[3];
            //         get_shape_grad_at_center(n, dN_dxi); 
            //         for (int d = 0; d < 3; d++) { 
            //             for (int i = 0; i < 3; i++) {
            //                 for (int j = 0; j < 3; j++) {
            //                     // Chain rule through the reference configuration
            //                     grad(n, d) += H(d, i) * inverse_reference_matrix(elem_gid, j, i) * dN_dxi[j];
            //                 }
            //             }
            //             sum_grad_sq += grad(n, d) * grad(n, d);
            //         }
            //     }

            //     // XPBD update logic for elements
            //     double old_lambda = volume_lambda(elem_gid);
            //     double delta_lambda = (-C - alpha_tilde_vol * old_lambda) / (sum_grad_sq + alpha_tilde_vol);
            //     volume_lambda(elem_gid) += delta_lambda;

            //     for (int n = 0; n < 8; n++) {
            //         int corner_gid = mesh.corners_in_elem(elem_gid, n);
            //         int node_gid = mesh.nodes_in_elem(elem_gid, n);

            //         for (int d = 0; d < 3; d++) {
            //             corner_delta(corner_gid, d) = (float)(delta_lambda * grad(n, d));
            //             Kokkos::atomic_add(&node.coords(node_gid, d), corner_delta(corner_gid, d));
            //         }
            //     }

                

            // });
            //}

            // Compute the volume of the element
            FOR_ALL(elem_gid, 0, mesh.num_elems, {

                const size_t num_nodes = 8;

                double x_array[8];
                double y_array[8];
                double z_array[8];

                // x, y, z coordinates of elem vertices
                auto x = ViewCArrayKokkos<double>(x_array, num_nodes);
                auto y = ViewCArrayKokkos<double>(y_array, num_nodes);
                auto z = ViewCArrayKokkos<double>(z_array, num_nodes);

                // get the coordinates of the nodes(rk,elem,node) in this element
                for (int node_lid = 0; node_lid < num_nodes; node_lid++) {
                    size_t node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
                    x(node_lid) = node_coords(node_gid, 0);
                    y(node_lid) = node_coords(node_gid, 1);
                    z(node_lid) = node_coords(node_gid, 2);
                }     // end for

                double twelth = 1. / 12.;

                // element volume
                double volume =
                    (x(1) * (y(2) * (-z(0) + z(3)) + y(4) * (z(0) - z(5)) + y(0) * (z(3) + z(2) - z(4) - z(5)) + y(7) * (-z(3) + z(5)) + y(5) * (z(0) - z(3) + z(4) - z(7)) + y(3) * (-z(0) - z(2) + z(5) + z(7))) +
                    x(6) * (y(0) * (-z(2) + z(4)) + y(7) * (z(3) + z(2) - z(4) - z(5)) + y(3) * (z(2) - z(7)) + y(2) * (z(0) - z(3) + z(4) - z(7)) + y(5) * (-z(4) + z(7)) + y(4) * (-z(0) - z(2) + z(5) + z(7))) +
                    x(2) * (y(1) * (z(0) - z(3)) + y(6) * (-z(0) + z(3) - z(4) + z(7)) + y(7) * (z(3) - z(6)) + y(3) * (z(0) + z(1) - z(7) - z(6)) + y(4) * (-z(0) + z(6)) + y(0) * (-z(1) - z(3) + z(4) + z(6))) +
                    x(5) * (y(0) * (z(1) - z(4)) + y(6) * (z(4) - z(7)) + y(3) * (-z(1) + z(7)) + y(1) * (-z(0) + z(3) - z(4) + z(7)) + y(4) * (z(0) + z(1) - z(7) - z(6)) + y(7) * (-z(1) - z(3) + z(4) + z(6))) +
                    x(7) * (y(1) * (z(3) - z(5)) + y(6) * (-z(3) - z(2) + z(4) + z(5)) + y(5) * (z(1) + z(3) - z(4) - z(6)) + y(4) * (z(5) - z(6)) + y(2) * (-z(3) + z(6)) + y(3) * (-z(1) + z(2) - z(5) + z(6))) +
                    x(0) * (y(3) * (z(1) - z(2)) + y(6) * (z(2) - z(4)) + y(5) * (-z(1) + z(4)) + y(1) * (-z(3) - z(2) + z(4) + z(5)) + y(2) * (z(1) + z(3) - z(4) - z(6)) + y(4) * (-z(1) + z(2) - z(5) + z(6))) +
                    x(3) * (y(0) * (-z(1) + z(2)) + y(5) * (z(1) - z(7)) + y(1) * (z(0) + z(2) - z(5) - z(7)) + y(6) * (-z(2) + z(7)) + y(7) * (z(1) - z(2) + z(5) - z(6)) + y(2) * (-z(0) - z(1) + z(7) + z(6))) +
                    x(4) *
                    (y(1) * (-z(0) + z(5)) + y(6) * (z(0) + z(2) - z(5) - z(7)) + y(2) * (z(0) - z(6)) + y(0) * (z(1) - z(2) + z(5) - z(6)) + y(7) * (-z(5) + z(6)) + y(5) * (-z(0) - z(1) + z(7) + z(6)))) *
                    twelth;

                gauss_point.fields(elem_gid) = static_cast<float>(volume);
            });
            

        }

        std::cout<<"Writing VTU file for step "<<step<<std::endl;
        int rank = 0;
        double time_value = static_cast<double>(step) * dt;
        write_vtu(mesh, node, gauss_point, rank, MPI_COMM_WORLD, step, time_value);
    }


    // Compute the shortest distance between any node in the mesh
    std::cout<<"Computing the shortest distance between any node in the mesh"<<std::endl;
    distance_lcl = 1.0e32;
    min_distance_calc = 1.0e32;
    FOR_REDUCE_MIN(elem_gid, 0, mesh.num_elems, distance_lcl, {

        double coords0[24];  // element coords
        ViewCArrayKokkos<double> coords(coords0, 8, 3);

        double distance0[28];  // array for holding distances between each node
        ViewCArrayKokkos<double> dist(distance0, 28);

        // Getting the coordinates of the element
        for (size_t node_lid = 0; node_lid < 8; node_lid++) {
            for (size_t dim = 0; dim < mesh.num_dims; dim++) {
                coords(node_lid, dim) = node_coords(mesh.nodes_in_elem(elem_gid, node_lid), dim);
            } // end for dim
        } // end for loop over node_lid

        // loop conditions needed for distance calculation
        size_t countA = 0;
        size_t countB = 1;
        size_t a;
        size_t b;
        size_t loop = 0;

        // Only works for 3D
        // Solving for the magnitude of distance between each node
        for (size_t i = 0; i < 28; i++) {
            a = countA;
            b = countB;

            // returns magnitude of distance between each node, 28 total options
            dist(i) = fabs(sqrt((pow((coords(b, 0) - coords(a, 0)), 2.0)
            + pow((coords(b, 1) - coords(a, 1)), 2.0)
            + pow((coords(b, 2) - coords(a, 2)), 2.0))));

            countB++;
            countA++;

            // tricky indexing
            if (countB > 7) {
                loop++;
                countB = 1 + loop;
                countA = 0;
            }
        } // end for i

        double dist_min = dist(0);

        for (int i = 0; i < 28; ++i) {
            dist_min = fmin(dist(i), dist_min);
        }

        if (dist_min < distance_lcl) {
            distance_lcl = dist_min;
        }
    }, min_distance_calc);  // end parallel reduction
    Kokkos::fence();

    std::cout<<"Minimum distance: "<<min_distance_calc<<std::endl;



    for(step = num_steps/repeats; step < num_steps; step++) {


        boundary_compliance = 0.0;
        // Linearly drop volume_compliance to zero after cutoff% of steps
        // volume_compliance = volume_compliance / 1.1;
        // Linearly drop edge_compliance to zero after cutoff% of steps
        // edge_compliance = edge_compliance / 1.1;


        double sub_step = dt/num_steps;
        
        double alpha_tilde_bdy = boundary_compliance / (sub_step * sub_step);
        double alpha_tilde_vol =  volume_compliance / (sub_step * sub_step);
        double alpha_tilde_edge =  edge_compliance / (sub_step * sub_step);
        
        volume_lambda.set_values(0.0);
        boundary_lambda.set_values(0.0);
        element_lambda_iso.set_values(0.0);
        element_lambda_diag.set_values(0.0);
        
        for(int iter = 0; iter < num_iterations; iter++) {
        
           
            // --- 2. Element Volumetric Constraint ---
            // FOR_ALL(elem_gid, 0, mesh.num_elems, {
            //     // for(int elem_gid=0; elem_gid<mesh.num_elems; elem_gid++) {
            //     double x1[8][3];
            //     for (int n = 0; n < 8; n++) {
            //         size_t node_gid = mesh.nodes_in_elem(elem_gid, n);
            //         for (int d = 0; d < 3; d++) x1[n][d] = node.coords(node_gid, d);
            //     }
            
            //     // Shape Matrix A (using i,j,k mapping)
            //     double A_mem[9];
            //     ViewCArrayKokkos<double> A(A_mem, 3, 3);
            //     A.set_values(0.0);

            //     for (int axis = 0; axis < 3; axis++) {
            //         for (int i = 0; i < 4; i++) {
            //             for (int d = 0; d < 3; d++) {
            //                 A(d, axis) += 0.25 * (x1[pos_nodes[axis][i]][d] - x1[neg_nodes[axis][i]][d]);
            //             }
            //         }
            //     }
            
            //     // Deformation Gradient F = A * invB
            //     double F_mem[9];
            //     ViewCArrayKokkos<double> F(F_mem, 3, 3);
            //     F.set_values(0.0);
                
            //     for(int i=0; i<3; i++)
            //         for(int j=0; j<3; j++)
            //             for(int k=0; k<3; k++)
            //                 F(i, j) += A(i, k) * inverse_reference_matrix(elem_gid, k, j);
            
            //     // C = det(F) - 1
            //     double detF = F(0,0)*(F(1,1)*F(2,2) - F(1,2)*F(2,1))
            //                 - F(0,1)*(F(1,0)*F(2,2) - F(1,2)*F(2,0))
            //                 + F(0,2)*(F(1,0)*F(2,1) - F(1,1)*F(2,0));
                
            //     double C = detF - 1.0;
            
            //     // Gradient of C (Cofactor matrix H)
            //     double H_mem[9];
            //     ViewCArrayKokkos<double> H(H_mem, 3, 3);
            //     H(0, 0) = F(1, 1)*F(2, 2) - F(1, 2)*F(2, 1); H(0, 1) = F(1, 2)*F(2, 0) - F(1, 0)*F(2, 2); H(0, 2) = F(1, 0)*F(2, 1) - F(1, 1)*F(2, 0);
            //     H(1, 0) = F(0, 2)*F(2, 1) - F(0, 1)*F(2, 2); H(1, 1) = F(0, 0)*F(2, 2) - F(0, 2)*F(2, 0); H(1, 2) = F(0, 1)*F(2, 0) - F(0, 0)*F(2, 1);
            //     H(2, 0) = F(0, 1)*F(1, 2) - F(0, 2)*F(1, 1); H(2, 1) = F(0, 2)*F(1, 0) - F(0, 0)*F(1, 2); H(2, 2) = F(0, 0)*F(1, 1) - F(0, 1)*F(1, 0);

            //     double grad_mem[24];
            //     ViewCArrayKokkos<double> grad(grad_mem, 8, 3);
            //     grad.set_values(0.0);
            //     double sum_grad_sq = 0;

            //     for (int n = 0; n < 8; n++) {
            //         double dN_dxi[3];
            //         get_shape_grad_at_center(n, dN_dxi); 
            //         for (int d = 0; d < 3; d++) { 
            //             for (int i = 0; i < 3; i++) {
            //                 for (int j = 0; j < 3; j++) {
            //                     // Chain rule through the reference configuration
            //                     grad(n, d) += H(d, i) * inverse_reference_matrix(elem_gid, j, i) * dN_dxi[j];
            //                 }
            //             }
            //             sum_grad_sq += grad(n, d) * grad(n, d);
            //         }
            //     }

            //     // XPBD update logic for elements
            //     double old_lambda = volume_lambda(elem_gid);
            //     double delta_lambda = (-C - alpha_tilde_vol * old_lambda) / (sum_grad_sq + alpha_tilde_vol);
            //     volume_lambda(elem_gid) += delta_lambda;

            //     for (int n = 0; n < 8; n++) {
            //         int corner_gid = mesh.corners_in_elem(elem_gid, n);
            //         int node_gid = mesh.nodes_in_elem(elem_gid, n);

            //         for (int d = 0; d < 3; d++) {
            //             corner_delta(corner_gid, d) = (float)(delta_lambda * grad(n, d));
            //             Kokkos::atomic_add(&node.coords(node_gid, d), corner_delta(corner_gid, d));
            //         }
            //     }

                
    
            // });

            // --- 4. Edge Length Constraint (Relaxation) ---
            FOR_ALL(node_gid, 0, mesh.num_nodes, {

                size_t node_a = node_gid;
                for(int node_lid = 0; node_lid < mesh.num_nodes_in_node(node_gid); node_lid++){


                    size_t node_b = mesh.nodes_in_node(node_gid, node_lid);
                    double L0 = min_distance_calc;

                    double dx[3];
                    double dist_sq = 0.0;
                    for(int d=0; d<3; d++) {
                        dx[d] = node.coords(node_a, d) - node.coords(node_b, d);
                        dist_sq += dx[d] * dx[d];
                    }
                    double dist = sqrt(dist_sq);
                    
                    if (dist < 1e-12) return; // Prevent division by zero

                    double C = dist - L0;
                    
                    // XPBD Denominator: sum(w*grad^2) + alpha_tilde
                    // w_a = 1.0, w_b = 1.0 -> denominator = 2.0 + alpha_tilde
                    double denominator = 2.0 + alpha_tilde_edge;

                    double old_lambda = edge_lambda(node_gid);
                    double delta_lambda = (-C - alpha_tilde_edge * old_lambda) / denominator;
                    edge_lambda(node_gid) += delta_lambda/mesh.num_nodes_in_node(node_gid);

                    // Unit gradient vector n
                    double n[3];
                    for(int d=0; d<3; d++) n[d] = dx[d] / dist;

                    // Update positions
                    for(int d=0; d<3; d++) {
                        double corr = delta_lambda * n[d];
                        Kokkos::atomic_add(&node.coords(node_a, d), corr);
                        Kokkos::atomic_add(&node.coords(node_b, d), -corr);
                    }
                }
            });

            // Element internal Diagonal Constraints
            FOR_ALL(elem_gid, 0, mesh.num_elems, {

                int num_sets = 4;

                int set_[num_sets*2];
                auto set = ViewCArrayKokkos<int>(set_, num_sets, 2);


                // Daigonals within the elelent
                set(0,0) = 0;
                set(0,1) = 7;
                set(1,0) = 1;
                set(1,1) = 6;
                set(2,0) = 2;
                set(2,1) = 5;
                set(3,0) = 3;
                set(3,1) = 4;

                double target = min_distance_calc*min_distance_calc + min_distance_calc*min_distance_calc + min_distance_calc*min_distance_calc;
                target = sqrt(target);

                for(int i = 0; i < 4; i++){

                    size_t node_a = mesh.nodes_in_elem(elem_gid, set(i,0));
                    size_t node_b = mesh.nodes_in_elem(elem_gid, set(i,1));

                    
                    double L0 = target;

                    double dx[3];
                    double dist_sq = 0.0;
                    for(int d=0; d<3; d++) {
                        dx[d] = node.coords(node_a, d) - node.coords(node_b, d);
                        dist_sq += dx[d] * dx[d];
                    }
                    double dist = sqrt(dist_sq);
                    
                    if (dist < 1e-12) return; // Prevent division by zero

                    double C = dist - L0;
                    
                    // XPBD Denominator: sum(w*grad^2) + alpha_tilde
                    // w_a = 1.0, w_b = 1.0 -> denominator = 2.0 + alpha_tilde
                    double denominator = 2.0 + alpha_tilde_edge;

                    double old_lambda = element_lambda_diag(elem_gid);
                    double delta_lambda = (-C - alpha_tilde_edge * old_lambda) / denominator;
                    element_lambda_diag(elem_gid) += delta_lambda / 4;

                    // Unit gradient vector n
                    double n[3];
                    for(int d=0; d<3; d++) n[d] = dx[d] / dist;

                    // Update positions
                    for(int d=0; d<3; d++) {
                        double corr = delta_lambda * n[d];
                        Kokkos::atomic_add(&node.coords(node_a, d), corr);
                        Kokkos::atomic_add(&node.coords(node_b, d), -corr);
                    }
                }
            
            });

            // Element face Diagonal Constraints
            FOR_ALL(elem_gid, 0, mesh.num_elems, {

                int num_sets = 8;

                int set_[num_sets*2];
                auto set = ViewCArrayKokkos<int>(set_, num_sets, 2);

                // Diagonals on the element face
                // -i face face
                set(0,0) = 0;
                set(0,1) = 6;
                set(1,0) = 3;
                set(1,1) = 4;
                // +i face face
                set(2,0) = 1;
                set(2,1) = 7;
                set(3,0) = 5;
                set(3,1) = 3;
                // -j face face
                set(4,0) = 0;
                set(4,1) = 5;
                set(5,0) = 1;
                set(5,1) = 4;
                // +j face face
                set(6,0) = 2;
                set(6,1) = 7;
                set(7,0) = 3;
                set(7,1) = 6;


                double target = min_distance_calc*min_distance_calc + min_distance_calc*min_distance_calc;
                target = sqrt(target);

                for(int i = 0; i < num_sets; i++){

                    size_t node_a = mesh.nodes_in_elem(elem_gid, set(i,0));
                    size_t node_b = mesh.nodes_in_elem(elem_gid, set(i,1));

                    
                    double L0 = target;

                    double dx[3];
                    double dist_sq = 0.0;
                    for(int d=0; d<3; d++) {
                        dx[d] = node.coords(node_a, d) - node.coords(node_b, d);
                        dist_sq += dx[d] * dx[d];
                    }
                    double dist = sqrt(dist_sq);
                    
                    if (dist < 1e-12) return; // Prevent division by zero

                    double C = dist - L0;
                    
                    // XPBD Denominator: sum(w*grad^2) + alpha_tilde
                    // w_a = 1.0, w_b = 1.0 -> denominator = 2.0 + alpha_tilde
                    double denominator = 2.0 + alpha_tilde_edge;

                    double old_lambda = element_lambda_diag(elem_gid);
                    double delta_lambda = (-C - alpha_tilde_edge * old_lambda) / denominator;
                    element_lambda_diag(elem_gid) += delta_lambda / 4;

                    // Unit gradient vector n
                    double n[3];
                    for(int d=0; d<3; d++) n[d] = dx[d] / dist;

                    // Update positions
                    for(int d=0; d<3; d++) {
                        double corr = delta_lambda * n[d];
                        Kokkos::atomic_add(&node.coords(node_a, d), corr);
                        Kokkos::atomic_add(&node.coords(node_b, d), -corr);
                    }
                }
            
            });


            


            // --- 3. Boundary Snap Constraint ---
            FOR_ALL(node_lid, 0, mesh.num_bdy_nodes, {
                size_t node_gid = mesh.bdy_nodes(node_lid);
                
                float p[3];
                for(int i=0; i<3; i++) p[i] = static_cast<float>(node.coords(node_gid, i));
                
                float closest[3];
                float dist_sq;
                // int facet_idx = stl_data.query_nearest_facet(p, dist_sq, closest);
                
                // float dist = sqrtf(dist_sq);

                // Calculate the distance using the previously saved nearest facet point
                dist_sq = 0.0f;
                for(int d = 0; d < 3; d++) {
                    dist_sq += (node.coords(node_gid, d) - nearest_facet_point(node_lid, d))*(node.coords(node_gid, d) - nearest_facet_point(node_lid, d));
                }

                double dist = sqrt(dist_sq);


                node.scalar_field(node_gid) = sqrt(dist_sq);

                // node.closest_facet_point(node_gid, 0) = closest[0] - node.coords(node_gid, 0);
                // node.closest_facet_point(node_gid, 1) = closest[1] - node.coords(node_gid, 1);
                // node.closest_facet_point(node_gid, 2) = closest[2] - node.coords(node_gid, 2);


                node.closest_facet_point(node_gid, 0) = nearest_facet_point(node_lid, 0) - node.coords(node_gid, 0);
                node.closest_facet_point(node_gid, 1) = nearest_facet_point(node_lid, 1) - node.coords(node_gid, 1);
                node.closest_facet_point(node_gid, 2) = nearest_facet_point(node_lid, 2) - node.coords(node_gid, 2);

                closest[0] = nearest_facet_point(node_lid, 0);
                closest[1] = nearest_facet_point(node_lid, 1);
                closest[2] = nearest_facet_point(node_lid, 2);

                if (dist < 1e-9) return; 
            
                float grad[3];
                for(int i=0; i<3; i++) grad[i] = (p[i] - closest[i]) / dist;
            
                // XPBD update for boundary
                double old_lambda = boundary_lambda(node_lid);
                double delta_lambda = (-dist - alpha_tilde_bdy * old_lambda) / (1.0 + alpha_tilde_bdy);
                boundary_lambda(node_lid) += delta_lambda;
                
                for(int d=0; d<3; d++) {
                    node.coords(node_gid, d) += delta_lambda * grad[d];
                }
            });

            // Compute the volume of the element
            FOR_ALL(elem_gid, 0, mesh.num_elems, {

                const size_t num_nodes = 8;

                double x_array[8];
                double y_array[8];
                double z_array[8];

                // x, y, z coordinates of elem vertices
                auto x = ViewCArrayKokkos<double>(x_array, num_nodes);
                auto y = ViewCArrayKokkos<double>(y_array, num_nodes);
                auto z = ViewCArrayKokkos<double>(z_array, num_nodes);

                // get the coordinates of the nodes(rk,elem,node) in this element
                for (int node_lid = 0; node_lid < num_nodes; node_lid++) {
                    size_t node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
                    x(node_lid) = node_coords(node_gid, 0);
                    y(node_lid) = node_coords(node_gid, 1);
                    z(node_lid) = node_coords(node_gid, 2);
                }     // end for

                double twelth = 1. / 12.;

                // element volume
                double volume =
                    (x(1) * (y(2) * (-z(0) + z(3)) + y(4) * (z(0) - z(5)) + y(0) * (z(3) + z(2) - z(4) - z(5)) + y(7) * (-z(3) + z(5)) + y(5) * (z(0) - z(3) + z(4) - z(7)) + y(3) * (-z(0) - z(2) + z(5) + z(7))) +
                    x(6) * (y(0) * (-z(2) + z(4)) + y(7) * (z(3) + z(2) - z(4) - z(5)) + y(3) * (z(2) - z(7)) + y(2) * (z(0) - z(3) + z(4) - z(7)) + y(5) * (-z(4) + z(7)) + y(4) * (-z(0) - z(2) + z(5) + z(7))) +
                    x(2) * (y(1) * (z(0) - z(3)) + y(6) * (-z(0) + z(3) - z(4) + z(7)) + y(7) * (z(3) - z(6)) + y(3) * (z(0) + z(1) - z(7) - z(6)) + y(4) * (-z(0) + z(6)) + y(0) * (-z(1) - z(3) + z(4) + z(6))) +
                    x(5) * (y(0) * (z(1) - z(4)) + y(6) * (z(4) - z(7)) + y(3) * (-z(1) + z(7)) + y(1) * (-z(0) + z(3) - z(4) + z(7)) + y(4) * (z(0) + z(1) - z(7) - z(6)) + y(7) * (-z(1) - z(3) + z(4) + z(6))) +
                    x(7) * (y(1) * (z(3) - z(5)) + y(6) * (-z(3) - z(2) + z(4) + z(5)) + y(5) * (z(1) + z(3) - z(4) - z(6)) + y(4) * (z(5) - z(6)) + y(2) * (-z(3) + z(6)) + y(3) * (-z(1) + z(2) - z(5) + z(6))) +
                    x(0) * (y(3) * (z(1) - z(2)) + y(6) * (z(2) - z(4)) + y(5) * (-z(1) + z(4)) + y(1) * (-z(3) - z(2) + z(4) + z(5)) + y(2) * (z(1) + z(3) - z(4) - z(6)) + y(4) * (-z(1) + z(2) - z(5) + z(6))) +
                    x(3) * (y(0) * (-z(1) + z(2)) + y(5) * (z(1) - z(7)) + y(1) * (z(0) + z(2) - z(5) - z(7)) + y(6) * (-z(2) + z(7)) + y(7) * (z(1) - z(2) + z(5) - z(6)) + y(2) * (-z(0) - z(1) + z(7) + z(6))) +
                    x(4) *
                    (y(1) * (-z(0) + z(5)) + y(6) * (z(0) + z(2) - z(5) - z(7)) + y(2) * (z(0) - z(6)) + y(0) * (z(1) - z(2) + z(5) - z(6)) + y(7) * (-z(5) + z(6)) + y(5) * (-z(0) - z(1) + z(7) + z(6)))) *
                    twelth;

                gauss_point.fields(elem_gid) = static_cast<float>(volume);
            });
            

        }

        std::cout<<"Writing VTU file for step "<<step<<std::endl;
        int rank = 0;
        double time_value = static_cast<double>(step) * dt;
        write_vtu(mesh, node, gauss_point, rank, MPI_COMM_WORLD, step, time_value);

    }
  

    

    } // end MATAR scope
    MATAR_FINALIZE();

    return 0;
}