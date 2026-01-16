//Scratch

       // --- 4. Strain Energy Constraint (St. Venant-Kirchhoff) ---
       for(int elem_gid=0; elem_gid < mesh.num_elems; elem_gid++) {

        // Compute the deformation graduent using the current positions and the B matrix
        double pos_array[8*3];
        
        // x, y, z coordinates of elem vertices
        auto pos = ViewCArrayKokkos<double>(pos_array, 8,3);

        // get the coordinates of the nodes(rk,elem,node) in this element
        for (int node_lid = 0; node_lid < 8; node_lid++) {

            pos(node_lid, 0) = node_coords(mesh.nodes_in_elem(elem_gid, node_lid), 0);
            pos(node_lid, 1) = node_coords(mesh.nodes_in_elem(elem_gid, node_lid), 1);
            pos(node_lid, 2) = node_coords(mesh.nodes_in_elem(elem_gid, node_lid), 2);
        } // end for

        for(int i=0; i<3; i++) {
            for(int j=0; j<3; j++) {
                F_matrix(elem_gid, i, j) = 0.0;
            }
        }

        for(int i=0; i<3; i++) {
            for(int j=0; j<3; j++) {
                for(int k=0; k<8; k++) {
                    F_matrix(elem_gid, i, j) += pos(k, i) * B_matrix(elem_gid, k, j);
                }
            }
        }

        // Compute the Green-Lagrange strain tensor
        double strain_[9];
        auto strain = ViewCArrayKokkos<double>(strain_, 3,3);
        for(int i=0; i<3; i++) {
            for(int j=0; j<3; j++) {
                strain(i, j) = 0.5 * (F_matrix(elem_gid, j, i) * F_matrix(elem_gid, i, j));
                if(i==j) {
                    strain(i, j) -= 0.5;
                }
            }
        }

        // Compute the first Piola-Kirchhoff tensor
        // P(F) = F*(2*mu*E + lambda*tr(E)*I)
        double P_[9];
        auto P = ViewCArrayKokkos<double>(P_, 3,3);
        for(int i=0; i<3; i++) {
            for(int j=0; j<3; j++) {
                P(i, j) = 0.0;
            }
        }

        double Youngs_modulus = 1.0;
        double poissons_ratio = 0.45;

        double lambda = Youngs_modulus * poissons_ratio / ((1.0 + poissons_ratio) * (1.0 - 2.0 * poissons_ratio));
        double mu = Youngs_modulus / (2.0 * (1.0 + poissons_ratio));  

        for(int i=0; i<3; i++) {
            for(int j=0; j<3; j++) {
                P(i, j) = F_matrix(elem_gid, i, j) * (2.0 * Youngs_modulus * strain(i, j) + lambda * (strain(0, 0) + strain(1, 1) + strain(2, 2)) * (i==j ? 1.0 : 0.0));
            }
        }

        // Compute the strain energy density C
        double tr_E = strain(0, 0) + strain(1, 1) + strain(2, 2);
        
        double E_sq_[9];
        auto E_sq = ViewCArrayKokkos<double>(E_sq_, 3,3);

        for(int i=0; i<3; i++) {
            for(int j=0; j<3; j++) {
                for(int k=0; k<3; k++) {
                    E_sq(i, j) += strain(i, k) * strain(k, j);
                }
            }
        }

        double tr_E_sq = E_sq(0, 0) + E_sq(1, 1) + E_sq(2, 2);

        double C = mu * tr_E_sq + 0.5 * lambda * tr_E * tr_E;

        // Compute the gradient of the strain energy density
        double grad_C_[9];
        auto grad_C = ViewCArrayKokkos<double>(grad_C_, 3,3);
        for(int i=0; i<3; i++) {
            for(int j=0; j<3; j++) {
                grad_C(i, j) = 0.0;
            }
        }
        
       
    }




    DCArrayKokkos<float> B_matrix(mesh.num_elems, 8, 3, "B_matrix");
    B_matrix.set_values(0.0);
    DCArrayKokkos<float> F_matrix(mesh.num_elems, 3, 3, "F_matrix");
    F_matrix.set_values(0.0);
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
            x(node_lid) = node_coords(mesh.nodes_in_elem(elem_gid, node_lid), 0);
            y(node_lid) = node_coords(mesh.nodes_in_elem(elem_gid, node_lid), 1);
            z(node_lid) = node_coords(mesh.nodes_in_elem(elem_gid, node_lid), 2);
        }     // end for


        double twelth = 1. / 12.;

        B_matrix(elem_gid, 0, 0) = (+y(1) * (-z(3) - z(2) + z(4) + z(5) )
            + y(3) * (+z(1) - z(2) )
            + y(2) * (+z(1) + z(3) - z(4) - z(6) )
            + y(4) * (-z(1) + z(2) - z(5) + z(6) )
            + y(5) * (-z(1) + z(4) )
            + y(6) * (+z(2) - z(4) ) ) * twelth;

        B_matrix(elem_gid, 1, 0) = (+y(0) * (+z(3) + z(2) - z(4) - z(5) )
            + y(3) * (-z(0) - z(2) + z(5) + z(7) )
            + y(2) * (-z(0) + z(3) )
            + y(4) * (+z(0) - z(5) )
            + y(5) * (+z(0) - z(3) + z(4) - z(7) )
            + y(7) * (-z(3) + z(5) ) ) * twelth;

        B_matrix(elem_gid, 2, 0) = (+y(0) * (-z(1) - z(3) + z(4) + z(6) )
            + y(1) * (+z(0) - z(3) )
            + y(3) * (+z(0) + z(1) - z(7) - z(6) )
            + y(4) * (-z(0) + z(6) )
            + y(7) * (+z(3) - z(6) )
            + y(6) * (-z(0) + z(3) - z(4) + z(7) ) ) * twelth;

        B_matrix(elem_gid, 3, 0) = (+y(0) * (-z(1) + z(2) )
            + y(1) * (+z(0) + z(2) - z(5) - z(7) )
            + y(2) * (-z(0) - z(1) + z(7) + z(6) )
            + y(5) * (+z(1) - z(7) )
            + y(7) * (+z(1) - z(2) + z(5) - z(6) )
            + y(6) * (-z(2) + z(7) ) ) * twelth;

        B_matrix(elem_gid, 4, 0) = (+y(0) * (+z(1) - z(2) + z(5) - z(6) )
            + y(1) * (-z(0) + z(5) )
            + y(2) * (+z(0) - z(6) )
            + y(5) * (-z(0) - z(1) + z(7) + z(6) )
            + y(7) * (-z(5) + z(6) )
            + y(6) * (+z(0) + z(2) - z(5) - z(7) ) ) * twelth;

        B_matrix(elem_gid, 5, 0) = (+y(0) * (+z(1) - z(4) )
            + y(1) * (-z(0) + z(3) - z(4) + z(7) )
            + y(3) * (-z(1) + z(7) )
            + y(4) * (+z(0) + z(1) - z(7) - z(6) )
            + y(7) * (-z(1) - z(3) + z(4) + z(6) )
            + y(6) * (+z(4) - z(7) ) ) * twelth;

        B_matrix(elem_gid, 6, 0) = (+y(0) * (-z(2) + z(4) )
            + y(3) * (+z(2) - z(7) )
            + y(2) * (+z(0) - z(3) + z(4) - z(7) )
            + y(4) * (-z(0) - z(2) + z(5) + z(7) )
            + y(5) * (-z(4) + z(7) )
            + y(7) * (+z(3) + z(2) - z(4) - z(5) ) ) * twelth;

        B_matrix(elem_gid, 7, 0) = (+y(1) * (+z(3) - z(5) )
            + y(3) * (-z(1) + z(2) - z(5) + z(6) )
            + y(2) * (-z(3) + z(6) )
            + y(4) * (+z(5) - z(6) )
            + y(5) * (+z(1) + z(3) - z(4) - z(6) )
            + y(6) * (-z(3) - z(2) + z(4) + z(5) ) ) * twelth;

        B_matrix(elem_gid, 0, 1) = (+z(1) * (-x(3) - x(2) + x(4) + x(5) )
            + z(3) * (+x(1) - x(2) )
            + z(2) * (+x(1) + x(3) - x(4) - x(6) )
            + z(4) * (-x(1) + x(2) - x(5) + x(6) )
            + z(5) * (-x(1) + x(4) )
            + z(6) * (+x(2) - x(4) ) ) * twelth;

        B_matrix(elem_gid, 1, 1) = (+z(0) * (+x(3) + x(2) - x(4) - x(5) )
            + z(3) * (-x(0) - x(2) + x(5) + x(7) )
            + z(2) * (-x(0) + x(3) )
            + z(4) * (+x(0) - x(5) )
            + z(5) * (+x(0) - x(3) + x(4) - x(7) )
            + z(7) * (-x(3) + x(5) ) ) * twelth;

        B_matrix(elem_gid, 2, 1) = (+z(0) * (-x(1) - x(3) + x(4) + x(6) )
            + z(1) * (+x(0) - x(3) )
            + z(3) * (+x(0) + x(1) - x(7) - x(6) )
            + z(4) * (-x(0) + x(6) )
            + z(7) * (+x(3) - x(6) )
            + z(6) * (-x(0) + x(3) - x(4) + x(7) ) ) * twelth;

        B_matrix(elem_gid, 3, 1) = (+z(0) * (-x(1) + x(2) )
            + z(1) * (+x(0) + x(2) - x(5) - x(7) )
            + z(2) * (-x(0) - x(1) + x(7) + x(6) )
            + z(5) * (+x(1) - x(7) )
            + z(7) * (+x(1) - x(2) + x(5) - x(6) )
            + z(6) * (-x(2) + x(7) ) ) * twelth;

        B_matrix(elem_gid, 4, 1) = (+z(0) * (+x(1) - x(2) + x(5) - x(6) )
            + z(1) * (-x(0) + x(5) )
            + z(2) * (+x(0) - x(6) )
            + z(5) * (-x(0) - x(1) + x(7) + x(6) )
            + z(7) * (-x(5) + x(6) )
            + z(6) * (+x(0) + x(2) - x(5) - x(7) ) ) * twelth;

        B_matrix(elem_gid, 5, 1) = (+z(0) * (+x(1) - x(4) )
            + z(1) * (-x(0) + x(3) - x(4) + x(7) )
            + z(3) * (-x(1) + x(7) )
            + z(4) * (+x(0) + x(1) - x(7) - x(6) )
            + z(7) * (-x(1) - x(3) + x(4) + x(6) )
            + z(6) * (+x(4) - x(7) ) ) * twelth;

        B_matrix(elem_gid, 6, 1) = (+z(0) * (-x(2) + x(4) )
            + z(3) * (+x(2) - x(7) )
            + z(2) * (+x(0) - x(3) + x(4) - x(7) )
            + z(4) * (-x(0) - x(2) + x(5) + x(7) )
            + z(5) * (-x(4) + x(7) )
            + z(7) * (+x(3) + x(2) - x(4) - x(5) ) ) * twelth;

        B_matrix(elem_gid, 7, 1) = (+z(1) * (+x(3) - x(5) )
            + z(3) * (-x(1) + x(2) - x(5) + x(6) )
            + z(2) * (-x(3) + x(6) )
            + z(4) * (+x(5) - x(6) )
            + z(5) * (+x(1) + x(3) - x(4) - x(6) )
            + z(6) * (-x(3) - x(2) + x(4) + x(5) ) ) * twelth;

        B_matrix(elem_gid, 0, 2) = (+x(1) * (-y(3) - y(2) + y(4) + y(5) )
            + x(3) * (+y(1) - y(2) )
            + x(2) * (+y(1) + y(3) - y(4) - y(6) )
            + x(4) * (-y(1) + y(2) - y(5) + y(6) )
            + x(5) * (-y(1) + y(4) )
            + x(6) * (+y(2) - y(4) ) ) * twelth;

        B_matrix(elem_gid, 1, 2) = (+x(0) * (+y(3) + y(2) - y(4) - y(5) )
            + x(3) * (-y(0) - y(2) + y(5) + y(7) )
            + x(2) * (-y(0) + y(3) )
            + x(4) * (+y(0) - y(5) )
            + x(5) * (+y(0) - y(3) + y(4) - y(7) )
            + x(7) * (-y(3) + y(5) ) ) * twelth;

        B_matrix(elem_gid, 2, 2) = (+x(0) * (-y(1) - y(3) + y(4) + y(6) )
            + x(1) * (+y(0) - y(3) )
            + x(3) * (+y(0) + y(1) - y(7) - y(6) )
            + x(4) * (-y(0) + y(6) )
            + x(7) * (+y(3) - y(6) )
            + x(6) * (-y(0) + y(3) - y(4) + y(7) ) ) * twelth;

        B_matrix(elem_gid, 3, 2) = (+x(0) * (-y(1) + y(2) )
            + x(1) * (+y(0) + y(2) - y(5) - y(7) )
            + x(2) * (-y(0) - y(1) + y(7) + y(6) )
            + x(5) * (+y(1) - y(7) )
            + x(7) * (+y(1) - y(2) + y(5) - y(6) )
            + x(6) * (-y(2) + y(7) ) ) * twelth;

        B_matrix(elem_gid, 4, 2) = (+x(0) * (+y(1) - y(2) + y(5) - y(6) )
            + x(1) * (-y(0) + y(5) )
            + x(2) * (+y(0) - y(6) )
            + x(5) * (-y(0) - y(1) + y(7) + y(6) )
            + x(7) * (-y(5) + y(6) )
            + x(6) * (+y(0) + y(2) - y(5) - y(7) ) ) * twelth;

        B_matrix(elem_gid, 5, 2) = (+x(0) * (+y(1) - y(4) )
            + x(1) * (-y(0) + y(3) - y(4) + y(7) )
            + x(3) * (-y(1) + y(7) )
            + x(4) * (+y(0) + y(1) - y(7) - y(6) )
            + x(7) * (-y(1) - y(3) + y(4) + y(6) )
            + x(6) * (+y(4) - y(7) ) ) * twelth;

        B_matrix(elem_gid, 6, 2) = (+x(0) * (-y(2) + y(4) )
            + x(3) * (+y(2) - y(7) )
            + x(2) * (+y(0) - y(3) + y(4) - y(7) )
            + x(4) * (-y(0) - y(2) + y(5) + y(7) )
            + x(5) * (-y(4) + y(7) )
            + x(7) * (+y(3) + y(2) - y(4) - y(5) ) ) * twelth;

        B_matrix(elem_gid, 7, 2) = (+x(1) * (+y(3) - y(5) )
            + x(3) * (-y(1) + y(2) - y(5) + y(6) )
            + x(2) * (-y(3) + y(6) )
            + x(4) * (+y(5) - y(6) )
            + x(5) * (+y(1) + y(3) - y(4) - y(6) )
            + x(6) * (-y(3) - y(2) + y(4) + y(5) ) ) * twelth;

    });