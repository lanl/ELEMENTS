/**********************************************************************************************
© 2020. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and
to permit others to do so.
This program is open source under the BSD-3 License.
Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:
1.  Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer.
2.  Redistributions in binary form must reproduce the above copyright notice, this list of
conditions and the following disclaimer in the documentation and/or other materials
provided with the distribution.
3.  Neither the name of the copyright holder nor the names of its contributors may be used
to endorse or promote products derived from this software without specific prior
written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**********************************************************************************************/
#include "matar.h"
#include <cassert>
#include <set>
#include <cstdio>

// This pulls in kokkos, matar, mesh, hash, ref_elem stuff, and PT-Scotch
#include "ELEMENTS.h"

//#undef NDEBUG     // Ensures NDEBUG is turned off

using namespace mtr;
using namespace swage;

int main(int argc, char* argv[]) {

    Kokkos::initialize(argc, argv);
    {

        // -----------------------------------------------------------------------
        // The idea: simulate two hexahedral elements sharing a face.
        //
        //   Element 0 occupies [0,1] x [0,1] x [0,1]
        //   Element 1 occupies [1,2] x [0,1] x [0,1]
        //
        //   Each element has 8 corner points.  The two elements share 4 points
        //   on the face x=1, so 16 corner points map to 12 unique nodes:
        //
        //       Element 0 corners (points 0-7):
        //         0: (0,0,0)   1: (1,0,0)   2: (1,1,0)   3: (0,1,0)
        //         4: (0,0,1)   5: (1,0,1)   6: (1,1,1)   7: (0,1,1)
        //
        //       Element 1 corners (points 8-15):
        //         8: (1,0,0)   9: (2,0,0)  10: (2,1,0)  11: (1,1,0)
        //        12: (1,0,1)  13: (2,0,1)  14: (2,1,1)  15: (1,1,1)
        //
        //   Coincident pairs (same (x,y,z), should share a node gid):
        //     point  1 == point  8  at (1,0,0)
        //     point  2 == point 11  at (1,1,0)
        //     point  5 == point 12  at (1,0,1)
        //     point  6 == point 15  at (1,1,1)
        // -----------------------------------------------------------------------

        const size_t num_points = 16;
        CArrayKokkos <double> corner_pts(num_points, 3, "corner_pts");

        RUN({
            // element 0 — unit cube [0,1]^3
            corner_pts(0,  0) = 0.0;  corner_pts(0,  1) = 0.0;  corner_pts(0,  2) = 0.0;
            corner_pts(1,  0) = 1.0;  corner_pts(1,  1) = 0.0;  corner_pts(1,  2) = 0.0;
            corner_pts(2,  0) = 1.0;  corner_pts(2,  1) = 1.0;  corner_pts(2,  2) = 0.0;
            corner_pts(3,  0) = 0.0;  corner_pts(3,  1) = 1.0;  corner_pts(3,  2) = 0.0;
            corner_pts(4,  0) = 0.0;  corner_pts(4,  1) = 0.0;  corner_pts(4,  2) = 1.0;
            corner_pts(5,  0) = 1.0;  corner_pts(5,  1) = 0.0;  corner_pts(5,  2) = 1.0;
            corner_pts(6,  0) = 1.0;  corner_pts(6,  1) = 1.0;  corner_pts(6,  2) = 1.0;
            corner_pts(7,  0) = 0.0;  corner_pts(7,  1) = 1.0;  corner_pts(7,  2) = 1.0;

            // element 1 — cube [1,2] x [0,1] x [0,1]
            // points 8,11,12,15 are coincident with points 1,2,5,6 from element 0
            corner_pts(8,  0) = 1.0;  corner_pts(8,  1) = 0.0;  corner_pts(8,  2) = 0.0;
            corner_pts(9,  0) = 2.0;  corner_pts(9,  1) = 0.0;  corner_pts(9,  2) = 0.0;
            corner_pts(10, 0) = 2.0;  corner_pts(10, 1) = 1.0;  corner_pts(10, 2) = 0.0;
            corner_pts(11, 0) = 1.0;  corner_pts(11, 1) = 1.0;  corner_pts(11, 2) = 0.0;
            corner_pts(12, 0) = 1.0;  corner_pts(12, 1) = 0.0;  corner_pts(12, 2) = 1.0;
            corner_pts(13, 0) = 2.0;  corner_pts(13, 1) = 0.0;  corner_pts(13, 2) = 1.0;
            corner_pts(14, 0) = 2.0;  corner_pts(14, 1) = 1.0;  corner_pts(14, 2) = 1.0;
            corner_pts(15, 0) = 1.0;  corner_pts(15, 1) = 1.0;  corner_pts(15, 2) = 1.0;
        });

        
        SpatialConnectivity_t sc;


        // build the bin mesh with a small padding around the domain [0,2]^3
        sc.build_bin_mesh(-0.1, -0.1, -0.1,
                           2.1,  2.1,  2.1,
                           10, 10, 10);

        CArrayKokkos <size_t> node_in_corner_point;
        size_t num_nodes = 0;
        sc.build_multi_node_connectivity(corner_pts, node_in_corner_point, num_nodes);


        printf("\n--- build_multi_node_connectivity test ---\n");
        printf("num_points = %zu,  num_nodes = %zu  (expected 12)\n\n",
               num_points, num_nodes);

        // --- check 1: correct number of unique nodes ---
        assert(num_nodes == 12 && "wrong number of unique nodes");

        // --- check 2: all node gids are in range [0, num_nodes) ---
        RUN({
            for (size_t pt = 0; pt < num_points; pt++) {
                assert(node_in_corner_point(pt) < num_nodes &&
                    "node gid out of range");
            }
        });

        // --- check 3: coincident point pairs share the same node gid ---
        // pairs: (1,8), (2,11), (5,12), (6,15)
        const size_t coincident_pairs[4][2] = {{1,8}, {2,11}, {5,12}, {6,15}};
        for (auto &pair : coincident_pairs) {
            size_t pt_a = pair[0];
            size_t pt_b = pair[1];
            RUN({
                assert(node_in_corner_point(pt_a) ==
                    node_in_corner_point(pt_b) &&
                    "coincident points have different node gids");
            });
        }

        // --- check 4: non-coincident points have different node gids ---
        // collect all 16 node gid assignments; exactly 12 must be unique
        RUN({
            size_t num_unique_nodes = 0;
            size_t unique_nodes[16];

            for (size_t pt = 0; pt < num_points; pt++){
                bool found = false;

                // loop over all nodes to see if the node exists
                for(size_t node_lid = 0; node_lid < num_unique_nodes; node_lid++){
                    if (node_in_corner_point(pt)==unique_nodes[node_lid]) found = true;
                } // end for

                if(found==false){
                    unique_nodes[num_unique_nodes]=node_in_corner_point(pt);
                    num_unique_nodes++;
                }
            } // end for
        
            assert(num_unique_nodes == 12 &&
                "wrong number of distinct node gids across all points");
        });

        // --- print the mapping for visual inspection ---
        printf("point -> node mapping:\n");
        RUN({
            for (size_t pt = 0; pt < num_points; pt++) {
                printf("  point %2zu  (%.1f, %.1f, %.1f)  ->  node %zu\n",
                    pt,
                    corner_pts(pt, 0),
                    corner_pts(pt, 1),
                    corner_pts(pt, 2),
                    node_in_corner_point(pt));
            }
        });

        printf("\nAll build_multi_node_connectivity checks passed.\n");

    } // end Kokkos scope
    Kokkos::finalize();
    return 0;
}