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
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cassert>

// This pulls in kokkos, matar, mesh, hash, ref_elem stuff, and PT-Scotch
#include "ELEMENTS.h"

using namespace mtr;
using namespace swage; // unstructured mesh and hash

int main(int argc, char** argv) {

MATAR_INITIALIZE(argc, argv);
{ // MATAR scope

    printf("\n--- ben key test ---\n");

    SpatialConnectivity_t sc;
    
    sc.build_bin_mesh(0.0, 0.0, 0.0,   // xmin, ymin, zmin
                      1.0, 1.0, 1.0,   // xmax, ymax, zmax
                      10, 10, 10);     // num bins per direction

    // bin spacing should be 0.1
    assert(fabs(sc.bin_dx - 0.1) < 1.e-14);

    // a point at (0.15, 0.25, 0.35) should land in bin (1, 2, 3)
    bin_keys_t bk = sc.get_bin_keys(0.15, 0.25, 0.35);
    assert(bk.i == 1 && bk.j == 2 && bk.k == 3 &&
                    "bin keys i,j,k are wrong");

    // get_bin_gid should give the same result as get_id_of_ijk
    size_t gid = sc.get_bin_gid(0.15, 0.25, 0.35);
    assert(gid == get_id_of_ijk(1, 2, 3, 10, 10) &&
                    "gid of i,j,k is wrong");

    // points at domain boundaries should clamp, not go out of bounds
    bk = sc.get_bin_keys(-1.0, -1.0, -1.0);  // below domain
    assert(bk.i == 0 && bk.j == 0 && bk.k == 0 &&
                    "bin keys i,j,k below domain are wrong");

    bk = sc.get_bin_keys(2.0, 2.0, 2.0);     // above domain
    assert(bk.i == 9 && bk.j == 9 && bk.k == 9 &&
                    "bin keys i,j,k above domain are wrong");

    printf("\nAll bin key checks passed.\n");

}
MATAR_FINALIZE();

} // end main

