/**********************************************************************************************
� 2020. Triad National Security, LLC. All rights reserved.
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

#ifndef STL_UTILS_H
#define STL_UTILS_H

#include "matar.h"

#include <Kokkos_Sort.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>

using namespace mtr;

// AABB Node structure optimized for GPU/Parallel access
struct AABBNode {
    float min_ext[3];
    float max_ext[3];
    int left = -1;   
    int right = -1;  
    int parent = -1; 
    int atomic_count = 0; // <--- Add this for the refit logic
};


struct stl_data{
    
    size_t num_facets;
    DCArrayKokkos<float> normal;   // (N, 3)
    DCArrayKokkos<float> vertices; // (N, 3, 3)
    DCArrayKokkos<float> center;   // (N, 3)

    // Global bounding box
    float min_x = 1e30f, max_x = -1e30f;
    float min_y = 1e30f, max_y = -1e30f;
    float min_z = 1e30f, max_z = -1e30f;

    // AABB tree data
    DCArrayKokkos<AABBNode> tree_nodes;
    DCArrayKokkos<int> sorted_facet_indices;
    int root_idx = 0;

    // Comparator for Kokkos::sort — must be at struct scope, not defined
    // inside buildAABBTree(), because CUDA's stub generator cannot handle
    // locally-defined types in kernel instantiations (ZZ-mangled names).
    struct MortonComparator {
        Kokkos::View<uint32_t*> codes;
        KOKKOS_INLINE_FUNCTION
        bool operator()(const int& i, const int& j) const {
            return codes(i) < codes(j);
        }
    };


    void initialize(size_t num_facets_in){
        this->num_facets = num_facets_in;
        normal = DCArrayKokkos<float>(num_facets, 3);
        vertices = DCArrayKokkos<float>(num_facets, 3, 3);
        center = DCArrayKokkos<float>(num_facets, 3);
    }


    // --- Bit Manipulation Helpers ---

    // Portable CLZ: pure integer arithmetic, no architecture guards needed.
    // Modern compilers (nvcc, hipcc, gcc, clang) recognize this pattern and
    // emit the hardware instruction (__clz on CUDA/HIP, lzcnt on x86).
    // Behavior is undefined for x == 0, matching __builtin_clz semantics.
    KOKKOS_INLINE_FUNCTION
    int clz32(uint32_t x) const {
        int n = 0;
        if ((x & 0xFFFF0000u) == 0) { n += 16; x <<= 16; }
        if ((x & 0xFF000000u) == 0) { n +=  8; x <<=  8; }
        if ((x & 0xF0000000u) == 0) { n +=  4; x <<=  4; }
        if ((x & 0xC0000000u) == 0) { n +=  2; x <<=  2; }
        if ((x & 0x80000000u) == 0) { n +=  1;           }
        return n;
    }

    KOKKOS_INLINE_FUNCTION
    uint32_t expandBits(uint32_t v) const {
        v = (v * 0x00010001u) & 0xFF0000FFu;
        v = (v * 0x00000101u) & 0x0F00F00Fu;
        v = (v * 0x00000011u) & 0xC30C30C3u;
        v = (v * 0x00000005u) & 0x49249249u;
        return v;
    }

    KOKKOS_INLINE_FUNCTION 
    uint32_t morton3D(float x, float y, float z) const {
        return (expandBits((uint32_t)x) << 2) | (expandBits((uint32_t)y) << 1) | expandBits((uint32_t)z);
    }


    KOKKOS_INLINE_FUNCTION 
    int common_prefix(const Kokkos::View<const uint32_t*>& codes, int n, int i, int j) {
        if (j < 0 || j > n - 1) return -1;
        const uint32_t a = codes(i);
        const uint32_t b = codes(j);
        if (a != b) {
            return clz32(a ^ b);
        }
        // Tie-break identical Morton codes with the indices so ordering is total
        return 32 + clz32(static_cast<uint32_t>(i ^ j));
    }

    

    void computeGlobalBoundingBox() {


        // --- X Dimension ---
        float loc_min_x = 1e30f;
        FOR_REDUCE_MIN_CLASS(i, 0, num_facets, loc_min_x, {
            for (int v = 0; v < 3; ++v) {
                if (vertices(i, v, 0) < loc_min_x) loc_min_x = vertices(i, v, 0);
            }
        }, min_x);

        float loc_max_x = -1e30f;
        FOR_REDUCE_MAX_CLASS(i, 0, num_facets, loc_max_x, {
            for (int v = 0; v < 3; ++v) {
                if (vertices(i, v, 0) > loc_max_x) loc_max_x = vertices(i, v, 0);
            }
        }, max_x);

        // --- Y Dimension ---
        float loc_min_y = 1e30f;
        FOR_REDUCE_MIN_CLASS(i, 0, num_facets, loc_min_y, {
            for (int v = 0; v < 3; ++v) {
                if (vertices(i, v, 1) < loc_min_y) loc_min_y = vertices(i, v, 1);
            }
        }, min_y);

        float loc_max_y = -1e30f;
        FOR_REDUCE_MAX_CLASS(i, 0, num_facets, loc_max_y, {
            for (int v = 0; v < 3; ++v) {
                if (vertices(i, v, 1) > loc_max_y) loc_max_y = vertices(i, v, 1);
            }
        }, max_y);

        // --- Z Dimension ---
        float loc_min_z = 1e30f;
        FOR_REDUCE_MIN_CLASS(i, 0, num_facets, loc_min_z, {
            for (int v = 0; v < 3; ++v) {
                if (vertices(i, v, 2) < loc_min_z) loc_min_z = vertices(i, v, 2);
            }
        }, min_z);

        float loc_max_z = -1e30f;
        FOR_REDUCE_MAX_CLASS(i, 0, num_facets, loc_max_z, {
            for (int v = 0; v < 3; ++v) {
                if (vertices(i, v, 2) > loc_max_z) loc_max_z = vertices(i, v, 2);
            }
        }, max_z);

        std::cout<<"Global bounding box: "<<min_x<<", "<<max_x<<", "<<min_y<<", "<<max_y<<", "<<min_z<<", "<<max_z<<std::endl;
    }


    // --- The Main Build Function ---
    void buildAABBTree() {

        // 1. Calculate Global Bounding Box for Morton Scaling
        computeGlobalBoundingBox();

        // 2. Generate Morton Codes
        DCArrayKokkos<uint32_t> morton_codes(num_facets, "morton_codes");
        this->sorted_facet_indices = DCArrayKokkos<int>(num_facets, "sorted_facet_indices");

        std::cout<<"Generating Morton codes"<<std::endl;
        FOR_ALL_CLASS(i, 0, num_facets, {
        
            float x = (center(i,0) - min_x) / (max_x - min_x) * 1023.0f;
            float y = (center(i,1) - min_y) / (max_y - min_y) * 1023.0f;
            float z = (center(i,2) - min_z) / (max_z - min_z) * 1023.0f;
            morton_codes(i) = morton3D(x, y, z);
            sorted_facet_indices(i) = i; // write into member view
        });
        MATAR_FENCE();

        std::cout<<"Morton codes generated"<<std::endl;

        // 3. Sort Indices by Morton Code
        auto morton_view = morton_codes.get_kokkos_dual_view().view_device();
        auto index_view = sorted_facet_indices.get_kokkos_dual_view().view_device();

        std::cout<<"Sorting indices by Morton code"<<std::endl;

        // 3. Call Kokkos::sort using the struct-scope MortonComparator
        MortonComparator comp = { morton_view };
        Kokkos::sort(index_view, comp);
        sorted_facet_indices.update_host();

        // Build a sorted Morton-code view so the hierarchy uses contiguous Morton ordering
        DCArrayKokkos<uint32_t> sorted_codes(num_facets, "sorted_codes");
        auto sorted_codes_v = sorted_codes.get_kokkos_dual_view().view_device();
        FOR_ALL(i, 0, num_facets, {
            sorted_codes_v(i) = morton_view(index_view(i));
        });
        sorted_codes.update_host();

        
        // 4. Initialize and Resize the Tree (host-side deterministic build)
        int num_leaves = num_facets;
        int num_internal = num_leaves - 1;
        int total_nodes = num_leaves + num_internal;

        std::cout<<"Number of internal nodes: "<<num_internal<<std::endl;
        std::cout<<"Total nodes: "<<total_nodes<<std::endl;

        this->tree_nodes = DCArrayKokkos<AABBNode>(total_nodes);

        // Host mirrors for deterministic serial build (works in serial/device-default too)
        sorted_codes.update_host(); // ensure host codes up to date
        sorted_facet_indices.update_host();
        tree_nodes.update_host();   // make sure host mirror exists

        auto codes_h = sorted_codes.host;
        auto idx_h   = sorted_facet_indices.host;

        // Initialize nodes on host
        for (int i = 0; i < total_nodes; ++i) {
            tree_nodes.host(i).atomic_count = 0;
            tree_nodes.host(i).parent = -1;
            tree_nodes.host(i).left = -1;
            tree_nodes.host(i).right = -1;
            for (int d = 0; d < 3; ++d) {
                tree_nodes.host(i).min_ext[d] = 1e30f;
                tree_nodes.host(i).max_ext[d] = -1e30f;
            }
        }

        auto clz32 = [this](uint32_t x) { return this->clz32(x); };

        auto prefix_host = [&](int i, int j) {
            if (j < 0 || j > num_leaves - 1) return -1;
            uint32_t a = codes_h(i);
            uint32_t b = codes_h(j);
            if (a != b) return clz32(a ^ b);
            return 32 + clz32(static_cast<uint32_t>(i ^ j));
        };

        auto find_split = [&](int first, int last) {
            // last is inclusive
            uint32_t first_code = codes_h(first);
            uint32_t last_code  = codes_h(last);
            if (first_code == last_code) {
                return (first + last) >> 1;
            }
            int common = clz32(first_code ^ last_code);
            int split = first;
            int step = last - first;
            do {
                step = (step + 1) >> 1;
                int new_split = split + step;
                if (new_split < last) {
                    uint32_t split_code = codes_h(new_split);
                    int split_common = clz32(first_code ^ split_code);
                    if (split_common > common) {
                        split = new_split;
                    }
                }
            } while (step > 1);
            return split;
        };

        std::cout<<"Initializing tree (host)"<<std::endl;
        int next_internal = 0;
        std::function<int(int,int)> build = [&](int first, int last_exclusive) -> int {
            // returns node index
            int range = last_exclusive - first;
            if (range == 1) {
                int leaf_idx = num_internal + first;
                return leaf_idx;
            }
            int last = last_exclusive - 1;
            int split = find_split(first, last);

            int my = next_internal++;
            int left_child  = build(first, split + 1);
            int right_child = build(split + 1, last_exclusive);

            tree_nodes.host(my).left = left_child;
            tree_nodes.host(my).right = right_child;
            tree_nodes.host(left_child).parent = my;
            tree_nodes.host(right_child).parent = my;
            return my;
        };

        int root = build(0, num_leaves);
        root_idx = root;
        tree_nodes.host(root).parent = -1;
        vertices.update_host();
        // Leaf bounds on host
        for (int i = 0; i < num_facets; ++i) {
            int leaf_idx = num_internal + i;
            int facet_idx = idx_h(i);
            for (int d = 0; d < 3; ++d) {
                float v0 = vertices.host(facet_idx, 0, d);
                float v1 = vertices.host(facet_idx, 1, d);
                float v2 = vertices.host(facet_idx, 2, d);
                tree_nodes.host(leaf_idx).min_ext[d] = std::fmin(std::fmin(v0, v1), v2);
                tree_nodes.host(leaf_idx).max_ext[d] = std::fmax(std::fmax(v0, v1), v2);
            }
        }

        // Bottom-up refit on host (reverse pre-order numbering ensures children > parent)
        for (int i = next_internal - 1; i >= 0; --i) {
            int l = tree_nodes.host(i).left;
            int r = tree_nodes.host(i).right;
            for (int d = 0; d < 3; ++d) {
                tree_nodes.host(i).min_ext[d] = std::min(tree_nodes.host(l).min_ext[d], tree_nodes.host(r).min_ext[d]);
                tree_nodes.host(i).max_ext[d] = std::max(tree_nodes.host(l).max_ext[d], tree_nodes.host(r).max_ext[d]);
            }
        }

        // Sync to device for any downstream GPU use
        tree_nodes.update_device();
        
    }

    void verifyAABBTree() {
        // 1. Sync data to host
        tree_nodes.update_host();
        sorted_facet_indices.update_host();
        
        int num_internal = num_facets - 1;
        int total_nodes = num_facets + num_internal;
        
        int errors = 0;
        int parent_mismatch = 0;
        int containment_errors = 0;

        // 2. Check Root Bounding Box against Global Bounding Box
        // Root is always index 0
        bool root_bounds_ok = true;
        float root_min[3] = {tree_nodes.host(0).min_ext[0], tree_nodes.host(0).min_ext[1], tree_nodes.host(0).min_ext[2]};
        float root_max[3] = {tree_nodes.host(0).max_ext[0], tree_nodes.host(0).max_ext[1], tree_nodes.host(0).max_ext[2]};
        
        if (std::abs(root_min[0] - min_x) > 1e-4 || std::abs(root_max[0] - max_x) > 1e-4) root_bounds_ok = false;
        
        if (!root_bounds_ok) {
            std::cout << "  [!] Warning: Root AABB does not match Global BBox." << std::endl;
            std::cout << "      Root Min: " << root_min[0] << " Global Min: " << min_x << std::endl;
        }

        // 3. Structural and Geometric Invariants
        for (int i = 0; i < total_nodes; i++) {
            // Check Parent-Child Linkage
            int p = tree_nodes.host(i).parent;
            if (i == 0) {
                if (p != -1) { 
                    std::cout << "  Error: Root has a parent index: " << p << std::endl;
                    errors++;
                }
            } else {
                if (p == -1) {
                    std::cout << "  Error: Orphan node detected at index: " << i << std::endl;
                    errors++;
                } else {
                    // Reciprocity check: If p is my parent, I must be p's left or right child
                    if (tree_nodes.host(p).left != i && tree_nodes.host(p).right != i) {
                        parent_mismatch++;
                    }
                }
            }

            // Check Bounding Box Containment (Internal Nodes Only)
            if (i < num_internal) {
                int l = tree_nodes.host(i).left;
                int r = tree_nodes.host(i).right;

                if (l != -1 && r != -1) {
                    for (int d = 0; d < 3; d++) {
                        // Left child containment
                        if (tree_nodes.host(l).min_ext[d] < tree_nodes.host(i).min_ext[d] - 1e-6 || 
                            tree_nodes(l).max_ext[d] > tree_nodes(i).max_ext[d] + 1e-6) {
                            containment_errors++;
                        }
                        // Right child containment
                        if (tree_nodes.host(r).min_ext[d] < tree_nodes.host(i).min_ext[d] - 1e-6 || 
                            tree_nodes.host(r).max_ext[d] > tree_nodes.host(i).max_ext[d] + 1e-6) {
                            containment_errors++;
                        }
                    }
                }
            }
        }

        std::cout << "  Total Nodes:          " << total_nodes << std::endl;
        std::cout << "  Structure Errors:     " << errors << std::endl;
        std::cout << "  Parent/Child Mismatch: " << parent_mismatch << std::endl;
        std::cout << "  Containment Errors:   " << containment_errors << std::endl;

        if (errors + parent_mismatch + containment_errors == 0) {
            std::cout << "  SUCCESS: AABB Tree Invariants Verified." << std::endl;
        } else {
            std::cout << "  FAILURE: Tree is corrupted." << std::endl;
        }
        std::cout << "---------------------------" << std::endl;
    }

    KOKKOS_INLINE_FUNCTION
    float point_aabb_dist_sq(const float p[3], const AABBNode& node) const {
        float dist_sq = 0.0f;
        for (int i = 0; i < 3; i++) {
            if (p[i] < node.min_ext[i]) {
                float diff = node.min_ext[i] - p[i];
                dist_sq += diff * diff;
            }
            else if (p[i] > node.max_ext[i]) {
                float diff = p[i] - node.max_ext[i];
                dist_sq += diff * diff;
            }
        }
        return dist_sq;
    }


    KOKKOS_INLINE_FUNCTION
    int query_nearest_facet(const float p[3], 
                            float& out_dist_sq, 
                            float out_closest_pt[3]) const {
        
        // Local views for high-performance device access

        int num_internal = (int)num_facets - 1;
        int stack[64];
        int stack_ptr = 0;
        
        stack[stack_ptr++] = 0; // Start at root node

        float min_dist_sq = 1e30f;
        int best_facet = -1;

        while (stack_ptr > 0) {
            int curr = stack[--stack_ptr];

            // 1. Prune: If the point is already further from this AABB than our best triangle, skip
            if (point_aabb_dist_sq(p, tree_nodes(curr)) >= min_dist_sq) continue;

            if (curr >= num_internal) {
                // 2. LEAF NODE: Perform exact triangle distance calculation
                int facet_idx = sorted_facet_indices(curr - num_internal);
                float cp[3];
                
                // Call your refactored get_closest_point_on_triangle
                get_closest_point_on_triangle(p, facet_idx, vertices, cp);
                
                float d2 = 0.0f;
                for(int d=0; d<3; ++d) {
                    float diff = p[d] - cp[d];
                    d2 += diff * diff;
                }
                
                if (d2 < min_dist_sq) {
                    min_dist_sq = d2;
                    best_facet = facet_idx;
                    for(int d=0; d<3; ++d) out_closest_pt[d] = cp[d];
                }
            } 
            else {
                // 3. INTERNAL NODE: Decide which child to visit first
                int l = tree_nodes(curr).left;
                int r = tree_nodes(curr).right;
                
                float dist_l = point_aabb_dist_sq(p, tree_nodes(l));
                float dist_r = point_aabb_dist_sq(p, tree_nodes(r));

                // Push the further child first so the closer child is popped next
                if (dist_l < dist_r) {
                    stack[stack_ptr++] = r;
                    stack[stack_ptr++] = l;
                } else {
                    stack[stack_ptr++] = l;
                    stack[stack_ptr++] = r;
                }
            }
        }
        
        out_dist_sq = min_dist_sq;
        return best_facet;
    }

    KOKKOS_INLINE_FUNCTION
    void get_closest_point_on_triangle(
        const float p[3], 
        const size_t f_idx, 
        const DCArrayKokkos<float>& verts, 
        float closest[3]) const
    {
        // Pull vertices from the STL view
        float a[3] = {verts(f_idx, 0, 0), verts(f_idx, 0, 1), verts(f_idx, 0, 2)};
        float b[3] = {verts(f_idx, 1, 0), verts(f_idx, 1, 1), verts(f_idx, 1, 2)};
        float c[3] = {verts(f_idx, 2, 0), verts(f_idx, 2, 1), verts(f_idx, 2, 2)};

        // Edge vectors
        float ab[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        float ac[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
        float ap[3] = {p[0] - a[0], p[1] - a[1], p[2] - a[2]};
        
        // Dot products for vertex region A
        float d1 = ab[0]*ap[0] + ab[1]*ap[1] + ab[2]*ap[2];
        float d2 = ac[0]*ap[0] + ac[1]*ap[1] + ac[2]*ap[2];
        if (d1 <= 0.0f && d2 <= 0.0f) {
            for(int i=0; i<3; ++i) closest[i] = a[i];
            return;
        }

        // Vertex region B
        float bp[3] = {p[0] - b[0], p[1] - b[1], p[2] - b[2]};
        float d3 = ab[0]*bp[0] + ab[1]*bp[1] + ab[2]*bp[2];
        float d4 = ac[0]*bp[0] + ac[1]*bp[1] + ac[2]*bp[2];
        if (d3 >= 0.0f && d4 <= d3) {
            for(int i=0; i<3; ++i) closest[i] = b[i];
            return;
        }

        // Edge region AB
        float vc = d1*d4 - d3*d2;
        if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) {
            float v = d1 / (d1 - d3);
            for(int i=0; i<3; ++i) closest[i] = a[i] + v * ab[i];
            return;
        }

        // Vertex region C
        float cp[3] = {p[0] - c[0], p[1] - c[1], p[2] - c[2]};
        float d5 = ab[0]*cp[0] + ab[1]*cp[1] + ab[2]*cp[2];
        float d6 = ac[0]*cp[0] + ac[1]*cp[1] + ac[2]*cp[2];
        if (d6 >= 0.0f && d5 <= d6) {
            for(int i=0; i<3; ++i) closest[i] = c[i];
            return;
        }

        // Edge region AC
        float vb = d5*d2 - d1*d6;
        if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
            float w = d2 / (d2 - d6);
            for(int i=0; i<3; ++i) closest[i] = a[i] + w * ac[i];
            return;
        }

        // Edge region BC
        float va = d3*d6 - d5*d4;
        if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) {
            float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
            for(int i=0; i<3; ++i) closest[i] = b[i] + w * (c[i] - b[i]);
            return;
        }

        // Face region
        float denom = 1.0f / (va + vb + vc);
        float v = vb * denom;
        float w = vc * denom;
        for(int i=0; i<3; ++i) closest[i] = a[i] + ab[i] * v + ac[i] * w;
    }

};





// ==============================================================
// ------------------- Binary STL Reader-------------------
// ==============================================================

// BINARY STL READER - (Note: it can ONLY read binary stl files)
void binary_stl_reader(std::string stl_file_path, stl_data& stl_data){
    // Open the binary file
    std::string filename = stl_file_path;
    auto input = std::ifstream{filename, std::ifstream::in | std::ifstream::binary};
    
    // Check that the file is actually open
    if (input.is_open()) {
        std::cout << "Opening .stl file... \n";
    }
    else {
        std::cout << "WARNING: .stl file is not open \n";
        std::cout << "File path: " << filename << "\n";
        throw std::runtime_error("**** Error in STL File Path ****");
    }
    
    // Initialize variables
    char heading[81];
    unsigned int n_facets;
    
    // Read the file's heading
    char* ptr1 = heading;
    input.read(ptr1, 80);
    heading[80] = '\0'; // to ensure a proper C string
    std::cout << "File heading: " << heading << "\n";
    
    // Read the number of facets in the file
    unsigned int* ptr2 = &n_facets;
    input.read(reinterpret_cast<char*>(ptr2), sizeof(unsigned int));
    std::cout << "Reading " << n_facets << " facets..." << "\n";


    stl_data.initialize(n_facets);
    
    // Initialize storage variables

    float normalp[3];
    float v1p[3];
    float v2p[3];
    float v3p[3];
    
    // Pull data from file
    float* ptr3 = normalp;
    float* ptr4 = v1p;
    float* ptr5 = v2p;
    float* ptr6 = v3p;
    for (size_t i = 0; i < n_facets; ++i) {
        size_t n_bytes = 3 * sizeof(float);
        input.read(reinterpret_cast<char*>(ptr3), n_bytes);
        input.read(reinterpret_cast<char*>(ptr4), n_bytes);
        input.read(reinterpret_cast<char*>(ptr5), n_bytes);
        input.read(reinterpret_cast<char*>(ptr6), n_bytes);
        input.seekg(2,input.cur);

        // read the normal
        for (int j=0; j<3; j++) {
            stl_data.normal.host(i,j) = normalp[j];
        }

        // read the vertices
        for (int j=0; j<3; j++) {
            stl_data.vertices.host(i,0,j) = v1p[j];
        }
        for (int j=0; j<3; j++) {
            stl_data.vertices.host(i,1,j) = v2p[j];
        }
        for (int j=0; j<3; j++) {
            stl_data.vertices.host(i,2,j) = v3p[j];
        }

        // center is I/O-bound arithmetic — compute here on host
        for (int j=0; j<3; j++) {
            stl_data.center.host(i, j) = (v1p[j] + v2p[j] + v3p[j]) / 3.0f;
        }
    }
    input.close();

    stl_data.normal.update_device();
    stl_data.vertices.update_device();
    stl_data.center.update_device();
}


#endif // STL_UTILS_H