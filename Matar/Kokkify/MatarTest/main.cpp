#include <iostream>
#include <chrono>         // To access timing calipers 
#include "pseudo_mesh.hpp"
#include "inherited_inits.hpp"

int main() {

    int size_i = 5, size_j = 4, size_k = 3;
    int big_i = 5000, big_j = 2000, big_k = 3000;


    // Kokkos GPU test
   
    Kokkos::initialize();
    {
    auto begin = std::chrono::high_resolution_clock::now();
//#ifdef NOTRAGGED
    using policy2D = Kokkos::MDRangePolicy< Kokkos::Rank<2> >;
    policy2D array_type = policy2D({0,0}, {size_i, size_j});
    policy2D matrix_type = policy2D({1,1}, {size_i+1, size_j+1});
    policy2D splice_type = policy2D({0,0}, {size_j, size_k});
    using policy3D = Kokkos::MDRangePolicy< Kokkos::Rank<3> >;
    policy3D array_type3 = policy3D({0,0,0}, {size_i, size_j, size_k});
    policy3D matrix_type3 = policy3D({1,1,1}, {size_i+1, size_j+1, size_k+1});
    

    printf("~~~~~~~~~~~~~~~~~GPU TEST~~~~~~~~~~~~~~~\n"); 

    int num_parent = 2; // number of materials
    Parent1D parent("parent", num_parent); // Initialize Kokkos View on the GPU of type material, size num_parent
    auto h_parent = Kokkos::create_mirror_view(parent); // Create a host view of the Kokkos View


    AllocateHost(h_parent, 0, BABY2_SIZE); // Function performed on Host to do raw Kokkos allocation of baby2 GPU space inside of Host data structure
    AllocateHost(h_parent, 1, BABY1_SIZE); // Function performed on Host to do raw Kokkos allocation of baby1 GPU space inside of Host data structure

    Kokkos::deep_copy(parent, h_parent); // deep copy Host data (allocated above) to the GPU Kokkos View. GPU View now has the class space allocated

    InitChildModels(parent, 0, baby2{}); // Kokkos Function to create new instances of the baby2 model on the GPU
    InitChildModels(parent, 1, baby1{1.4,1.0}); // Kokkos Function to create new instances of the baby1 models on the GPU

    // Model test, also shows a Kokkos reduction
    double value_1;
    Kokkos::parallel_reduce(
        "CheckValues",                             
        num_parent,                   
        KOKKOS_LAMBDA(const int idx, real_t &lsum)
        {
             lsum += parent(idx).child->math(2.0, 4.0);
        }
        , value_1);                                

    printf("value %f\n", value_1);



    ClearDeviceModels(parent); // Kokkos Function to call deconstructors of objects on the GPU

    FreeHost(h_parent); // Function performed on Host to free the allocated GPU classes inside of the Host mirror


    /*
    auto partest = InheritedArray2L <parent_models> (num_parent);
    partest.AllocateHost<child_models>(BABY2_SIZE, partest(0, 0).child);
    partest.AllocateHost<child_models>(BABY1_SIZE, partest(1, 0).child);
    partest.AllocateGPU();
    partest.InitModels(partest(0, 1).child, baby2{});
    partest.InitModels(partest(1, 1).child, baby1{1.4,1.0});


    Kokkos::parallel_for("Initialize", num_parent, KOKKOS_LAMBDA(const int i) {
            });
    Kokkos::fence();
    */   


    printf("Pseudo Mesh Kokkos\n");
    pseudo_mesh pmesh;
    pmesh.init(size_i, size_j);

    Kokkos::parallel_for("PseudoMesh", 1, KOKKOS_LAMBDA(const int&) {
            //pmesh->init(size_i, size_j);
        });

    //using TeamPolicy = Kokkos::TeamPolicy<Kokkos::Cuda>;
    //Kokkos::parallel_for ("PseudoMesh", TeamPolicy(1, Kokkos::AUTO), KOKKOS_LAMBDA (const TeamPolicy::member_type &teamMember) {
        //pmesh->init(size_i, size_j);
    //    });

    auto carray_view = ViewCArrayKokkos <real_t> (pmesh.carray.pointer(), size_i, size_j);
    auto cmatrix_view = ViewCMatrixKokkos <real_t> (pmesh.cmatrix.pointer(), size_i, size_j);
    auto farray_view = ViewFArrayKokkos <real_t> (pmesh.farray.pointer(), size_i, size_j);
    auto fmatrix_view = ViewFMatrixKokkos <real_t> (pmesh.fmatrix.pointer(), size_i, size_j);

    printf("\nCArray\n");
    Kokkos::parallel_for("CArrayValues", array_type, KOKKOS_LAMBDA(const int i, const int j) {
            pmesh.carray(i, j) = (real_t) j + i * size_j;
            printf("(%d, %d) %lf %lf\n", i, j, pmesh.carray(i, j), carray_view(i, j));
        });
    Kokkos::fence();
    printf("\nCMatrix\n");
    Kokkos::parallel_for("CMatrixValues", matrix_type, KOKKOS_LAMBDA(const int i, const int j) {
            pmesh.cmatrix(i, j) = (real_t) (j - 1) + (i - 1) * size_j;
            printf("(%d, %d) %lf %lf\n", i, j, pmesh.cmatrix(i, j), cmatrix_view(i, j));
        });
    Kokkos::fence();
    printf("\nFArray\n");
    Kokkos::parallel_for("FArrayValues", array_type, KOKKOS_LAMBDA(const int i, const int j) {
            pmesh.farray(i, j) = (real_t) i + j * size_i;
            printf("(%d, %d) %lf %lf\n", i, j, pmesh.farray(i, j), farray_view(i, j));
        });
    Kokkos::fence();
    printf("\nFMatrix\n");
    Kokkos::parallel_for("FMatrixValues", matrix_type, KOKKOS_LAMBDA(const int i, const int j) {
            pmesh.fmatrix(i, j) = (real_t) (i - 1) + (j - 1) * size_i;
            printf("(%d, %d) %lf %lf\n", i, j, pmesh.fmatrix(i, j), fmatrix_view(i, j));
        });
    Kokkos::fence();

    printf("\nRaggedRight\n");
    Kokkos::parallel_for ("RaggedRightTeam", TeamPolicy(size_i, Kokkos::AUTO), KOKKOS_LAMBDA (const TeamPolicy::member_type &teamMember) {
            const int i = teamMember.league_rank();
            const int i_stride = pmesh.raggedright.stride(i);
            Kokkos::parallel_for (Kokkos::TeamThreadRange (teamMember, i_stride), [=] (const int j) {
                pmesh.raggedright(i, j) = (double) i + j * i_stride; // not the exact placement, needs start index
                printf("%d %d %lf\n", i, j, pmesh.raggedright(i, j));
            });
        });
    Kokkos::fence();
    
    printf("\nRaggedDown\n");
    Kokkos::parallel_for ("RaggedDownTeam", TeamPolicy(size_j, Kokkos::AUTO), KOKKOS_LAMBDA (const TeamPolicy::member_type &teamMember) {
            const int j = teamMember.league_rank();
            const int j_stride = pmesh.raggeddown.stride(j);
            Kokkos::parallel_for (Kokkos::TeamThreadRange (teamMember, j_stride), [=] (const int i) {
                pmesh.raggeddown(i, j) = (double) j + i * j_stride; // not the exact placement, needs start index
                printf("%d %d %lf\n", i, j, pmesh.raggeddown(i, j));
            });
        });
    Kokkos::fence();

    /*------------------------- Stream Benchmarks ----------------------------------------*/
#ifdef STREAM
    // Kokkos stream benchmark test (for FMatrixKokkos and CArrayKokkos objects)
    printf("\n Kokkos stream benchmark (FMatrixKokkos and CArrayKokkos)\n");
    
    double scalar = 3.0;
    size_t nsize = 256 * 256 * 256;
    pseudo_mesh bench;
    bench.init(nsize);

    printf("Stream benchmark with %d elements.\n", nsize);
    
    Kokkos::parallel_for("Initialize", nsize, KOKKOS_LAMBDA(const int i) {
            bench.arr1(i) = 1.0;
            bench.arr2(i) = 2.0;
            bench.arr3(i) = 0.0;
            });
    Kokkos::fence();
   
    // Perform copy benchmark
    begin = std::chrono::high_resolution_clock::now();
     
    Kokkos::parallel_for("Copy", nsize, KOKKOS_LAMBDA(const int i) {
            bench.arr3(i) = bench.arr1(i);
            });
    Kokkos::fence();
    

    std::chrono::duration<double> copy_time = std::chrono::high_resolution_clock::now() - begin;

    std::cout << "Copy time: " << copy_time.count() << " s." << std::endl;

    // Perform scaling benchmark
    begin = std::chrono::high_resolution_clock::now();
    
    Kokkos::parallel_for("Scale", nsize, KOKKOS_LAMBDA(const int i) {
            bench.arr2(i) = scalar * bench.arr3(i);
            });
    Kokkos::fence();
    
    std::chrono::duration<double> scale_time = std::chrono::high_resolution_clock::now() - begin;

    std::cout << "Scale time: " << scale_time.count() << " s." << std::endl;

    // Perform sum benchmark
    begin = std::chrono::high_resolution_clock::now();
    
    Kokkos::parallel_for("Sum", nsize, KOKKOS_LAMBDA(const int i) {
            bench.arr3(i) = bench.arr1(i) + bench.arr2(i);
            });
    Kokkos::fence();
    
    std::chrono::duration<double> sum_time = std::chrono::high_resolution_clock::now() - begin; 

    std::cout << "Sum time: " << sum_time.count() << " s." << std::endl;

    // Perform triad benchmark
    begin = std::chrono::high_resolution_clock::now();
    
    Kokkos::parallel_for("Triad", nsize, KOKKOS_LAMBDA(const int i) {
            bench.arr1(i) = bench.arr2(i) + 
                                  (scalar * bench.arr3(i));
            });
    Kokkos::fence();
    
    std::chrono::duration<double> triad_time = std::chrono::high_resolution_clock::now() - begin;

    std::cout << "Triad time: " << triad_time.count() << " s." << std::endl;
//#endif

    // Extra test
    size_t nsize_3D = 256;
    int scalar = 5;
    using policy3D = Kokkos::MDRangePolicy< Kokkos::Rank<3> >;

    // Create 3D CArrayKokkos objects
    auto arr1_3D = CArrayKokkos <real_t> (nsize_3D, nsize_3D, nsize_3D);
    auto arr2_3D = CArrayKokkos <real_t> (nsize_3D, nsize_3D, nsize_3D);
    auto arr3_3D = CArrayKokkos <real_t> (nsize_3D, nsize_3D, nsize_3D);

    // Create 3D Kokkos View objects
    RMatrix3D kv_arr1_3D("kv_arr1_3D", nsize_3D, nsize_3D, nsize_3D);
    RMatrix3D kv_arr2_3D("kv_arr2_3D", nsize_3D, nsize_3D, nsize_3D);
    RMatrix3D kv_arr3_3D("kv_arr3_3D", nsize_3D, nsize_3D, nsize_3D);

    // Initialize 3D CArrayKokkos and 3D Kokkos View objects
    policy3D array_type_STREAM = policy3D({0, 0, 0},
                                          {nsize_3D, nsize_3D, nsize_3D});

    Kokkos::parallel_for("Initialize (3D)", array_type_STREAM, 
                         KOKKOS_LAMBDA(const int i, const int j, const int k) {
                arr1_3D(i, j, k) = 1.0;
                arr2_3D(i, j, k) = 2.0;
                arr3_3D(i, j, k) = 0.0;

                kv_arr1_3D(i, j, k) = 1.0;
                kv_arr2_3D(i, j, k) = 2.0;
                kv_arr3_3D(i, j, k) = 0.0;
                });
    Kokkos::fence();

        begin = std::chrono::high_resolution_clock::now();
        Kokkos::parallel_for("Triad (3D)", array_type_STREAM, 
                         KOKKOS_LAMBDA(const int i, const int j, const int k) {
                arr1_3D(i, j, k) = arr2_3D(i, j, k) + (scalar * arr3_3D(i, j, k));
                });
        std::chrono::duration<double> triad_time_3D = std::chrono::high_resolution_clock::now() - begin;
        //Kokkos::fence();
        
    std::cout << "Triad time: " << triad_time_3D.count() << " s." << std::endl;

        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_for("Triad (3D KV)", array_type_STREAM, 
                         KOKKOS_LAMBDA(const int i, const int j, const int k) {
                kv_arr1_3D(i, j, k) = kv_arr2_3D(i, j, k) + (scalar * kv_arr3_3D(i, j, k));
                });
        std::chrono::duration<double> triad_time_kv_3D = std::chrono::high_resolution_clock::now() - begin;
        //Kokkos::fence();



    std::cout << "View Triad time: " << triad_time_kv_3D.count() << " s." << std::endl;
#endif

#ifdef NOTSTREAM
    size_t matrix_size = 64 * 64;
    auto stride_test = CArrayKokkos <size_t> (matrix_size);

    Kokkos::parallel_for("StrideValuesTime", matrix_size, KOKKOS_LAMBDA(const int i) {
            stride_test(i) = i + 1;
        });

    begin = std::chrono::high_resolution_clock::now();

    ProfileRegionStart("Matar");
    RaggedRightArrayKokkos <real_t> raggedright_test;
    raggedright_test = RaggedRightArrayKokkos <real_t> (stride_test);
    Kokkos::parallel_for ("RaggedRightTeamTime", TeamPolicy(matrix_size, Kokkos::AUTO), KOKKOS_LAMBDA (const TeamPolicy::member_type &teamMember) {
            const int i = teamMember.league_rank();
            const int i_stride = raggedright_test.stride(i);
            Kokkos::parallel_for (Kokkos::TeamThreadRange (teamMember, i_stride), [=] (const int j) {
                raggedright_test(i, j) = (double) i + j * i_stride; // not the exact placement, needs start index
                real_t tempval = raggedright_test(i, j);
            });
        });
    Kokkos::fence();
    ProfileRegionEnd();
    std::chrono::duration<double> raggedright_time = std::chrono::high_resolution_clock::now() - begin; 

    std::cout << "raggedright time: " << raggedright_time.count() << " s." << std::endl;

    begin = std::chrono::high_resolution_clock::now();

    ProfileRegionStart("KV");;
    RMatrix2D kv_raggedright("kv_raggedright", matrix_size, matrix_size); 
    Kokkos::parallel_for ("RaggedRightTeamTime2", TeamPolicy(matrix_size, Kokkos::AUTO), KOKKOS_LAMBDA (const TeamPolicy::member_type &teamMember) {
            const int i = teamMember.league_rank();
            const int i_stride = stride_test(i);
            Kokkos::parallel_for (Kokkos::TeamThreadRange (teamMember, i_stride), [=] (const int j) {
                kv_raggedright(i, j) = (double) i + j * i_stride; // not the exact placement, needs start index
                real_t tempval = raggedright_test(i, j);
            });
        });
    Kokkos::fence();
    ProfileRegionEnd();
    std::chrono::duration<double> kv_raggedright_time = std::chrono::high_resolution_clock::now() - begin; 

    std::cout << "kv_raggedright time: " << kv_raggedright_time.count() << " s." << std::endl;
#endif 

    }
    Kokkos::finalize();


    printf("--- finished ---\n");

    return 0;
}
