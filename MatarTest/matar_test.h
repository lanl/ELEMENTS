#ifndef MATAR_TEST_H
#define MATAR_TEST_H

#include <iostream>
#include <stdio.h>
//#include <Kokkos_Core.hpp>
#include "kokkos_alias.h"


template <typename T>
class c_array_t {
    
private:
    size_t dim1, dim2, dim3, dim4, dim5, dim6, length;
    T * this_array;
    
public:
    // default constructor
    c_array_t ();
    
    //--- 1D array ---
   
    // overloaded constructor
    c_array_t (size_t some_dim1);

    c_array_t (size_t some_dim1, size_t some_dim2);

    // overload opperator() to access data as array(i) where i=[0:N-1]
    T& operator() (size_t i);

    T& operator() (size_t i, size_t j);

    // overload copy assignment operator
    c_array_t& operator= (const c_array_t& temp);

    // deconstructor
    ~c_array_t ();
};


template <typename T>
c_array_t<T>::c_array_t() {}

template <typename T>
c_array_t<T>::c_array_t(size_t some_dim1) {
    dim1 = some_dim1;
    length = dim1;
    this_array = new T[length];
}

template <typename T>
T& c_array_t<T>::operator() (size_t i) {
    return this_array[i];
}

template <typename T>
c_array_t<T>::c_array_t(size_t some_dim1, size_t some_dim2) {
    dim1 = some_dim1;
    dim2 = some_dim2;
    length = dim1 * dim2;
    this_array = new T[length];
}

template <typename T>
T& c_array_t<T>::operator() (size_t i, size_t j) {
    return this_array[j + i*dim2];
}

// for the object assignment THIS = c_array_t <> TEMP(n,m,etc.) 
template <typename T>
c_array_t<T>& c_array_t<T>::operator= (const c_array_t& temp) {
    // do nothing if the assignment is of form x = x
    if (this != &temp) { 
        // Assign old elements equal to the new ones
        dim1 = temp.dim1;
        dim2 = temp.dim2;
        dim3 = temp.dim3;
        dim4 = temp.dim4;
        dim5 = temp.dim5;
        dim6 = temp.dim6;
        length = temp.length;
        // Allocate new array
        this_array = new T[length];
    }

}


template <typename T>
c_array_t<T>::~c_array_t () {
    delete[] this_array;
}


//Kokkos
template <typename T>
class k_array_t {
    
private:
    size_t dim1, dim2, dim3, dim4, dim5, dim6, length;
    T * this_array;
    
public:
    // default constructor
    //KOKKOS_FUNCTION
    k_array_t ();
    
    //--- 1D array ---
   
    // overloaded constructor
    //KOKKOS_FUNCTION
    k_array_t (size_t some_dim1);

    //KOKKOS_FUNCTION
    k_array_t (size_t some_dim1, size_t some_dim2);

    // overload opperator() to access data as array(i) where i=[0:N-1]
    KOKKOS_FUNCTION
    T& operator() (size_t i);

    KOKKOS_FUNCTION
    T& operator() (size_t i, size_t j);

    // overload copy assignment operator
    //KOKKOS_FUNCTION
    k_array_t& operator= (const k_array_t& temp);

    // destructor
    KOKKOS_FUNCTION
    ~k_array_t ();
};


template <typename T>
//KOKKOS_FUNCTION
k_array_t<T>::k_array_t() {}

template <typename T>
//KOKKOS_FUNCTION
k_array_t<T>::k_array_t(size_t some_dim1) {
    //Kokkos::parallel_for("PseudoMesh", 1, KOKKOS_LAMBDA(const int&) {
    dim1 = some_dim1;
    length = dim1;
    this_array = new T[length];
    //});
}

template <typename T>
KOKKOS_FUNCTION
T& k_array_t<T>::operator() (size_t i) {
    return this_array[i];
}

template <typename T>
//KOKKOS_FUNCTION
k_array_t<T>::k_array_t(size_t some_dim1, size_t some_dim2) {
    //Kokkos::parallel_for("PseudoMesh", 1, KOKKOS_LAMBDA(const int&) {
    dim1 = some_dim1;
    dim2 = some_dim2;
    length = dim1 * dim2;
    this_array = new T[length];
    //});
}

template <typename T>
KOKKOS_FUNCTION
T& k_array_t<T>::operator() (size_t i, size_t j) {
    return this_array[i + j*dim1];
}

// for the object assignment THIS = k_array_t <> TEMP(n,m,etc.) 
template <typename T>
//KOKKOS_FUNCTION
k_array_t<T>& k_array_t<T>::operator= (const k_array_t& temp) {
    // do nothing if the assignment is of form x = x
    Kokkos::parallel_for("PseudoMesh", 1, KOKKOS_LAMBDA(const int&) {
    if (this != &temp) { 
        // Assign old elements equal to the new ones
        dim1 = temp.dim1;
        dim2 = temp.dim2;
        dim3 = temp.dim3;
        dim4 = temp.dim4;
        dim5 = temp.dim5;
        dim6 = temp.dim6;
        length = temp.length;
        // Allocate new array
        this_array = new T[length];
    }
    });
    // Return *this
    return *this;
}

template <typename T>
KOKKOS_FUNCTION
k_array_t<T>::~k_array_t () {
    delete[] this_array;
}

//template <typename T>

template <typename T>
class k_test_t {

using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
private:
    size_t dim1, dim2, dim3, dim4, dim5, dim6, length;
    //T * this_array;
    //Kokkos::View<T *,Layout,ExecSpace> this_array; 
    TArray1D this_array;
public:
    // default constructor
    //KOKKOS_FUNCTION
    k_test_t ();
    
    //--- 1D array ---
   
    // overloaded constructor
    //KOKKOS_FUNCTION
    k_test_t (size_t some_dim1);

    //KOKKOS_FUNCTION
    k_test_t (size_t some_dim1, size_t some_dim2);

    // overload opperator() to access data as array(i) where i=[0:N-1]
    KOKKOS_FUNCTION
    T& operator() (size_t i) const;

    KOKKOS_FUNCTION
    T& operator() (size_t i, size_t j) const;

    // overload copy assignment operator
    //KOKKOS_FUNCTION
    k_test_t& operator= (const k_test_t& temp);

    // destructor
    //KOKKOS_FUNCTION
    ~k_test_t ();
};


template <typename T>
//KOKKOS_FUNCTION
k_test_t<T>::k_test_t() {}

template <typename T>
//KOKKOS_FUNCTION
k_test_t<T>::k_test_t(size_t some_dim1) {

    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;

    dim1 = some_dim1;
    length = dim1;
    //this_array = new T[length];
    //this_array = (T *) Kokkos::kokkos_malloc<Kokkos::CudaSpace>(length*sizeof(T));
    this_array = TArray1D("this_array", length);
}

template <typename T>
KOKKOS_FUNCTION
T& k_test_t<T>::operator() (size_t i) const {
    return this_array(i);
}

template <typename T>
//KOKKOS_FUNCTION
k_test_t<T>::k_test_t(size_t some_dim1, size_t some_dim2) {

    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;

    dim1 = some_dim1;
    dim2 = some_dim2;
    length = dim1 * dim2;
    //this_array = new T[length];
    //this_array = (T *) Kokkos::kokkos_malloc<Kokkos::CudaSpace>(length*sizeof(T));
    //this_array = Kokkos::View<T *,Layout,ExecSpace>("this_array", length);
    this_array = TArray1D("this_array", length);
}

template <typename T>
KOKKOS_FUNCTION
T& k_test_t<T>::operator() (size_t i, size_t j) const {
    return this_array(j + i*dim2);
}

// for the object assignment THIS = k_array_t <> TEMP(n,m,etc.) 
template <typename T>
//KOKKOS_FUNCTION
k_test_t<T>& k_test_t<T>::operator= (const k_test_t& temp) {

    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;

    // do nothing if the assignment is of form x = x
    if (this != &temp) { 
        // Assign old elements equal to the new ones
        dim1 = temp.dim1;
        dim2 = temp.dim2;
        dim3 = temp.dim3;
        dim4 = temp.dim4;
        dim5 = temp.dim5;
        dim6 = temp.dim6;
        length = temp.length;
        // Allocate new array
        //this_array = new T[length];
        //this_array = (T *) Kokkos::kokkos_malloc<Kokkos::CudaSpace>(length*sizeof(T));
        //this_array = Kokkos::View<T *,Layout,ExecSpace>("this_array", length);
        this_array = TArray1D("this_array", length);
    }
    // Return *this
    return *this;
}

template <typename T>
//KOKKOS_FUNCTION
k_test_t<T>::~k_test_t () {
    //delete[] this_array;
    //Kokkos::kokkos_free(this_array);
}

#endif
