/*****************************************************************************
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


#ifndef MATAR_H
#define MATAR_H 


//==============================================================================
//   MATAR > 'MAT'rix & 'AR'rays
//
//     This code is for creating and viewing matricies and arrays using
//     standard notation e.g., A(i,j,k).  The indicies for a matrix
//     start at 1 for consistancy with mathematics.  The indicies for
//     an array start at 0 for consistancy with the C++ langauge.  The
//     entries in the arrays and matrices are stored in a 1D array that
//     is contiguous in memory, but are accessed using e.g., A(i,j,k).
//     The classes that have f_array in the name stride through the 1D  
//     array following the Fortran convention (i.e., first index varies
//     the quickest).  The classes with c_array in the name stride through
//     the 1D array following the C convention (i.e., last index varies
//     the quickest).
//
//   view_f_array(array_in, dim1, ...)
//   f_array_t(dim1,...)
//   view_f_matrix(matrix_in, dim1, ...)
//   f_matrix_t(dim1,...)
//
//   view_c_array(array_in, dim1, ...)
//   c_array_t(dim1,...)
//   view_c_matrix(matrix_in, dim1, ...)
//   c_matrix_t(dim1,...)
//
//
//==============================================================================
//
//  matar.cpp
//
//     g++ -02 --std=c++14 matar.cpp
//
//  Created by Nathaniel Morgan on 5/18/19.
//
//==============================================================================

//#include <iostream>
#include <stdio.h>


//==============================================================================
//   Fortran stride array related classes  (first index varies the quickest)
//==============================================================================

// view a 1D vector as an array(i,...,n), where i=[0:N-1]
template <typename T>
class view_f_array {

private:
    size_t dim1, dim2, dim3, dim4, dim5, dim6;
    T * this_array;
    
public:
    
    // default constructor
    view_f_array ();
    
    
    //--- 2D array ---
    
    // overloaded constructor
    view_f_array (T *some_array,
                size_t some_dim1,
                size_t some_dim2)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        this_array = some_array;
    }
    
    // overload operator() to access data as array(i,j)
    //  where i=[0:N-1], j=[0:N-1]
    inline T& operator()(size_t i, size_t j) const
    {
        return this_array[i + j*dim1];
    }
    
    
    //--- 3D array ---
    
    // overloaded constructor
    view_f_array (T *some_array,
                size_t some_dim1,
                size_t some_dim2,
                size_t some_dim3)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        this_array = some_array;
    }
    
    // overload operator() to access data as array(i,j,k)
    //  where i=[0:N-1], j=[0:N-1], k=[0:N-1]
    inline T& operator()(size_t i, size_t j, size_t k) const
    {
        return this_array[i + j*dim1 + k*dim1*dim2];
    }

    
    //--- 4D array ---
    
    // overloaded constructor
    view_f_array (T *some_array,
                size_t some_dim1,
                size_t some_dim2,
                size_t some_dim3,
                size_t some_dim4)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        dim4 = some_dim4;
        this_array = some_array;
    }
    
    // overload operator() to access data as array(i,j,k,l)
    //  where i=[0:n-1], j=[0:N-1], k=[0:N-1], l=[0:N-1]
    inline T& operator()(size_t i, size_t j, size_t k, size_t l) const
    {
        return this_array[i + j*dim1 + k*dim1*dim2 + l*dim1*dim2*dim3];
    }

    
    //--- 5D array ---
    
    // overloaded constructor
    view_f_array (T *some_array,
                size_t some_dim1,
                size_t some_dim2,
                size_t some_dim3,
                size_t some_dim4,
                size_t some_dim5)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        dim4 = some_dim4;
        dim5 = some_dim5;
        this_array = some_array;
    }
    
    // overload operator() to access data as array(i,j,k,l,m)
    //  where i=[0:N-1], j=[0:N-1], k=[0:N-1], l=[0:N-1], m=[0:N-1]
    inline T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const
    {
        return this_array[i + j*dim1 + k*dim1*dim2 + l*dim1*dim2*dim3
                        + m*dim1*dim2*dim3*dim4];
    }

    
    //--- 6D array ---
    
    // overloaded constructor
    view_f_array (T *some_array,
                size_t some_dim1,
                size_t some_dim2,
                size_t some_dim3,
                size_t some_dim4,
                size_t some_dim5,
                size_t some_dim6)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        dim4 = some_dim4;
        dim5 = some_dim5;
        dim6 = some_dim6;
        this_array = some_array;
    }
    
    // overload operator() to access data as array(i,j,k,l,m,n)
    //  where i=[0:N-1], j=[0:N-1], k=[0:N-1], l=[0:N-1], m=[0:N-1]
    inline T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const
    {
        return this_array[i + j*dim1 + k*dim1*dim2 + l*dim1*dim2*dim3
                        + m*dim1*dim2*dim3*dim4 + n*dim1*dim2*dim3*dim4*dim5];
    }

    
}; // end of view_f_array



// create a 1D vector that is accessed as an array(i,...,n), where i=[0:N-1]
template <typename T>
class f_array_t {
    
private:
    size_t dim1, dim2, dim3, dim4, dim5, dim6;
    T * this_array;
    
public:
    
    // default constructor
    f_array_t ();
    
   //--- 1D array ---
   
   // overloaded constructor
   f_array_t (size_t some_dim1)
   {
      dim1 = some_dim1;
      this_array = new T[dim1];
   }
   
   // overload operator() to access data as array(i) where i=[0:N-1]
   inline T& operator()(size_t i) const
   {
      return this_array[i];
   }
   
    //--- 2D array ---
    
    // overloaded constructor
    f_array_t (size_t some_dim1, size_t some_dim2)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        this_array = new T[dim1*dim2];
    }
    
    // overload operator() to access data as array(i,j)
    //  where i=[0:N-1], j=[0:N-1]
    inline T& operator()(size_t i, size_t j) const
    {
        return this_array[i + j*dim1];
    }
    
    
    //--- 3D array ---
    
    // overloaded constructor
    f_array_t (size_t some_dim1, size_t some_dim2, size_t some_dim3)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        this_array = new T[some_dim1*some_dim2*some_dim3];
    }
    
    // overload operator() to access data as array(i,j,k)
    //  where i=[0:N-1], j=[0:N-1], k=[0:N-1]
    inline T& operator()(size_t i, size_t j, size_t k) const
    {
        return this_array[i + j*dim1 + k*dim1*dim2];
    }
    
    
    //--- 4D array ---
    
    // overloaded constructor
    f_array_t (size_t some_dim1,
               size_t some_dim2,
               size_t some_dim3,
               size_t some_dim4)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        dim4 = some_dim4;
        this_array = new T[some_dim1*some_dim2*some_dim3*some_dim4];
    }
    
    // overload operator() to access data as array(i,j,k,l)
    //  where i=[0:n-1], j=[0:N-1], k=[0:N-1], l=[0:N-1]
    inline T& operator()(size_t i, size_t j, size_t k, size_t l) const
    {
        return this_array[i + j*dim1 + k*dim1*dim2 + l*dim1*dim2*dim3];
    }
    
    
    //--- 5D array ---
    
    // overloaded constructor
    f_array_t (size_t some_dim1,
               size_t some_dim2,
               size_t some_dim3,
               size_t some_dim4,
               size_t some_dim5)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        dim4 = some_dim4;
        dim5 = some_dim5;
        this_array = new T[some_dim1*some_dim2*some_dim3*some_dim4*some_dim5];
    }
    
    // overload operator() to access data as array(i,j,k,l,m)
    //  where i=[0:N-1], j=[0:N-1], k=[0:N-1], l=[0:N-1], m=[0:N-1]
    inline T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const
    {
        return this_array[i + j*dim1 + k*dim1*dim2 + l*dim1*dim2*dim3
                        + m*dim1*dim2*dim3*dim4];
    }
    
    
    //--- 6D array ---
    
    // overloaded constructor
    f_array_t (size_t some_dim1,
               size_t some_dim2,
               size_t some_dim3,
               size_t some_dim4,
               size_t some_dim5,
               size_t some_dim6)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        dim4 = some_dim4;
        dim5 = some_dim5;
        dim6 = some_dim6;
        this_array = new T[some_dim1*some_dim2*some_dim3
                          *some_dim4*some_dim5*some_dim6];
    }
    
    // overload operator() to access data as array(i,j,k,l,m,n)
    //  where i=[0:N-1], j=[0:N-1], k=[0:N-1], l=[0:N-1], m=[0:N-1]
    inline T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const
    {
        return this_array[i + j*dim1 + k*dim1*dim2 + l*dim1*dim2*dim3
                        + m*dim1*dim2*dim3*dim4 + n*dim1*dim2*dim3*dim4*dim5];
    }

    
    
    // deconstructor
    ~f_array_t ( )
    {
        delete[] this_array;
    }
    
}; // end of f_array_t





//==============================================================================
//   C stride array related classes  (last index varies the quickest)
//==============================================================================

// view a 1D vector as an array(i,...,n), where i=[0:N-1]
template <typename T>
class view_c_array {

private:
    size_t dim1, dim2, dim3, dim4, dim5, dim6;
    T * this_array;
    
public:
    
    // default constructor
    view_c_array ();
    
    
    //--- 1D array ---
    
    // overloaded constructor
    view_c_array (T *some_array,
                size_t some_dim1)
    {
        dim1 = some_dim1;
        this_array = some_array;
    }
    
    // overload operator() to access data as array(i,j)
    //  where i=[0:N-1], j=[0:N-1]
    inline T& operator()(size_t i) const
    {
        return this_array[i];
    }
    
    
    //--- 2D array ---
    
    // overloaded constructor
    view_c_array (T *some_array,
                size_t some_dim1,
                size_t some_dim2)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        this_array = some_array;
    }
    
    // overload operator() to access data as array(i,j)
    //  where i=[0:N-1], j=[0:N-1]
    inline T& operator()(size_t i, size_t j) const
    {
        return this_array[j + i*dim2];
    }
    
    
    //--- 3D array ---
    
    // overloaded constructor
    view_c_array (T *some_array,
                size_t some_dim1,
                size_t some_dim2,
                size_t some_dim3)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        this_array = some_array;
    }
    
    // overload operator() to access data as array(i,j,k)
    //  where i=[0:N-1], j=[0:N-1], k=[0:N-1]
    inline T& operator()(size_t i, size_t j, size_t k) const
    {
        return this_array[k + j*dim3 + i*dim3*dim2];
    }

    
    //--- 4D array ---
    
    // overloaded constructor
    view_c_array (T *some_array,
                size_t some_dim1,
                size_t some_dim2,
                size_t some_dim3,
                size_t some_dim4)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        dim4 = some_dim4;
        this_array = some_array;
    }
    
    // overload operator() to access data as array(i,j,k,l)
    //  where i=[0:n-1], j=[0:N-1], k=[0:N-1], l=[0:N-1]
    inline T& operator()(size_t i, size_t j, size_t k, size_t l) const
    {
        return this_array[l + k*dim4 + j*dim4*dim3 + i*dim4*dim3*dim2];
    }

    
    //--- 5D array ---
    
    // overloaded constructor
    view_c_array (T *some_array,
                size_t some_dim1,
                size_t some_dim2,
                size_t some_dim3,
                size_t some_dim4,
                size_t some_dim5)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        dim4 = some_dim4;
        dim5 = some_dim5;
        this_array = some_array;
    }
    
    // overload operator() to access data as array(i,j,k,l,m)
    //  where i=[0:N-1], j=[0:N-1], k=[0:N-1], l=[0:N-1], m=[0:N-1]
    inline T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const
    {
        return this_array[m + l*dim5 + k*dim5*dim4 + j*dim5*dim4*dim3
                        + i*dim5*dim4*dim3*dim2];
    }

    
    //--- 6D array ---
    
    // overloaded constructor
    view_c_array (T *some_array,
                size_t some_dim1,
                size_t some_dim2,
                size_t some_dim3,
                size_t some_dim4,
                size_t some_dim5,
                size_t some_dim6)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        dim4 = some_dim4;
        dim5 = some_dim5;
        dim6 = some_dim6;
        this_array = some_array;
    }
    
    // overload operator() to access data as array(i,j,k,l,m,n)
    //  where i=[0:N-1], j=[0:N-1], k=[0:N-1], l=[0:N-1], m=[0:N-1]
    inline T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const
    {
        return this_array[n + m*dim6 + l*dim6*dim5 + k*dim6*dim5*dim4 
                        + j*dim6*dim5*dim4*dim3 + i*dim6*dim5*dim4*dim3*dim2];
    }

    
}; // end of view_c_array


// create a 1D vector that is accessed as an array(i,...,n), where i=[0:N-1]
template <typename T>
class c_array_t {
    
private:
    size_t dim1, dim2, dim3, dim4, dim5, dim6;
    T * this_array;
    
public:
    
    // default constructor
    c_array_t ();
    
   //--- 1D array ---
   
   // overloaded constructor
   c_array_t (size_t some_dim1)
   {
      dim1 = some_dim1;
      this_array = new T[dim1];
   }
   
   // overload operator() to access data as array(i) where i=[0:N-1]
   inline T& operator()(size_t i) const
   {
      return this_array[i];
   }
   
    //--- 2D array ---
    
    // overloaded constructor
    c_array_t (size_t some_dim1, size_t some_dim2)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        this_array = new T[dim1*dim2];
    }
    
    // overload operator() to access data as array(i,j)
    //  where i=[0:N-1], j=[0:N-1]
    inline T& operator()(size_t i, size_t j) const
    {
        return this_array[j + i*dim2];
    }
    
    
    //--- 3D array ---
    
    // overloaded constructor
    c_array_t (size_t some_dim1, size_t some_dim2, size_t some_dim3)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        this_array = new T[some_dim1*some_dim2*some_dim3];
    }
    
    // overload operator() to access data as array(i,j,k)
    //  where i=[0:N-1], j=[0:N-1], k=[0:N-1]
    inline T& operator()(size_t i, size_t j, size_t k) const
    {
        return this_array[k + j*dim3 + i*dim3*dim2];
    }
    
    
    //--- 4D array ---
    
    // overloaded constructor
    c_array_t (size_t some_dim1,
               size_t some_dim2,
               size_t some_dim3,
               size_t some_dim4)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        dim4 = some_dim4;
        this_array = new T[some_dim1*some_dim2*some_dim3*some_dim4];
    }
    
    // overload operator() to access data as array(i,j,k,l)
    //  where i=[0:n-1], j=[0:N-1], k=[0:N-1], l=[0:N-1]
    inline T& operator()(size_t i, size_t j, size_t k, size_t l) const
    {
        return this_array[l + k*dim4 + j*dim4*dim3 + i*dim4*dim3*dim2];
    }
    
    
    //--- 5D array ---
    
    // overloaded constructor
    c_array_t (size_t some_dim1,
               size_t some_dim2,
               size_t some_dim3,
               size_t some_dim4,
               size_t some_dim5)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        dim4 = some_dim4;
        dim5 = some_dim5;
        this_array = new T[some_dim1*some_dim2*some_dim3*some_dim4*some_dim5];
    }
    
    // overload operator() to access data as array(i,j,k,l,m)
    //  where i=[0:N-1], j=[0:N-1], k=[0:N-1], l=[0:N-1], m=[0:N-1]
    inline T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const
    {
        return this_array[m + l*dim5 + k*dim5*dim4 + j*dim5*dim4*dim3
                        + i*dim5*dim4*dim3*dim2];
    }
    
    
    //--- 6D array ---
    
    // overloaded constructor
    c_array_t (size_t some_dim1,
               size_t some_dim2,
               size_t some_dim3,
               size_t some_dim4,
               size_t some_dim5,
               size_t some_dim6)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        dim4 = some_dim4;
        dim5 = some_dim5;
        dim6 = some_dim6;
        this_array = new T[some_dim1*some_dim2*some_dim3
                          *some_dim4*some_dim5*some_dim6];
    }
    
    // overload operator() to access data as array(i,j,k,l,m,n)
    //  where i=[0:N-1], j=[0:N-1], k=[0:N-1], l=[0:N-1], m=[0:N-1]
    inline T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const
    {
        return this_array[n + m*dim6 + l*dim6*dim5 + k*dim6*dim5*dim4 
                        + j*dim6*dim5*dim4*dim3 + i*dim6*dim5*dim4*dim3*dim2];
    }

    
    
    // deconstructor
    ~c_array_t ( )
    {
        delete[] this_array;
    }
    
}; // end of c_array_t





//==============================================================================
//   Fortran stride matrix related classes
//==============================================================================

// view a 1D vector as a Matrix(i,...,n), where i=[1:N]
template <typename T>
class view_f_matrix {
    
private:
    size_t dim1, dim2, dim3, dim4, dim5, dim6;
    T * this_matrix;
    
public:
    
    // default constructor
    view_f_matrix ();
    
    
    //--- 1D matrix ---
    
    // overloaded constructor
    view_f_matrix (T *some_matrix,
                size_t some_dim1)
    {
        dim1 = some_dim1;
        this_matrix = some_matrix;
    }
    
    // overload operator() to access data as array(i,j)
    //  where i=[0:N-1], j=[0:N-1]
    inline T& operator()(size_t i) const
    {
        return this_matrix[i-1];
    }
    
    
    //--- 2D matrix ---
    
    // overloaded constructor
    view_f_matrix (T *some_matrix, size_t some_dim1, size_t some_dim2)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        this_matrix = some_matrix;
    }
    
    // overload operator() to access data as matrix(i,j) where i=[1:N], j=[1:N]
    inline T& operator()(size_t i, size_t j) const
    {
        return this_matrix[(i-1) + (j-1)*dim1];
    }
    
    
    //--- 3D matrix ---
    
    // overloaded constructor
    view_f_matrix (T *some_matrix, size_t some_dim1, size_t some_dim2, size_t some_dim3)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        this_matrix = some_matrix;
    }
    
    // overload operator() to access data as matrix(i,j,k)
    //  where i=[1:N], j=[1:N], k=[1:N]
    inline T& operator()(size_t i, size_t j, size_t k) const
    {
        return this_matrix[(i-1) + (j-1)*dim1 + (k-1)*dim1*dim2];
    }
    
    
    //--- 4D matrix ---
    
    // overloaded constructor
    view_f_matrix (T *some_matrix,
                 size_t some_dim1,
                 size_t some_dim2,
                 size_t some_dim3,
                 size_t some_dim4)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        dim4 = some_dim4;
        this_matrix = some_matrix;
    }
    
    // overload operator() to access data as matrix(i,j,k,l)
    //  where i=[1:N], j=[1:N], k=[1:N], l=[1:N]
    inline T& operator()(size_t i, size_t j, size_t k, size_t l) const
    {
        return this_matrix[(i-1) + (j-1)*dim1 + (k-1)*dim1*dim2
                      + (l-1)*dim1*dim2*dim3];
    }
    
    
    //--- 5D matrix ---
    
    // overloaded constructor
    view_f_matrix (T *some_matrix,
                 size_t some_dim1,
                 size_t some_dim2,
                 size_t some_dim3,
                 size_t some_dim4,
                 size_t some_dim5)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        dim4 = some_dim4;
        dim5 = some_dim5;
        this_matrix = some_matrix;
    }
    
    // overload operator() to access data as matrix(i,j,k,l,m)
    //  where i=[1:N], j=[1:N], k=[1:N], l=[1:N], m=[1:N]
    inline T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const
    {
        return this_matrix[(i-1) + (j-1)*dim1 + (k-1)*dim1*dim2
                      + (l-1)*dim1*dim2*dim3 + (m-1)*dim1*dim2*dim3*dim4];
    }
    
    
    //--- 6D matrix ---
    
    // overloaded constructor
    view_f_matrix (T *some_matrix,
                 size_t some_dim1,
                 size_t some_dim2,
                 size_t some_dim3,
                 size_t some_dim4,
                 size_t some_dim5,
                 size_t some_dim6)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        dim4 = some_dim4;
        dim5 = some_dim5;
        dim6 = some_dim6;
        this_matrix = some_matrix;
    }
    
    // overload operator() to access data as matrix(i,j,k,l,m,n)
    //  where i=[1:N], j=[1:N], k=[1:N], l=[1:N], m=[1:N]
    inline T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const
    {
        return this_matrix[(i-1) + (j-1)*dim1 + (k-1)*dim1*dim2
                         + (l-1)*dim1*dim2*dim3 + (m-1)*dim1*dim2*dim3*dim4 +
                           (n-1)*dim1*dim2*dim3*dim4*dim5];
    }
    
    
}; // end of view_f_matrix



// create a 1D vector that is accessed as an matrix(i,...,n), where i=[0:N-1]
template <typename T>
class f_matrix_t {
    
private:
    size_t dim1, dim2, dim3, dim4, dim5, dim6;
    T * this_matrix;
    
public:
    
   // default constructor
   f_matrix_t ();

   //--- 1D matrix ---
   
   // overloaded constructor
   f_matrix_t (size_t some_dim1)
   {
      dim1 = some_dim1;
      this_matrix = new T[dim1];
   }
   
   // overload operator() to access data as matrix(i) where i=[1:N]
   inline T& operator()(size_t i) const
   {
     return this_matrix[i-1];
   }
   
    
    //--- 2D matrix ---
    
    // overloaded constructor
    f_matrix_t (size_t some_dim1, size_t some_dim2)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        this_matrix = new T[dim1*dim2];
    }
    
    // overload operator() to access data as matrix(i,j) where i=[1:N], j=[1:N]
    inline T& operator()(size_t i, size_t j) const
    {
        return this_matrix[(i-1) + (j-1)*dim1];
    }
    
    
    //--- 3D array ---
    
    // overloaded constructor
    f_matrix_t (size_t some_dim1, size_t some_dim2, size_t some_dim3)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        this_matrix = new T[some_dim1*some_dim2*some_dim3];
    }
    
    // overload operator() to access data as matrix(i,j,k)
    //  where i=[1:N], j=[1:N], k=[1:N]
    inline T& operator()(size_t i, size_t j, size_t k) const
    {
        return this_matrix[(i-1) + (j-1)*dim1 + (k-1)*dim1*dim2];
    }
    
    
    //--- 4D array ---
    
    // overloaded constructor
    f_matrix_t (size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        dim4 = some_dim4;
        this_matrix = new T[some_dim1*some_dim2*some_dim3*some_dim4];
    }
    
    // overload operator() to access data as matrix(i,j,k,l)
    //  where i=[1:N], j=[1:N], k=[1:N], l=[1:N]
    inline T& operator()(size_t i, size_t j, size_t k, size_t l) const
    {
        return this_matrix[(i-1) + (j-1)*dim1 + (k-1)*dim1*dim2
                      + (l-1)*dim1*dim2*dim3];
    }
    
    
    //--- 5D array ---
    
    // overloaded constructor
    f_matrix_t (size_t some_dim1,
                size_t some_dim2,
                size_t some_dim3,
                size_t some_dim4,
                size_t some_dim5)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        dim4 = some_dim4;
        dim5 = some_dim5;
        this_matrix = new T[some_dim1*some_dim2*some_dim3*some_dim4*some_dim5];
    }
    
    // overload operator() to access data as matrix(i,j,k,l,m)
    //  where i=[1:N], j=[1:N], k=[1:N], l=[1:N], m=[1:N]
    inline T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const
    {
        return this_matrix[(i-1) + (j-1)*dim1 + (k-1)*dim1*dim2
                         + (l-1)*dim1*dim2*dim3 + (m-1)*dim1*dim2*dim3*dim4];
    }
    
    
    //--- 6D array ---
    
    // overloaded constructor
    f_matrix_t (size_t some_dim1,
                size_t some_dim2,
                size_t some_dim3,
                size_t some_dim4,
                size_t some_dim5,
                size_t some_dim6)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        dim4 = some_dim4;
        dim5 = some_dim5;
        dim6 = some_dim6;
        this_matrix = new T[some_dim1*some_dim2*some_dim3
                       *some_dim4*some_dim5*some_dim6];
    }
    
    // overload operator() to access data as matrix(i,j,k,l,m,n)
    //  where i=[1:N], j=[1:N], k=[1:N], l=[1:N], m=[1:N]
    inline T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const
    {
        return this_matrix[(i-1) + (j-1)*dim1 + (k-1)*dim1*dim2
                         + (l-1)*dim1*dim2*dim3 + (m-1)*dim1*dim2*dim3*dim4 +
                           (n-1)*dim1*dim2*dim3*dim4*dim5];
    }
    
    
    
    // deconstructor
    ~f_matrix_t ( )
    {
        delete[] this_matrix;
    }
    
}; // end of f_matrix_t






//==============================================================================
//   C stride matrix related classes  (last index varies the quickest)
//==============================================================================

// view a 1D vector as an array(i,...,n), where i=[0:N-1]
template <typename T>
class view_c_matrix {

private:
    size_t dim1, dim2, dim3, dim4, dim5, dim6;
    T * this_matrix;
    
public:
    
    // default constructor
    view_c_matrix ();
    
    
    //--- 1D array ---
    
    // overloaded constructor
    view_c_matrix (T *some_matrix,
                size_t some_dim1)
    {
        dim1 = some_dim1;
        this_matrix = some_matrix;
    }
    
    // overload operator() to access data as array(i,j)
    //  where i=[0:N-1], j=[0:N-1]
    inline T& operator()(size_t i) const
    {
        return this_matrix[i-1];
    }
    
    
    //--- 2D array ---
    
    // overloaded constructor
    view_c_matrix (T *some_matrix,
                size_t some_dim1,
                size_t some_dim2)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        this_matrix = some_matrix;
    }
    
    // overload operator() to access data as array(i,j)
    //  where i=[0:N-1], j=[0:N-1]
    inline T& operator()(size_t i, size_t j) const
    {
        return this_matrix[(j-1) + (i-1)*dim2];
    }
    
    
    //--- 3D array ---
    
    // overloaded constructor
    view_c_matrix (T *some_matrix,
                size_t some_dim1,
                size_t some_dim2,
                size_t some_dim3)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        this_matrix = some_matrix;
    }
    
    // overload operator() to access data as array(i,j,k)
    //  where i=[0:N-1], j=[0:N-1], k=[0:N-1]
    inline T& operator()(size_t i, size_t j, size_t k) const
    {
        return this_matrix[(k-1) + (j-1)*dim3 + (i-1)*dim3*dim2];
    }

    
    //--- 4D array ---
    
    // overloaded constructor
    view_c_matrix (T *some_matrix,
                size_t some_dim1,
                size_t some_dim2,
                size_t some_dim3,
                size_t some_dim4)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        dim4 = some_dim4;
        this_matrix = some_matrix;
    }
    
    // overload operator() to access data as array(i,j,k,l)
    //  where i=[0:n-1], j=[0:N-1], k=[0:N-1], l=[0:N-1]
    inline T& operator()(size_t i, size_t j, size_t k, size_t l) const
    {
        return this_matrix[(l-1) + (k-1)*dim4 + (j-1)*dim4*dim3 + (i-1)*dim4*dim3*dim2];
    }

    
    //--- 5D array ---
    
    // overloaded constructor
    view_c_matrix (T *some_matrix,
                size_t some_dim1,
                size_t some_dim2,
                size_t some_dim3,
                size_t some_dim4,
                size_t some_dim5)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        dim4 = some_dim4;
        dim5 = some_dim5;
        this_matrix = some_matrix;
    }
    
    // overload operator() to access data as array(i,j,k,l,m)
    //  where i=[0:N-1], j=[0:N-1], k=[0:N-1], l=[0:N-1], m=[0:N-1]
    inline T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const
    {
        return this_matrix[(m-1) + (l-1)*dim5 + (k-1)*dim5*dim4 + (j-1)*dim5*dim4*dim3
                         + (i-1)*dim5*dim4*dim3*dim2];
    }

    
    //--- 6D array ---
    
    // overloaded constructor
    view_c_matrix (T *some_matrix,
                   size_t some_dim1,
                   size_t some_dim2,
                   size_t some_dim3,
                   size_t some_dim4,
                   size_t some_dim5,
                   size_t some_dim6)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        dim4 = some_dim4;
        dim5 = some_dim5;
        dim6 = some_dim6;
        this_matrix = some_matrix;
    }
    
    // overload operator() to access data as array(i,j,k,l,m,n)
    //  where i=[0:N-1], j=[0:N-1], k=[0:N-1], l=[0:N-1], m=[0:N-1]
    inline T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const
    {
        return this_matrix[(n-1) + (m-1)*dim6 + (l-1)*dim6*dim5 + (k-1)*dim6*dim5*dim4 
                         + (j-1)*dim6*dim5*dim4*dim3 + (i-1)*dim6*dim5*dim4*dim3*dim2];
    }

    
}; // end of view_c_matrix


// create a 1D vector that is accessed as an array(i,...,n), where i=[0:N-1]
template <typename T>
class c_matrix_t {
    
private:
    size_t dim1, dim2, dim3, dim4, dim5, dim6;
    T * this_matrix;
    
public:
    
    // default constructor
    c_matrix_t ();
    
    
   //--- 1D matrix ---
   
   // overloaded constructor
   c_matrix_t (size_t some_dim1)
   {
      dim1 = some_dim1;
      this_matrix = new T[dim1];
   }
   
   // overload operator() to access data as array(i) where i=[0:N-1]
   inline T& operator()(size_t i) const
   {
      return this_matrix[i-1];
   }
   
    //--- 2D array ---
    
    // overloaded constructor
    c_matrix_t (size_t some_dim1, size_t some_dim2)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        this_matrix = new T[dim1*dim2];
    }
    
    // overload operator() to access data as array(i,j)
    //  where i=[0:N-1], j=[0:N-1]
    inline T& operator()(size_t i, size_t j) const
    {
        return this_matrix[(j-1) + (i-1)*dim2];
    }
    
    
    //--- 3D array ---
    
    // overloaded constructor
    c_matrix_t (size_t some_dim1, size_t some_dim2, size_t some_dim3)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        this_matrix = new T[some_dim1*some_dim2*some_dim3];
    }
    
    // overload operator() to access data as array(i,j,k)
    //  where i=[0:N-1], j=[0:N-1], k=[0:N-1]
    inline T& operator()(size_t i, size_t j, size_t k) const
    {
        return this_matrix[(k-1) + (j-1)*dim3 + (i-1)*dim3*dim2];
    }
    
    
    //--- 4D array ---
    
    // overloaded constructor
    c_matrix_t (size_t some_dim1,
                size_t some_dim2,
                size_t some_dim3,
                size_t some_dim4)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        dim4 = some_dim4;
        this_matrix = new T[some_dim1*some_dim2*some_dim3*some_dim4];
    }
    
    // overload operator() to access data as array(i,j,k,l)
    //  where i=[0:n-1], j=[0:N-1], k=[0:N-1], l=[0:N-1]
    inline T& operator()(size_t i, size_t j, size_t k, size_t l) const
    {
        return this_matrix[(l-1) + (k-1)*dim4 + (j-1)*dim4*dim3 + (i-1)*dim4*dim3*dim2];
    }
    
    
    //--- 5D array ---
    
    // overloaded constructor
    c_matrix_t (size_t some_dim1,
                size_t some_dim2,
                size_t some_dim3,
                size_t some_dim4,
                size_t some_dim5)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        dim4 = some_dim4;
        dim5 = some_dim5;
        this_matrix = new T[some_dim1*some_dim2*some_dim3*some_dim4*some_dim5];
    }
    
    // overload operator() to access data as array(i,j,k,l,m)
    //  where i=[0:N-1], j=[0:N-1], k=[0:N-1], l=[0:N-1], m=[0:N-1]
    inline T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const
    {
        return this_matrix[(m-1) + (l-1)*dim5 + (k-1)*dim5*dim4 + (j-1)*dim5*dim4*dim3
                         + (i-1)*dim5*dim4*dim3*dim2];
    }
    
    
    //--- 6D array ---
    
    // overloaded constructor
    c_matrix_t (size_t some_dim1,
                size_t some_dim2,
                size_t some_dim3,
                size_t some_dim4,
                size_t some_dim5,
                size_t some_dim6)
    {
        dim1 = some_dim1;
        dim2 = some_dim2;
        dim3 = some_dim3;
        dim4 = some_dim4;
        dim5 = some_dim5;
        dim6 = some_dim6;
        this_matrix = new T[some_dim1*some_dim2*some_dim3
                          *some_dim4*some_dim5*some_dim6];
    }
    
    // overload operator() to access data as array(i,j,k,l,m,n)
    //  where i=[0:N-1], j=[0:N-1], k=[0:N-1], l=[0:N-1], m=[0:N-1]
    inline T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const
    {
        return this_matrix[(n-1) + (m-1)*dim6 + (l-1)*dim6*dim5 + (k-1)*dim6*dim5*dim4 
                         + (j-1)*dim6*dim5*dim4*dim3 + (i-1)*dim6*dim5*dim4*dim3*dim2];
    }

    
    
    // deconstructor
    ~c_matrix_t ( )
    {
        delete[] this_matrix;
    }
    
}; // end of c_matrix_t

#endif //MATAR_H