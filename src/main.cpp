#include <cmath>
#include <iostream>
#include <initializer_list>
#include <memory>
#include <map>
#include <stdexcept>
#include <utility>

// DISCLAIMER
//
// The preliminary Spec Tests that you are able to run only test that your code compilesand provides an indication about how well your code can be tested with the complete Spec Test.The preliminary Spec Tests DO NOT TEST the correct functioning of your code.
// The maximum score you can obtain with the preliminary Spec Test is 10 / 100. If you obtain this maximum score it indicates that your code is likely compatible with the complete Spec Test.
// It is your responsibility to test that your code is correct.
//
// After the deadline, your code will be tested with the complete Spec Test which does test the correct functioning of your code.
// The maximum score you can obtain with the complete Spec Test is 100 / 100.
// The complete Spec Test will determine your grade for this assignment.

///////////////////////////////////////////////////////////////////////////////

// General problem description

// Consider the following initial boundary value problem:
// ∂u∂t−αΔuu(x,t)u(x,0)=0=0=∏k=0n−1sin(πxk)on Ω,∀x∈∂Ω,∀x∈Ω

// where n
// is the number of dimensions of the domain, Ω=[0,1]n is the domain (the unit square in 2D or the unit cube in 3D), ∂Ω is the boundary of the domain, u is the temperature as a function of space x=(x0,…,xn−1)⊤∈Rn and time t and Δ is the Laplace operator, in an n

// -dimensional space given by

// Δ:=∑k=0n−1∂2∂x2k.

// The exact solution to above initial boundary value problem is given by

// u(x,t)=e−nπ2αt∏k=0n−1sin(πxk)=e−nπ2αtu(x,0).

// We are going to apply a finite difference discretisation. We create a mesh of mn
// equidistant interior nodes, m nodes per dimension, hence the distance between two nodes is Δx=1/(m+1). All nodes are labeled from 0 to mn (excluding) such that (with a few possible exceptions at the boundary) node j+1 is the neighbour of node j in dimension 0, node j+m is the neighbour of node j in dimension 1, and in general, node j+mk is the neighbour of node j in dimension k. The following illustrates a two-dimensional mesh with m=3

// :

// x
//  1
//  ↑
// 1    ┏━━━━━┯━━━━━━┯━━━━━┯━━━━━┓
//      ┃     │      │     │     ┃
//      ┃     │      │     │     ┃
// 0.75 ┠─────6──────7─────8─────┨
//      ┃     │      │     │     ┃
//      ┃     │      │     │     ┃
// 0.5  ┠─────3──────4─────5─────┨
//      ┃     │      │     │     ┃
//      ┃     │      │     │     ┃
// 0.25 ┠─────0──────1─────2─────┨
//      ┃     │      │     │     ┃
//      ┃     │      │     │     ┃
// 0    ┗━━━━━┷━━━━━━┷━━━━━┷━━━━━┛
//      0   0.25   0.5   0.75   1  → x
//                                0

// Approximating the Laplace operator with a second-order central finite difference discretization and the remaining semi-discrete differential equation with the Backward Euler method gives the following system of equations

// Mwl+1=wl

// where vector wl∈Rmn
// is the discrete approximation of the heat u at time t=lΔt, with vector component i refering to node i. The matrix M

// is given by

// Mij=Iij−αΔtΔx2∑k=0n−1Dkij∀i,j∈0,1,…,mn,

// where I
// is the identity matrix and Dk is the discrete approximation of the second-order derivative in dimension k∈0,1,…,n−1

// :

// DkijDkijDkij=−2=1=1if j=i,if j is a left neighbour of i in dimension k,if j is a right neighbour of i in dimension k.

// Note: Not all nodes have left or right neighbours in every dimension. In the example mesh given above node 0
// has only right neighbours: node 1 in dimension 0 and node 3 in dimension 1. Node 5 has two left neighbours (node 4 in dimension 0 and node 2 in dimension 1), but only one right neighbour: node 8 in dimension 1.

///////////////////////////////////////////////////////////////////////////////

// 1. VECTOR

// Create the class

// template <typename T>
// class Vector
// {...};

// whereby the Vector’s elements are of type T. 
// The Vector class must provide the following functionality:

// Constructors and destructor

//     A default constructor that sets the length to zero.
//     A copy constructor and a move constructor that creates a Vector from another Vector.
//     A constructor that takes a length as an argument and allocates the internal 
//       data structures accordingly.
//     A constructor that takes an initialiser list representing the contents of this Vector, 
//       allocates the internal data structures and initialises the Vector’s content accordingly.
//     A destructor.

// Operators and functions

//     A copy assignment operator and a move assignment operator from another Vector.
//     An operator[](int i) that returns a reference to the i-th entry of the vector. 
//        Implement an overload of this operator that returns a constant reference. 
//        Both operators can be used to access the entries in functions that are implemented outside the Vector class.
//     Arithmetic operators operator+ and operator- to add and subtract two Vectors. 
//         These operators must support Vectors of different types, whereby the resulting 
//         Vector has to be of the type that dominates (e.g., double dominates float). If the Vectors have different lengths, all operators must throw an exception.
//     Arithmetic operators operator* between a scalar and a Vector (w=s⋅v

// ) and a Vector and a scalar (w=v⋅s

//     ), whereby the scalar and the Vector can have different types and the resulting Vector must be of the dominating type.
//     A function len that returns the length of the Vector. This function can be used to retrieve the length in functions that are implemented outside the Vector class.

// Create a function dot that computes the standard inner product of two Vectors. The function must have the following signature:

// template<typename T, typename U>
// typename std::common_type<T,U>::type
// dot(const Vector<T>& lhs, 
//     const Vector<U>& rhs)
// {...}

// If the Vectors have different lengths, the dot function must throw an exception.


template <typename T>
class Vector
{
    // Your implementation of the Vector class starts here
    Vector 

};

template<typename T, typename U>
typename std::common_type<T, U>::type
dot(const Vector<T>& lhs,
    const Vector<U>& rhs)
{

    // Your implementation of the dot function starts here
}

// [JF] Vector Unit test 
void TestVector() 
{

}

///////////////////////////////////////////////////////////////////////////////

// 2. MATRIX

// Create the class

// template <typename T> 
// class Matrix
// {...};

// that represents a sparse matrix, whereby the Matrix’s elements are of type T. A sparse 
// matrix is a matrix in which most elements are zero so that it makes sense to store only 
// the non-zero entries to save memory.

// Consider for example the matrix

// A=⎛⎝⎜⎜⎜1001020003501001⎞⎠⎟⎟⎟

// with only 7 (out of 16) non-zero entries. If we store the non-zero entries in the form

// (0,0)(0,3)(1,1)(1,2)(2,2)(3,0)(3,3)→1→1→2→3→5→1→1

// then the matrix-vector multiplication

// y=Ax

// can be realised in a loop over all non-zero entries ij
// of the matrix A

// with the update formula

// yi:=yi+Aijxj

// Hint: Use the std::map class from the C++ standard library to store the non-zero matrix entries. 
// Use std::pair<int, int> as the key type to store the row and column 
// positions and T as the data type to store the value at the particular matrix position.

// Even though std::map allows to search for a key it is not a good idea to use this feature for 
// the sparse matrix-vector multiplication since this is very inefficient. Instead, use an iterator 
// over all entries of the map and extract row and column positions and the data value via

// int i   = it->first[0];
// int j   = it>-first[1];
// T value = it->second;

// The Matrix class must provide the following functionality:

//     A constructor that accepts two integer values (number of rows and columns) 
// and initialises the internal data structures.

//     Matrix<double> M(10, 20); // initialise M with 10 rows and 20 columns

// It is not allowed to change the dimensions of the matrix afterwards.

// A destructor.

// An operator[](const std::pair<int, int>& ij) that returns the matrix entry ij (i.e. the entry at row i
// and column j

// ) by reference.

// M[{0,0}] = 1.0; // set value at row 0, column 0 to 1.0
// M[{1,2}] = 2.0; // set value at row 1, column 2 to 2.0

// Implement an overload of this operators that returns a constant reference and 
// throws an exception if the entry is not present

// std::cout << M[{0,0}] << std::endl; // prints 1.0
// std::cout << M[{3,3}] << std::endl; // throws an exception

// An operator* between a Matrix and a Vector object that implements the 
// sparse matrix-vector product. The operator must have the following signature:

// template<typename T, typename U>
// Vector<typename std::common_type<T,U>::type>
// operator*(const Matrix<T>& lhs, 
//           const Vector<U>& rhs)
// { ... }

// If the dimension of the Matrix and the Vector are not compatible, the operator must throw an exception.

template <typename T>
class Matrix
{
    // Start your implementation of the matrix class here
};

// [JF] Matrix unit test
void TestMatrix() 
{

}

///////////////////////////////////////////////////////////////////////////////

// 3. Conjugate gradient method

//     Create a function named cg that solves a linear system using the 
// Conjugate Gradient method (CG), given in pseudocode by....

// p_0 = r_0 = b - A x_0
// for k = 0, 1, ..., maxiter-1
//     alpha_k = dot(r_k, r_k) / dot(A p_k, p_k)
//     x_(k+1) = x_k + alpha_k p
//     r_(k+1) = r_k - alpha_k A p

//     if dot(r_(k+1), r_(k+1)) < tol*tol
//        stop

//     beta_k  = dot(r_(k+1), r_(k+1)) / dot(r_k, r_k)
//     p_(k+1) = r_(k+1) + beta_k p_k


// Here, A is a symmetric positive definite matrix, 
// x_k, r_k and p_k are vectors, 
// dot is the standard l2
// -inner product, 
// tol is an absolute tolerance for the residual 
// maxiter is the maximum allowed number of iterations.

// The cg function must have the following signature: DONE 

// The third argument serves both as the initial guess 
// (x_0 in the pseudocode) and as the result 
// (x_k+1 in the pseudocode, where k is the last iteration).  

// The function must return the number of iterations used to 
// achieve the desired tolerance if the Conjugate Gradient
//  method converged within maxiter iterations and -1 otherwise, 
// that is, if maxiter iterations have been reached without 
// reaching convergence. -> DONE 

template<typename T>
int cg(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x, T tol = (T)1e-8, int maxiter = 100)
{


// p_0 = r_0 = b - A x_0
// for k = 0, 1, ..., maxiter-1
//     alpha_k = dot(r_k, r_k) / dot(A p_k, p_k)
//     x_(k+1) = x_k + alpha_k p;
//     r_(k+1) = r_k - alpha_k A p;

//     if dot(r_(k+1), r_(k+1)) < tol*tol
//        return k;

//     beta_k  = dot(r_(k+1), r_(k+1)) / dot(r_k, r_k)
//     p_(k+1) = r_(k+1) + beta_k p_k

    for(int k = 0 ; k < limit; k++)
    {
        bool convergence = false;
        if (convergence)
        {
            return k;
        }
    }


    // Your implementation of the cg function starts here
    return -1;
}

///////////////////////////////////////////////////////////////////////////////

// 4. Finite difference discretization of the heat equation

//     Create the class

//     template <int n, typename T> Heat
//     {...};

// that represents the n-dimensional initial boundary value problem given in the introduction and has the following functionality:
// * Create a constructor that accepts three arguments: the diffusion coefficient alpha (of type T), the number of points per dimension m (of type int) and the time-step size dt (of type T). The constructor must create the iteration matrix M given in the introduction and store the result as an attribute with type Matrix<T>.
// * Create a method Vector<T> exact(T t) const that returns the exact solution at time t evaluated at all interior grid points (see introduction).
// * Create a method Vector<T> solve(T t) const that solves the initial boundary value problem given in the introduction using the procedure given at the end of the introduction until time t and returns the numerical solution at the last time step.

template <int n, typename T>
class Heat
{
    // Your implementation of the heat class starts here
};

///////////////////////////////////////////////////////////////////////////////

// 5. Testing your implementation

//     Test your implementation by creating a solver for the one-dimensional initial value problem with diffusion coefficient alpha=0.3125, time step dt=0.001 and m=99, solve the problem until time t=1 and compare the numerical solution with the exact solution.

//     To verify that matrix M is assembled correctly, set alpha=0.3125, dt=0.1, m=3 and compare it with

// All zero entries have been left out.

// Test your implementation by creating a solver for the two-dimensional initial value problem with diffusion coefficient alpha=0.3125, time step dt=0.001 and m=99, solve the problem until time t=0.5 and compare the numerical solution with the exact solution.

// To verify that matrix M is assembled correctly, set alpha=0.3125, dt=0.1, m=3 and compare it with

int main(int argc, char* argv[])
{
    std::cout << "kaas" << std::endl;
    std::cout << "henk" << std::endl;
    std::cout << "fristi" << std::endl;
    return 0;
}