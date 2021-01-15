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

// General problem description
//
// Consider the following initial boundary value problem :
//


template <typename T>
class Vector
{
    // Your implementation of the Vector class starts here
};

template<typename T, typename U>
typename std::common_type<T, U>::type
dot(const Vector<T>& lhs,
    const Vector<U>& rhs)
{
    // Your implementation of the dot function starts here
}




template <typename T>
class Matrix
{
    // Start your implementation of the matrix class here
};

template<typename T>
int cg(const Matrix<T>& A,
    const Vector<T>& b,
    Vector<T>& x,
    T                tol = (T)1e-8,
    int              maxiter = 100)
{
    // Your implementation of the cg function starts here
}

template <int n, typename T>
class Heat
{
    // Your implementation of the heat class starts here
};

int main(int argc, char* argv[])
{

    return 0;
}