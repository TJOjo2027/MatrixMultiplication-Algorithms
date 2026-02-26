/*
* Author: Oluwatosin Ojo
* Assignment Title: Program 2 - SquareMatrix
* Assignment Description: This assignment explores the different ways to
*       computer matrix multiplication using different algorithm types.
* Due Date: 2/22/2026
* Date Created: 2/10/2026
* Date Last Modified: 2/22/2026
*/

#ifndef SQUAREMATRIX_H
#define SQUAREMATRIX_H


#include <ostream>
#include <iomanip>
using namespace std;

// early struct declarations
struct ThreadArguments;
struct PARAM;

struct SquareMatrix{
    size_t dim;
    int** data; // points to a [dim x dim] square matrix

    /*
    * description: constructor for SquareMatrix object
    * return: none, creates SquareMatrix Object
    * precondition: assume dim is the dimension of the matrix
    * postcondition: creates a square matrix of dimensions dim x dim
    *
    */
    SquareMatrix(size_t dim) : dim(dim) {
        // allocate memory for the first data pointer
        data = new int*[dim];
        // allocate memory for each block inside the data pointer
        for (size_t i = 0; i < dim; i++) {
            data[i] = new int[dim];
            for (size_t j = 0; j < dim; j++) {
                // fill with zeros
                data[i][j] = 0;
            }
        }
    }

    /*
    * description: destructor for SquareMatrix object
    * return: none, destorys SquareMatrix Object
    * precondition: assume there is a SquareMatrix to be destroyed
    * postcondition: deallocates the memory being used
    *
    */
    ~SquareMatrix() {
        // delete the data being pointed to
        for (size_t i = 0 ; i < dim; i++) {
            delete[] data[i];
        }
        // delete the data pointer
        delete[] data;
    }

    /*
    * description: displays the matrix object to output stream
    * return: none, creates SquareMatrix Object
    * precondition: the matrix calling this has been constructed
    * postcondition: displays the matrix to the output stream
    *
    */
    void display(ostream& os) {
        for (size_t r = 0; r < dim; r++) {
            for (size_t c = 0; c < dim; c++) {
                os << setw(4) << data[r][c];
            }
            os << endl;
        }
    }
};

struct ThreadArguments {
    const SquareMatrix* A;
    const SquareMatrix* B;
    SquareMatrix* C;
    size_t A_subMatrix_1;
    size_t  B_subMatrix_1;
    size_t A_subMatrix_2;
    size_t B_subMatrix_2;
    size_t C_subMatrix;

    /*
    * description: creates a struct instance of ThreadArguments
    * return: none, creates ThreadArguments Object
    * precondition: the parameters needed are declared already
    * postcondition: creates the ThreadArguments object
    *
    */
    ThreadArguments(const SquareMatrix* A, const SquareMatrix* B,
        SquareMatrix* C, size_t a_sub_1, size_t b_sub_1,
        size_t a_sub_2, size_t b_sub_2, size_t c_sub) {
        this -> A = A;
        this -> B = B;
        this -> C = C;
        this -> A_subMatrix_1 = a_sub_1;
        this -> B_subMatrix_1 = b_sub_1;
        this -> A_subMatrix_2 = a_sub_2;
        this -> B_subMatrix_2 = b_sub_2;
        this -> C_subMatrix = c_sub;
    }
};


struct PARAM {
    size_t dimensions;
    size_t A_row_offset, A_col_offset;
    size_t B_row_offset, B_col_offset;
    size_t C_row_offset, C_col_offset;
    const SquareMatrix *A, *B;
    SquareMatrix *C;

    /*
    * description: creates a struct instance of PARAM
    * return: none, creates PARAM Object
    * precondition: the parameters needed are declared already
    * postcondition: creates the PARAM object
    *
    */
    PARAM() {
        dimensions = 0;
        A_row_offset = 0;
        A_col_offset = 0;
        B_row_offset = 0;
        B_col_offset = 0;
        C_row_offset = 0;
        C_col_offset = 0;
        A = nullptr;
        B = nullptr;
        C = nullptr;
    }
};

/*
* description: multiplies matrices with brute force
* return: pointer to multiplied matrix
* precondition: the two matrices to be multiplied
* postcondition: gives you the address to allocated memory holding the result
*
*/
SquareMatrix* BruteForce(const SquareMatrix& A, const SquareMatrix& B);

/*
* description: multiplies matrices with naive block recursion
* return: pointer to multiplied matrix
* precondition: the two matrices to be multiplied
* postcondition: gives you the address to allocated memory holding the result
*
*/
SquareMatrix* DivideAndConquer(const SquareMatrix& A, const SquareMatrix& B);

/*
* description: multiplies matrices with divide and conquer techniques and
* parallel programming (threading)
* return: pointer to multiplied matrix
* precondition: the two matrices to be multiplied
* postcondition: gives you the address to allocated memory holding the result
*
*/
SquareMatrix* ThreadedDivideAndConquer(const SquareMatrix& A, const SquareMatrix& B);

/*
* description: multiplies matrices by following Strassen's algorithm
* return: pointer to multiplied matrix
* precondition: the two matrices to be multiplied
* postcondition: gives you the address to allocated memory holding the result
*
*/
SquareMatrix* Strassen(const SquareMatrix& A, const SquareMatrix& B);


// helper function declarations

/*
* description: finds starting point for iteration through submatrix
* return: nothing, does work internally
* precondition: have parameters available, especially the quadrant string
* postcondition: gives row and column start variable their values
*
*/
void getQuadrantBounds(size_t dim, size_t quadrant,
    size_t& row_start, size_t& col_start);

/*
* description: multiplies the submatrix of two matrices into another matrix's
* submatrix.
* return: nothing, does work internally by reference
* precondition: matrixes are of valid sizes and quadrants are of correct value.
* postcondition: populates the sub array with the calculated values
*
*/
void SubMatrixMultiplication(const SquareMatrix& A, const SquareMatrix& B,
    SquareMatrix& C, size_t matrixA_quadrant, size_t matrixB_quadrant,
    size_t matrixC_quadrant);

/*
* description: void* function for the work done by the thread
* return: nothing, handled by the thread class
* precondition: variable p (will be cast to thread arguments in function)
* postcondition: allows the thread to do work in parallel wih others.
*
*/
void* threadWork(void* p);

/*
* description: adds two matrices together
* return: a SquareMatrix pointer pointing to the memory containing the sum of
* the matrices.
* precondition: valid matrices A and B
* postcondition: allocated memory containing the sum of two matrices.
*
*/
SquareMatrix* addMatrix(const SquareMatrix& A, const SquareMatrix& B);

/*
* description: subtracts two matrices
* return: a SquareMatrix pointer pointing to the memory containing the result
* of the matrices.
* precondition: valid matrices A and B
* postcondition: allocated memory containing the result of two matrices.
*
*/
SquareMatrix* subtractMatrix(const SquareMatrix& A, const SquareMatrix& B);

/*
* description: writes data from bigger matrix to a matrix 1/4 of its size
* return: nothing, everything is internal
* precondition: correctly sized matrices with the correct quadrant string
* postcondition: copies data from one SquareMatrix object to another smaller
* one.
*
*/
void populateToSubMatrix(const SquareMatrix& whole_matrix,
    SquareMatrix& sub_matrix, size_t quadrant);

/*
* description: writes data from a smaller matrix to another 4x its size.
* return: nothing, does everything internally
* precondition: correctly sized matrices with the correct quadrant string
* postcondition: copies data from one SquareMatrix object to another larger
* one.
*
*/
void writeFromSubMatrix(SquareMatrix& whole_matrix,
    const SquareMatrix& sub_matrix, string quadrant);

/*
* description: configures the PARAM object being used during the recursion.
* return: nothing, but plays the important role of setting the correct values
* in the PARAM object.
* precondition: correct and valid arguments and a functioning PARAM object
* postcondition: sets the PARAM object up for success during its recursive
* calls.
*
*/
void recursionBounds(PARAM& param, size_t a_quadrant, size_t b_quadrant,
    size_t c_quadrant, size_t A_row_init, size_t A_col_init,
    size_t B_row_init, size_t B_col_init,
    size_t C_row_init, size_t C_col_init, size_t init_dim);

/*
* description: implements the recursion behavior for the divide and conquer
* function.
* return: nothing, does all its work internally and recursively.
* precondition: have a configured PARAM object
* postcondition: calculates the correct values for the data of the result
* matrix.
*
*/
void divideAndConquerHelper(PARAM& params);

#endif //SQUAREMATRIX_H
