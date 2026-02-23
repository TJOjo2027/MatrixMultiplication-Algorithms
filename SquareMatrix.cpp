/*
* Author: Oluwatosin Ojo
* Assignment Title: Program 2 - SquareMatrix
* Assignment Description: This assignment explores the different ways to
*       computer matrix multiplication using different algorithm types.
* Due Date: 2/22/2026
* Date Created: 2/10/2026
* Date Last Modified: 2/22/2026
*/

#include <iostream>
#include "SquareMatrix.h"

void getQuadrantBounds(size_t dim, string quadrant,
    size_t& row_start, size_t& col_start) {

    if (quadrant == "11") {
        // row: row: 0 -> (dim/2) - 1
        row_start = 0;
        // column: 0 -> (dim/2) - 1
        col_start = 0;
    }
    else if (quadrant == "12") {
        // row: 0 -> (dim/2) - 1
        row_start = 0;
        // column: (dim/2) -> dim -1
        col_start = dim/2;
    }
    else if (quadrant == "21") {
        // row: (dim/2) -> dim
        row_start = dim/2;
        // column: 0 -> (dim/2) - 1
        col_start = 0;
    }
    else if (quadrant == "22") {
        // row: (dim/2) -> dim
        row_start = dim/2;
        // column: (dim/2) -> dim - 1
        col_start = dim/2;
    }
    else {
        // for debug
        cout << "Wrong quadrant input" << endl;
    }
}

void SubMatrixMultiplication(const SquareMatrix& A, const SquareMatrix& B,
    SquareMatrix& C, string matrixA_quadrant, string matrixB_quadrant,
    string matrixC_quadrant) {

    size_t matrix_dim = C.dim;

    // find the quadrants dimensions of the matrices assuming this order:

    /*
     * |-----------|
     * | A12 | A12 |
     * |-----|-----|
     * | A21 | A22 |
     * |_____|_____|
     */

    PARAM multParam;
    multParam.dimensions = matrix_dim;
    multParam.A = &A;
    multParam.B = &B;
    multParam.C = &C;

    // get Matrix A matrix offsets
    getQuadrantBounds(matrix_dim, matrixA_quadrant,
        multParam.A_row_offset,
        multParam.A_col_offset);

    // get Matrix B matrix offsets
    getQuadrantBounds(matrix_dim, matrixB_quadrant,
        multParam.B_row_offset,
        multParam.B_col_offset);

    // get Matrix C matrix offsets
    getQuadrantBounds(matrix_dim, matrixC_quadrant,
        multParam.C_row_offset,
        multParam.C_col_offset);

    // loops in the submatrix
    for (size_t row = 0; row < matrix_dim / 2; row++) {
        for (size_t col = 0; col < matrix_dim / 2; col++) {
            int value = 0;
            for (size_t k = 0; k < matrix_dim / 2; k++) {
                value += A.data[row + multParam.A_row_offset][k + multParam.A_col_offset] *
                    B.data[k + multParam.B_row_offset][col + multParam.B_col_offset];
            }
            C.data[row + multParam.C_row_offset][col + multParam.C_col_offset] += value;
        }
    }
}

void* threadWork(void* p) {
    ThreadArguments* arg = static_cast<ThreadArguments*>(p);
    // do the thread's work

    // first term
    SubMatrixMultiplication(*arg->A, *arg->B, *arg->C,
        arg->A_subMatrix_1, arg->B_subMatrix_1,
        arg->C_subMatrix);

    // second term
    SubMatrixMultiplication(*arg->A, *arg->B, *arg->C,
        arg->A_subMatrix_2, arg->B_subMatrix_2,
        arg->C_subMatrix);
    return nullptr;
}

SquareMatrix* addMatrix(const SquareMatrix& A, const SquareMatrix& B) {
    // iterates throught the matrices and adds
    SquareMatrix* result = new SquareMatrix(A.dim);
    for (size_t r = 0; r < result -> dim; r++) {
        for (size_t c = 0; c < result -> dim; c++) {
            result->data[r][c] = A.data[r][c] + B.data[r][c];
        }
    }
    return result;
}

SquareMatrix* subtractMatrix(const SquareMatrix& A, const SquareMatrix& B) {
    // iterates through the matrices and performs subtraction of B from A
    SquareMatrix* result = new SquareMatrix(A.dim);
    for (size_t r = 0; r < result -> dim; r++) {
        for (size_t c = 0; c < result -> dim; c++) {
            result->data[r][c] = A.data[r][c] - B.data[r][c];
        }
    }
    return result;
}

void populateToSubMatrix(const SquareMatrix& whole_matrix,
        SquareMatrix& sub_matrix, string quadrant) {
    size_t row_start, col_start;
    getQuadrantBounds(whole_matrix.dim, quadrant, row_start, col_start);

    // fill sub_matrix with the quadrant chosen form whole_matrix
    for (size_t row = 0; row < whole_matrix.dim / 2; row++) {
        for (size_t col = 0; col < whole_matrix.dim / 2; col++) {
            sub_matrix.data[row][col] =
                whole_matrix.data[row + row_start][col + col_start];
        }
    }
}

void writeFromSubMatrix(SquareMatrix& whole_matrix,
        const SquareMatrix& sub_matrix, string quadrant) {
    size_t row_start, col_start;
    getQuadrantBounds(whole_matrix.dim, quadrant, row_start, col_start);

    // fill sub_matrix with the quadrant chosen form whole_matrix
    for (size_t row = 0; row < whole_matrix.dim / 2; row++) {
        for (size_t col = 0; col < whole_matrix.dim / 2; col++) {
            whole_matrix.data[row + row_start][col + col_start] =
                sub_matrix.data[row][col];
        }
    }
}

void recursionBounds(PARAM& param, string a_quadrant, string b_quadrant,
        string c_quadrant, size_t A_row_init, size_t A_col_init,
        size_t B_row_init, size_t B_col_init,
        size_t C_row_init, size_t C_col_init, size_t init_dim) {

    // load the param function with the values from the getQuadrantBounds()
    // function
    getQuadrantBounds(init_dim, a_quadrant, param.A_row_offset, param.A_col_offset);
    getQuadrantBounds(init_dim, b_quadrant, param.B_row_offset, param.B_col_offset);
    getQuadrantBounds(init_dim, c_quadrant, param.C_row_offset, param.C_col_offset);

    // builds on the prior values of the offsets
    param.A_row_offset += A_row_init;
    param.A_col_offset += A_col_init;
    param.B_row_offset += B_row_init;
    param.B_col_offset += B_col_init;
    param.C_row_offset += C_row_init;
    param.C_col_offset += C_col_init;
}

void divideAndConquerHelper(PARAM& params) {
    // base case
    if (params.dimensions == 1) {
        params.C->data[params.C_row_offset][params.C_col_offset] +=
            params.A->data[params.A_row_offset][params.A_col_offset] *
                params.B->data[params.B_row_offset][params.B_col_offset];
        return;
    }

    // recursive case

    // split into the 4 quadrants
    size_t initial_dim = params.dimensions;
    params.dimensions /= 2;

    /*
     * |-----------|
     * | A12 | A12 |
     * |-----|-----|
     * | A21 | A22 |
     * |_____|_____|
     */

    const string A11 = "11";
    const string A12 = "12";
    const string A21 = "21";
    const string A22 = "22";

    size_t A_row_init = params.A_row_offset;
    size_t A_col_init = params.A_col_offset;
    size_t B_row_init = params.B_row_offset;
    size_t B_col_init = params.B_col_offset;
    size_t C_row_init = params.C_row_offset;
    size_t C_col_init = params.C_col_offset;

    // Matrix C - QUADRANT 1

    // section 1
    recursionBounds(params, A11, A12, A12,
        A_row_init, A_col_init,
        B_row_init, B_col_init,
        C_row_init, C_col_init,
        initial_dim);

    // recursive call
    divideAndConquerHelper(params);
    // restore param.dimensions for recursion's alterations
    params.dimensions = initial_dim / 2;

    // section 2
    recursionBounds(params, A12, A22, A12,
        A_row_init, A_col_init,
        B_row_init, B_col_init,
        C_row_init, C_col_init,
        initial_dim);

    // recursive call
    divideAndConquerHelper(params);
    // restore param.dimensions for recursion's alterations
    params.dimensions = initial_dim / 2;

    // Matrix C - QUADRANT 2

    // section 1
    recursionBounds(params, A11, A11, A11,
        A_row_init, A_col_init,
        B_row_init, B_col_init,
        C_row_init, C_col_init,
        initial_dim);

    // recursive call
    divideAndConquerHelper(params);
    // restore param.dimensions for recursion's alterations
    params.dimensions = initial_dim / 2;

    // section 2
    recursionBounds(params, A12, A21, A11,
        A_row_init, A_col_init,
        B_row_init, B_col_init,
        C_row_init, C_col_init,
        initial_dim);

    // recursive call
    divideAndConquerHelper(params);
    // restore param.dimensions for recursion's alterations
    params.dimensions = initial_dim / 2;

    // Matrix C - QUADRANT 3

    // section 1
    recursionBounds(params, A21, A11, A21,
        A_row_init, A_col_init,
        B_row_init, B_col_init,
        C_row_init, C_col_init,
        initial_dim);

    // recursive call
    divideAndConquerHelper(params);
    // restore param.dimensions for recursion's alterations
    params.dimensions = initial_dim / 2;

    // section 2
    recursionBounds(params, A22, A21, A21,
        A_row_init, A_col_init,
        B_row_init, B_col_init,
        C_row_init, C_col_init,
        initial_dim);

    // recursive call
    divideAndConquerHelper(params);
    // restore param.dimensions for recursion's alterations
    params.dimensions = initial_dim / 2;

    // Matrix C - QUADRANT 4

    // section 1
    recursionBounds(params, A21, A12, A22,
        A_row_init, A_col_init,
        B_row_init, B_col_init,
        C_row_init, C_col_init,
        initial_dim);

    // recursive call
    divideAndConquerHelper(params);
    // restore param.dimensions for recursion's alterations
    params.dimensions = initial_dim / 2;

    // section 2
    recursionBounds(params, A22, A22, A22,
        A_row_init, A_col_init,
        B_row_init, B_col_init,
        C_row_init, C_col_init,
        initial_dim);

    // recursive call
    divideAndConquerHelper(params);
    // no need to restore because function is done
}

SquareMatrix* BruteForce(const SquareMatrix& A, const SquareMatrix& B) {
    SquareMatrix* C;
    C = new SquareMatrix(A.dim);
    // iterates through the matrix with a triple-nested loop
    for (size_t r = 0; r < C -> dim; r++) {
        for (size_t c = 0; c < C -> dim; c++) {
            int val = 0;
            // implementing general matrix multiplication formula
            for (size_t k = 0; k < C -> dim; k++) {
                val += A.data[r][k] * B.data[k][c];
            }
            C -> data[r][c] = val;
        }
    }
    return C;
}

SquareMatrix* DivideAndConquer(const SquareMatrix& A, const SquareMatrix& B) {
    SquareMatrix* C = new SquareMatrix(A.dim);

    // initialize PARAM object
    PARAM matrixParams;
    matrixParams.dimensions = A.dim;
    matrixParams.A = &A;
    matrixParams.B = &B;
    matrixParams.C = C;
    matrixParams.A_row_offset = 0;
    matrixParams.A_col_offset = 0;
    matrixParams.B_row_offset = 0;
    matrixParams.B_col_offset = 0;
    matrixParams.C_row_offset = 0;
    matrixParams.C_col_offset = 0;

    // recursion lies in here
    divideAndConquerHelper(matrixParams);

    return C;
}

SquareMatrix* ThreadedDivideAndConquer(const SquareMatrix& A, const SquareMatrix& B) {
    // case where the matrix has a dimension of 1
    if (A.dim == 1) {
        SquareMatrix *C = new SquareMatrix(A.dim);
        C -> data[0][0] = A.data[0][0] * B.data[0][0];
        return C;
    }

    SquareMatrix* C = new SquareMatrix(A.dim);

    // constants to be used
    const string A11 = "11";
    const string A12 = "12";
    const string A21 = "21";
    const string A22 = "22";

    // create the threads
    pthread_t t1, t2, t3, t4;

    // create the thread argument objects
    ThreadArguments arg1(&A, &B, C,
        A11, A11,
        A12, A21,
        A11);

    ThreadArguments arg2(&A, &B, C,
        A11, A12,
        A12, A22,
        A12);

    ThreadArguments arg3(&A, &B, C,
        A21, A11,
        A22, A21,
        A21);

    ThreadArguments arg4(&A, &B, C,
        A21, A12,
        A22, A22,
        A22);

    // SUBMATRIX - C11
    pthread_create(&t1, nullptr, threadWork, &arg1);

    // SUBMATRIX - C12
    pthread_create(&t2, nullptr, threadWork, &arg2);

    // SUBMATRIX - C21
    pthread_create(&t3, nullptr, threadWork, &arg3);

    // SUBMATRIX - C22
    pthread_create(&t4, nullptr, threadWork, &arg4);

    // join the threads
    pthread_join(t1, nullptr);
    pthread_join(t2, nullptr);
    pthread_join(t3, nullptr);
    pthread_join(t4, nullptr);

    return C;
}

SquareMatrix* Strassen(const SquareMatrix& A, const SquareMatrix& B) {
    // case where the matrix has a dimension of 1
    if (A.dim == 1) {
        SquareMatrix *C = new SquareMatrix(A.dim);
        C -> data[0][0] = A.data[0][0] * B.data[0][0];
        return C;
    }

    // S Matrices creation
    SquareMatrix* S1;
    SquareMatrix* S2;
    SquareMatrix* S3;
    SquareMatrix* S4;
    SquareMatrix* S5;
    SquareMatrix* S6;
    SquareMatrix* S7;
    SquareMatrix* S8;
    SquareMatrix* S9;
    SquareMatrix* S10;

    // P Matrices creation
    SquareMatrix* P1;
    SquareMatrix* P2;
    SquareMatrix* P3;
    SquareMatrix* P4;
    SquareMatrix* P5;
    SquareMatrix* P6;
    SquareMatrix* P7;

    // QUADRANT STRING CONSTANTS
    const string Q11 = "11";
    const string Q12 = "12";
    const string Q21 = "21";
    const string Q22 = "22";

    // A Submatrices
    SquareMatrix* A11 = new SquareMatrix(A.dim / 2);
    populateToSubMatrix(A, *A11, Q11);
    SquareMatrix* A12 = new SquareMatrix(A.dim / 2);
    populateToSubMatrix(A, *A12, Q12);
    SquareMatrix* A21 = new SquareMatrix(A.dim / 2);
    populateToSubMatrix(A, *A21, Q21);
    SquareMatrix* A22 = new SquareMatrix(A.dim / 2);
    populateToSubMatrix(A, *A22, Q22);


    // B Submatrices
    SquareMatrix* B11 = new SquareMatrix(A.dim / 2);
    populateToSubMatrix(B, *B11, Q11);
    SquareMatrix* B12 = new SquareMatrix(A.dim / 2);
    populateToSubMatrix(B, *B12, Q12);
    SquareMatrix* B21 = new SquareMatrix(A.dim / 2);
    populateToSubMatrix(B, *B21, Q21);
    SquareMatrix* B22 = new SquareMatrix(A.dim / 2);
    populateToSubMatrix(B, *B22, Q22);

    // creating C matrix
    SquareMatrix* C = new SquareMatrix(A.dim);

    // S1
    S1 = subtractMatrix(*B12, *B22);

    // S2
    S2 = addMatrix(*A11, *A12);

    // S3
    S3 = addMatrix(*A21, *A22);

    // S4
    S4 = subtractMatrix(*B21, *B11);

    // S522222222222222222222222
    S5 = addMatrix(*A11, *A22);

    // S6
    S6 = addMatrix(*B11, *B22);

    // S7
    S7 = subtractMatrix(*A12, *A22);

    // S8
    S8 = addMatrix(*B21, *B22);

    // S9
    S9 = subtractMatrix(*A11, *A21);

    // S10
    S10 = addMatrix(*B11, *B12);

    // P1
    P1 = BruteForce(*A11, *S1);

    // P2
    P2 = BruteForce(*S2, *B22);

    // P3
    P3 = BruteForce(*S3, *B11);

    // P4
    P4 = BruteForce(*A22, *S4);

    // P5
    P5 = BruteForce(*S5, *S6);

    // P6
    P6 = BruteForce(*S7, *S8);

    // P7
    P7 = BruteForce(*S9, *S10);

    // populating C submatrices

    SquareMatrix* tempMatrix1;
    SquareMatrix* tempMatrix2;
    SquareMatrix* C11;
    SquareMatrix* C12;
    SquareMatrix* C21;
    SquareMatrix* C22;

    // C11
    tempMatrix1 = addMatrix(*P5, *P4);
    tempMatrix2 = subtractMatrix(*tempMatrix1, *P2);

    // hold the value of P5 + P4 - P2 + P6
    C11 = addMatrix(*tempMatrix2, *P6);
    writeFromSubMatrix(*C, *C11, Q11);
    // deleting to prevent dangling pointers
    delete tempMatrix1;
    delete tempMatrix2;
    delete C11;

    // C12
    C12 = addMatrix(*P1, *P2);
    writeFromSubMatrix(*C, *C12, Q12);

    delete C12;

    // C21
    C21 = addMatrix(*P3, *P4);
    writeFromSubMatrix(*C, *C21, Q21);

    delete C21;

    // C22
    tempMatrix1 = addMatrix(*P5, *P1);
    tempMatrix2 = subtractMatrix(*tempMatrix1, *P3);
    C22 = subtractMatrix(*tempMatrix2, *P7);
    // tempMatrix1 now holds P5 + P1 - P3 - P7
    writeFromSubMatrix(*C, *C22, Q22);

    delete tempMatrix1;
    delete tempMatrix2;
    delete C22;

    // delete allocated memory from function calls + allocations
    delete A11; delete A12; delete A21; delete A22;
    delete B11; delete B12; delete B21; delete B22;
    delete S1; delete S2; delete S3; delete S4; delete S5;
    delete S6; delete S7; delete S8; delete S9; delete S10;
    delete P1; delete P2; delete P3; delete P4; delete P5;
    delete P6; delete P7;

    return C;
}