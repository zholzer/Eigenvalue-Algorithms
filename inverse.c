#include <stdio.h>
#include <complex.h>
#include <math.h>
// function to find the inverse of a matrix
// functions included:
// ---> displayMatrix (displays a given matrix in the terminal)
// ---> GetCofactor (finds the cofactor of a matrix at a given element)
// ---> det2by2 (finds the determinant of a 2x2 matrix)
// ---> GetTranspose (computes the transpose of a given matrix)
// ---> Get Determinant (computes the determinant for any NxN matrix- 2x2 determinant function is included in this)
// 
// there will also be a main script
// process in the main script for finding the inverse of a matrix:
// step 1. determinant of entire matrix (if determinant = 0, there is no inverse)
// step 2. cofactor matrix 
//      // find the cofactor of each element from original matrix
//      // take the determinant of each cofactor matrix and multiply by the element we are at
//      // include changing +/- signs
//      // this new value becomes the associated element of the "cofactor matrix"
// step 3. adjoint = (cofactor matrix)T (transpose)
// step 4. inverse = adjoint/determinant

// "display matrix function"
// takes in a declared matrix A of size n x n and size n
void displayMatrix(int n, double complex A[n][n]){
    int i, j; 
    for (i=0; i<n; i++){
      for (j=0; j<n; j++){
            
            printf("%.2f %+.2fi ", creal(A[i][j]),cimag(A[i][j])); 
        
      }
      printf("\n");
    }
}  

// cofactor function (only works for higher than 2x2 matrices, cofactor for 2x2 is written directly in main )
// GetCofactor is a function to compute the cofactor of a matrix A[p][q]
//
// inputs: int size of matrix N, given matrix A[N][N], empty matrix cof[N][N] to hold output,
//         ints m/n to specify the element we're at in A
void GetCofactor(int N, double complex A[N][N], double complex cof[N-1][N-1], int m, int n)
{
    int i, j, r, c;

    // initialize the index
    i = 0; j = 0;

    // for loop thorugh each element
    for (r = 0; r < N; r++)
    {
        for (c = 0; c < N; c++)
        {
            // identify the rows/cols that we are NOT at
            if (r != m && c != n)
            {
                cof[i][j] = A[r][c];
                j++;

                // once the row is filled, we reset for next row
                if(j == N-1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }

}

// function for the determinant of a 2x2 matrix
// inputs: n size matrix, A[n][n] matrix
double complex det2by2(int n, double complex A[n][n])
{

    // if statement if the matrix is not a 2x2
    if (n != 2)
    {
        printf("Please input a 2x2 matrix.\n");
        return 1;
    }
    
    double complex a, b, c, d, det;

    a = A[0][0]; b = A[0][1];
    c = A[1][0]; d = A[1][1];

    det = a*d - b*c;

    return det;

}


// take and NxN matrix and computes the transpose
// inputs: matrix size N, given matrix, empty matrix to hold the transpose matrix
void GetTranspose(int N, double complex matrix[N][N], double complex matrixT[N][N])
{
    int row, col;

    for (row = 0; row < N; row++)
    {
        for (col = 0; col < N; col++)
        {
            matrixT[col][row] = matrix[row][col];
        }
    }

}

// determinant function (for any NxN matrix)
// inputs: size N, matrix
// output: determinant
double complex GetDeterminant(int N, double complex matrix[N][N])
{
    int sign, c;
    double complex determinant, intmatrix[N-1][N-1];
    
    // initialize
    determinant = 0 + 0*I;

    // if the matrix is 2x2
    if (N == 2)
        {
            determinant = det2by2(N, matrix);
        }
    
    // else, break down the matrix
    else{
        // for each column in the first row:
        for (c = 0; c < N; c++){

            // call the cofactor of the row/col we are at
            GetCofactor(N, matrix, intmatrix, 0, c);

            // define sign, save the element we are at, call determinant function again
            sign = pow(-1, c);
            determinant = determinant + (sign * matrix[0][c]*GetDeterminant( N-1, intmatrix ));
        }
    }

    return determinant;
}


// main function
int main (void)
{
    int i, j, N, sign;
    
    // define some matrix 
    N = 4;
    double complex A[4][4] = { { 1, 2-3*I, -9, 4},
                    { 2+3*I, 5, 4, -12},
                    {2, 4, 5+I, 0},
                    {2, 4, 5+I, 9}
                    };

    // initialize for cofactor function
    double complex determinant, det, cofM[N][N];    // determinant, and cofactor matrix (det is the determinant in the cofactor loop)

        // get the determinant of the matrix
        determinant = GetDeterminant(N, A);
        printf("Determinant= %.2f %+.2fi\n", creal(determinant),cimag(determinant));

        // if the determinant = 0, there is no inverse!!
        if (determinant == 0) {
            printf("Determinant = 0. There is no inverse. Please input another matrix.\n");
            return 1;
        }

    // find the Adjoint- first find the cofactor matrix
    // find the cofactor for each element, find the determinant of that matrix, determinant = that element
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            // if statement if the matrix is a 2x2
            if (N == 2) { 
                double complex cof; // empty matrix for cofactor 
                int row,col;
                for (row = 0; row < N; row++){
                    for (col = 0; col < N; col++){
                        if (row != i && col != j)
                        {
                            cof = A[row][col];
                            det = cof;
                        }
                    }
                }
            }

            // for all other N sized matrices:
            else {
                double complex cof[N-1][N-1]; // empty matrix for cofactor
                // call cofactor function
                GetCofactor(N, A, cof, i, j);

                // call the determinant function 
                det = GetDeterminant(N-1, cof);
           }


            // save the determinant in the corresponding element
            // this is the cofactor matrix cofM
            sign = pow(-1, i+j);
            cofM[i][j] = det*sign;

        }
    }

    // display the cofactor matrix
    printf("The cofactor matrix:\n");
    displayMatrix(N,cofM);

    // now that we have the cofactor matrix, the transpose = Adjugate
    // adj will hold the adjugate matrix
    double complex adj[N][N];
    GetTranspose(N, cofM, adj);

    printf("The adjugate:\n");
    displayMatrix(N,adj);

    // inverse = Adjugate/determinant
    double complex inverse[N][N];

    // divide the determinant from each element
    for (i = 0; i < N; i++) 
    {
        for (j = 0; j < N; j++) 
        {
            inverse[i][j] = adj[i][j] / determinant;
        }
    }

    // print
    printf("The Inverse: \n");
    displayMatrix(N,inverse);

    return 0;
}
