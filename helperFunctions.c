#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "functions.h"

void conjugateTranspose(int n, int m, double complex A[n][m], double complex AH[m][n]){
    int i,j; 
    for (i=0; i<n; i++){
        for (j=0; j<m; j++){
            AH[j][i] = conj(A[i][j]); // get conjugate of values and then allocate as the transpose to AH
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

void displayMatrix(int n, int m, double complex A[n][m]){
    int i, j; 
    for (i=0; i<n; i++){
      for (j=0; j<m; j++){
          printf("%.2lf %+.2lfi ", creal(A[i][j]), cimag(A[i][j])); 
      }
      printf("\n");
    }
}

void fillIdentityN(int n, double complex id[n][n]){ 
  int i,j; 

  for (i=0; i<n; i++){
    for (j=0; j<n; j++){
      
      if (i == j){
        id[i][j] = 1.0; // sets diagonal as 1
      }
  
      else{
        id[i][j] = 0.0; // sets non-diagonal as 0
      }
        	}
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
            // method is too intensive with larger matrices
            // define sign, save the element we are at, call determinant function again
            sign = pow(-1, c);
            determinant = determinant + (sign * matrix[0][c]*GetDeterminant( N-1, intmatrix ));
        }
    }

    return determinant;
}

int GetInverse(int N, double complex A[N][N], double complex inverse[N][N])
{
    int i, j, sign;
    
    // initialize for cofactor function
    double complex determinant, det, cofM[N][N];    // determinant, and cofactor matrix (det is the determinant in the cofactor loop)

        // get the determinant of the matrix
        determinant = GetDeterminant(N, A);
        //printf("Determinant= %.2f %+.2fi\n", creal(determinant),cimag(determinant));

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
        //printf("The cofactor matrix:\n");
        //displayMatrix(N,cofM);

    // now that we have the cofactor matrix, the transpose = Adjugate
    // adj will hold the adjugate matrix
    double complex adj[N][N];
    GetTranspose(N, cofM, adj);

        //printf("The adjugate:\n");
        //displayMatrix(N,adj);

    // inverse = Adjugate/determinant
    // divide the determinant from each element
    for (i = 0; i < N; i++) 
    {
        for (j = 0; j < N; j++) 
        {
            inverse[i][j] = adj[i][j] / determinant;
        }
    }

    // print
    //printf("The Inverse: \n");
    //displayMatrix(N,N,inverse);

    return 0;
}

void GetTranspose(int N, double complex matrix[N][N], double complex matrixT[N][N])
{
    int row, col;

    for (row = 0; row < N; row++)
    {
        for (col = 0; col < N; col++)
        {
            matrixT[col][row] = matrix[row][col]; //gets the transpose
        }
    }

}

void matrix_addition(int n, int m, double complex A[n][m], double complex B[n][m], double complex APlusB[n][m]){
    int i; int j; 

    for(i=0; i < n; ++i){
        for(j=0; j < m; ++j){
            APlusB[i][j] = A[i][j] + B[i][j]; // adds matrix elements
        }
    }
}

// multiplies two matrices of any compatible size
void matrix_multiplication(int a, int b, double complex matrix1[a][b], int m, int n, double complex matrix2[m][n], double complex matAXmatB[a][n]){
    int i; int j; int k;
    // error messages used to find source of segfault ****

    if(b != m){ //number of columns of matrix1 is not the same as number of rows of matrix2
        printf("Cannot do the multiplication. \n");
    } 
    else{   

        for(i = 0; i < a; i++){
            for(j = 0; j < n; j++){
                matAXmatB[i][j] = 0.0 + 0.0*I;
            }
        }  

        for(i = 0; i < a; i++){
            for(j = 0; j < n; j++){
                for(k = 0; k < m; k++){
                    matAXmatB[i][j] += matrix1[i][k]*matrix2[k][j];
                }
            }
        }
    }
}

double my_cabs(double complex x){
    return sqrt(creal(x)*creal(x) + cimag(x)*cimag(x)); // complex absolute value
}

// takes the complex norm of a vector
double norm (int n, double complex vec[n][1]){
    double term, sum, sol;  
    int i;
    // error code for an empty vector
    if (n == 0){
        printf("Error: Empty Vector. Retry.\n");
        return 1;
    }

    sum = 0;              // initialize the sum
    for (i = 0; i < n; ++i)
    {
        term = pow(creal(vec[i][0]),2) + pow(cimag(vec[i][0]),2);
        sum = sum + term;
    }
    // after the loop, square root the sum
    sol = sqrt(sum);
    return sol;
}

void scalarByMatrixMultiplication(double complex scalar, int n, int m, double complex A[n][m], double complex scalarXA[n][m]){
    int i; int j; 

    for(i=0; i < n; ++i){
        for(j=0; j < m; ++j){
            scalarXA[i][j] = scalar*A[i][j]; // multiplies matrix entries by a scalar
        }
    }
}

// sets a vector with normalized values (1/sqrt(dim))
void setNormalVec(int n, double complex v[n][1]){
    int i; 
    for (i=0; i<n; i++){
        v[i][0] = 1/sqrt(n);
    }
}

// sets the T matix for Lanczos, puts alphas on the diagonal and betas one above and below the diagonal
void setT(int n, double complex alpha[1][n], double complex beta[1][n], double complex T[n][n]){
    int i, j; 
    
    for (i=0; i<n; i++){
        for (j=0; j<n; j++){
            if (i == j){
                T[i][j] = alpha[0][j];
            }
            else if ((i-1) == j){
                T[i][j] = beta[0][j];
            }
            else if ((i+1) == j){
                T[i][j] = beta[0][j-1];
            }
            else{T[i][j] = 0;}
        }
    }
}
