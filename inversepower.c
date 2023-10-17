#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

void displayMatrix(int n, int m, double complex A[n][m]){
    int i, j; 
    for (i=0; i<n; i++){
      for (j=0; j<m; j++){
            //printf("%f ", (float) A[i][j]);
          printf("%.2lf %+.2lfi ", creal(A[i][j]), cimag(A[i][j])); 
      }
      printf("\n");
    }
}

void matrix_multiplication(int a, int b, double complex matrix1[a][b], int m, int n, double complex matrix2[m][n], double complex matAXmatB[a][n]){
    int i; int j; int k;
    fflush(stdout);
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
    return sqrt(creal(x)*creal(x) + cimag(x)*cimag(x));
}

void conjugateTranspose(int n, int m, double complex A[n][m], double complex AH[m][n]){
    int i,j; 
    for (i=0; i<n; i++){
        for (j=0; j<m; j++){
            AH[j][i] = conj(A[i][j]); // get conjugate of values and then allocate as the transpose to AH
	    }
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

void scalarByMatrixMultiplication(double complex scalar, int n, int m, double complex A[n][m], double complex scalarXA[n][m]){
    int i; int j; 

    for(i=0; i < n; ++i){
        for(j=0; j < m; ++j){
            scalarXA[i][j] = scalar*A[i][j]; 
        }
    }
}

void matrix_addition(int n, int m, double complex A[n][m], double complex B[n][m], double complex APlusB[n][m]){
    int i; int j; 

    for(i=0; i < n; ++i){
        for(j=0; j < m; ++j){
            APlusB[i][j] = A[i][j] + B[i][j]; 
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

int GetInverse(int N, double complex A[N][N], double complex inverse[N][N])
{

    int i, j, sign;
    

    // initialize for cofactor function
    double complex determinant, det, cofM[N][N];    // determinant, and cofactor matrix (det is the determinant in the cofactor loop)

        // get the determinant of the matrix
        determinant = GetDeterminant(N, A);
       // printf("Determinant= %.2f %+.2fi\n", creal(determinant),cimag(determinant));

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
   // printf("The Inverse: \n");
   //displayMatrix(N,inverse);

    return 0;
}

void setNormalVec(int n, double complex v[n][1]){
    int i; 
    for (i=0; i<n; i++){
        v[i][0] = 1/sqrt(n);
    }
}

void vectorTranspose(int N, double complex matrix[N][1], double complex matrixT[1][N])
{
    int row;

    for (row = 0; row < N; row++)
    {
       
        matrixT[0][row] = matrix[row][0];
    }

}

void inverseIteration(int n, double complex matrixA[n][n]){
    double complex Id[n][n];
    fillIdentityN(n, Id);
    double complex muId[n][n];
    double complex Ainv[n][n];
    int max_iterations = 100; int j = 0;
    //double complex bH[n][1]; 
    double complex Ab[n][1];
    double complex eig[1][1];
    double complex AminMuI[n][n];
    double complex newEigVecT[1][n];
    // normalized vector as starting guess
    double complex guess_eigenvector[n][1];
    setNormalVec(n, guess_eigenvector);
    double complex eigenvalue0 = -1.5 + 0.0*I;
    double complex eigenvalue = eigenvalue0;

    // create identity function
    // - mew * I
    scalarByMatrixMultiplication(-1*(eigenvalue), n, n, Id, muId);
    // A + (- mew * I)
   matrix_addition(n, n, matrixA, muId, AminMuI);
    GetInverse(n,AminMuI,Ainv);

    while (j < max_iterations){
        double complex new_eigenvector[n][1];

        //Multiply matrix A and guess eigenvector
        matrix_multiplication(n,n, Ainv,n,1,guess_eigenvector,new_eigenvector); //new_eigenvector is matrixA*guessvectorb
        
        //Find the maximum element of new_eigenvector 
        double max_val = my_cabs(new_eigenvector[0][0]); //cabs computes the complex absolute value
        int i;
        for (i = 1; i < n; i++){
            if(my_cabs(new_eigenvector[i][0]) > max_val){
                max_val = my_cabs(new_eigenvector[i][0]); 
            } 
        } 
        // Normalize the new_eigenvector by dividing it by the maximum element 
        for (i = 0; i < n; i++) {
                new_eigenvector[i][0] = new_eigenvector[i][0] / max_val; 
            }  
        printf("Eigenvalue: %f %+f*i \n",creal(eigenvalue),cimag(eigenvalue));
        
        conjugateTranspose(n,1,new_eigenvector,newEigVecT);
        matrix_multiplication(n,n, matrixA, n,1,new_eigenvector,Ab);
        matrix_multiplication(1,n,newEigVecT,n,1,Ab,eig);
        
        eigenvalue = eig[0][0]; 
        
        //Keep looping until eigvenvector of nth interation is equal to eigenvector of (n-1)th iteration
        double check = 0.0;
        for (i = 0; i < n; i++) {
            check += my_cabs(new_eigenvector[i][0] - guess_eigenvector[i][0]);
        }

        // Update a new eigenvector as guess_eigenvector for the next iteration
        for (i = 0; i < n; i++) {
            guess_eigenvector[i][0] = new_eigenvector[i][0];
        }
        j++;


        //break out of the loop if they are equal 
        if (my_cabs(check) < 1e-6) {
            printf("Dominant Eigenvalue: %f %+f*i \n",creal(eigenvalue),cimag(eigenvalue));
            break; 
        }
        if (j >= max_iterations){
            printf("Reached maximum iterations. CANNOT FIND REAL EIGENVALUE\n");
            printf("Ignore the rest!\n"); 
        }
    } 
    printf("Number of iterations: %d \n",j);
}



int main() {
    double complex A[3][3] = {{1, 2, 3+I}, {2, 2, 3}, {3-I, 3, 3}};
    //double complex guess_eigenvector[2][1] = {{1.0 + 0.0*I},{1.0 + 0.0*I}};
    //double complex eigenvalue = 1.0 + 0.0*I; 
    int n=3; 
    printf("Matrix A:\n");
    displayMatrix(n,n,A);

    inverseIteration(n, A);

    //printf("Dominant Eigenvalue: %f %+f*i \n",creal(eigenvalue),cimag(eigenvalue));
    //printf("Eigenvector: \n");
   // for (i = 0; i < n; i++){
       // printf("%f %+fi\n", creal(guess_eigenvector[i][0]), cimag(guess_eigenvector[i][0]));
    //}
   return 0;
}
