#include <stdio.h>
#include <math.h>

void display_matrix(int nRows, int nCols, double matrix[nRows][nCols]){
    int row; int column; 
    for(row=0; row < nRows; ++row){
        for(column=0; column < nCols; ++column)
            printf("%2f ",matrix[row][column]); 
        printf("\n");
    }
}  

void matrix_multiplication(int n, double matrix1[n][n], double matrix2[n][1],double result[n][1]){
    int row; int column; int k;
    for(row = 0; row < n; row++){
        for(column = 0; column < 1; column++){
            result[row][column]=0;
        }
    }

    for(row = 0; row < n; row++){
        for(column = 0; column < 1; column++){
            for(k = 0; k < n; k++){
                result[row][column] += matrix1[row][k]*matrix2[k][column];
            }
        }
    }
    
    printf("The Product of two matrix: \n"); 
    display_matrix(n,1,result); 
}

int norm (int n, double vec[]){ 
    int i;
    double sum = 0.0;              // initialize the sum

    for (i = 0; i < n; ++i){
        sum += vec[i]*vec[i];
    }
    // after the loop, square root the sum
    return sqrt(sum);
}

void scalarToMatrixMultiplication(int n, int m, double scalar, double A[n][m], double scalar_result[n][m]){
    int i; int j; 
    for(i=0; i < n; ++i){
        for(j=0; j < m; ++j){
            scalar_result[i][j] = A[i][j]*scalar; 
        }
    }
    printf("b_{k+1}: \n");
    display_matrix(n,m,scalar_result);
}

int areMatricesEqual(int n, int m, double matrix1[n][m], double matrix2[n][m]) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            if (fabs(matrix1[i][j] - matrix2[i][j]) > 1e-6) {
                return 0; // Matrices are not equal
            }
        }
    }
    return 1; // Matrices are equal
}

int main(void){
    double matrixA[2][2] = 
    {
        {0.0,  1.0},
        {1.0,  1.0}
    };
    double matrixb[2][1] = 
    {
        {1.0},
        {1.0}
    };
    printf("Matrix A:\n");
    display_matrix(2,2,matrixA); 

    printf("Matrix b:\n");
    display_matrix(2,1,matrixb); 
    int n=2; int m=1;
    double eigenvalue[n][1]; //e A*b
    double result_multi[n][1];
    double prev_eigenvalue[n][1]; 
    double eigenvalue_norm; 
    
    
    double result_scalar[n][1]; //eigenvector 
    double prev_result_scalar[n][1];
    //scalarToMatrixMultiplication(n,m,(1/result_norm),result,scalar_result);
    
    int max_iterations = 100; 
    int iterations = 0; 
    while(iterations < max_iterations){
        matrix_multiplication(n, matrixA, matrixb, eigenvalue);
        eigenvalue_norm = norm(n, eigenvalue);

        printf("Iteration %d - the norm of A*b %f\n", iterations + 1, eigenvalue_norm);
        scalarToMatrixMultiplication(n, m, (1 / eigenvalue_norm), eigenvalue, result_multi);

        if (areMatricesEqual(n,1,eigenvalue,prev_eigenvalue)){
            printf("Converged after %d iterations.\n", iterations);
            printf("Dominant Eigenvalue: %f\n", 1 / eigenvalue_norm);
            printf("Eigenvector:\n");
            display_matrix(n, 1, result_multi);
            break;
        }

        // Update ||A*b|| for the next iteration
        for (int i = 0; i < n; i++) {
            prev_eigenvalue[i][0] = eigenvalue[i][0]; 
        }
        // Update matrixb for the next iteration
        for (int i = 0; i < n; i++){
            matrixb[i][0] = result_multi[i][0];
        }

        iterations++;
        
    }

    if (iterations >= max_iterations) {
        printf("Reached maximum iterations without convergence.\n");
    }
    return 0; 
}

