#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

void display_matrix(int nRows, int nCols, double complex matrix[nRows][nCols]){
    int row; int column; 
    for(row=0; row < nRows; ++row){
        for(column=0; column < nCols; ++column)
            printf("%lf %+lfi ", creal(matrix[row][column]), cimag(matrix[row][column])); 
        printf("\n");
    }
}  

void matrix_multiplication(int n, double complex matrix1[n][n], double complex matrix2[n],double complex result[n]){
    int row; int column; 
    for(row = 0; row < n; row++){
        result[row] = 0;
        for(column = 0; column < n; column++){
            result[row] += matrix1[row][column]*matrix2[column]; 
        }
    }
}

void power_iteration(int n, double complex matrixA[n][n], double complex guess_eigenvector[n], double complex *eigenvalue){
    int max_iterations = 100;
    int j;
    while (j < max_iterations){
        double complex new_eigenvector[n];

        //Multiply matrix A and guess eigenvector
        matrix_multiplication(n,matrixA,guess_eigenvector,new_eigenvector); //new_eigenvector is matrixA*guessvectorb
        
        //Find the maximum element of new_eigenvector 
        double max_val = cabs(new_eigenvector[0]); //cabs computes the complex absolute value
        int i;
        for (i = 1; i < n; i++){
            if(cabs(new_eigenvector[i]) > max_val){
                max_val = cabs(new_eigenvector[i]); 
            } 
        } 
        // Normalize the new_eigenvector by dividing it by the maximum element 
        for (i = 0; i < n; i++) {
                new_eigenvector[i] = new_eigenvector[i] / max_val; 
            }  

        *eigenvalue = max_val; 
        
        //Keep looping until eigvenvector of nth interation is equal to eigenvector of (n-1)th iteration
        double complex check = 0.0 + 0.0*I;
        for (i = 0; i < n; i++) {
            check += cabs(new_eigenvector[i] - guess_eigenvector[i]);
        }

        // Update a new eigenvector as guess_eigenvector for the next iteration
        for (i = 0; i < n; i++) {
            guess_eigenvector[i] = new_eigenvector[i];
        }
        j++;

        //break out of the loop if they are equal 
        if (cabs(check) < 1e-6) {
            break; 
        }
        if (j >= max_iterations){
            printf("Reached maximum iterations. CANNOT FIND REAL EIGENVECTOR AND EIGENVALUE\n");
        }
    } 
    printf("Number of iterations: %d \n",j);
}

int main() {
    double complex A[2][2] = {
        {0.0  + 0.0*I,  1.0 + 0.0*I},
        {-1.0 + 0.0*I, 0.0 + 0.0*I}};
    double complex guess_eigenvector[2] = {1.0 + 0.0*I,1.0 + 0.0*I};
    double complex eigenvalue = 0.0 + 0.0*I; 
    int iterations = 0;
    int n=2; 
    int i;
    printf("Matrix A:\n");
    display_matrix(2,2,A); 

    power_iteration(2, A, guess_eigenvector, &eigenvalue);

    printf("Dominant Eigenvalue: %lf %+lfi \n",creal(eigenvalue),cimag(eigenvalue));
    printf("Eigenvector: \n");
    for (i = 0; i < n; i++){
        printf("%lf %+lfi \n", creal(guess_eigenvector[i]), cimag(guess_eigenvector[i]));
    }
    return 0; 
}
