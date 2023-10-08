#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void display_matrix(int nRows, int nCols, double matrix[nRows][nCols]){
    int row; int column; 
    for(row=0; row < nRows; ++row){
        for(column=0; column < nCols; ++column)
            printf("%2f ",matrix[row][column]); 
        printf("\n");
    }
}  

void matrix_multiplication(int n, double matrix1[n][n], double matrix2[n],double result[n]){
    int row; int column; 
    for(row = 0; row < n; row++){
        result[row] = 0;
        for(column = 0; column < n; column++){
            result[row] += matrix1[row][column]*matrix2[column]; 
        }
    }
}

void power_iteration(int n, double matrixA[n][n], double guess_eigenvector[n], double *eigenvalue, int *iterations){
    double prev_eigenvalue = 0.0;
    while (1){ //run continually until break
        double new_eigenvalue = 0.0; 
        double new_eigenvector[n];
        double prev_eigenvector[n];
        //Multiply matrix A and guess eigenvector
        matrix_multiplication(n,matrixA,guess_eigenvector,new_eigenvector); //new_eigenvector is matrixA*guessvectorb
        
        //Find the maximum element of new_eigenvector 
        double max_val = new_eigenvector[0];
        for (int i = 1; i < n; i++){
            if(new_eigenvector[i] > max_val){
                max_val = new_eigenvector[i]; 
            } 
        } 
        // Normalize the new_eigenvector by dividing it by the maximum element 
        for (int i = 0; i < n; i++) {
                new_eigenvector[i] = new_eigenvector[i] / max_val; 
            }  

        *eigenvalue = max_val; 
        
        //Keep looping until eigvenvector of nth interation is equal to eigenvector of (n-1)th iteration
        double check = 0.0;
        for (int i = 0; i < n; i++) {
            check += fabs(new_eigenvector[i] - guess_eigenvector[i]);
        }

        // Update a new eigenvector as guess_eigenvector for the next iteration
        for (int i = 0; i < n; i++) {
            guess_eigenvector[i] = new_eigenvector[i];
        }
        (*iterations)++; //count the number of iterations needed
        
        //break out of the loop if they are equal 
        if (check < 1e-6) {
            break; 
        }
    } 
}

int main() {
    double A[2][2] = {
        {0.0,  1.0},
        {1.0,  1.0}};
    double guess_eigenvector[2] = {1.0,1.0};
    double eigenvalue = 0.0; 
    int iterations = 0;
    int n=2; 
    
    printf("Matrix A:\n");
    display_matrix(2,2,A); 

    power_iteration(2, A, guess_eigenvector, &eigenvalue,&iterations);

    printf("Dominant Eigenvalue: %lf \n",eigenvalue);
    printf("Eigenvector: \n");
    for (int i = 0; i < n; i++){
        printf("%lf\n", guess_eigenvector[i]);
    }
    printf("Number of iterations: %d \n",iterations);
    return 0; 
}
