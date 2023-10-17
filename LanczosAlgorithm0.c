#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include "functions.h"

/*int main(void){
    int iter = 2;
    int dim = 3;
    double complex A[3][3] = {{1, 2, 3+I}, {2, 2, 3}, {3-I, 3, 3}};
    LanczosAlgorithm(dim, A, iter);
    return 0;
}*/

void LanczosAlgorithm(int d, double complex A[d][d], int itr){
   // initialize variables
    int dim = d;
    int iter = itr;
    double complex q[dim][1];
    double complex r[dim][1];
    double complex alpha[1][iter];
    double complex beta[1][iter];

    // store intermediate calculations
    double complex qh[1][dim];
    double complex alphaJ[1][1];
    double complex minusAlphaq[1][1];
    double complex betaJ = 0.0 + 0.0*I;
    double complex Aq[dim][1];
    //double complex (*Aq)[1] = malloc(sizeof(double complex[dim][1]));
    double complex betajMin1XV[dim][1];
    int i, j;
    
    // q = normalized vector
    setNormalVec(dim, q); 

    // r = aq, size: dimxdim * dimx1 = dimx1
    matrix_multiplication(dim, dim, A, dim, 1, q, r); 

    // alpha1 = q*r, size: 1xdim * dimx1 = 1 x 1
    conjugateTranspose(dim, 1, q, qh); // conjugate of q
    matrix_multiplication(dim, dim, qh, dim, 1, r, alphaJ);  // set alpha1
    alpha[0][0] = alphaJ[0][0]; // store alpha1
   
    // r = r - alpha1q
    scalarByMatrixMultiplication(-1*alphaJ[0][0], dim, 1, q, minusAlphaq); // get - alpha1q
    matrix_addition(dim, 1, r, minusAlphaq, r);
    
    // beta1 = ||r||
    betaJ = norm(dim, r);
   
   // this breaks? :(
    beta[0][0] = betaJ; // store beta1
    
    // start loop at j = 2 until betaj = 0
    for(j = 1; j<iter; j++) { // n = m by default

        // v = q
        complex double v[dim][1];
        for (i = 0; i<dim; i++){
            v[i][0] = q[i][0];
        }
        // q = r/beta[j-1]
        scalarByMatrixMultiplication(1/betaJ, dim, 1, r, q);
        
        // r = Aq - betaj-1v
        // breaks here :(
        printf("end");
        fflush(stdout);
        matrix_multiplication(dim, dim, A, dim, 1, q, Aq);

        scalarByMatrixMultiplication(-betaJ, dim, 1, v, betajMin1XV);
        matrix_addition(dim, 1, Aq, betajMin1XV, r);
        
        // aj = q*r
        conjugateTranspose(dim, 1, q, qh);
        matrix_multiplication(1, dim, qh, dim, 1, r, alphaJ);

        // r = r - alphajq
        scalarByMatrixMultiplication(-1*alphaJ[0][0], dim, 1, q, minusAlphaq); // get - alpha1q
        matrix_addition(dim, 1, r, minusAlphaq, r);
       
        // betaj = ||r||
        betaJ = norm(dim, r);

        // store alpha and beta
        alpha[0][j] = alphaJ[0][0];
        beta[0][j] = betaJ;
        
        if (betaJ == 0){break;}
    }

    // get the eigenvalue of T
    double complex T[iter][iter];
    double complex eigenValue = 0.0 + 0.0*I;
    double complex vecGuess[iter][1];
    setT(iter, alpha, beta, T);
    setNormalVec(iter, vecGuess);
    // Next step: Call power iteration to get eigenvalue of T.
    
    power_iteration(iter, T, vecGuess, &eigenValue);
    
    printf("The dominant eigenvalue of A produced by the Lanczos Algorithm is: %lf %+lf*i \n",creal(eigenValue),cimag(eigenValue));
    printf("This answer took %d iterations in the Lanczos Algorithm.\n",iter);
    
}

// reference pg. 178 https://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter10.pdf
// https://www.youtube.com/watch?v=2Y1ZDQw_2zw


//functions used
// sets v1 as arbritrary vector with norm = 1
void setNormalVec(int n, double complex v[n][1]){
    int i; 
    for (i=0; i<n; i++){
        v[i][0] = 1/sqrt(n);
    }
}

void matrix_multiplication(int a, int b, double complex matrix1[a][b], int m, int n, double complex matrix2[m][n], double complex matAXmatB[a][n]){
    int i; int j; int k;
    printf("In Multiply\n");
    printf("%d %d %d %d\n", a, b, m, n);
    fflush(stdout);
    if(b != m){ //number of columns of matrix1 is not the same as number of rows of matrix2
        printf("Cannot do the multiplication. \n");
    } 
    else{   
        printf("Setting to Zero\n");
        fflush(stdout);
        for(i = 0; i < a; i++){
            for(j = 0; j < n; j++){
                printf("%d, %d\n", i, j);
                fflush(stdout);
                matAXmatB[i][j] = 0.0 + 0.0*I;
            }
        }  
        printf("Multiplying\n");
        fflush(stdout);
        for(i = 0; i < a; i++){
            for(j = 0; j < n; j++){
                for(k = 0; k < m; k++){
                    printf("%d, %d, %d\n", i, j, k);
                    fflush(stdout);
                    matAXmatB[i][j] += matrix1[i][k]*matrix2[k][j];
                }
            }
        }
    }
}

void conjugateTranspose(int n, int m, double complex A[n][m], double complex AH[m][n]){
    int i,j; 
    for (i=0; i<n; i++){
        for (j=0; j<m; j++){
            AH[j][i] = conj(A[i][j]); // get conjugate of values and then allocate as the transpose to AH
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

// check size
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

double my_cabs(double complex x){
    return sqrt(creal(x)*creal(x) + cimag(x)*cimag(x));
}

void power_iteration(int n, double complex matrixA[n][n], double complex guess_eigenvector[n][1], double complex *eigenvalue){
    int max_iterations = 100; int j = 0;
    while (j < max_iterations){
        double complex new_eigenvector[n][1];

        //Multiply matrix A and guess eigenvector
        matrix_multiplication(n, n, matrixA, n, 1, guess_eigenvector,new_eigenvector); //new_eigenvector is matrixA*guessvectorb
        
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

        *eigenvalue = max_val; 
        
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
            break; 
        }
        if (j >= max_iterations){
            printf("Reached maximum iterations. CANNOT FIND REAL EIGENVECTOR AND EIGENVALUE\n");
            printf("Ignore the rest!\n"); 
        }
    } 
    printf("Number of iterations inside the Power Iteration: %d \n",j);
}
    }
}
