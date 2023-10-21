#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include "functions.h"

// Luna's power iteration
void power_iteration(int n, double complex matrixA[n][n]){
    int max_iterations = 100; int j = 0;
    // normalized vector as starting guess
    double complex guess_eigenvector[n][1];
    setNormalVec(n, guess_eigenvector);
    double complex eigenvalue = 0.0 + 0.0*I;
    while (j < max_iterations){
        double complex new_eigenvector[n][1];

        //Multiply matrix A and guess eigenvector
        matrix_multiplication(n,n,matrixA,n,1,guess_eigenvector,new_eigenvector); //new_eigenvector is matrixA*guessvectorb
        
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

        eigenvalue = max_val; 
        
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
            printf("Dominant Eigenvalue: %lf %+lf*i \n",creal(eigenvalue),cimag(eigenvalue));
            break; 
        }
        if (j >= max_iterations){
            printf("Reached maximum iterations. CANNOT FIND REAL EIGENVECTOR AND EIGENVALUE\n");
            printf("Ignore the rest!\n"); 
        }
    } 
    printf("Number of iterations: %d \n",j);
}

void inverseIteration(int n, double complex matrixA[n][n]){
    double complex Id[n][n];
    double complex muId[n][n];
    double complex Ainv[n][n];
    int max_iterations = 100; int i, j = 0;
    //double complex bH[n][1]; 
    //double complex Ab[n][1];
    //double complex eig[1][1];
    double complex AminMuI[n][n];
    //double complex newEigVecT[1][n];
    // normalized vector as starting guess
    double complex guess_eigenvector[n][1];
    double complex new_eigenvector[n][1];
    setNormalVec(n, guess_eigenvector);
    double complex eigenvalue = 0.0 + 0.0*I;
    double complex eigenvalue0;

    // run a couple iterations of Power Algorithm to get guess for eigenvalue
    for (j = 0; j < 3; j++){
        matrix_multiplication(n,n, matrixA,n,1,guess_eigenvector,new_eigenvector); //new_eigenvector is matrixA*guessvectorb
        //Find the maximum element of new_eigenvector 
        double max_val = my_cabs(new_eigenvector[0][0]); //cabs computes the complex absolute value
        for (i = 1; i < n; i++){
            if(my_cabs(new_eigenvector[i][0]) > max_val){
                max_val = my_cabs(new_eigenvector[i][0]); 
            } 
        }
    
        // Normalize the new_eigenvector by dividing it by the maximum element 
        for (i = 0; i < n; i++) {
                new_eigenvector[i][0] = new_eigenvector[i][0] / max_val; 
        }  
        eigenvalue0 = max_val; // initial guess for the eigenvalue
        
        for (i = 0; i < n; i++) {
            guess_eigenvector[i][0] = new_eigenvector[i][0];
        }
    }
    eigenvalue = eigenvalue0;
    
    // create identity function
    fillIdentityN(n, Id);
    // - mew * I
    scalarByMatrixMultiplication(-1*(eigenvalue), n, n, Id, muId);
    // A + (- mew * I)
    matrix_addition(n, n, matrixA, muId, AminMuI);
    GetInverse(n,AminMuI,Ainv);

    while (j < max_iterations){
        //double complex new_eigenvector[n][1];

        //Multiply matrix A and guess eigenvector
        matrix_multiplication(n,n, Ainv,n,1,guess_eigenvector,new_eigenvector); 
        
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
        eigenvalue = max_val; 
       
        //Keep looping until eigvenvector of nth interation is equal to eigenvector of (n-1)th iteration
        double check = 0.0;
        for (i = 0; i < n; i++) {
            check += my_cabs(abs(new_eigenvector[i][0]) - abs(guess_eigenvector[i][0]));
        }
        
        // Update a new eigenvector as guess_eigenvector for the next iteration
        for (i = 0; i < n; i++) {
            guess_eigenvector[i][0] = new_eigenvector[i][0];
        }
        j++;


        //break out of the loop if they are equal 
        if (my_cabs(check) < 1e-6) {
            eigenvalue = 1/eigenvalue + eigenvalue0;
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


void rayleighIteration(int n, double complex matrixA[n][n]){
    double complex Id[n][n];
    fillIdentityN(n, Id);
    double complex muId[n][n];
    double complex Ainv[n][n];
    int max_iterations = 100; int j = 0;
    double complex bH[n][1]; 
    double complex Ab[n][1];
    double complex rayNum[1][1];
    double complex rayDen[1][1];
    double complex AminMuI[n][n];
    // normalized vector as starting guess
    double complex guess_eigenvector[n][1];
    setNormalVec(n, guess_eigenvector);
    double complex eigenvalue = 1.0 + 0.0*I;
    while (j < max_iterations){
        double complex new_eigenvector[n][1];

        // create identity function
        // - mew * I
        scalarByMatrixMultiplication(-1*(eigenvalue), n, n, Id, muId);
        // A + (- mew * I)
        matrix_addition(n, n, matrixA, muId, AminMuI);
        GetInverse(n,AminMuI,Ainv);

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

        // update mu
        // calculate numerator bH(A)(b)
        conjugateTranspose(n, 1, guess_eigenvector, bH);
        matrix_multiplication(n,n,matrixA,n,1,guess_eigenvector,Ab);
        matrix_multiplication(1,n,bH,n,1,Ab,rayNum);
        // calculate denominator bH(b)
        matrix_multiplication(1,n,bH,n,1,guess_eigenvector,rayDen);
        eigenvalue = rayNum[0][0] / rayDen[0][0];

        //Keep looping until eigvenvector of nth interation is equal to eigenvector of (n-1)th iteration
        double check = 0.0;
        for (i = 0; i < n; i++) {
            check += my_cabs(abs(new_eigenvector[i][0]) - abs(guess_eigenvector[i][0]));
        }

        // Update a new eigenvector as guess_eigenvector for the next iteration
        for (i = 0; i < n; i++) {
            guess_eigenvector[i][0] = new_eigenvector[i][0];
        }
        j++;

        //break out of the loop if they are equal 
        if (my_cabs(check) <= 1e-6) {
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

void LanczosAlgorithm(int dim, double complex A[dim][dim], int iter){
   // was experimenting with calloc/malloc
   // initialize variables
    //double complex q[dim][1];
    double complex (*q)[1];
    q = calloc(dim, sizeof(*q));
    double complex r[dim][1];
    double complex alpha[1][iter];
    double complex beta[1][iter];

    // store intermediate calculations
    double complex qh[1][dim];
    double complex alphaJ[1][1];
    double complex minusAlphaq[dim][1];
    double complex betaJ = 0.0 + 0.0*I;
    //double complex Aq[dim][1];
    double complex (*Aq)[1];
    Aq = calloc(dim, sizeof(*Aq));
    //double complex (*Aq)[1] = malloc(sizeof(double complex[dim][1]));
    double complex betajMin1XV[dim][1];
    int i, j;
    
    // q = normalized vector
    setNormalVec(dim, q); 

    // r = aq, size: dimxdim * dimx1 = dimx1
    matrix_multiplication(dim, dim, A, dim, 1, q, r); 

    // alpha1 = q*r, size: 1xdim * dimx1 = 1 x 1
    conjugateTranspose(dim, 1, q, qh); // conjugate of q
    //displayMatrix(dim, 1, r);
    matrix_multiplication(1, dim, qh, dim, 1, r, alphaJ);  // set alpha1
    alpha[0][0] = alphaJ[0][0]; // store alpha1

    // r = r - alpha1q
    scalarByMatrixMultiplication(-1*alphaJ[0][0], dim, 1, q, minusAlphaq); // get - alpha1q
    matrix_addition(dim, 1, r, minusAlphaq, r);

    // beta1 = ||r||
    betaJ = norm(dim, r);
   //printf("%lf", betaJ);
   
   // this breaks sometimes? :( ********
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
        // breaks here :(  *********
       
        matrix_multiplication(dim, dim, A, dim, 1, q, Aq);
        // error messages used to find location of segfault
       
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
    // free pointers
    free(q);
    free(Aq);

    // get the eigenvalue of T
    double complex T[iter][iter];
    double complex vecGuess[iter][1];
    setT(iter, alpha, beta, T);
    setNormalVec(iter, vecGuess);
    // Next step: Call power iteration to get eigenvalue of T.
    
    power_iteration(iter, T);
    
    printf("This answer took %d iterations in the Lanczos Algorithm.\n",iter);
    
}

// references pg. 178 https://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter10.pdf
// https://www.youtube.com/watch?v=2Y1ZDQw_2zw
