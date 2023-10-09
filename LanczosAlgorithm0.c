#include <stdio.h>
#include <math.h>
#include <complex.h>
// declare functions
void setq(int n, double complex v[n][1]);
void matrix_multiplication(int a, int b, complex double matrix1[a][b], int m, int n, complex double matrix2[m][n], complex double matAXmatB[a][n]);
void conjugateTranspose(int n, int m, double complex A[n][m], double complex AH[m][n]);
void scalarByMatrixMultiplication(double complex scalar, int n, int m, double complex A[n][m], double complex scalarXA[n][m]);
void matrix_addition(int n, int m, double complex A[n][m], double complex B[n][m], double complex APlusB[n][m]);
double norm (int n, double complex vec[n][1]);
void displayMatrix(int n, int m, double complex A[n][m]);
void setT(int n, double complex alpha[1][n], double complex beta[1][n], double complex T[n][n]);

int main(void){
    // initialize matrix

    int dim = 3;
    double complex A[3][3] = {{1,1,1}, {1,2,1}, {1,1,3}};
    int iter = 2;

   // initialize variables
    double complex q[2][1] = {{0},{1}};
    double complex Q[dim][iter];
    double complex r[dim][1];
    double complex alpha[1][iter];
    double complex beta[1][iter];
    double complex T[dim][iter];

    // store intermediate calculations
    int i, j;
    double complex qh[1][dim];
    double complex alphaJ[1][1];
    double complex minusAlphaq[1][1];
    double complex betaJ;
    double complex Aq[dim][1];
    double complex betajMin1XV[dim][1];

    // q = normalized vector
    setq(dim, q); 
    // Q1 = q
    for (i = 0; i<dim; i++){
        Q[i][0] = q[i][0];
    }

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
        // Qj = [Qj-1, q]
        for (i = 0; i<dim; i++){
            Q[i][0] = q[i][0];
        }

        // r = Aq - betaj-1v
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

    setT(dim, alpha, beta, T);
    displayMatrix(iter, iter, T);
    return 0;
// eigenpair = (value, vector) = (lambda, vector)

}

// reference pg. 178 https://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter10.pdf
// https://www.youtube.com/watch?v=2Y1ZDQw_2zw

// sets v1 as arbritrary vector with norm = 1
void setq(int n, double complex v[n][1]){
    int i; 
    for (i=0; i<n; i++){
        if(i == 1){
            v[i][0] = 1;
        }
        else{
            v[i][0] = 0;
        }
    }
}

void matrix_multiplication(int a, int b, complex double matrix1[a][b], int m, int n, complex double matrix2[m][n], complex double matAXmatB[a][n]){
    int i; int j; int k;
    if(b != m){ //number of columns of matrix1 is not the same as number of rows of matrix2
        printf("Can not do the multiplication. \n");
    } 
    else{
        for(i = 0; i < a; i++){
            for(j = 0; j < n; j++){
                for(k = 0; k < m; k++){
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
          printf("%.2f %+.2fi ", creal(A[i][j]), cimag(A[i][j])); 
      }
      printf("\n");
    }
}

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
