#include <stdio.h>
#include <complex.h>

// takes in a declared complex double matrix A of size n x n, complex double matrix AH of size n x n, and size n
// fills AH with the conjugate transpose of A
void conjugateTranspose(int n, double complex A[n][n], double complex AH[n][n]){
    int i,j; 
    for (i=0; i<n; i++){
        for (j=0; j<n; j++){
            AH[j][i] = conj(A[i][j]); // get conjugate of values and then allocate as the transpose to AH
	    }
    }
}

// takes in a declared complex double matrix A of size n x n and size n
// displays matrix, altered from Luna's code
void displayMatrix(int n, double complex A[n][n]){
    int i, j; 
    for (i=0; i<n; i++){
      for (j=0; j<n; j++){
            printf("%.2f %+.2fi ", creal(A[i][j]), cimag(A[i][j])); 
      }
      printf("\n");
    }
}  

int main(void){
// tests a size 3 complex matrix
  int size = 3;
  double complex testA[3][3] = {{1.0, 2.0 + I, 3.0}, {4.0, 5.0, 6.0 + I}, {7.0 + I, 8.0, 9.0}};
  double complex testAH[size][size];
  conjugateTranspose(size, testA, testAH);
  printf("A\n");
  displayMatrix(size, testA);
  printf("\n");
  printf("AH\n");
  displayMatrix(size, testAH);
  return 0;
}
