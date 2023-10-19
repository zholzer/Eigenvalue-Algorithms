#include <stdio.h>
#include <complex.h>

void scalarToMatrixMultiplication(int n, double complex scalar, double complex A[n][n], double complex scalarXA[n][n]){
    int i; int j; 

    for(i=0; i < n; ++i){
        for(j=0; j < n; ++j){
            scalarXA[i][j] = scalar*A[i][j]; 
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
  int size = 2;
  double complex scalar = 3;
  double complex testA[2][2] = {{1 , 2}, {2, 1}};
  double complex testScalarXA[size][size];
  scalarToMatrixMultiplication(size, scalar, testA, testScalarXA);
  printf("A\n");
  displayMatrix(size, testA);
  printf("\n");
  printf("scalarXA\n");
  displayMatrix(size, testScalarXA);
  return 0;
}
