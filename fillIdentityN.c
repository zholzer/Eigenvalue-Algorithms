#include <stdio.h>
#include <complex.h>

// takes in a declared complex double matrix id of size n x n and size n
// fills with the identity matrix at that size
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
// tests for size 4 identity matrix
  int size = 4;
  double complex idTest[size][size];
  fillIdentityN(size, idTest);
  displayMatrix(size, idTest);
  return 0;
}
