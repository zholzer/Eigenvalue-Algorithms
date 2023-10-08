#include <stdio.h>
#include <complex.h>

// display a matrix function 
void displayMatrix(int n, double complex A[n][n])
{
    int i, j; 
    for (i=0; i<n; i++){
      for (j=0; j<n; j++){
            
            printf("%.2f %+.2fi ", creal(A[i][j]),cimag(A[i][j])); 
        
      }
      printf("\n");
    }
}  


// determinant of a 2x2 matrix
// inputs: n size matrix, A[n][n] matrix
int det2by2(int n, double complex A[n][n])
{

    // if statement if the matrix is not a 2x2
    if (n != 2)
    {
        printf("Please input a 2x2 matrix.\n");
        return 1;
    }
    
    complex double a, b, c, d, det;

    a = A[0][0]; b = A[0][1];
    c = A[1][0]; d = A[1][1];


    det = a*d - b*c;

    displayMatrix(n, A);
    printf("The determinant is %.2f %+.2fi ", creal(det),cimag(det));

    return 0;

}


// main code to call 
int main(void)
{

    int n;
    //double complex det;


    // define some matrix
    n = 2;
    double complex A[2][2] = { { 5, -2 },
                               { 1,  0}                  
                               };

    det2by2(n,A);

    return 0;
}
