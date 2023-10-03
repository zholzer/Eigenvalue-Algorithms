#include <stdio.h>
// function to find the cofactor of a matrix
// what it will consist:
//
// GetCofactor is a function to compute the cofactor of a matrix A[p][q]
// within this function:
//
// inputs: int size of matrix N, given matrix A[N][N], empty matrix cof[N][N] to hold output,
//         ints m/n to specify the element we're at in A
//
// outputs: gets the cofactor of A[m][n]
//
// there will also be a main script
// within the main script, we define the matrix
// and run a for loop to call the cofactor at each element A[p][q]

// "display matrix function"
// takes in a declared matrix A of size n x n and size n
void displayMatrix(int n, int A[n][n]){
    int i, j; 
    for (i=0; i<n; i++){
      for (j=0; j<n; j++){
            
            printf("%i ", A[i][j]); 
        
      }
      printf("\n");
    }
}  

// cofactor function
void GetCofactor(int N, int A[N][N], int cof[N][N], int m, int n)
{
    int i, j, r, c;

    // initialize the index
    i = 0; j = 0;

    // for loop thorugh each element
    for (r = 0; r < N; r++)
    {
        for (c = 0; c < N; c++)
        {
            // identify the rows/cols that we are NOT at
            if (r != n && c != m)
            {
                cof[i][j] = A[r][c];
                j++;

                // once the row is filled, we reset for next row
                if(j == N-1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
    
    
    printf("The cofactor for A[%i][%i] =\n", m, n);
    displayMatrix(N, cof);
}


// main function
int main (void)
{
    int i, j, N;
    // some matrix

    N = 4;
    int A[4][4] = { { 5, -2, 2, 7 },
                    { 1, 0, 0, 3 },
                    { -3, 1, 5, 0 },
                    { 3, -1, -9, 4 } };

    int cof[N][N];  // empty matrix for cofactor

    // the way we will use this to find inverse:
    // int adj[N][N]; float inv[N][N];
    // int adjoint(A,adj); inverse(A,inv)
    // lines above call the adjoint function and inverse function
    // within the adjoint function, the cofactor function is called as so:

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            // call cofactor function
            GetCofactor(N, A, cof, i, j);
        }
    }
    return 0;
}
