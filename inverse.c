#include <stdio.h>
#include <complex.h>
#include <math.h>
// function to find the inverse of a matrix
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
void displayMatrix(int n, double complex A[n][n]){
    int i, j; 
    for (i=0; i<n; i++){
      for (j=0; j<n; j++){
            
            printf("%.2f %+.2fi ", creal(A[i][j]),cimag(A[i][j])); 
        
      }
      printf("\n");
    }
}  

// cofactor function (only works for higher than 2x2 matrices)
void GetCofactor(int N, double complex A[N][N], double complex cof[N][N], int m, int n)
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
            if (r != m && c != n)
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

}

// function to rewrite the matrices from an NxN to an (N-1)x(N-1)
// this is for after calling cofactor function
void rewrite(int N, double complex cof[N][N], double complex cofrewrite[N-1][N-1])
    {
        int r, c;
         for (r = 0; r < N-1; r++)
         {
              for (c = 0; c < N-1; c++)
               {
              cofrewrite[r][c] = cof[r][c];
            }
        }
    
    }


// function for the determinant of a 2x2 matrix
// inputs: n size matrix, A[n][n] matrix
double complex det2by2(int n, double complex A[n][n])
{

    // if statement if the matrix is not a 2x2
    if (n != 2)
    {
        printf("Please input a 2x2 matrix.\n");
        return 1;
    }
    
    double complex a, b, c, d, det;

    a = A[0][0]; b = A[0][1];
    c = A[1][0]; d = A[1][1];


    det = a*d - b*c;

    //displayMatrix(n, A);
    //printf("The determinant is %.2f %+.2fi ", creal(det),cimag(det));

    return det;

}


// take and NxN matrix and computes the transpose
void GetTranspose(int N, double complex matrix[N][N], double complex matrixT[N][N])
{
    int row, col;

    for (row = 0; row < N; row++)
    {
        for (col = 0; col < N; col++)
        {
            matrixT[col][row] = matrix[row][col];
        }
    }

}

// determinant function (for 3x3 + sized matrices)
double complex GetDeterminant(int N, double complex matrix[N][N])
{
    int row, col, sign;
    double complex det, determinant, cof[N][N], cofrewrite[N-1][N-1];
    
    determinant = 0 + 0*I;
    
    // expand the first row:
    row = 0;
    for ( col = 0; col < N; col++)
        {
            // call cofactor function
            GetCofactor(N, matrix, cof, row, col);
            rewrite(N, cof, cofrewrite);

            //if ()
            // call the determinant function 
            det = det2by2(N-1, cofrewrite);

            // save the determinant in a sum
            sign = pow(-1, row+col);
            determinant = determinant + sign*det*matrix[row][col];

            //printf("Determinant at [%i][%i] = %.2f %+.2fi\n", row, col, creal(determinant),cimag(determinant));
        }
    return determinant;
}


// main function
int main (void)
{
    int i, j, N, sign;
    
    // define some matrix 
    N = 2;
    double complex A[2][2] = { { 2, 3},
                    { 4, 2}
                    };

    // initialize for cofactor function
    double complex determinant, det, cofM[N][N];    // determinant, and cofactor matrix (det is the determinant in the cofactor loop)

        // get the determinant of the matrix
        if (N == 2)
        {
            determinant = det2by2(N, A);
            printf("Determinant= %.2f %+.2fi\n", creal(determinant),cimag(determinant));
        }
        else 
        {
            determinant = GetDeterminant(N, A);
            printf("Determinant= %.2f %+.2fi\n", creal(determinant),cimag(determinant));
        }

        // if the determinant = 0, there is no inverse!!
        if (determinant == 0) {
            printf("Determinant = 0. There is no inverse. Please input another matrix.\n");
            return 1;
        }

    // find the Adjoint- first find the cofactor matrix
    // find the cofactor for each element, find the determinant of that matrix, determinant = that element
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            // if statement if the matrix is a 2x2
            if (N == 2) { 
                double complex cof; // empty matrix for cofactor 
                int row,col;
                for (row = 0; row < N; row++){
                    for (col = 0; col < N; col++){
                        if (row != i && col != j)
                        {
                            cof = A[row][col];
                            det = cof;
                        }
                    }
                }
            }

            else {
                double complex cof[N][N], cofrewrite[N-1][N-1]; // empty matrix for cofactor
                // for all other N sized matrices:
                // call cofactor function
                GetCofactor(N, A, cof, i, j);
                rewrite(N, cof, cofrewrite);
                displayMatrix(N,cofrewrite);

                // call the determinant function 
                det = det2by2(N-1, cofrewrite);
           }


            // save the determinant in the corresponding element
            // this is the cofactor matrix cofM
            sign = pow(-1, i+j);
            cofM[i][j] = det*sign;

        }
    }

    // display the cofactor matrix
    printf("The cofactor matrix:\n");
    displayMatrix(N,cofM);

    // now that we have the cofactor matrix, the transpose = Adjugate
    // adj will hold the adjugate matrix
    double complex adj[N][N];
    GetTranspose(N, cofM, adj);

    printf("The adjugate:\n");
    displayMatrix(N,adj);

    // inverse = Adjugate/determinant
    double complex inverse[N][N];

    // divide the determinant from each element
    for (i = 0; i < N; i++) 
    {
        for (j = 0; j < N; j++) 
        {
            inverse[i][j] = adj[i][j] / determinant;
        }
    }

    // print
    printf("The Inverse: \n");
    displayMatrix(N,inverse);

    return 0;
}
