// Add shared functions to this file.

#include <stdio.h>

void display_matrix(int nRows, int nCols, int matrix[nRows][nCols]){
    int row; int column; 
    for(row=0; row < nRows; ++row){
        for(column=0; column < nCols; ++column)
            printf("%5i ",matrix[row][column]); 
        printf("\n");
    }
} 

void matrix_addition(int nRows, int nCols, int A[nRows][nCols], int B[nRows][nCols]){
     int row; int column; 
    int C[nRows][nCols];

    for(row=0; row < nRows; ++row){
        for(column=0; column < nCols; ++column){
            C[row][column] = A[row][column] + B[row][column]; 
        }
    }
    printf("The Addition of two matrix: \n"); 
    display_matrix(nRows, nCols,C);
}

void matrix_multiplication(int a, int b, int matrix1[][100], int m, int n, int matrix2[][100]){
    int row; int column; int k;
    int C[a][n];
    if(b != m) //number of columns of matrix1 is not the same as number of rows of matrix2
        printf("Can not do the multiplication. \n"); 
    else
        for(row = 0; row < a; row++){
            for(column = 0; column < n; column++){
                C[row][column]=0;
            }
        }

        for(row = 0; row < a; row++){
            for(column = 0; column < n; column++){
                for(k = 0; k < m; k++){
                    C[row][column] += matrix1[row][k]*matrix2[k][column];
                }
            }
        }
        printf("The Product of two matrix: \n"); 
        display_matrix(a,n,C); 
}

int main() {
    int row; int column; int n; 
    int A[100][100];
    int b[100][100];

    printf("Enter the order of your matrix A: ");
    scanf("%d",&n);

    printf("Enter the elements of your square matrix A: \n");
    for (row=0; row < n; ++row){
        for (column=0; column < n; ++column){
            printf("A[%d][%d]= " ,row,column); 
            scanf("%d", &A[row][column]); 
            }  
        }

    printf("Enter the elements of your b vector: \n");
    for (row=0; row < n; ++row){
        for(column=0; column < 1; ++column){
            printf("b[%d]= ",row);
            scanf("%d", &b[row][column]);
        }  
    }
    matrix_multiplication(n,n,A,n,1,b);
    return 0;       
}
