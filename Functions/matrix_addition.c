#include <stdio.h>

int main(void){
    void matrix_addition(int nRows, int nCols, int A[nRows][nCols],int B[nRows][nCols]);
    void display_matrix(int nRows, int nCols, int matrix[nRows][nCols]); 
    int matrix1[3][3] = 
    {
        {3,  23, 55},
        {12, 42, 10},
        {-2, 12,  9},
    };

    int matrix2[3][3] = 
    {
        {7,  20, 12},
        {16, 18, -10},
        {-2, -9,  7},
    };
    printf("Matrix 1:\n");
    display_matrix(3,3,matrix1); 

    printf("Matrix 2:\n");
    display_matrix(3,3,matrix2);

    printf("Addition of Matrix 1 and Matrix 2:\n");
    matrix_addition(3,3,matrix1,matrix2);


}
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

    for(row=0; row < 3; ++row){
        for(column=0; column < 3; ++column){
            C[row][column] = A[row][column] + B[row][column]; 
            //printf("%5i ",C[row][column]);
        }
    }
    
    display_matrix(nRows, nCols,C);
}
