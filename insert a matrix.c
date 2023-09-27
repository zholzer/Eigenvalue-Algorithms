// Add shared functions to this file.

int main(){
    int row; int column; int n; 
    float A[100][100];
    float b[100];

    printf("Enter the order of your matrix A: ");
    scanf("%d",&n);

    printf("Enter the elements of your square matrix A: \n");
    for (row=0; row < n; ++row){
        for (column=0; column < n; ++column){
            printf("A[%d][%d]= " ,row,column); 
            scanf("%f", &A[row][column]); 
            }
        }
        
            
    printf("Enter the elements of your b vector: \n");
    for (row=0; row < n; ++row){
        printf("b[%d]= ",row);
        scanf("%f", &b[row]);
    }
    return 0;
        
}
    cross_product(A,b);
    return 0;
}
