// Add shared functions to this file.
void cross_product(float vector_A[10][10], float vector_b[10]){
    float P[10];
    P[0] = vector_A[0][0]*vector_b[0] + vector_A[0][1]*vector_b[1];
    P[1] = vector_A[1][0]*vector_b[0] + vector_A[1][1]*vector_b[1];
    float Ab[] = {P[1],P[0]};
    printf("A*b = %f %f \n", Ab[0],Ab[1]);
}

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
