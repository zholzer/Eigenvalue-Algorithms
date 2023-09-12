// Add shared functions to this file.
void cross_product(float vector_A[10][10], float vector_b[10]){
    float P[10];
    P[0] = vector_A[0][0]*vector_b[0] + vector_A[0][1]*vector_b[1];
    P[1] = vector_A[1][0]*vector_b[0] + vector_A[1][1]*vector_b[1];
    float Ab[] = {P[1],P[0]};
    printf("A*b = %f %f \n", Ab[0],Ab[1]);
}

int main(){
    int i,j,n; 
    float A[10][10],b[10],P[10];
    printf("\nEnter the elements of your 2x2 matrix: \n"); 
    for (i = 0; i < 2; i++) {
        for (j = 0; j < 2; j++) {
            printf("A[%d][%d]=", i,j);
            scanf("%f",& A[i][j]);
            }    
        }
    
    printf("\nEnter the elements of your 2-element vector: \n");
    for (i = 0; i < 2; i++){
        printf("b[%d]=",i);
        scanf("%f",& b[i]);
    }
    cross_product(A,b);
    return 0;
}
