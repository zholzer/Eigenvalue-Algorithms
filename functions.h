// Although header files usually have function descriptions, you will find them in the functions themselves. 
// This is to highlight the thought process of how we implemented them.
// Declarations:
// Main Algorithms
void power_iteration(int n, double complex matrixA[n][n]);
void inverseIteration(int n, double complex matrixA[n][n]);
void rayleighIteration(int n, double complex matrixA[n][n]);
void LanczosAlgorithm(int dim, double complex A[dim][dim], int iter);

// Functions
void conjugateTranspose(int n, int m, double complex A[n][m], double complex AH[m][n]);
double complex det2by2(int n, double complex A[n][n]);
void displayMatrix(int n, int m, double complex A[n][m]);
void fillIdentityN(int n, double complex id[n][n]);
void GetCofactor(int N, double complex A[N][N], double complex cof[N-1][N-1], int m, int n);
double complex GetDeterminant(int N, double complex matrix[N][N]);
int GetInverse(int N, double complex A[N][N], double complex inverse[N][N]);
void GetTranspose(int N, double complex matrix[N][N], double complex matrixT[N][N]);
void matrix_addition(int n, int m, double complex A[n][m], double complex B[n][m], double complex APlusB[n][m]);
void matrix_multiplication(int a, int b, double complex matrix1[a][b], int m, int n, double complex matrix2[m][n], double complex matAXmatB[a][n]);
double my_cabs(double complex x);
double norm (int n, double complex vec[n][1]);
void scalarByMatrixMultiplication(double complex scalar, int n, int m, double complex A[n][m], double complex scalarXA[n][m]);
void setNormalVec(int n, double complex v[n][1]);
void setT(int n, double complex alpha[1][n], double complex beta[1][n], double complex T[n][n]);
