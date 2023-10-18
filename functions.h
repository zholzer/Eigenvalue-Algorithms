// header file
void setNormalVec(int n, double complex v[n][1]);
void matrix_multiplication(int a, int b, double complex matrix1[a][b], int m, int n, double complex matrix2[m][n], double complex matAXmatB[a][n]);
void conjugateTranspose(int n, int m, double complex A[n][m], double complex AH[m][n]);
void scalarByMatrixMultiplication(double complex scalar, int n, int m, double complex A[n][m], double complex scalarXA[n][m]);
void matrix_addition(int n, int m, double complex A[n][m], double complex B[n][m], double complex APlusB[n][m]);
double norm (int n, double complex vec[n][1]);
void displayMatrix(int n, int m, double complex A[n][m]);
void setT(int n, double complex alpha[1][n], double complex beta[1][n], double complex T[n][n]);
void LanczosAlgorithm(int d, double complex A[d][d], int itr);
double my_cabs(double complex x);
void rayleighIteration(int n, double complex matrixA[n][n]);
void GetCofactor(int N, double complex A[N][N], double complex cof[N-1][N-1], int m, int n);
double complex det2by2(int n, double complex A[n][n]);
void GetTranspose(int N, double complex matrix[N][N], double complex matrixT[N][N]);
double complex GetDeterminant(int N, double complex matrix[N][N]);
int GetInverse(int N, double complex A[N][N], double complex inverse[N][N]);
void power_iteration(int n, double complex matrixA[n][n]);
