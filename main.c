#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include "eigenvalueFunctions.h"

int main(void){
    printf("Welcome to the Power Iteration program.\n");

    // Allows user to switch between our test files.
    printf("Please choose which test file to run by entering an integer between 1 and 5.\n");
    int fileNum;
    char fileName[15];
    // gets case number from user, currently not in any particular order
    scanf("%d", &fileNum);
    switch (fileNum){
    case 1:
    // test file one breaks in Lanczos ******
        strcpy(fileName, "testCase1.txt");
        printf("You chose test file one with a 3 x 3 matrix.\n");
        break;

    case 2:
        strcpy(fileName, "testCase2.txt");
        printf("You chose test file two with a 5 x 5 matrix.\n");
        break;

    case 3:
        strcpy(fileName, "testCase3.txt");
        printf("You chose test file three with a 10 x 10 real matrix.\n");
        break;

    case 4:
        strcpy(fileName, "testCase4.txt");
        printf("You chose test file four with a 20 x 20 real matrix.\n");
        break;

    case 5:
        strcpy(fileName, "testCase5.txt");
        printf("You chose test file five with a 50 x 50 real matrix.\n");
        break;
    
    default:
        printf("Input must be an integer between 1 and 5.\n");
        break;
    }

    // Get size of test file matrix. 
    int i, j;
    double c, n;
    int dim;
    // Opens file
    FILE *testFile;
    testFile = fopen(fileName, "r");
    // increases count for each long float in test file
    int count = 0;
    while (1){
        fscanf(testFile, "%lf", &c);
        if(feof(testFile)){break;}
        count++;
    }
    // file consists of two equally sized n by n matrices; 
    // the first matrix is read in as the real values (a of a+bi)
    // the second is read in as the imaginary values (b of a+bi)
    // the size of the matrix is the # digits / 2 square rooted since n by n
    dim = sqrt(count/2);
    //double complex A[dim][dim];
    //double complex (*A)[dim] = malloc(sizeof(double complex[dim][dim]));
    double complex (*A)[dim];
    A = calloc(dim, sizeof(*A));
    // rewind file pointers to store matrix values
    rewind(testFile);
    // real matrix
    for (i=0; i<dim; i++){
        for (j=0; j<dim; j++){
            fscanf(testFile, "%lf", &n);
            A[i][j] = n;
            if(feof(testFile)){break;}
        }
    }
    // imaginary matrix
    for (i=0; i<dim; i++){
        for (j=0; j<dim; j++){
            fscanf(testFile, "%lf", &n);
            A[i][j] = A[i][j] + n*I;
            if(feof(testFile)){break;}
        }
    }
    // done with the file
    fclose(testFile);
    
    char disp;
    printf("Would you like to display the selected matrix? Default is no. (y = yes / n = no) \n");
    scanf(" %c", &disp);
    if (disp == 'y'){
        printf("This test file contains the matrix: \n");
        displayMatrix(dim, dim, A);
    }

    // allows user to input how many iterations they want Lanczos to run
    int iter;
    int booL = 0;
    while(booL == 0){
        printf("Input number of iterations of the Lanczos Algorithm as an integer (must be less then %d).\n", dim);
        scanf("%d", &iter);
        if ((iter > dim) || (iter < 1)){
            printf("Must be less then %d. Retry. \n", dim);
        }
        else{booL = 1;}
    }

    // calls algorithms
    printf("The Power Iteration produces ");
    power_iteration(dim, A); // runs with 10 by 10
    if (dim < 20){
        printf("The Inverse Iteration produces ");
        inverseIteration(dim, A);
    }
    else{printf("Too big for the Inverse Algorithm (the inverse is too expensive). \n");}
    if (dim < 10){
        printf("The Rayleigh Iteration produces ");
        rayleighIteration(dim, A); // very slow with 10 by 10
    }
    else{printf("Too big for the Rayleigh Algorithm (the inverse is too expensive). \n");}
    printf("The Lanczos Algorithm produces ");
    LanczosAlgorithm(dim, A, iter); // segfault 10 by 10
    // free pointer
    free(A);
    return 0;
}

// references
// https://stackoverflow.com/questions/36890624/malloc-a-2d-array-in-c
//stackoverflow.com/questions/29977084/calloc-a-two-dimensional-array
