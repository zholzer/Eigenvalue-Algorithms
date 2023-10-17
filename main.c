#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include "functions.h"
// function declarations

int main(void){
    printf("Welcome to the Power Iteration program.\n");

    // Allows user to switch between our test files.
    printf("Please choose which test file to run by entering an integer between 1 and 5.\n");
    int fileNum;
    char fileName[15];
    scanf("%d", &fileNum);
    switch (fileNum){
    case 1:
        strcpy(fileName, "testCase1.txt");
        printf("You chose test file one with a 10 x 10 matrix.\n");
        break;

    case 2:
        strcpy(fileName, "testCase2.txt");
        printf("You chose test file two with a 3 x 3 matrix.\n");
        break;

    case 3:
        strcpy(fileName, "testCase3.txt");
        printf("You chose test file three with a 10 x 10 real matrix.\n");
        break;

    case 4:
        strcpy(fileName, "testCase4.txt");
        printf("You chose test file four with a 10 x 10 real matrix.\n");
        break;

    case 5:
        strcpy(fileName, "testCase5.txt");
        printf("You chose test file five with a 10 x 10 real matrix.\n");
        break;
    
    default:
        printf("Input must be an integer between 1 and 5.\n");
        break;
    }

    // Opens and reads indicated test file to store the values in matrix A
    int i, j;
    double c, n;
    int count, dim;
    //char contents[500];
    // initialize matrix
    FILE *testFile;
    testFile = fopen(fileName, "r");
    while (1){
        fscanf(testFile, "%lf", &c);
        if(feof(testFile)){break;}
        count++;
    }
    
    dim = sqrt(count/2);
    //double complex A[dim][dim];
    double complex (*A)[dim] = malloc(sizeof(double complex[dim][dim]));

    rewind(testFile);
    for (i=0; i<dim; i++){
        for (j=0; j<dim; j++){
            fscanf(testFile, "%lf", &n);
            A[i][j] = n;
            if(feof(testFile)){break;}
        }
    }
    for (i=0; i<dim; i++){
        for (j=0; j<dim; j++){
            fscanf(testFile, "%lf", &n);
            A[i][j] = A[i][j] + n*I;
            if(feof(testFile)){break;}
        }
    }
    fclose(testFile);
    
    printf("This test file contains the matrix: \n");
    displayMatrix(dim, dim, A);

    int iter;
    int booL = 0;
    while(booL == 0){
        printf("Input number of iterations as an integer (must be less then %d).\n", dim);
        scanf("%d", &iter);
        if (iter > dim){
            printf("Must be less then %d. Retry. \n", dim);
        }
        else{booL = 1;}
    }
    //printf("The Power Iteration produces ");
    //power_iteration(dim, A); // runs with 10 by 10
    //printf("The Power Iteration produces ");
    //rayleighIteration(dim, A); // very slow with 10 by 10
    printf("The Lanczos Algorithm produces ");
    LanczosAlgorithm(dim, A, iter); // segfault 10 by 10
    //free(A);
    return 0;
}

// https://stackoverflow.com/questions/36890624/malloc-a-2d-array-in-c
// https://stackoverflow.com/questions/36890624/malloc-a-2d-array-in-c
