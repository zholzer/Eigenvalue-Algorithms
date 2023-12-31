#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include "functions.h"

int main(void){
    printf("Welcome to the Power Iteration program.\n");

    // Allows user to switch between our test files.
    printf("Please choose which test file to run by entering an integer between 1 and 5.\nThe matrices are of size 3 x 3, 5 x 5, 10 x 10, 20 x 20, and 50 x 50 respectivly.\n");
    char fileNum;
    char fileName[15];
    // gets case number from user, currently not in any particular order
    int booL0 = 0;
    while(booL0 == 0){
        scanf(" %c", &fileNum);
        int temp = fileNum - '0';
        if (temp >= 1 && temp <= 5){
            booL0 = 1;
        }
        else{printf("Invalid test case selection. Retry. \n");}
    }

    // switches the file to be read based on input
    switch (fileNum){
    case '1':
    // test file one breaks in Lanczos ******
        strcpy(fileName, "testCase1.txt");
        printf("You chose test file one with a 3 x 3 matrix.\n");
        break;

    case '2':
        strcpy(fileName, "testCase2.txt");
        printf("You chose test file two with a 5 x 5 matrix.\n");
        break;

    case '3':
        strcpy(fileName, "testCase3.txt");
        printf("You chose test file three with a 10 x 10 real matrix.\n");
        break;

    case '4':
        strcpy(fileName, "testCase4.txt");
        printf("You chose test file four with a 20 x 20 real matrix.\n");
        break;

    case '5':
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
    // does the user want to see the matrix? some are very large
    char disp;
    printf("Would you like to display the selected matrix? (y = yes / n = no) \n");
    int booL1 = 1;
    while(booL1 == 1){
         scanf(" %c", &disp);
        if (disp == 'y'){
            printf("This test file contains the matrix: \n");
            displayMatrix(dim, dim, A);
            booL1 = 0;
        }
        else if (disp == 'n'){
            booL1 = 0;
        }
        else{printf("Please input y for yes or n for no. Retry. \n");}
    }
   
    // allows user to input how many iterations they want Lanczos to run
    printf("Input number of iterations of the Lanczos Algorithm as an integer (must be less then %d).\nThis will determine the size of T, a smaller matrix with the leading eigenvalues of A.\n", dim);
    int iter, temp;
    int booL2 = 0;
    while(booL2 == 0){
        scanf("%d", &iter);
        if ((iter < dim) && (iter > 1)){
            booL2 = 1;
        }
        else{
            printf("Must be less then %d. Retry. \nInput number of iterations of the Lanczos Algorithm as an integer.\n", dim);
            // https://stackoverflow.com/questions/34165592/scanf-skipped-after-reading-in-integer-in-c-in-while-loop
            // clears scanf so it asks for input again
            while((temp = getchar()) != '\n');
        }
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
// https://stackoverflow.com/questions/36890624/malloc-a-2d-array-in-c
//stackoverflow.com/questions/29977084/calloc-a-two-dimensional-array
