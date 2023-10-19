// this program asks the user to input a vector, and prints the norm
#include <stdio.h>
#include <math.h>
#include <complex.h>

double norm (double complex vec[], int n)
{
    double term, sum;  
    int ii;

    // error code for an empty vector
    if (n == 0){
        printf("Error: Empty Vector. Retry.\n");
        return 1;
    }


    // verify what they entered
    printf("You entered: [ ");

    for (ii = 0; ii < n; ++ii){
    printf("%.2f %+.2fi ", creal(vec[ii]),cimag(vec[ii]));
    }

    printf("] with %i elements, ",n);


    sum = 0;              // initialize the sum

    for (ii = 0; ii < n; ++ii){
        term = pow(creal(vec[ii]),2) + pow(cimag(vec[ii]),2);
        sum = sum + term;
    }

    // after the loop, square root the sum
    double sol = sqrt(sum);

    // print the norm
    printf("and the norm is %f\n",sol);

    return sol;
}

int main (void)
{
    // define a vector 
    double complex vec[] = {1 + 2*I, 7 +  I};
    int n = sizeof(vec) / sizeof(vec[0]);

    // run the program
    norm(vec,n);

    return 0;
}
