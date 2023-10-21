This C project was conducted for COMP 526 Computational Methods for Scientists taught by Professor Miguel Dumett
Topic: Eigenvalue algorithms-Power method, Shifted inverse power method, Rayleigh quotient method, Arnoldi/Lanczos method 
Group members: Zoe Holzer, Luna Huynh, Emma Topolcsik 

Since we restrict to work with a symmetric or Hermitian matrix, we implemented four algorithms: power method, inverse power method, Rayleigh quotient method and Lanczos method to find eigenvalues and their corresponding eigenvectors.
This cproject folder contains functions.h, helperFunctions.c, main.c, PowerAlgorithms.c, testCase1.txt, testCase2.txt, testCase3.txt, testCase4.txt, testCase5.txt, and README.txt.
    functions.h is the header file that contains all functions needed to run the program. It helps reduce the complexity and number of lines of code. 
    helperFunctions,c is the c file that contains source code for essential tool functions that support the four iteration methods such as inverse of a matrix, transpose of a matrix, identity matrix, etc. 
    PowerAlgorithms.c is the c file that contains source code for power iteration function, inverse iteration function, Rayleigh iteration function and Lanzcos algorithm function. 
    Five test cases are text files that contain a 3x3 matrix, a 5x5 matrix, a 10x10 matrix, a 20x20 matrix, and a 50x50 matrix, respectively. 
    README.txt (this file) introduces the project and explains how to run this iteration program. 
    
Please follow the following steps to run our power iteration program:
1. To compile the program, in the terminal type this command: gcc -Wall helperFunctions.c main.c PowerAlgorithms.c -lm 
2. To run the program, in the terminal type this command: ./a.out
3. Follow steps by steps what the program is asking the user to input 
    - First, the program asks the user to choose test file (total of 5 test cases), the user need to enter an integer between 1 and 5
    - The program asks if the user wants to display the matrix, the user need to enter y for yes or n for n
    - Lastly, the program asks the user to enter the number of iterations for the Lanczos iteration 
