#include "Matrix.h"
#include <iostream>
#include "Utility.h"

#define PRINTL(x) std::cout<< x;
#define LOG(x) std::cout<< x <<std::endl;

void choleskyExample()
{
	/*    
	 * 	Linear System: Ax = b
	 * 
	 * 	   25  15  -5   1
	 * 
	 *     15  18   0   1
	 * 
	 *     -5   0  11   1
	 * 
	 * 
	 */
	std::vector<double> values  = {25,15,-5,1,15,18,0,1,-5,0,11,1}; //Values
	int rows = 3; //Numbers of rows
	int cols = 4; //We have 4 cols because the last one contains the known terms
	Matrix A = Matrix(rows,cols,values); //Inizialize the matrix
	
	LOG("Matrix A")
	LOG("")
	A.Print();
	LOG("")
	
	Matrix B = A.choleskyDecomposition(); //It performs the decomposition and store it in B
	
	LOG("Matrix A decomposed")
	LOG("")
	B.Print();
	LOG("")
	
	/*
	 * Linear system B after cholesky decomposition:
	 * 
	 *      5   0   0   1
	 * 
	 *      3   3   0   1
	 * 
	 *     -1   1   3   1
	 * 
	 */
	
	 Matrix B_t = B.transpose(); // Transpose the matrix B
	
	 std::vector<double> y = B.forwardSubstitution(); //It performs the forward substitution
	
	/*
	 * Now i'm going to substitute the last column (known terms)
	 * With the partial solution y i just found
	 */
	 
	 B_t.substituteCol(cols - 1, y);
	
	 std::vector<double> x = B_t.backwardSubstitution(); // It performs the backward substitution
	 LOG("Solutions: ")
	 printVector(x);
	
}


void gaussianExample()
{
	
	std::vector<double> values  = {1,-3,1,4,2,-8,8,-2,-6,3,-15,9}; //Values
	int rows = 3; //Numbers of rows
	int cols = 4; //We have 4 cols because the last one contains the known terms
	Matrix A_WP = Matrix(rows,cols,values); //Inizialize the matrix
	Matrix A_PP = A_WP;
	Matrix A_CP = A_WP;
	
	A_WP.gaussianElimination("no_pivot");
	A_PP.gaussianElimination("partial_pivot");
	A_CP.gaussianElimination("complete_pivot");
	
	std::vector<double> x_wp = A_WP.backwardSubstitution();
	std::vector<double> x_pp = A_PP.backwardSubstitution();
	std::vector<double> x_cp = A_CP.backwardSubstitution();
	
	LOG("Solution without pivoting :")
	
	printVector(x_wp);
	
	LOG("Solution with partial pivoting :")

	printVector(x_pp);
	
	LOG("Solution with complete pivoting :")

	printVector(x_cp);
	
	
	
}