#pragma once
#include <vector>
#include <string>

class Matrix
{
	
private:

	int m_rows;
	int m_cols;
	std::vector < std::vector < double > > m_matrix;
	void initZeros(int, int);
	void initOnes(int, int);
	
public:

	Matrix(int, int, std::vector<double>);
	Matrix(int, int, std::string);
	Matrix operator+( const Matrix&);
	Matrix operator*( const Matrix&);
	Matrix operator=( const Matrix&);
	Matrix operator-( const Matrix&);
	Matrix transpose();
	Matrix choleskyDecomposition();
	void gaussianElimination(std::string);
	std::vector<double> backwardSubstitution();
	std::vector<double> forwardSubstitution();
	int getHeight();
	int getLength();
	void setValue(int,int,double);
	double getValue(int,int);
	void sumValue(int,int,double);
	void swapRows(int, int);
	void swapCols(int,int);
	void substituteCol(int, std::vector<double>);
	void Print() const;
	
};