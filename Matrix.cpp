#include "Matrix.h"
#include <iostream>
#include <exception>
#include <random>
#include <iomanip>
#include <algorithm>
#include <limits>
#include <cmath>
//#include <math>

#define PRINTL(x) std::cout<< x;
#define LOG(x) std::cout<< x <<std::endl;


void Matrix::initOnes(int rows, int cols)
{
	m_rows = rows;
	m_cols = cols;
	m_matrix.resize(m_rows);
	for( int i = 0; i < m_cols; i++)
	{
		m_matrix[i].resize(m_cols,1);
	}
}

void Matrix::initZeros(int rows, int cols)
{
	m_rows = rows;
	m_cols = cols;
	m_matrix.resize(m_rows);
	for( int i = 0; i < m_cols; i++)
	{
		m_matrix[i].resize(m_cols,0);
	}
}

///<summary> Constructor, it inizializes the matrix with a vector<double> of values</summary>
///<param name="rows"> number of rows </param>
///<param name="cols"> number of cols </param>
///<param name="values"> vector<double> of values </param>
///<returns> A matrix rows x cols which contains values </returns>
Matrix::Matrix(int rows, int cols, std::vector<double> values)
{
	initZeros(rows, cols);
	for( int i = 0; i < m_rows; i++)
	{
		for( int j = 0; j < m_cols; j++)
		{
			m_matrix[i][j] = values[m_cols * i + j];
		}
	}
}


///<summary> Constructor, it inizializes the matrix with zeros or ones or at random</summary>
///<param name="rows"> number of rows </param>
///<param name="cols"> number of cols </param>
///<param name="type"> string, it can be "zeros", "ones" or "random"</param>
///<returns> A matrix rows x cols which is inizialized as type </returns>
Matrix::Matrix(int rows, int cols, std::string type)
{
	if( type == "zeros")
	{
		initZeros(rows, cols);
	}
	else if(type == "ones")
	{
		initOnes(rows, cols);
	}
	else if(type == "random")
	{
		initZeros(rows,cols);
		constexpr int MIN = -100;
		constexpr int MAX = 100;
		std::srand(std::time(nullptr));
		for( int i = 0; i < m_rows; i++)
		{
			for( int j = 0; j < m_cols; j++)
			{
				m_matrix[i][j] = MIN + (double)(rand()) / ((double)(RAND_MAX/(MAX - MIN)));
			}
		}
	}
	else
	{
		std::cerr<<" Matrix can be random, zeros or ones";
	}
}


Matrix Matrix::operator=(const Matrix& b)
{
	Matrix copy = Matrix(m_rows,m_cols,"zeros");
	for(int i = 0; i < m_rows; i++)
	{
		for( int j = 0; j < m_cols; j++)
		{
			copy.setValue(i,j,m_matrix[i][j]);
		}
	}
}


Matrix Matrix::operator*(const Matrix& b)
{
	if(m_cols != b.m_rows)
		std::cerr<<"Matrix sizes are not compatible";
	Matrix result = Matrix(m_rows, m_cols, "zeros");
	for( int i = 0; i < m_rows; i++)
	{
		for( int j = 0; j < m_cols; j++)
		{
			for( int k = 0; k < m_cols; k++)
			{
				result.sumValue(i,j, m_matrix[i][k] * b.m_matrix[k][j]);
			}
		}
	}
	return result;
}


Matrix Matrix::operator-(const Matrix& b)
{
	if(m_cols != b.m_cols ||  m_rows != b.m_rows)
		std::cerr<<"Matrix sizes are not compatible";
	Matrix result = Matrix(m_rows, m_cols, "zeros");
	for( int i = 0; i < m_rows; i++)
	{
		for( int j = 0; j < m_cols; j++)
		{
			result.setValue(i,j, m_matrix[i][j] - b.m_matrix[i][j]);
		}
	}
	return result;
}


Matrix Matrix::operator+(const Matrix& b)
{
	if(m_rows != b.m_rows || m_cols != b.m_cols)
		std::cerr<<"Matrix sizes are different!";
	Matrix result = Matrix(m_rows, m_cols, "zeros");
	for( int i = 0; i < m_rows; i++)
	{
		for( int j = 0; j < m_cols; j++)
		{
			result.setValue(i,j, m_matrix[i][j] + b.m_matrix[i][j]);
		}
	}
	return result;
}


Matrix Matrix::transpose()
{
	Matrix transposed = Matrix(m_rows, m_cols, "zeros");
	for( int i = 0; i < m_rows; i++)
	{
		transposed.setValue(i,m_cols-1, m_matrix[i][m_cols-1]);
		transposed.setValue(i,i,m_matrix[i][i]); 
	}
	for( int i = 1; i < m_rows; i++)
	{
		for(int j = 0; j < i ; j++)
		{
			transposed.setValue(i,j,m_matrix[j][i]);
			transposed.setValue(j,i,m_matrix[i][j]);
		}
	}
	return transposed;
}


Matrix Matrix::choleskyDecomposition()
{
	double sqrt_term;
	Matrix decomposed = Matrix(m_rows, m_cols, "zeros");
	for( int i = 0; i < m_rows; i++)
	{
		sqrt_term = std::sqrt(m_matrix[i][i]);
		decomposed.setValue(i,i,sqrt_term);
		decomposed.setValue(i,m_cols-1, m_matrix[i][m_cols-1]);
		if( i < m_rows)
		{
			for( int c = i + 1; c < m_rows; c++)
			{
				decomposed.setValue(c,i,(m_matrix[c][i]/sqrt_term));
			}
			for( int x = i; x < m_rows; x++)
			{
				for( int y = i; y < m_rows; y++)
				{
					double vector_val_a = decomposed.getValue(x,i);
					double vector_val_b = decomposed.getValue(y,i);

					m_matrix[x][y] = m_matrix[x][y] - ( vector_val_a * vector_val_b);
				}
			}
		}
	}
	return decomposed;
}



void Matrix::gaussianElimination(std::string type)
{
	if( type == "no_pivot"){
		for(int i = 0; i < m_rows - 1; i++)
		{
			for( int j = i + 1; j < m_rows; j++)
			{
				double multiplier = m_matrix[j][i]/m_matrix[i][i];
				for( int k = 0; k < m_cols; k++)
				{
					m_matrix[j][k] = m_matrix[j][k] - ( multiplier * m_matrix[i][k] );
				}
				
			}
		}
	}
	else if(type == "partial_pivot")
	{
		constexpr double DOUBLE_MIN = std::numeric_limits<double>::lowest();

		for(int i = 0; i < m_rows - 1; i++)
		{
			int max_index = i;
			double max_pivot = DOUBLE_MIN;
			for(int k = i; k < m_rows; k++)
			{
				if(std::abs(m_matrix[k][i]) > max_pivot )
				{
					max_pivot = std::abs(m_matrix[k][i]);
					max_index = k;
				}
			}
			if(max_index != i)
			swapRows(i, max_index);
			for( int j = i + 1; j < m_rows; j++)
			{
				double multiplier = m_matrix[j][i]/m_matrix[i][i];
				for( int k = 0; k < m_cols; k++)
				{
					m_matrix[j][k] = m_matrix[j][k] - ( multiplier * m_matrix[i][k] );
				}
				
			}
		}
	}
	else if( type == "complete_pivot")
	{
		constexpr double DOUBLE_MIN = std::numeric_limits<double>::lowest();

		for(int i = 0; i < m_rows - 1; i++)
		{
			int max_row_index = i;
			int max_col_index = i;
			double max_pivot = DOUBLE_MIN;
			for(int k = i; k < m_rows; k++)
			{
				for( int j = i; j < m_rows; j++)
				{
					if(std::abs(m_matrix[k][i]) > max_pivot )
					{
						max_pivot = std::abs(m_matrix[k][i]);
						max_row_index = k;
						max_col_index = j;
					}
				}
			}
			if(max_row_index != i)
				swapRows(i, max_row_index);
			if(max_col_index != i)
				swapCols(i, max_col_index);
				
				
			for( int j = i + 1; j < m_rows; j++)
			{
				double multiplier = m_matrix[j][i]/m_matrix[i][i];
				for( int k = 0; k < m_cols; k++)
				{
					m_matrix[j][k] = m_matrix[j][k] - ( multiplier * m_matrix[i][k] );
				}
				
			}
		}		
	}
}


std::vector<double> Matrix::forwardSubstitution()
{
		
	if( m_cols != m_rows + 1)
		std::cerr<<" This matrix is not a linear system!";
		
	std::vector<double> solutions(m_rows,0);

	double b_1 = m_matrix[0][m_cols - 1];
	
	solutions[0] = b_1 / m_matrix[0][0];
	
	for(int i = 1; i < m_rows; i++)
	{
		double multiplier_i = 1 / m_matrix[i][i];
		double b_i = m_matrix[i][m_cols - 1];
		double sum = 0;
		for( int k = 0 ; k < i  ; k++)
		{
			sum += m_matrix[i][k] * solutions[k];
		}
		solutions[i] = multiplier_i * (b_i - sum);
	}
	return solutions;
	
}

std::vector<double> Matrix::backwardSubstitution()
{
	
	if( m_cols != m_rows + 1)
		std::cerr<<" This matrix is not a linear system!";
		
	std::vector<double> solutions(m_rows,0);

	double b_n = m_matrix[m_rows - 1][m_cols - 1];
	
	solutions[m_rows - 1] = b_n / m_matrix[m_rows - 1][m_rows - 1];
	
	for(int i = m_rows - 2; i >= 0; i--)
	{
		double multiplier_i = 1 / m_matrix[i][i];
		double b_i = m_matrix[i][m_cols - 1];
		double sum = 0;
		for( int k = i + 1 ; k < m_cols - 1 ; k++)
		{
			sum += m_matrix[i][k] * solutions[k];
		}
		solutions[i] = multiplier_i * (b_i - sum);
	}
	return solutions;
}



int Matrix::getHeight()
{
	return m_rows;
}

int Matrix::getLength()
{
	return m_cols;
}

void Matrix::setValue(int row, int col, double value)
{
	m_matrix[row][col] = value;
}


double Matrix::getValue(int row, int col)
{
	return m_matrix[row][col];
}


void Matrix::sumValue(int row, int col, double value)
{
	m_matrix[row][col] += value;
}





void Matrix::swapRows(int i, int j)
{
	std::vector<double> temp = m_matrix[i];
	m_matrix[i] = m_matrix[j];
	m_matrix[j] = temp;
}

void Matrix::swapCols(int i, int j)
{
	double temp;
	for( int k = 0; k < m_rows; k++)
	{
		temp = m_matrix[k][i];
		m_matrix[k][i] = m_matrix[k][j];
		m_matrix[k][j] = temp;
	}
}

void Matrix::substituteCol(int c, std::vector<double> col)
{
	for( int i = 0; i < m_rows; i++ )
	{
		m_matrix[i][c] = col[i];
	}
}



void Matrix::Print() const
{	
	for( int i = 0; i < m_rows; i++)
	{
		for( int j = 0; j < m_cols; j++)
		{
			PRINTL(std::setprecision(5)<<""<< m_matrix[i][j])
			PRINTL("     ");
		}
		LOG("")
	}
}