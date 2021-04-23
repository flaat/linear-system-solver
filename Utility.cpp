#include "Utility.h"
#include <iostream>
#include <iomanip>


#define LOG(x) std::cout<< x <<std::endl;

void printVector(std::vector<double> y)
{
	int len = y.size();
	
	for( int i = 0; i < len; i++)
	{
		LOG(std::setprecision(20)<<"x_"<<i<<" = "<<y[i])
	}
}