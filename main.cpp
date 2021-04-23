#include <iostream>
#include "Matrix.h"
#include <iomanip>
#include "Example.h"


#define LOG(x) std::cout<< x <<std::endl;


int main(int argc, char **argv)
{
	LOG("-------------------------------------")
	LOG("|         CHOLESKY EXAMPLE          |")
	LOG("-------------------------------------")


	choleskyExample();
	
	LOG("-------------------------------------")
	LOG("|         GAUSSIAN EXAMPLE          |")
	LOG("-------------------------------------")
	
	gaussianExample();

}
