
#ifndef __Ising_hpp__
#define __Ising_hpp__

#include <armadillo>
#include <iostream>
#include <cmath>
#include <fstream>


using namespace arma;
using namespace std;

class Ising {

public:

	//First we declare the member variables

	double T_ , L_;

	//Then, we declare also the constructor

	Ising() {}

	Ising(double T_in, double L_in);


	//Finally, we declare some methods for calculating some important things

	void Ising::matrix(bool random);
	
};


#endif