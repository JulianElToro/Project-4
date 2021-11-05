
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

	double T_;

	int L_;

	//Then, we declare also the constructor

	Ising() {}

	Ising(double T_in, int L_in);


	//Finally, we declare some methods for calculating some important things

	void create_matrix(imat& A, bool random);
	void flip_spin(mat S);
	void energy_spin(mat S);
	void magnetization_spin(mat S);

};


#endif