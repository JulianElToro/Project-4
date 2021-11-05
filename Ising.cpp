#include "Ising.hpp"


Ising::Ising(double T_in, double L_in) {

	//We assign the introduced values to the member variables

	T_ = T_in;
	L_ = L_in;

}



//Method that adds a particle to the trap by copying an input Particle

void Ising::matrix(bool random) {

	if (random) {

		mat S = randi(L_, L_, distr_param(0, 1));

		S.replace(0, -1);
	}
	
	else{

		mat S = ones(L_, L_);
	}

}

void Ising::flip_spin(mat  S) {

	

}