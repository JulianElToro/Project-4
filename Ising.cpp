#include "Ising.hpp"


Ising::Ising(double T_in, int L_in) {

	//We assign the introduced values to the member variables

	T_ = T_in;
	L_ = L_in;

}



//Method that adds a particle to the trap by copying an input Particle

void Ising::create_matrix(imat& A ,bool random) {

	if (random) {

		imat S = randi(L_, L_, distr_param(0, 1));

		S.replace(0, -1);
	}
	
	else{

		imat S = ones(L_, L_);
	}

}

void Ising::flip_spin(mat S) {

	S( rand() % L_ , rand() % L_) = - S(rand() % L_, rand() % L_);	

}

void Ising::energy_spin(mat S_) {

	double E = 0.0;
	/*

	for (int j = 1; j < L_ + 2; j++) {
		
		for (int i = 1; i < L_ + 2; i++) {

			E = -(S_(i, j) * S_(i, j + 1) + S_(i, j) * S_(i, j - 1) + S_(i, j) * S_(i + 1, j) + S_(i, j) * S_(i - 1, j));

		}

		E+=
	}

	*/

	//double LS_ = S_.n_rows;

	for (int j = 1; j < L_ + 2; j++){
		
		for (int i = 0; i < L_ + 2; i++) {

			E = -( S_(i, j) * S_(i + 1, j) );

		}

		E +=
	}

	for (int j = 0; j < L_ + 2; j++) {

		for (int i = 1; i < L_ + 2; i++) {

			E = -(S_(i, j) * S_(i, j + 1));

		}

		E +=
	}


	double e = E / (L_ * L_);

}

void Ising::magnetization_spin(mat S) {

	double m = abs( accu(S) ) / (L_ * L_);

}

void Ising::Cv(mat S) {

	double e2 ;

	for (int j = 1; j < L_ + 2; j++) {

		for (int i = 1; i < L_ + 2; i++) {

			E = -(S_(i, j) * S_(i, j + 1) + S_(i, j) * S_(i, j - 1) + S_(i, j) * S_(i + 1, j) + S_(i, j) * S_(i - 1, j));

		}

	}

	double Cv = energy_spin(S)*energy_spin(S) / T_;


}