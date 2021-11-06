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

void Ising::boundary_conditions(mat S) {

	//First we extract the vectors in the borders

	rowvec up = S.row(0);
	rowvec down = S.row(L_ - 1);
	vec left = S.col(0);
	vec right = S.col(L_ - 1);


	up.insert_rows(0, 1);
	up.insert_rows(L_ + 2, 1);
	down.insert_rows(0, 1);
	down.insert_rows(L_2 + 2, 1);


	//Now we create a new matrix S_  inserting these vectors in t>

	mat S_ = S;

	S_.insert.cols(0, left);
	S_.insert.cols(L_ + 2, right);
	S_.insert_rows(0, up);
	S_.insert_rows(L_ + 2, down);

}

double Ising::energy_spin(mat S) {

	double E = 0.0;

	//meter aqui la función de Elena

	mat S_;
	
	for (int j = 1; j < L_ + 2; j++){
		
		for (int i = 0; i < L_ + 2; i++) {

			E += -( S_(i, j) * S_(i + 1, j) );

		}

	}

	for (int j = 0; j < L_ + 2; j++) {

		for (int i = 1; i < L_ + 2; i++) {

			E += -(S_(i, j) * S_(i, j + 1));

		}

	}


	double e = E / (L_ * L_);

}

double Ising::magnetization_spin(mat S) {

	double m = abs( accu(S) ) / (L_ * L_);

}

double Ising::Cv(mat S) {

	double E2;

	mat S_;

	for (int j = 1; j < L_ + 2; j++) {

		for (int i = 0; i < L_ + 2; i++) {

			E2 += (S_(i, j) * S_(i + 1, j))* (S_(i, j) * S_(i + 1, j));

		}

	}

	for (int j = 0; j < L_ + 2; j++) {

		for (int i = 1; i < L_ + 2; i++) {

			E2 += (S_(i, j) * S_(i, j + 1))*(S_(i, j) * S_(i, j + 1));

		}

	}


	double e2 = E2 / (L_ * L_);

	
	double Cv = ( e2 - energy_spin(S)*energy_spin(S) ) / T_;


}