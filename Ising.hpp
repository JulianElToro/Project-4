
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

	void create_matrix(imat& S, bool random);
	void flip_spin(imat S, int& k, int& l);
	void boundary_conditions(imat& S_, imat S);
	double energy_spin(imat S_);
	double magnetization_spin(imat S_);
	double Cv(imat S_, double mean_e, double mean_e2);
	double X(double mean_m, double mean_m2);
	double acceptance(imat S_, int k, int l);
	void MCMC(imat S, imat& S_, int k, int l);

};


#endif