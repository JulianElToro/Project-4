#include "Ising.hpp"


int main() {

	//First of all, we create a file to store the data we are going to obtain

	ofstream ofile;
	ofile.open("Ising_L20_T1_histogram.txt");  //That is an example name, but with this program we have generated four different files
	ofile << scientific;




	//Then, we introduce the characteristics of the system

	int L = 20;  //Length of the lattice
	double N_size = 1.0 * (L * L);  //Number of spins in the lattice
	double T = 1.0;  //Temperature. This is something you can change for generating different data


	//Set up for generating random numbers

	random_device rd;
	mt19937 mt(rd());


	//Finally we create our system

	Ising my_system(T, L, N_size, mt);




	//Now, we create some elements that we will fill up later

	mat S(L, L);  //Matrix representing the lattice

	int k, l;

	double e_value, q, dE, dM;


	//We fill S with a random configuration (if true) or with all the spins being up (if false)

	my_system.create_matrix(S, true);




	//Once we have done this previous settings, let's start with the MCMC method. For that, we first set the number of MCMC cycles that we want to perform.

	int MC_cycles = 1000000;
	int n = 100000;  //Burn-in time

	for (int i = 0; i < MC_cycles; i++) {

		//We perform a cycle of MCMC

		my_system.MCMC(S, k, l, q, dE, dM);


		//We start storing the data only after the burn_in time

		if (i == n){

			//Here we calculate the energy of the first step after the burn-in time and store it

			e_value = my_system.energy_spin(S);

			ofile << e_value << endl;

		}

		if (i > n){

			//Here we just update the energy value and we store it

			e_value += q * (dE / N_size);

			ofile << e_value << endl;

		}

	}

	ofile.close();

	return 0;

}
