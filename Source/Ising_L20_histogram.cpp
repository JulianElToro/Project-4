#include "Ising.hpp"


int main() {

	//First of all, we create a file to store the data we are going to obtain

	ofstream ofile1;
	ofile1.open("Ising_L20_T1_histogram.txt");
	ofile1 << scientific;




	//Then, we introduce the characteristics of the system

	int L = 20;  //Length of the lattice
	double T = 1.0;  //Temperature. This is something you can change for generating different data


	//Set up for generating random numbers

	random_device rd;
	mt19937 mt(rd());


	//Finally we create our system

	Ising my_system(T, L, mt);




	//Now, we create some elements that we will fill up later

	mat S(L, L);  //Matrix representing the lattice

	int k, l;


	//We fill S with a random configuration (if true) or with all the spins being up (if false)

	my_system.create_matrix(S, true);




	//Once we have done this previous settings, let's start with the MCMC method. For that, we first set the number of MCMC cycles that we want to perform.

	int MC_cycles = 100000;

	for (int i = 0; i < MC_cycles; i++) {

		//We perform a cycle of MCMC

		my_system.MCMC(S, k, l);


		if (i >= 15000){

			//We store the energy and the magnetization per spin of the current state

			ofile1 << my_system.energy_spin(S) << endl;

		}

	}


	//Finally, we obtain the mean of both quantities and store them in another file

	ofile1.close();

	return 0;

}
