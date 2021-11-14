#include "Ising.hpp"


int main() {

	//First of all, we create a file to store the data we are going to obtain

	ofstream ofile1;
	ofile1.open("Ising_L20_T1_histogram.txt");
	ofile1 << scientific;

	ofstream ofile2;
        ofile2.open("Ising_L20_T1_means.txt");
        ofile2 << scientific;





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

	double e_sum, e_mean, m_sum, m_mean;
        e_sum = 0;
	e_mean = 0;
	m_sum = 0;
	m_mean = 0;


	//We fill S with a random configuration (if true) or with all the spins being up (if false)

	my_system.create_matrix(S, false);




	//Once we have done this previous settings, let's start with the MCMC method. For that, we first set the number of MCMC cycles that we want to perform.

	int MC_cycles = 1000;

	for (int i = 0; i < MC_cycles; i++) {

		//We perform a cycle of MCMC

		my_system.MCMC(S, k, l);


		//Then we sum the energy and magnetization per spin of the new state to the one from the states we've sampled before that to do the mean later

		e_sum += my_system.energy_spin(S);

		m_sum += abs(my_system.magnetization_spin(S));


		//We store the energy and the magnetization per spin of the current state

		ofile1 << my_system.energy_spin(S) << "   " << my_system.magnetization_spin(S) << endl;

	}


	//Finally, we obtain the mean of both quantities and store them in another file

	e_mean = e_sum / MC_cycles;

	m_mean = m_sum / MC_cycles;

	ofile2 << e_mean << "   " << m_mean << endl;


	ofile1.close();

	ofile2.close();

	return 0;

}
