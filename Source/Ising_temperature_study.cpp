#include "Ising.hpp"


int main() {

	//First of all, we create a file to store the data we are going to obtain

	ofstream ofile1;
	ofile1.open("Ising_L40_temperature_study.txt");  //That is an example name, but with this file we have generated four different files
	ofile1 << scientific;


	//Some integers that we'll need several times

	int MC_cycles = 100000;
	int N = 10;
        int n = 15000;


	//Then, we introduce the characteristics of the system

	int L = 40;  //Length of the lattice
	vec T = linspace(2.1, 2.4, N);  //Vector of temperatures


	//Set up for generating random numbers

	random_device rd;
	mt19937 mt(rd());


	//Now, we create some elements that we will fill up later

        mat S(L, L);  //Matrix representing the lattice

        int k, l;

        double e_sum, e2_sum, e_mean, e2_mean, m_sum, m2_sum, m_mean, m2_mean, Cv, X;




	for(int j = 0; j < N; j++){

		e_sum = 0;
        	e2_sum = 0;
		e_mean = 0;
		e2_mean = 0;
        	m_sum = 0;
        	m2_sum = 0;
		m_mean = 0;
		m2_mean = 0;
		Cv = 0;
		X = 0;


		//Finally we create our system

		Ising my_system(T(j), L, mt);


		//We fill S with a random configuration (if true) or with all the spins being up (if false)

		my_system.create_matrix(S, true);


		//Once we have done this previous settings, let's start with the MCMC method. For that, we first set the number of MCMC cycles that we want to perform.

		for (int i = 0; i < MC_cycles; i++) {

			//We perform a cycle of MCMC

			my_system.MCMC(S, k, l);


			if (i >= n){


				//Then we sum the energy and magnetization per spin of the new state to the one from the states we've sampled before that and do the mean

				e_sum += my_system.energy_spin(S);
				e2_sum += (my_system.energy_spin(S)) * (my_system.energy_spin(S));

				m_sum += abs(my_system.magnetization_spin(S));
				m2_sum += (my_system.magnetization_spin(S)) * (my_system.magnetization_spin(S));


			}

		}


		e_mean = e_sum / (MC_cycles - n);
		e2_mean = e2_sum / (MC_cycles - n);

		m_mean = m_sum / (MC_cycles - n);
                m2_mean = m2_sum / (MC_cycles - n);

		Cv = my_system.Cv(e_mean, e2_mean);

        	X = my_system.X(m_mean, m2_mean);


		//Finally, we introduce all this results to the file opened before

        	ofile1 << T(j) << "   " << e_mean << "   " << m_mean << "   " << Cv << "   " << X << endl;

	}

	ofile1.close();

	return 0;

}
