#include "Ising.hpp"


int main() {

	//First of all, we create a file to store the data we are going to obtain

	ofstream ofile;
	ofile.open("Ising_L60_ts.txt");  //That is an example name, but you can change it if you want to generate different data with this program
	ofile << scientific;


	//Some integers that we'll need several times

	int MC_cycles = 2000000;  //Number of MCMC cycles
	int N = 100;  //Number of temperatures that we will study
        int n = 100000;  //Burn in time(we will take into account the energies and magnetizations from this step on)


	//Then, we introduce the characteristics of the system

	int L = 60;  //Length of the lattice
	double N_size = 1.0 * (L * L);  //Number of elements inside the matrix that represents the lattice
	vec T = linspace(2.1, 2.4, N);  //Vector of temperatures


	//Set up for generating random numbers

	random_device rd;
	mt19937 mt(rd());


	//Now, we create some elements that we will fill up later

        mat S(L, L);  //Matrix representing the lattice

        int k, l;

        double e_value, m_value, e_mean, e2_mean, m_mean, m2_mean, Cv_value, X_value, q, dE, dM;


	//Then, we start the loop over the temperatures

	for(int j = 0; j < N; j++){

		e_mean = 0.0;
		e2_mean = 0.0;
		m_mean = 0.0;
		m2_mean = 0.0;
		Cv_value = 0.0;
		X_value = 0.0;
		dM = 0.0;


		//First, we create our system for the temperature of this step

		Ising my_system(T(j), L, N_size, mt);


		//We fill S with a random configuration (if true) or with all the spins being up (if false)

		my_system.create_matrix(S, true);


		//Now, let's start with the MCMC method

		for (int i = 0; i < MC_cycles; i++) {

			//We perform a cycle of MCMC

			my_system.MCMC(S, k, l, q, dE, dM);


			//We calulate here the energy for the first step after the burn-in time

			if (i == n){

				e_value = my_system.energy_spin(S);

				m_value = my_system.magnetization_spin(S);

				e_mean = e_value;

				m_mean = abs(m_value);

				e2_mean = (e_value * e_value);

				m2_mean = (m_value * m_value);

			}



			//After that, we update the energy and the magnetization in every step and we accumulate them to do the mean at the end

			if (i > n){

				e_value += q * (dE / N_size);  //e of this step
				m_value += q * (dM / N_size);  //m of this step

				e_mean += e_value;  //Sum of all the energies from the step n
				m_mean += abs(m_value);  //Sum of all the magnetizations from the step n

				e2_mean += (e_value * e_value);  //Same for e²
				m2_mean += (m_value * m_value);  //Same for m²

			}

		}


		//Finally, we calculate the means, Cv, X and store everything

		e_mean = e_mean / (MC_cycles - n);
		m_mean = m_mean / (MC_cycles - n);

		e2_mean = e2_mean / (MC_cycles - n);
                m2_mean = m2_mean / (MC_cycles - n);

		Cv_value = my_system.Cv(e_mean, e2_mean);

        	X_value = my_system.X(m_mean, m2_mean);

        	ofile << T(j) << "   " << e_mean << "   " << m_mean << "   " << Cv_value << "   " << X_value << endl;

	}

	ofile.close();

	return 0;

}
