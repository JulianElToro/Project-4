#include "Ising.hpp"


int main() {

	//First of all, we create a file to store the data we are going to obtain

	ofstream ofile;
	ofile.open("Ising_L20_T1_ordered.txt");  //That is an example name, but with this file we have generated four different files
	ofile << scientific;




	//Then, we introduce the characteristics of the system

	int L = 20;  //Length of the lattice
	double N_size = 1.0 * (L*L);
	double T = 1.0;  //Temperature. This is something you can change for generating different data


	//Set up for generating random numbers

	random_device rd;
	mt19937 mt(rd());


	//Finally we create our system

	Ising my_system(T, L, N_size, mt);




	//Now, we create some elements that we will fill up later

	mat S(L, L);

        int k, l;

        double e_value, m_value, e_sum, m_sum, e_mean, e2_mean, m_mean, m2_mean, Cv_value, X_value, q, dE, dM;
	e_sum = 0.0;
	m_sum = 0.0;
	e_mean = 0.0;
        m_mean = 0.0;
        e2_mean = 0.0;
        m2_mean = 0.0;
        Cv_value = 0.0;
        X_value = 0.0;
        dM = 0.0;



	//We fill S with a random configuration (if true) or with all the spins being up (if false)

	my_system.create_matrix(S, false);


	//We calculate the energy and the magnetization of this first matrix

	e_value = my_system.energy_spin(S);

        m_value = my_system.magnetization_spin(S);




	//Once we have done this previous settings, let's start with the MCMC method. For that, we first set the number of MCMC cycles that we want to perform

	int MC_cycles = 500000;

	for (int i = 0; i < MC_cycles; i++) {

		//We do one time the MCMC, obtaining a new state (or keeping the one from before)

                my_system.MCMC(S, k, l, q, dE, dM);


                //Then, we calculate the energy and the magnetization per spin of this step and we add those values to the previous ones in order to calculate the mean afterwards

                e_value += q * (dE / N_size);  //e of this step

                m_value += q * (dM / N_size);  //m of this step

                e_sum += e_value;  //Sum of all the energies from this step an the ones before

                m_sum += abs(m_value);  //Sum of all the magnetizations from this step an the ones before


		//Finally we calculate the means after the i+1 steps that we already did

		e_mean = e_sum / ( (1.0*i) + 1.0 );

		m_mean = m_sum / ( (1.0*i) + 1.0 );



		//We introduce the new mean of each thing to the file opened before

		ofile << e_mean << "   " << m_mean << endl;

	}

	ofile.close();

	return 0;

}
