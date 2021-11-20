#include "Ising.hpp"


int main() {

	//First of all, we create a file to store the data we are going to obtain

	ofstream ofile1;
	ofile1.open("Ising_L20_T1_ordered.txt");  //That is an example name, but with this file we have generated four different files
	ofile1 << scientific;




	//Then, we introduce the characteristics of the system

	int L = 20;  //Length of the lattice
	double N = static_cast<double>(L) * static_cast<double>(L);
	double T = 1.0;  //Temperature. This is something you can change for generating different data


	//Set up for generating random numbers

	random_device rd;
	mt19937 mt(rd());


	//Finally we create our system

	Ising my_system(T, L, N, mt);




	//Now, we create some elements that we will fill up later

	mat S(L, L);

        int k, l;

        double e_, m_, e_sum, m_sum, e_mean, e2_mean, m_mean, m2_mean, Cv_value, X_value, q, dE, dM;
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

	e_ = my_system.energy_spin(S);

        m_ = my_system.magnetization_spin(S);




	//Once we have done this previous settings, let's start with the MCMC method. For that, we first set the number of MCMC cycles that we want to perform.

	int MC_cycles = 500000;

	for (int i = 0; i < MC_cycles; i++) {

		//We do one time the MCMC, obtaining a new state (or keeping the one from before)

                my_system.MCMC(S, k, l, q, dE, dM);


                //Then, we calculate the energy and the magnetization per spin and their squares. We storage them in a file and we also sum t>

                e_ += q * (dE / N);

                m_ += q * ( dM / N );

                e_sum += e_;

                m_sum += abs(m_);

		e_mean = e_sum / (1.0*i+1);

		m_mean = m_sum / (1.0*i+1);



		//Finally, we introduce the new mean of each thing to the file opened before

		ofile1 << e_mean << "   " << m_mean << endl;

	}

	ofile1.close();

	return 0;

}
