#include "Ising.hpp"
#include "omp.h"
#include <iomanip>
#include <assert.h>
#include <string>

using namespace std;


int main() {

	//Some integers that we'll need several times

	int MC_cycles = 1000000;  //Number of MCMC cycles
	int N = 100;  //Number of temperatures in which we divide the interval of temperature
	int n = 100000;  //Number of data that we neglect from the beggining(due to burn-in time)


	//Then, we introduce the characteristics of the system

	int L = 40;  //Length of the lattice
	double N_size = 1.0 * L * L;
	vec T = linspace(2.1, 2.4, N);  //Vector of temperatures


	//Set up for generating random numbers

	random_device rd;
	mt19937 mt(rd());


	/*//Common name to all the files

	const string output_file_name = "Ising_tsp.txt";*/


	ofstream ofile;
        ofile.open("Ising_L40_tsp.txt");
        ofile << scientific;




	//Here we start the parallelize area

	#pragma omp parallel
	{

		//Now, we define some elements that we will fill out later

		mat S(L, L);  //Matrix representing the lattice

		int k, l;

		double e_, m_, e_mean, e2_mean, m_mean, m2_mean, Cv_value, X_value, q, dE, dM;


		/*// Prepare for file output

		static int print_prec = 10;  //Number of decimal positions.(???)

		// Each thread will get its own output file name

		const int my_thread = omp_get_thread_num();
		ofstream ofile;
		string my_output_file_name = output_file_name + ".thread_" + to_string(my_thread);
		ofile.open(my_output_file_name.c_str(), ofstream::trunc);  // ofstream::trunc makes sure old content is deleted
		ofile << setprecision(print_prec) << scientific;*/



		//Now we start the paralellized loop over j

		#pragma omp for ordered

		for (int j = 0; j < N; j++) {

			//We set to zero the quantities from the last loop

			e_mean = 0.0;
			m_mean = 0.0;
			e2_mean = 0.0;
			m2_mean = 0.0;
			Cv_value = 0.0;
			X_value = 0.0;
			dM = 0.0;

			//We create our system

			Ising my_system(T(j), L, N_size, mt);


			//We fill S with a random configuration (if true) or with all the spins being up (if false)

			my_system.create_matrix(S, false);

			//Once we have done this previous settings, let's start with the MCMC method

			for (int i = 0; i < MC_cycles; i++) {

				//We perform a cycle of MCMC

				my_system.MCMC(S, k, l, q, dE, dM);


				//We calculate the energy for the step n(the first after the burning time)

				if (i == n) {

					e_ = my_system.energy_spin(S);

					m_ = my_system.magnetization_spin(S);

					e_mean = e_;

					m_mean = abs(m_);

					e2_mean = (e_ * e_);

					m2_mean = (m_ * m_);

				}


				//After that, we sum in every step the diference between the energy of the step before and the current one

				if (i > n) {

					//Then we sum the energy and magnetization per spin of the new state to the one from the states we've sampled before that(and also their squares)

					e_ += q * (dE / N_size); //e of this step

					m_ += q * (dM / N_size); //m of this step

					e_mean += e_; // Sum of all the energy from the step n

					m_mean += abs(m_); //Sum of all the magnetizations from the step n

					e2_mean += (e_ * e_); // Same for e²

					m2_mean += (m_ * m_); //Same for m²

				}

			}


			//Finally, we calculate the means, Cv, X and store everything

			e_mean = e_mean / (MC_cycles - n);

			m_mean = m_mean / (MC_cycles - n);

			e2_mean = e2_mean / (MC_cycles - n);

			m2_mean = m2_mean / (MC_cycles - n);

			Cv_value = my_system.Cv(e_mean, e2_mean);

			X_value = my_system.X(m_mean, m2_mean);



			//Finally, we introduce all this results to the file opened before

			#pragma omp ordered

			ofile << T(j) << "   " << e_mean << "   " << m_mean << "   " << Cv_value << "   " << X_value << endl;

		}

	}

	ofile.close();

	return 0;

}
