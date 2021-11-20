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
        int n = 15000;  //Number of data that we neglect from the beggining(due to burn-in time)


	//Then, we introduce the characteristics of the system

	int L = 20;  //Length of the lattice
	vec T = linspace(2.1, 2.4, N);  //Vector of temperatures


	//Set up for generating random numbers

	random_device rd;
	mt19937 mt(rd());


	//Common name to all the files

	const string output_file_name = "Ising_tsp.txt";




	//Here we start the parallelize area

	#pragma omp parallel
	{

		//Now, we define some elements that we will fill out later

		mat S(L, L);  //Matrix representing the lattice

	        int k, l;

        	double e_, e2_, e_sum, e2_sum, e_mean, e2_mean, m_, m2_, m_sum, m2_sum, m_mean, m2_mean, Cv, X, q, dE;


		// Prepare for file output

    		static int print_prec = 10;  //Number of decimal positions.(???)

    		// Each thread will get its own output file name

    		const int my_thread = omp_get_thread_num();
    		ofstream ofile;
    		string my_output_file_name = output_file_name + ".thread_" + to_string(my_thread);
    		ofile.open(my_output_file_name.c_str(), ofstream::trunc);  // ofstream::trunc makes sure old content is deleted
		ofile << setprecision(print_prec) << scientific;


		//Now we start the paralellized loop over j

		#pragma omp for

			for(int j = 0; j < N; j++){

				//We set to zero the quantities from the last loop

				e_ = 0;
				e2_ = 0;
				e_sum = 0;
        			e2_sum = 0;
				e_mean = 0;
				e2_mean = 0;
				m_ = 0;
				m2_ = 0;
        			m_sum = 0;
        			m2_sum = 0;
				m_mean = 0;
				m2_mean = 0;
				Cv = 0;
				X = 0;


				//We create our system

				Ising my_system(T(j), L, mt);


				//We fill S with a random configuration (if true) or with all the spins being up (if false)

				my_system.create_matrix(S, true);


				//Once we have done this previous settings, let's start with the MCMC method

				for (int i = 0; i < MC_cycles; i++){

					//We perform a cycle of MCMC

					my_system.MCMC(S, k, l, q, dE);


					//We calculate the energy for the step n(the first after the burning time)

					if (i == n){

						e_ = my_system.energy_spin(S);
                                                e2_ = (my_system.energy_spin(S)) * (my_system.energy_spin(S));
						e_sum = e_;
						e2_sum = e2_;


						m_ = abs(my_system.magnetization_spin(S));
                                                m2_ = (my_system.magnetization_spin(S)) * (my_system.magnetization_spin(S));
						m_sum = m_;
						m2_sum = m2_;

					}


					//After that, we sum in every step the diference between the energy of the step before and the current one

					if (i > n){

        	                                double dM = S(k,l) * 2.0 / (L*L);

						//Then we sum the energy and magnetization per spin of the new state to the one from the states we've sampled before that(and also their squares)

						e_ += q * ( dE / (L*L) );  //e of this step

						e_sum += e_;  //Sum of all the energy from the step n
						e2_sum += (e_ * e_);  //Same for e²


                                                m_ += q * dM;  //m of this step

						m_sum += abs(m_);  //Sum of all the magnetizations from the step n
						m2_sum += (m_ * m_);  //Same for m²

					}

				}


				//Finally, we calculate the means, Cv, X and store everything

				e_mean = e_sum / (MC_cycles - n);
		                e2_mean = e2_sum / (MC_cycles - n);

                		m_mean = m_sum / (MC_cycles - n);
                		m2_mean = m2_sum / (MC_cycles - n);

                		Cv = my_system.Cv(e_mean, e2_mean);

                		X = my_system.X(m_mean, m2_mean);


                		//Finally, we introduce all this results to the file opened before

                		ofile << T(j) << "   " << e_mean << "   " << m_mean << "   " << Cv << "   " << X << endl;


			}

	}

	return 0;

}
