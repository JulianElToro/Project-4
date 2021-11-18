#include "Ising.hpp"
#include "omp.h"
#include <iomanip>
#include <assert.h>
#include <string>

using namespace std;


int main() {

	//Some integers that we'll need several times

	int MC_cycles = 750000;  //Number of MCMC cycles
	int N = 100;  //Number of temperatures in which we divide the interval of temperature
        int n = 15000;  //Number of data that we neglect from the beggining(due to burn-in time)


	//Then, we introduce the characteristics of the system

	int L = 80;  //Length of the lattice
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

        	double e_sum, e2_sum, e_mean, e2_mean, m_sum, m2_sum, m_mean, m2_mean, Cv, X;


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


				//We create our system

				Ising my_system(T(j), L, mt);


				//We fill S with a random configuration (if true) or with all the spins being up (if false)

				my_system.create_matrix(S, true);


				//Once we have done this previous settings, let's start with the MCMC method

				for (int i = 0; i < MC_cycles; i++){

					//We perform a cycle of MCMC

					my_system.MCMC(S, k, l);


					if (i >= n){

						//Then we sum the energy and magnetization per spin of the new state to the one from the states we've sampled before that(and also their squares)

						e_sum += my_system.energy_spin(S);
						e2_sum += (my_system.energy_spin(S)) * (my_system.energy_spin(S));

						m_sum += abs(my_system.magnetization_spin(S));
						m2_sum += (my_system.magnetization_spin(S)) * (my_system.magnetization_spin(S));

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
