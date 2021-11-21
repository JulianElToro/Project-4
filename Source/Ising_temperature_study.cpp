#include "Ising.hpp"
#include "omp.h"
#include <iomanip>
#include <assert.h>
#include <string>

using namespace std;


int main() {

	//First of all, we create a file to store the data we are going to obtain

	ofstream ofile;
        ofile.open("Ising_L20_ts.txt");  //That is an example name, but you can change it if you want to generate different data with this program
        ofile << scientific;


	//Some integers that we'll need several times

	int MC_cycles = 2000000;  //Number of MCMC cycles
	int N = 100;  //Number of temperatures in which we divide the interval of temperature
	int n = 100000;  //Number of data that we neglect from the beggining(due to burn-in time)


	//Then, we introduce the characteristics of the system

	int L = 20;  //Length of the lattice
	double N_size = 1.0 * L * L;
	vec T = linspace(2.1, 2.4, N);  //Vector of temperatures


	//Set up for generating random numbers

	random_device rd;
	mt19937 mt(rd());



	//We will pass thorugh here if "-fopenmp" is including while compiling

	#ifdef _OPENMP
	{

		#pragma omp parallel
		{

			//We define some elements that we will fill out later

			mat S(L, L);  //Matrix representing the lattice

			int k, l;

			double e_value, m_value, e_mean, e2_mean, m_mean, m2_mean, Cv_value, X_value, q, dE, dM;




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

						e_value = my_system.energy_spin(S);

						m_value = my_system.magnetization_spin(S);

						e_mean = e_value;

						m_mean = abs(m_value);

						e2_mean = (e_value * e_value);

						m2_mean = (m_value * m_value);

					}


					//After that, we update the energy and the magnetization in every step and we accumulate them to do the mean at the end

					if (i > n) {

						e_value += q * (dE / N_size); //e of this step

						m_value += q * (dM / N_size); //m of this step

						e_mean += e_value; //Sum of all the energies from the step n

						m_mean += abs(m_value); //Sum of all the magnetizations from the step n

						e2_mean += (e_value * e_value); // Same for e²

						m2_mean += (m_value * m_value); //Same for m²

					}

				}


				//Finally, we calculate the means, Cv, X and store everything

				e_mean = e_mean / (MC_cycles - n);

				m_mean = m_mean / (MC_cycles - n);

				e2_mean = e2_mean / (MC_cycles - n);

				m2_mean = m2_mean / (MC_cycles - n);

				Cv_value = my_system.Cv(e_mean, e2_mean);

				X_value = my_system.X(m_mean, m2_mean);


				#pragma omp ordered  //For the file to be in the correct order

				ofile << T(j) << "   " << e_mean << "   " << m_mean << "   " << Cv_value << "   " << X_value << endl;

			}

		}


		ofile.close();

	}




	//If we don't include "-fopenmp" in the compiling, we will go through this non-parallelized region, where the steps are more or less the same that in the parallelized one

	#else
	{

        	mat S(L, L);

        	int k, l;

        	double e_value, m_value, e_mean, e2_mean, m_mean, m2_mean, Cv_value, X_value, q, dE, dM;



        	for(int j = 0; j < N; j++){

        	        e_mean = 0.0;
                	e2_mean = 0.0;
                	m_mean = 0.0;
                	m2_mean = 0.0;
                	Cv_value = 0.0;
                	X_value = 0.0;
                	dM = 0.0;


                	Ising my_system(T(j), L, N_size, mt);


                	my_system.create_matrix(S, true);


                	for (int i = 0; i < MC_cycles; i++) {

                        	my_system.MCMC(S, k, l, q, dE, dM);


                        	if (i == n){

                                	e_value = my_system.energy_spin(S);

                                	m_value = my_system.magnetization_spin(S);

                                	e_mean = e_value;

                                	m_mean = abs(m_value);

                                	e2_mean = (e_value * e_value);

                                	m2_mean = (m_value * m_value);

                        	}


                        	if (i > n){

                                	e_value += q * (dE / N_size);
                                	m_value += q * (dM / N_size);

                                	e_mean += e_value;
                                	m_mean += abs(m_value);

                                	e2_mean += (e_value * e_value);
                                	m2_mean += (m_value * m_value);

                        	}

                	}



                	e_mean = e_mean / (MC_cycles - n);
                	m_mean = m_mean / (MC_cycles - n);

                	e2_mean = e2_mean / (MC_cycles - n);
                	m2_mean = m2_mean / (MC_cycles - n);

                	Cv_value = my_system.Cv(e_mean, e2_mean);

                	X_value = my_system.X(m_mean, m2_mean);

                	ofile << T(j) << "   " << e_mean << "   " << m_mean << "   " << Cv_value << "   " << X_value << endl;

		}

		ofile.close();

	}

	#endif

	return 0;

}
