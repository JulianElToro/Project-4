#include "Ising.hpp"

int main() {

        //First we create a file in order to take a look of our results

        ofstream ofile;
        ofile.open("Ising_L2_means_Cv_X.txt");
        ofile << scientific;



        //Then, we set up the characteristics of our system

        int L = 2;  //Lenght of the lattice
	double N_size = 1.0 * (L * L);  //Number of spins that are inside the lattice
        double T = 1.0;  //Temperature in J/kB


	//Set up for generating random numbers

        random_device rd;
        mt19937 mt(rd());


        //We create the system

        Ising my_system(T, L, N_size, mt);




        //After that, we create elements to storage important things later

        mat S(L, L);

        int k, l;

        double e_value, m_value, e_mean, e2_mean, m_mean, m2_mean, Cv_value, X_value, q, dE, dM;
        e_mean = 0.0;
        m_mean = 0.0;
	e2_mean = 0.0;
	m2_mean = 0.0;
        Cv_value = 0.0;
        X_value = 0.0;
	dM = 0.0;


        //We fill the system matrix with spins that cann be random (if true) or all of them pointing up (if false)

        my_system.create_matrix(S, true);


	//We calculate the energy and the magnetization of this first matrix

	e_value = my_system.energy_spin(S);

	m_value = my_system.magnetization_spin(S);




        //Now that we have done all the previous set up, we can start with the MCMC method. For that, we first set the number of cycles that  we will perform

        int MC_cycles = 1000000;

        for (int i = 0; i < MC_cycles; i++) {

                //We do one time the MCMC, obtaining a new state (or keeping the one from before)

                my_system.MCMC(S, k, l, q, dE, dM);


                //Then, we calculate the energy and the magnetization per spin of this step and we add those values to the previous ones in order to calculate the mean afterwards

		e_value += q * (dE / N_size);  //e of this step

                m_value += q * (dM / N_size);  //m of this step

		e_mean += e_value;  //Sum of all the energies from this step an the ones before

		m_mean += abs(m_value);  //Sum of all the magnetizations from this step and the ones before

                e2_mean += (e_value * e_value);  //Same with e²

                m2_mean += (m_value * m_value);  //Same with m²

        }




	//Finally we calculate the mean of the energy and the magnetization per spin, and of their squares. Finally, we obtain Cv and X and we storage everything in another file

        e_mean = e_mean / MC_cycles;

        m_mean = m_mean / MC_cycles;

	e2_mean = e2_mean / MC_cycles;

	m2_mean = m2_mean / MC_cycles;

        Cv_value = my_system.Cv(e_mean, e2_mean);

        X_value = my_system.X(m_mean, m2_mean);

        ofile << e_mean << "    " << m_mean << "    " << e2_mean << "    " << m2_mean << "    " << Cv_value << "    " << X_value << endl;
	
        ofile.close();


        return 0;

}
