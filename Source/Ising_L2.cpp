#include "Ising.hpp"

int main() {

        //First we create two files in order to take a look up our results

        ofstream ofile1;
        ofile1.open("MC_test_1000_MC_cycles.txt");
        ofile1 << scientific;

        ofstream ofile2;
        ofile2.open("MC_test_1000_MC_cycles2.txt");
        ofile2 << scientific;



        //Then, we set up the characteristics of our system

        int L = 2;  //Lenght of the lattice
        double T = 1.0;  //Temperature in J/kB


	//Set up for generating random numbers

        random_device rd;
        mt19937 mt(rd());


        //We create the system

        Ising my_system(T, L, mt);




        //After that, we create elements to storage important things later

        mat S(L, L);

        int k, l;

        double e_mean, m_mean, e2_mean, m2_mean, Cv_value, X_value;
        e_mean = 0.0;
        m_mean = 0.0;
        e2_mean = 0.0;
        m2_mean = 0.0;
        Cv_value = 0.0;
        X_value = 0.0;


        //We fill the system matrix with spins that cann be random (if true) or all of them pointing up (if false)

        my_system.create_matrix(S, true);




        //Now that we have done all the previous set up, we can start with the MCMC method. For that, we first set the number of cycles

        int MC_cycles = 1000;

        for (int i = 0; i < MC_cycles; i++) {

                //We do one time the MCMC, obtaining a new state (or keeping the one from before)

                my_system.MCMC(S, k, l);


                //Then, we calculate the energy and the magnetization per spin and their squares. We storage them in a file and we also sum t>

                e_mean += my_system.energy_spin(S);

                m_mean += abs(my_system.magnetization_spin(S));

                e2_mean += (my_system.energy_spin(S))*(my_system.energy_spin(S));

                m2_mean += (my_system.magnetization_spin(S))*(my_system.magnetization_spin(S));

                ofile1 << my_system.energy_spin(S) << "   " << my_system.magnetization_spin(S) << endl;

        }




	//After that, we calculate the mean of the energy and the magnetization per spin, and of their squares. Finally, we obtain Cv and X and we storage everything in another file

        e_mean = e_mean / MC_cycles;

        m_mean = m_mean / MC_cycles;

        e2_mean = e2_mean / MC_cycles;

        m2_mean = m2_mean / MC_cycles;

        Cv_value = my_system.Cv(e_mean, e2_mean);

        X_value = my_system.X(m_mean, m2_mean);

        ofile2 << e_mean << "    " << m_mean << "    " << e2_mean << "    " << m2_mean << endl;

        ofile1.close();

        ofile2.close();


        return 0;

}
